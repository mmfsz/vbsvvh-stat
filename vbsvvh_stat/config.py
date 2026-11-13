"""Configuration: the analysis-wide tables, a channel's channel.yaml and a result's YAML."""

import os
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parents[1]

#: Datacard era tags, chronological.
ERAS = ("2016preVFP", "2016postVFP", "2017", "2018",
        "2022", "2022EE", "2023", "2023BPix", "2024", "2025")
RUNS = {"Run2": ERAS[:4], "Run3": ERAS[4:]}
#: Centre-of-mass energy tag of a nuisance name -> the eras it is correlated across.
ENERGIES = {"13TeV": RUNS["Run2"], "13p6TeV": RUNS["Run3"]}
#: Prefix of the nuisances specific to this analysis (its CMS analysis number).
ANALYSIS_PREFIX = "CMS_HIG26013"
#: How a channel delivers a signal theory variation: ``raw`` (cross section included),
#: ``acceptance`` (divided by the generator-level ratio) or ``skim`` (divided by the
#: varied yield of a post-skim normalisation region, ``theory_norm``). See theory.py.
THEORY_MODES = ("raw", "acceptance", "skim")


class ConfigError(ValueError):
    pass


def load_yaml(path):
    with open(path) as f:
        return yaml.safe_load(f) or {}


def _str_keys(d):
    return {str(k): v for k, v in d.items()}


class Analysis:
    """The analysis-wide tables: nuisances.yaml, signals.yaml, lumi.yaml."""

    def __init__(self, root=REPO):
        self.nuisances = load_yaml(root / "nuisances.yaml")
        for name, entry in self.nuisances.items():
            unknown = set(entry or {}) - {"core", "description", "acceptance"}
            if unknown:
                raise ConfigError(f"nuisances.yaml: {name}: unknown key(s) {sorted(unknown)}")
        self.signals = load_yaml(root / "signals.yaml")
        lumi = load_yaml(root / "lumi.yaml")
        self.luminosity = _str_keys(lumi["luminosity"])
        self.year = {k: str(v) for k, v in _str_keys(lumi["year"]).items()}
        self.lumi = {name: {"kappa": _str_keys(entry["kappa"]),
                            "description": entry.get("description")}
                     for name, entry in lumi["uncertainties"].items()}
        self.master = root / "systematics" / "systematics_master.yml"
        from .theory import load_table
        self.theory = load_table(root / "theory_inclusive.yaml")

    def is_theory(self, name):
        """True for a signal theory nuisance whose normalisation the result decides."""
        entry = self.nuisance_entry(name)
        return bool(entry and (self.nuisances.get(entry) or {}).get("acceptance"))

    def nuisance_entry(self, name):
        """The nuisances.yaml key that ``name`` is an instance of, or None."""
        if name in self.nuisances:
            return name
        for key in self.nuisances:
            if "{era}" in key and any(key.format(era=e) == name for e in ERAS):
                return key
        return None

    def lumi_kappas(self, era):
        """``{nuisance: kappa}`` acting on ``era``."""
        year = self.year[era]
        return {n: e["kappa"][year] for n, e in self.lumi.items() if year in e["kappa"]}


def load_channel(path):
    """channel.yaml -> dict with ``tag``, ``shapes`` (absolute), ``abcd``, ``provenance``."""
    path = Path(path).resolve()
    raw = load_yaml(path)
    unknown = set(raw) - {"channel", "abcd", "provenance", "theory", "theory_norm"}
    if unknown:
        raise ConfigError(f"{path}: unknown key(s) {sorted(unknown)}")
    tag = raw.get("channel")
    if not tag or not str(tag).isidentifier() or "__" in str(tag):
        raise ConfigError(f"{path}: 'channel' must be a tag like MET1FJ (letters, digits, single underscores)")
    abcd = raw.get("abcd")
    if abcd is not None:
        missing = [k for k in ("background", "A", "B", "C", "D") if k not in abcd]
        if missing:
            raise ConfigError(f"{path}: abcd is missing {missing}")
        if len({abcd[r] for r in "ABCD"}) != 4:
            raise ConfigError(f"{path}: abcd regions must be four distinct regions")
    theory = raw.get("theory")
    if theory is not None and theory not in THEORY_MODES:
        raise ConfigError(f"{path}: 'theory' must be one of {list(THEORY_MODES)}, not {theory!r}")
    norm = raw.get("theory_norm")
    if theory == "skim" and not norm:
        raise ConfigError(f"{path}: 'theory: skim' needs 'theory_norm', the region holding the "
                          "post-skim yield the variations were divided by")
    if norm and theory != "skim":
        raise ConfigError(f"{path}: 'theory_norm' is only for 'theory: skim'")
    if norm and abcd and norm in abcd.values():
        raise ConfigError(f"{path}: 'theory_norm' {norm!r} is an ABCD region")
    shapes = path.parent / "shapes.root"
    if not shapes.exists():
        raise ConfigError(f"{path}: no shapes.root beside it")
    return {"tag": tag, "shapes": shapes, "abcd": abcd, "theory": theory, "theory_norm": norm,
            "provenance": raw.get("provenance")}


def load_result(path, analysis):
    """A result YAML, with ``runs``/``fits`` left for :func:`resolve_runs` once eras are known."""
    path = Path(path).resolve()
    raw = load_yaml(path)
    unknown = set(raw) - {"name", "channels", "signals", "runs", "fits", "impacts", "closure", "unblind",
                          "theory"}
    if unknown:
        raise ConfigError(f"{path}: unknown key(s) {sorted(unknown)}")
    for key in ("name", "channels", "signals"):
        if key not in raw:
            raise ConfigError(f"{path}: '{key}' is required")
    if raw.get("unblind", False) not in (True, False):
        raise ConfigError(f"{path}: 'unblind' must be true or false")
    theory = raw.get("theory", "acceptance")
    if theory not in ("acceptance", "total"):
        raise ConfigError(f"{path}: 'theory' must be acceptance or total, not {theory!r}")
    bad = [s for s in raw["signals"] if s not in analysis.signals]
    if bad:
        raise ConfigError(f"{path}: signals {bad} not in signals.yaml")
    channels = {}
    for tag, entry in raw["channels"].items():
        given = entry["input"] if isinstance(entry["input"], list) else [entry["input"]]
        inputs = [load_channel(Path(path.parent, i)) for i in given]
        for ch in inputs:
            if ch["tag"] != tag:
                raise ConfigError(f"{path}: channel {tag!r} points to a channel.yaml with tag {ch['tag']!r}")
            if ch["abcd"] != inputs[0]["abcd"]:
                raise ConfigError(f"{path}: the inputs of channel {tag!r} differ in their abcd block")
            if (ch["theory"], ch["theory_norm"]) != (inputs[0]["theory"], inputs[0]["theory_norm"]):
                raise ConfigError(f"{path}: the inputs of channel {tag!r} differ in their theory declaration")
        channels[tag] = {"tag": tag, "abcd": inputs[0]["abcd"], "theory": inputs[0]["theory"],
                         "theory_norm": inputs[0]["theory_norm"],
                         "inputs": inputs,
                         "shapes": [ch["shapes"] for ch in inputs]}
    name = raw["name"]
    results_dir = Path(os.environ.get("VBSVVH_STAT_RESULTS", REPO / "results")) / name
    return {
        "name": name,
        "path": path,
        "raw": raw,
        "channels": channels,
        "signals": list(raw["signals"]),
        "runs": raw.get("runs") or {},
        "fits": raw.get("fits"),
        "impacts": raw.get("impacts") or {},
        "closure": raw.get("closure") or {},
        "unblind": raw.get("unblind", False),
        "theory": theory,
        "results_dir": results_dir,
        "work_dir": results_dir / "work",
    }


def resolve_runs(overrides, eras_present):
    """Planes of the result: ``{plane: {"run", "eras", "substitute"}}``, and the run -> planes map.

    By default every run's eras that are present form one plane named after the run.
    An override entry without ``eras`` modifies a default run (``substitute``,
    ``split``); one with ``eras`` defines its own run, and those eras leave the
    default runs.
    """
    explicit = {str(e) for o in overrides.values() if "eras" in o for e in o["eras"]}
    runs = {}
    for run, eras in RUNS.items():
        o = overrides.get(run, {})
        if "eras" in o:
            continue
        subst = {str(k): str(v) for k, v in (o.get("substitute") or {}).items()}
        members = [e for e in eras if e not in explicit and (e in eras_present or e in subst)]
        if members:
            runs[run] = dict(o, eras=members, substitute=subst)
    for run, o in overrides.items():
        if "eras" in o:
            subst = {str(k): str(v) for k, v in (o.get("substitute") or {}).items()}
            runs[run] = dict(o, eras=[str(e) for e in o["eras"]], substitute=subst)
        elif run not in RUNS:
            raise ConfigError(f"runs: {run!r} is not a run (Run2, Run3) and has no 'eras'")
    for run, o in runs.items():
        for target, source in o["substitute"].items():
            if source not in eras_present:
                raise ConfigError(f"runs.{run}.substitute: source era {source} not in the inputs")
            if target not in o["eras"]:
                o["eras"].append(target)

    planes, members = {}, {}
    for run, o in runs.items():
        if o.get("split"):
            names = [f"{run}_{e}" for e in o["eras"]]
            for e, plane in zip(o["eras"], names):
                planes[plane] = {"run": run, "eras": [e], "substitute": o["substitute"]}
        else:
            names = [run]
            planes[run] = {"run": run, "eras": o["eras"], "substitute": o["substitute"]}
        members[run] = names
    return planes, members


def resolve_fits(fits, run_names):
    """``{fit: [run, ...]}``: the given ones, or each run plus all runs together."""
    if fits:
        for fit, runs in fits.items():
            bad = [r for r in runs if r not in run_names]
            if bad:
                raise ConfigError(f"fits.{fit}: unknown run(s) {bad}")
        return {f: list(r) for f, r in fits.items()}
    out = {r: [r] for r in run_names}
    if len(run_names) > 1:
        out["".join(run_names)] = list(run_names)
    return out
