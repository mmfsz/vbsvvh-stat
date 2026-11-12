"""Datacards: one per (channel, plane, signal).

A plane is the set of eras summed into one set of bins, by default one run
(see config.resolve_runs). For each channel and plane:

1. eras replaced by a stand-in (``substitute``) take the source era's histograms,
   scaled to the target era's luminosity, with nuisances tagged ``_<source>``
   renamed ``_<target>``;
2. luminosity nuisances (lumi.yaml) are attached per era to signal and MC
   processes, as Up/Down variations;
3. signal theory variations (``acceptance: true`` in nuisances.yaml) are brought
   from the normalisation the channel delivered to the one the result asks for
   (``theory: acceptance | total``), with the generator-level ratios of
   theory_inclusive.yaml (see theory.py);
4. the eras are summed: nominal and Sumw2 add, and each nuisance is summed
   varied where an era carries it and nominal where it does not. On a summed
   bin this is the diluted kappa of an era-tagged nuisance,
   ``1 + (Y_era / Y_tot) (k_era - 1)``, and the summed ratio of a correlated one.

A card is ``counting`` (numbers, lnN, one MC statistical lnN
``CMS_HIG26013_<tag>_mcstat_<region>_<plane>_<process>`` per process and bin) when
every region has one bin, else ``shape`` (templates in a ROOT file beside the
cards, shape nuisances, autoMCStats).

Names carry the channel tag so nothing is shared between channels by accident:
bins ``<tag>_<region>_<plane>``, backgrounds ``<tag>_<process>``. Signals keep
their analysis-wide names, so a hypothesis is one process in every channel.

ABCD (counting regions only): the data-driven background has rate 1 and
absolute-yield rateParams ``CMS_HIG26013_<tag>_<plane>_abcd_{B,C,D}`` seeded with the
observed count N in [N - 3 sqrt N, N + 3 sqrt N], A = B*C/D, and the
non-closure lnN ``CMS_HIG26013_<tag>_<plane>_abcd_closure`` on A: 1 + sqrt(delta^2 +
error^2) from the MC closure and its MC statistical error, or 1 + delta for a value
given in the result.
"""

import json
import math

import numpy as np

from .config import ANALYSIS_PREFIX, ConfigError, resolve_fits, resolve_runs
from .inputs import DATA, Hist, Shapes, add, background_kind, scale, validate
from .naming import write_declaration
from .theory import factors

#: A variation within this of 1 in both directions is not written.
NULL_KAPPA = 1e-4
#: Floor for a varied/nominal ratio, which lnN cannot take as zero.
MIN_RATIO = 1e-6


class Plane:
    """One channel's histograms summed over the eras of one plane."""

    def __init__(self, channel, shapes, name, cfg, analysis, theory="acceptance"):
        self.channel = channel
        self.theory = theory
        self.tag = channel["tag"]
        self.name = name
        self.analysis = analysis
        sources = {e: cfg["substitute"].get(e, e) for e in cfg["eras"]}
        self.eras = [e for e, s in sources.items() if s in shapes.eras]
        self.substitute = {e: s for e, s in sources.items() if e != s and e in self.eras}
        per_era = [self._era(shapes, e, sources[e]) for e in self.eras]

        self.nominal, self.varied = {}, {}
        for hists in per_era:
            for (r, p, n, d), h in hists.items():
                if n is None:
                    self.nominal[(r, p)] = add(self.nominal.get((r, p)), h)
        names = {(r, p, n) for hists in per_era for (r, p, n, d) in hists if n}
        for r, p, n in names:
            up = dn = None
            for hists in per_era:
                nom = hists.get((r, p, None, None))
                if nom is None:
                    continue
                up = add(up, hists.get((r, p, n, "Up"), nom))
                dn = add(dn, hists.get((r, p, n, "Down"), nom))
            self.varied[(r, p, n)] = (up, dn)

    def _era(self, shapes, era, source):
        factor = 1.0
        if source != era:
            lumi = self.analysis.luminosity
            factor = lumi[era] / lumi[source]
        out = {}
        for (e, r, p, n, d), h in shapes.hists.items():
            if e != source:
                continue
            if n and source != era and n.endswith("_" + source):
                n = n[: -len(source)] + era
            out[(r, p, n, d)] = scale(h, factor)
        for (r, p, n, d), h in list(out.items()):
            if n and self.kind(p) == "signal" and self.analysis.is_theory(n):
                f = factors(self.analysis.theory, self.analysis.nuisance_entry(n), p, source,
                            self.channel["theory"], self.theory)
                if f is not None:
                    out[(r, p, n, d)] = scale(h, f[0] if d == "Up" else f[1])
        for (r, p, n, d), h in list(out.items()):
            if n is None and self.kind(p) in ("signal", "mc"):
                for lumi_n, k in self.analysis.lumi_kappas(era).items():
                    out[(r, p, lumi_n, "Up")] = Hist(h.values * k, h.variances, h.edges)
                    out[(r, p, lumi_n, "Down")] = Hist(h.values / k, h.variances, h.edges)
        return out

    def kind(self, process):
        return background_kind(self.channel, process, self.analysis)

    @property
    def regions(self):
        regs = sorted({r for r, _ in self.nominal})
        abcd = self.channel["abcd"]
        first = [abcd[x] for x in "ABCD"] if abcd else []
        return first + [r for r in regs if r not in first]

    def bin(self, region):
        return f"{self.tag}_{region}_{self.name}"

    def process_name(self, process):
        return process if self.kind(process) == "signal" else f"{self.tag}_{process}"

    def integral(self, region, process):
        h = self.nominal.get((region, process))
        return (float(h.values.sum()), float(h.variances.sum())) if h else (0.0, 0.0)

    @property
    def shape_mode(self):
        return any(len(h.values) > 1 for h in self.nominal.values())

    def backgrounds(self):
        """MC backgrounds present, and the ABCD background, which needs no histogram."""
        procs = {p for _, p in self.nominal if self.kind(p) == "mc"}
        if self.channel["abcd"]:
            procs.add(self.channel["abcd"]["background"])
        return sorted(procs)

    def columns(self, signal):
        """``[(region, process)]`` written in the card for ``signal``."""
        abcd = self.channel["abcd"]
        abcd_regions = {abcd[x] for x in "ABCD"} if abcd else set()
        cols = []
        for r in self.regions:
            for p in [signal] + self.backgrounds():
                if self.kind(p) == "data_driven":
                    if r in abcd_regions:
                        cols.append((r, p))
                elif self.integral(r, p)[0] > 0:
                    cols.append((r, p))
        return cols

    def closure(self, override):
        """``{"delta", "error", "method", "kappa"}`` of the ABCD non-closure.

        The MC statistical error enters the kappa in quadrature, so a closure test
        limited by MC statistics cannot give a smaller uncertainty than it resolves.
        """
        if override is not None:
            if isinstance(override, bool) or not isinstance(override, (int, float)):
                raise ConfigError(f"closure.{self.tag}.{self.name}: the non-closure must be a number, "
                                  f"not {override!r}")
            return {"delta": float(override), "error": None, "method": "config", "kappa": 1.0 + override}
        abcd = self.channel["abcd"]
        n = {x: self.integral(abcd[x], abcd["background"]) for x in "ABCD"}
        if any(n[x][0] <= 0 for x in "ABCD"):
            raise ConfigError(f"{self.tag}/{self.name}: MC closure needs MC background in A, B, C, D; "
                              "give the closure in the result's 'closure'")
        c = n["A"][0] * n["D"][0] / (n["B"][0] * n["C"][0])
        err = c * math.sqrt(sum(n[x][1] / n[x][0] ** 2 for x in "ABCD"))
        delta = abs(1.0 - c)
        return {"delta": delta, "error": err, "method": "mc", "kappa": 1.0 + math.hypot(delta, err)}


def _fmt(x):
    return f"{x:.8g}"


def _row(label, cells, kind="", table="process"):
    """A table row, laid out by :func:`_layout` once every row of the card is known."""
    return (table, label, kind, [str(c) for c in cells])


def _layout(lines):
    """Card text: each table's columns as wide as their widest cell, two spaces apart."""
    rows = [line for line in lines if isinstance(line, tuple)]
    label_w = max(len(r[1]) for r in rows)
    kind_w = max(len(r[2]) for r in rows)
    widths = {}
    for table, _, _, cells in rows:
        w = widths.setdefault(table, [0] * len(cells))
        widths[table] = [max(a, len(c)) for a, c in zip(w, cells)]
    out = []
    for line in lines:
        if isinstance(line, str):
            out.append(line)
            continue
        table, label, kind, cells = line
        body = "  ".join(f"{c:<{w}}" for c, w in zip(cells, widths[table]))
        out.append(f"{label:<{label_w}}  {kind:<{kind_w}}  {body}".rstrip())
    return "\n".join(out) + "\n"


def _kappa(plane, region, process, nuisance):
    nom = plane.integral(region, process)[0]
    pair = plane.varied.get((region, process, nuisance))
    if pair is None or nom <= 0:
        return None
    up, dn = (max(h.values.sum() / nom, MIN_RATIO) for h in pair)
    if abs(up - 1) < NULL_KAPPA and abs(dn - 1) < NULL_KAPPA:
        return None
    return f"{_fmt(dn)}/{_fmt(up)}"


def card_text(plane, signal, closure, shapes_file=None, blind=True):
    """The datacard for ``signal`` on ``plane``. ``shapes_file`` is given in shape mode.

    A blinded counting card has no observation line: its workspace holds no data, so
    the -t -1 fits are unchanged and an observed fit stops for lack of data. A shape
    card keeps its data (``observation -1``, read from the shapes file).
    """
    cols = plane.columns(signal)
    bins = [plane.bin(r) for r in plane.regions]
    procs = [signal] + [p for p in plane.backgrounds() if any(p == q for _, q in cols)]
    index = {p: i for i, p in enumerate(procs)}
    abcd = plane.channel["abcd"]
    shape = shapes_file is not None

    lines = [f"# channel {plane.tag}, plane {plane.name}, eras {', '.join(plane.eras)}, signal {signal}"]
    for target, source in plane.substitute.items():
        lines.append(f"# era {target}: stand-in from {source} scaled by "
                     f"{_fmt(plane.analysis.luminosity[target] / plane.analysis.luminosity[source])}")
    if abcd:
        if closure["error"] is None:
            lines.append(f"# ABCD non-closure delta = {_fmt(closure['delta'])} ({closure['method']}) "
                         f"-> lnN {_fmt(closure['kappa'])}")
        else:
            lines.append(f"# ABCD non-closure delta = {_fmt(closure['delta'])}, MC stat error = "
                         f"{_fmt(closure['error'])} ({closure['method']}) -> lnN 1 + sqrt(delta^2 + error^2) "
                         f"= {_fmt(closure['kappa'])}")
    if plane.channel.get("theory"):
        lines.append(f"# signal theory variations: delivered {plane.channel['theory']}, "
                     f"written {plane.theory} (generator-level ratios, theory_inclusive.yaml)")
    blind = blind and not shape
    if blind:
        lines.append("# BLINDED: no observed data; only -t -1 fits are possible")
    lines += ["imax *", "jmax *", "kmax *", "-" * 80]
    if shape:
        lines += [f"shapes * * {shapes_file} $CHANNEL/$PROCESS $CHANNEL/$PROCESS__$SYSTEMATIC", "-" * 80]
    if not blind:
        obs = ["-1" if shape else _fmt(plane.integral(r, DATA)[0]) for r in plane.regions]
        lines += [_row("bin", bins, table="bin"), _row("observation", obs, table="bin"), "-" * 80]
    lines += [
        _row("bin", [plane.bin(r) for r, _ in cols]),
        _row("process", [plane.process_name(p) for _, p in cols]),
        _row("process", [str(index[p]) for _, p in cols]),
        _row("rate", ["-1" if shape else ("1" if plane.kind(p) == "data_driven" else _fmt(plane.integral(r, p)[0]))
                      for r, p in cols]),
        "-" * 80,
    ]

    written = set(cols)
    nuisances = sorted({n for (r, p, n) in plane.varied if (r, p) in written and plane.kind(p) != "data_driven"})
    for n in nuisances:
        cells = []
        for r, p in cols:
            k = None if plane.kind(p) == "data_driven" else _kappa(plane, r, p, n)
            cells.append("-" if k is None else ("1" if shape else k))
        if any(c != "-" for c in cells):
            lines.append(_row(n, cells, "shape" if shape else "lnN"))
    if not shape:
        for r, p in cols:
            y, var = plane.integral(r, p)
            if plane.kind(p) == "data_driven" or y <= 0 or var <= 0:
                continue
            cells = [_fmt(1 + math.sqrt(var) / y) if (rr, pp) == (r, p) else "-" for rr, pp in cols]
            lines.append(_row(f"{ANALYSIS_PREFIX}_{plane.tag}_mcstat_{r}_{plane.name}_{p}", cells, "lnN"))
    if abcd:
        bkg = abcd["background"]
        cells = [_fmt(closure["kappa"]) if (r, p) == (abcd["A"], bkg) else "-" for r, p in cols]
        lines.append(_row(f"{ANALYSIS_PREFIX}_{plane.tag}_{plane.name}_abcd_closure", cells, "lnN"))
        lines.append("-" * 80)
        params = {x: f"{ANALYSIS_PREFIX}_{plane.tag}_{plane.name}_abcd_{x}" for x in "ABCD"}
        bkg_name = plane.process_name(bkg)
        for x in "BCD":
            n = plane.integral(abcd[x], DATA)[0]
            lo, hi = max(0.0, n - 3 * math.sqrt(n)), n + 3 * math.sqrt(n)
            lines.append(f"{params[x]} rateParam {plane.bin(abcd[x])} {bkg_name} {_fmt(n)} [{_fmt(lo)},{_fmt(hi)}]")
        lines.append(f"{params['A']} rateParam {plane.bin(abcd['A'])} {bkg_name} (@0*@1/@2) "
                     f"{params['B']},{params['C']},{params['D']}")
    if shape:
        lines += [f"{b} autoMCStats 10 0 1" for b in bins]
    return _layout(lines)


def write_shapes(plane, path):
    """The card-local ROOT file of a shape-mode plane, with tagged names."""
    import hist
    import uproot

    def th1(h):
        out = hist.Hist(hist.axis.Variable(h.edges), storage=hist.storage.Weight())
        out.view()["value"], out.view()["variance"] = h.values, h.variances
        return out

    abcd = plane.channel["abcd"]
    with uproot.recreate(path) as f:
        if abcd:
            unit = Hist(np.ones(1), np.zeros(1), np.array([0.0, 1.0]))
            for x in "ABCD":
                f[f"{plane.bin(abcd[x])}/{plane.process_name(abcd['background'])}"] = th1(unit)
        for (r, p), h in plane.nominal.items():
            if plane.kind(p) == "data_driven":
                continue
            name = DATA if p == DATA else plane.process_name(p)
            f[f"{plane.bin(r)}/{name}"] = th1(h)
        for (r, p, n), (up, dn) in plane.varied.items():
            if plane.kind(p) in ("signal", "mc"):
                f[f"{plane.bin(r)}/{plane.process_name(p)}__{n}Up"] = th1(up)
                f[f"{plane.bin(r)}/{plane.process_name(p)}__{n}Down"] = th1(dn)


def _check_common_nuisances(shapes, channels, analysis):
    """A core nuisance carried by one channel must be carried by all."""
    carried = {}
    for tag, sh in shapes.items():
        ch = channels[tag]
        procs = [p for p in sh.processes() if background_kind(ch, p, analysis) in ("signal", "mc")]
        carried[tag] = {analysis.nuisance_entry(n) for p in procs for n in sh.nuisances(process=p)}
    core = {e for e in set().union(*carried.values()) - {None}
            if (analysis.nuisances.get(e) or {}).get("core")}
    return [f"{tag} does not carry {entry}, a core nuisance that other channels carry"
            for entry in sorted(core) for tag in carried if entry not in carried[tag]]


def write_cards(result, analysis, dictionary, log=print):
    """Validate the channels, write every card. Returns the card index for the fits."""
    shapes = {tag: Shapes(ch["shapes"]) for tag, ch in result["channels"].items()}
    errors = []
    for tag, ch in result["channels"].items():
        errs, warns, _ = validate(ch, analysis, shapes[tag])
        errors += [f"{tag}: {e}" for e in errs]
        for w in warns:
            log(f"warning: {tag}: {w}")
    if len(shapes) > 1:
        errors += _check_common_nuisances(shapes, result["channels"], analysis)
    if errors:
        raise ConfigError("inputs do not pass validation:\n  " + "\n  ".join(errors))

    eras_present = {e for sh in shapes.values() for e in sh.eras}
    planes, members = resolve_runs(result["runs"], eras_present)
    fits = resolve_fits(result["fits"], list(members))
    used = {p for runs in fits.values() for r in runs for p in members[r]}
    planes = {p: cfg for p, cfg in planes.items() if p in used}
    members = {r: ps for r, ps in members.items() if set(ps) & used}
    missing = sorted({e for p in planes.values() for e in p["eras"]} - set(analysis.luminosity))
    if missing:
        raise ConfigError(f"era(s) {missing} not implemented yet (no luminosity in lumi.yaml)")

    cards_dir = result["results_dir"] / "cards"
    cards_dir.mkdir(parents=True, exist_ok=True)
    index = {s: {} for s in result["signals"]}
    summary, nuisances = {}, set()
    for tag, ch in result["channels"].items():
        for name, cfg in planes.items():
            plane = Plane(ch, shapes[tag], name, cfg, analysis, result["theory"])
            if not plane.eras:
                continue
            closure = None
            if ch["abcd"]:
                override = result["closure"].get(tag, {})
                closure = plane.closure(override.get(name, override.get(cfg["run"])))
            root = None
            if plane.shape_mode:
                root = cards_dir / f"{tag}_{name}.root"
                write_shapes(plane, root)
            summary.setdefault(tag, {})[name] = {
                "eras": plane.eras, "substitute": plane.substitute,
                "mode": "shape" if root else "counting", "closure": closure,
                "theory": {"delivered": ch["theory"], "written": result["theory"]},
                "blinded": not result["unblind"] and not root,
            }
            for signal in result["signals"]:
                if not any(plane.integral(r, signal)[0] > 0 for r in plane.regions):
                    continue
                text = card_text(plane, signal, closure, str(root) if root else None,
                                 blind=not result["unblind"])
                path = cards_dir / signal / f"{tag}_{name}.txt"
                path.parent.mkdir(exist_ok=True)
                path.write_text(text)
                index[signal].setdefault(name, []).append(path)
                nuisances |= {line.split()[0] for line in text.splitlines()
                              if len(line.split()) > 1 and line.split()[1] in ("lnN", "shape")}
                log(f"card {path.relative_to(result['results_dir'])}")

    write_declaration(cards_dir / "systematics.yml", nuisances, analysis, dictionary)
    with open(cards_dir / "summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    return {"planes": planes, "members": members, "fits": fits, "cards": index, "summary": summary}
