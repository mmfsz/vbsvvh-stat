"""Generator-level (inclusive) ratios of the signal theory variations.

A theory variation of the signal (QCD scales, parton shower) moves both the cross
section and the acceptance. Which of the two the card carries is decided once, by the
result (``theory: acceptance | total``), not by each channel:

* a channel declares in channel.yaml whether its variations are ``raw`` (varied yield
  over nominal yield, cross section included) or already ``acceptance``-only
  (each divided by the same variation's generator-level ratio);
* the card writer divides or multiplies by the generator-level ratio of this table
  wherever the two differ (:func:`factors`).

The table, ``theory_inclusive.yaml``, holds per nuisance, signal hypothesis and era
the ratio of the varied to the nominal generator-level cross section, from the
``Runs`` tree of the NanoAOD skims (``LHEScaleSumw``, ``PSSumw``). A hypothesis is
the sum of the four VBS VVH processes, so its ratio is their cross-section-weighted
mean. Build it with ``vbsvvh-stat theory-table``.
"""

import re
from pathlib import Path

import yaml

from .config import ConfigError

#: The processes summed into one signal hypothesis.
SIGNAL_PROCESSES = ("VBSWWH_OS", "VBSWWH_SS", "VBSWZH", "VBSZZH")

#: nuisance -> (Runs branch, index of Up, index of Down, index of the nominal or None).
#: NanoAOD LHEScaleWeight: [muR, muF] = 0:[.5,.5] 1:[.5,1] 2:[.5,2] 3:[1,.5] 4:[1,1]
#: 5:[1,2] 6:[2,.5] 7:[2,1] 8:[2,2]. PSWeight: 0 ISR=2, 1 FSR=2, 2 ISR=0.5, 3 FSR=0.5.
SOURCES = {
    "QCDscale_fac_vbsvvh": ("LHEScaleSumw", 5, 3, 4),
    "QCDscale_ren_vbsvvh": ("LHEScaleSumw", 7, 1, 4),
    "ps_isr": ("PSSumw", 0, 2, None),
    "ps_fsr": ("PSSumw", 1, 3, None),
}

#: Skim directory suffix -> datacard era tag. Longest suffix first.
ERA_SUFFIXES = (
    ("_Run3Summer24Summer24for2025", "2025"),
    ("_Run3Summer24_ext1", "2024"),
    ("_Run3Summer24", "2024"),
    ("_UL16APV", "2016preVFP"),
    ("_UL16", "2016postVFP"),
    ("_UL17", "2017"),
    ("_UL18", "2018"),
)

_SAMPLE = re.compile(r"^(?P<proc>VBSWWH_OS|VBSWWH_SS|VBSWZH|VBSZZH)_c2v(?P<c2v>\w+?)_c3_(?P<c3>\w+?)$")


def signal_name(c2v, c3):
    return f"Signal_C2V_{c2v}_C3_{c3}"


def sample_point(signal):
    """``Signal_C2V_1p5_C3_1p0`` -> ``c2v1p5_c3_1p0``, the sample-name spelling."""
    c2v, c3 = re.fullmatch(r"Signal_C2V_(\w+?)_C3_(\w+)", signal).groups()
    return f"c2v{c2v}_c3_{c3}"


def parse_sample(dirname):
    """``VBSWZH_c2v1p5_c3_1p0_UL18`` -> ("VBSWZH", "Signal_C2V_1p5_C3_1p0", "2018"), else None."""
    for suffix, era in ERA_SUFFIXES:
        if dirname.endswith(suffix):
            m = _SAMPLE.match(dirname[: -len(suffix)])
            if m:
                return m["proc"], signal_name(m["c2v"], m["c3"]), era
            return None
    return None


def run_of(era):
    return "run2" if era in ("2016preVFP", "2016postVFP", "2017", "2018") else "run3"


def sample_ratios(files):
    """``{nuisance: (up, down)}`` of one sample, summed over the ``Runs`` trees of ``files``."""
    import numpy as np
    import uproot

    sumw, sums = 0.0, {}
    for path in files:
        runs = uproot.open(path)["Runs"]
        branches = sorted({src[0] for src in SOURCES.values()} & set(runs.keys()))
        arr = runs.arrays(["genEventSumw"] + branches, library="np")
        for i, w in enumerate(arr["genEventSumw"]):
            sumw += w
            for b in branches:
                v = np.asarray(arr[b][i], dtype=float) * w
                sums[b] = v if b not in sums else sums[b] + v
    out = {}
    for name, (branch, up, dn, nom) in SOURCES.items():
        v = sums.get(branch)
        if v is None or len(v) <= max(up, dn):
            continue
        v = v / sumw
        ref = v[nom] if nom is not None else 1.0
        out[name] = (float(v[up] / ref), float(v[dn] / ref))
    return out


def build_table(skim_dirs, xsec, log=print):
    """``{nuisance: {signal: {era: [up, down]}}}`` from the signal sample directories.

    ``xsec``: ``{"run2"|"run3": {"<process>_c2v.._c3..": xsec}}``. A hypothesis needs all
    four processes in an era; an incomplete one is reported and left out.
    """
    per = {}
    for base in skim_dirs:
        for d in sorted(Path(base).iterdir()):
            parsed = parse_sample(d.name)
            if parsed is None:
                continue
            files = sorted(d.glob("*.root"))
            if not files:
                continue
            proc, signal, era = parsed
            per.setdefault((signal, era), {}).setdefault(proc, []).extend(files)

    table = {}
    for (signal, era), procs in sorted(per.items()):
        missing = [p for p in SIGNAL_PROCESSES if p not in procs]
        if missing:
            log(f"skip {signal} {era}: no {missing}")
            continue
        point = sample_point(signal)
        weights, ratios = {}, {}
        for p in SIGNAL_PROCESSES:
            key = f"{p}_{point}"
            try:
                weights[p] = float(xsec[run_of(era)][key])
            except KeyError:
                raise ConfigError(f"no cross section for {key} ({run_of(era)})") from None
            ratios[p] = sample_ratios(procs[p])
        total = sum(weights.values())
        for name in SOURCES:
            if not all(name in ratios[p] for p in SIGNAL_PROCESSES):
                continue
            up = sum(weights[p] * ratios[p][name][0] for p in SIGNAL_PROCESSES) / total
            dn = sum(weights[p] * ratios[p][name][1] for p in SIGNAL_PROCESSES) / total
            table.setdefault(name, {}).setdefault(signal, {})[era] = [round(up, 6), round(dn, 6)]
        log(f"{signal} {era}: " + ", ".join(
            f"{n} {table[n][signal][era][0]:.4f}/{table[n][signal][era][1]:.4f}"
            for n in SOURCES if era in table.get(n, {}).get(signal, {})))
    return table


def write_table(path, table, provenance):
    header = ("# Generator-level ratio varied/nominal of the signal theory variations, [Up, Down],\n"
              "# per nuisance, signal hypothesis and era. Written by `vbsvvh-stat theory-table`; do\n"
              "# not edit by hand.\n")
    with open(path, "w") as f:
        f.write(header)
        yaml.safe_dump({"provenance": provenance, "ratios": table}, f, sort_keys=True, width=100)


def load_table(path):
    path = Path(path)
    if not path.exists():
        return {}
    with open(path) as f:
        return (yaml.safe_load(f) or {}).get("ratios", {})


#: (what the channel delivers, what the result wants) -> how the inclusive ratio r enters.
_CONVERSION = {
    ("raw", "acceptance"): -1,      # divide by r
    ("acceptance", "total"): +1,    # multiply by r
    ("raw", "total"): 0,
    ("acceptance", "acceptance"): 0,
}


def factors(table, nuisance, signal, era, delivered, wanted):
    """``(f_up, f_down)`` multiplying the delivered Up/Down histograms, or None if they stay."""
    power = _CONVERSION[(delivered, wanted)]
    if power == 0:
        return None
    try:
        up, dn = table[nuisance][signal][era]
    except KeyError:
        raise ConfigError(
            f"{nuisance} on {signal} in {era}: delivered {delivered!r}, result wants {wanted!r}, but "
            "theory_inclusive.yaml has no generator-level ratio for it (run `vbsvvh-stat theory-table`)"
        ) from None
    return up ** power, dn ** power
