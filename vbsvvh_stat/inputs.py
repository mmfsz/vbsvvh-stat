"""Read a channel's shapes.root and check it against the input contract (INPUTS.md)."""

import os
import re
from collections import namedtuple

import numpy as np
import uproot

from .config import ANALYSIS_PREFIX, ENERGIES, ERAS, ConfigError

Hist = namedtuple("Hist", "values variances edges")

DATA = "data_obs"
_NAME = re.compile(r"[A-Za-z][A-Za-z0-9_]*$")
_RESERVED = re.compile(rf"^(lumi_|{ANALYSIS_PREFIX}_)")
#: Nuisances named after the signal act on the signal processes only.
SIGNAL_ONLY_SUFFIX = "_vbsvvh"
#: A variation whose varied/nominal ratio leaves this band is reported.
_EXTREME = (0.5, 2.0)


def scale(h, s):
    return Hist(h.values * s, h.variances * s * s, h.edges)


def add(a, b):
    if a is None:
        return b
    return Hist(a.values + b.values, a.variances + b.variances, a.edges)


class Shapes:
    """Every histogram of a shapes.root, keyed ``(era, region, process, nuisance, direction)``.

    ``nuisance`` and ``direction`` are None for a nominal (and for ``data_obs``).
    Keys that do not parse are kept in ``bad_keys`` for :func:`validate`.
    Several files can be given, each holding different eras of the channel.
    """

    def __init__(self, paths):
        paths = [paths] if isinstance(paths, (str, os.PathLike)) else paths
        self.hists = {}
        self.no_sumw2 = []
        self.bad_keys = []
        for path in paths:
            hists = self._read(path)
            shared = sorted({k[0] for k in hists} & {k[0] for k in self.hists})
            if shared:
                raise ConfigError(f"era(s) {shared} of {path} are also in another input of the channel")
            self.hists.update(hists)

    def _read(self, path):
        hists = {}
        with uproot.open(path) as f:
            for key, cls in f.classnames(recursive=True).items():
                key = key.split(";")[0]
                if not cls.startswith("TH1"):
                    continue
                parts = key.split("/")
                if len(parts) != 3:
                    self.bad_keys.append(key)
                    continue
                era, region, name = parts
                parsed = _parse(name)
                if parsed is None:
                    self.bad_keys.append(key)
                    continue
                obj = f[key]
                values = np.asarray(obj.values(flow=False), dtype=float)
                sumw2 = obj.member("fSumw2")
                if len(sumw2) == 0:
                    if parsed[0] != DATA:
                        self.no_sumw2.append(key)
                    variances = np.zeros_like(values)
                else:
                    variances = np.asarray(sumw2, dtype=float)[1:-1]
                edges = np.asarray(obj.axis().edges(), dtype=float)
                hists[(era, region) + parsed] = Hist(values, variances, edges)
        return hists

    @property
    def eras(self):
        return [e for e in ERAS if any(k[0] == e for k in self.hists)]

    def regions(self, era=None):
        return sorted({k[1] for k in self.hists if era is None or k[0] == era})

    def processes(self, era=None, region=None):
        return sorted({k[2] for k in self.hists
                       if (era is None or k[0] == era) and (region is None or k[1] == region)})

    def nuisances(self, era=None, process=None):
        return sorted({k[3] for k in self.hists if k[3]
                       and (era is None or k[0] == era) and (process is None or k[2] == process)})

    def nbins(self, region):
        return {len(h.values) for k, h in self.hists.items() if k[1] == region}


def _parse(name):
    """``proc`` -> (proc, None, None); ``proc__nuisUp`` -> (proc, nuis, "Up"); else None."""
    if "__" not in name:
        return (name, None, None) if _NAME.match(name) else None
    proc, _, rest = name.partition("__")
    for d in ("Up", "Down"):
        if rest.endswith(d) and _NAME.match(proc) and _NAME.match(rest[:-len(d)] or "-"):
            return proc, rest[:-len(d)], d
    return None


def background_kind(channel, process, analysis):
    """``signal``, ``data``, ``data_driven`` or ``mc``."""
    if process == DATA:
        return "data"
    if process in analysis.signals:
        return "signal"
    if channel["abcd"] and process == channel["abcd"]["background"]:
        return "data_driven"
    return "mc"


def validate(channel, analysis, shapes=None):
    """``(errors, warnings, summary)`` for one channel. ``summary`` lists MC yields only.

    Data are never printed, so a signal region cannot be unblinded by validating it.
    """
    shapes = shapes or Shapes(channel["shapes"])
    errors, warnings = [], []
    errors += [f"key does not follow <era>/<region>/<process>[__<nuisance>Up|Down]: {k}"
               for k in shapes.bad_keys]
    errors += [f"no Sumw2 stored (MC statistical uncertainty lost): {k}" for k in shapes.no_sumw2]
    errors += [f"unknown era {e!r}" for e in sorted({k[0] for k in shapes.hists} - set(ERAS))]
    errors += [f"era {e} is not implemented yet (no luminosity in lumi.yaml)"
               for e in shapes.eras if e not in analysis.luminosity]
    if not shapes.hists:
        errors.append("no histograms")
        return errors, warnings, []

    abcd = channel["abcd"]
    if abcd:
        for r in "ABCD":
            if abcd[r] not in shapes.regions():
                errors.append(f"abcd region {r}={abcd[r]!r} not in the shapes")
            elif shapes.nbins(abcd[r]) != {1}:
                errors.append(f"abcd region {abcd[r]!r} must be a 1-bin (counting) region")

    for region in shapes.regions():
        if len(shapes.nbins(region)) > 1:
            errors.append(f"region {region!r}: binning differs between histograms")
        if _RESERVED.search(region) or "__" in region:
            errors.append(f"region name {region!r} is reserved or contains '__'")

    signals = [p for p in shapes.processes() if p in analysis.signals]
    if not signals:
        errors.append(f"no signal process; signals are named as in signals.yaml: {list(analysis.signals)}")
    for p in shapes.processes():
        if p.startswith("Signal") and p not in analysis.signals:
            errors.append(f"process {p!r} looks like a signal but is not in signals.yaml")
        if _RESERVED.search(p):
            errors.append(f"process name {p!r} is reserved")
        if p != DATA and p not in analysis.signals:
            bad = [n for n in shapes.nuisances(process=p) if n.endswith(SIGNAL_ONLY_SUFFIX)]
            if bad:
                errors.append(f"{p}: {bad} act on the signal only")

    for era in shapes.eras:
        for region in shapes.regions(era):
            procs = shapes.processes(era, region)
            if DATA not in procs:
                errors.append(f"{era}/{region}: no {DATA}")
            for p in procs:
                if p == DATA:
                    continue
                nom = shapes.hists.get((era, region, p, None, None))
                if nom is None:
                    errors.append(f"{era}/{region}/{p}: variations without a nominal")
                    continue
                if not np.all(np.isfinite(nom.values)):
                    errors.append(f"{era}/{region}/{p}: non-finite content")
                if (nom.values < 0).any():
                    warnings.append(f"{era}/{region}/{p}: negative bins")
                for n in shapes.nuisances(era, p):
                    up = shapes.hists.get((era, region, p, n, "Up"))
                    dn = shapes.hists.get((era, region, p, n, "Down"))
                    if (up is None) != (dn is None):
                        errors.append(f"{era}/{region}/{p}__{n}: only one of Up/Down")
                    if up is None or dn is None or nom.values.sum() <= 0:
                        continue
                    for d, h in (("Up", up), ("Down", dn)):
                        ratio = h.values.sum() / nom.values.sum()
                        if not _EXTREME[0] <= ratio <= _EXTREME[1]:
                            warnings.append(f"{era}/{region}/{p}__{n}{d}: varied/nominal = {ratio:.3g}")

    for n in shapes.nuisances():
        if analysis.nuisance_entry(n) is None:
            errors.append(f"nuisance {n!r} is not in nuisances.yaml")
    # A nuisance without an era tag is correlated across eras (across the eras of
    # one run, for a centre-of-mass energy tag), so an era that has the process
    # but not the nuisance would silently be read as nominal.
    for p in shapes.processes():
        if p == DATA or background_kind(channel, p, analysis) == "data_driven":
            continue
        carried = {n for n in shapes.nuisances(process=p)}
        for n in sorted(carried):
            if any(n.endswith("_" + e) for e in ERAS):
                continue
            scope = next((eras for tag, eras in ENERGIES.items() if n.endswith("_" + tag)), ERAS)
            missing = [e for e in shapes.eras if e in scope and p in shapes.processes(e)
                       and n not in shapes.nuisances(e, p)]
            if missing:
                errors.append(f"{p}: correlated nuisance {n} missing in era(s) {missing}")

    if abcd:
        for era in shapes.eras:
            for r in "BCD":
                h = shapes.hists.get((era, abcd[r], DATA, None, None))
                if h is not None and h.values.sum() <= 0:
                    errors.append(f"{era}/{abcd[r]}: no data in ABCD region {r}")

    summary = []
    for era in shapes.eras:
        for region in shapes.regions(era):
            mc = {p: shapes.hists[(era, region, p, None, None)].values.sum()
                  for p in shapes.processes(era, region)
                  if p != DATA and (era, region, p, None, None) in shapes.hists}
            nn = len(shapes.nuisances(era))
            summary.append(f"{era:12s} {region:12s} " + "  ".join(f"{p}={v:.4g}" for p, v in mc.items())
                           + f"  ({nn} nuisances)")
    return errors, warnings, summary
