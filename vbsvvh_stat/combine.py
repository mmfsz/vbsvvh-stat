"""Running Combine: one work directory per (fit, signal) holding the combined card,
the workspace, the limit and the impacts."""

import json
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor

import uproot

#: ``quantileExpected`` -> key; Combine stores the quantiles as floats.
QUANTILES = ((0.025, "exp-2"), (0.16, "exp-1"), (0.5, "exp0"), (0.84, "exp+1"), (0.975, "exp+2"))
_QTOL = 5e-3
#: The ABCD rateParams and non-closure, drawn apart from the other nuisances in the impacts.
_ABCD = re.compile(r"_abcd_(A|B|C|D|closure)$")


def run(cmd, cwd, log, stdout=None):
    """Run ``cmd`` in ``cwd``; output to ``log`` (stdout to ``stdout`` if given). Raises on failure."""
    with open(cwd / log, "w") as err:
        out = open(cwd / stdout, "w") if stdout else err
        try:
            proc = subprocess.run(cmd, cwd=cwd, stdout=out, stderr=err)
        finally:
            if stdout:
                out.close()
    if proc.returncode:
        tail = "\n".join((cwd / log).read_text().splitlines()[-20:])
        raise RuntimeError(f"{' '.join(map(str, cmd))} failed in {cwd}:\n{tail}")


def fit_cards(index):
    """``{(fit, signal): [card, ...]}`` from the card index of ``cards.write_cards``."""
    out = {}
    for fit, runs in index["fits"].items():
        planes = [p for r in runs for p in index["members"][r]]
        for signal, by_plane in index["cards"].items():
            cards = [c for p in planes for c in by_plane.get(p, [])]
            if cards:
                out[(fit, signal)] = cards
    return out


def workdir(result, fit, signal):
    return result["work_dir"] / fit / signal


def limit(result, fit, signal, cards):
    """Workspace and expected limit of one (fit, signal); returns the quantiles."""
    wdir = workdir(result, fit, signal)
    wdir.mkdir(parents=True, exist_ok=True)
    run(["combineCards.py", "-S"] + [f".={c}" for c in cards], wdir, "combineCards.log", stdout="combined.txt")
    with open(wdir / "validate.log", "w") as f:
        subprocess.run(["ValidateDatacards.py", "combined.txt", "--printLevel", "1"], cwd=wdir,
                       stdout=f, stderr=subprocess.STDOUT)
    run(["text2workspace.py", "combined.txt", "-o", "ws.root", "-m", "125"], wdir, "text2workspace.log")
    run(["combine", "-M", "AsymptoticLimits", "-d", "ws.root", "-t", "-1", "--run", "blind",
         "-m", "125", "-n", ".limit"], wdir, "limit.log")
    return read_limit(wdir / "higgsCombine.limit.AsymptoticLimits.mH125.root")


def read_limit(path):
    tree = uproot.open(path)["limit"]
    values = tree["limit"].array(library="np")
    quantiles = tree["quantileExpected"].array(library="np")
    out = {}
    for q, key in QUANTILES:
        hits = [float(v) for v, qe in zip(values, quantiles) if abs(qe - q) < _QTOL]
        if not hits:
            raise RuntimeError(f"{path}: no quantile {q}")
        out[key] = hits[0]
    return out


def limits(result, index, jobs=4, log=print):
    """Expected limits of every (fit, signal), written to ``limits.json``."""
    todo = fit_cards(index)
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        futures = {key: pool.submit(limit, result, *key, cards) for key, cards in todo.items()}
    out = {}
    for (fit, signal), fut in futures.items():
        out.setdefault(fit, {})[signal] = q = fut.result()
        log(f"limit {fit:12s} {signal:24s} {q['exp0']:.4g} [{q['exp-1']:.4g}, {q['exp+1']:.4g}]")
    path = result["results_dir"] / "limits.json"
    with open(path, "w") as f:
        json.dump(out, f, indent=2)
    return out


def impacts(result, fit, signal, limits_, jobs=4, plots_only=False):
    """Nuisance impacts of one (fit, signal) on the Asimov dataset, and their plots."""
    wdir = workdir(result, fit, signal)
    if plots_only:
        if not (wdir / "impacts.json").exists():
            raise RuntimeError(f"no impacts.json in {wdir}: run the impacts first")
    else:
        if not (wdir / "ws.root").exists():
            raise RuntimeError(f"no workspace in {wdir}: run the limits first")
        r_max = max(2 * limits_[fit][signal]["exp+2"], 5.0)
        common = ["combineTool.py", "-M", "Impacts", "-d", "ws.root", "-m", "125", "-t", "-1",
                  "--expectSignal", "1", "--robustFit", "1", "--setParameterRanges", f"r={-r_max:g},{r_max:g}"]
        run(common + ["--doInitialFit"], wdir, "impacts_initial.log")
        run(common + ["--doFits", "--parallel", str(jobs)], wdir, "impacts_fits.log")
        run(common + ["-o", "impacts.json"], wdir, "impacts_collect.log")
    out_dir = result["results_dir"] / "impacts"
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = out_dir / f"impacts_{fit}_{signal}"
    data = json.loads((wdir / "impacts.json").read_text())
    for suffix, params in split_impacts(data["params"]).items():
        name = f"impacts{suffix}.json"
        if suffix:
            (wdir / name).write_text(json.dumps(dict(data, params=params)))
        run(["plotImpacts.py", "-i", name, "-o", f"{stem}{suffix}"], wdir, f"impacts_plot{suffix}.log")
    (out_dir / f"impacts_{fit}_{signal}.json").write_text((wdir / "impacts.json").read_text())
    return stem.with_suffix(".pdf")


def split_impacts(params):
    """``{suffix: params}`` of the impact figures: every parameter (``""``), and, when
    there are ABCD terms, those alone (``"_abcd"``) and all the others (``"_noabcd"``),
    which the ABCD terms would otherwise squeeze on a common axis."""
    abcd = [p for p in params if _ABCD.search(p["name"])]
    rest = [p for p in params if not _ABCD.search(p["name"])]
    if not abcd or not rest:
        return {"": params}
    return {"": params, "_abcd": abcd, "_noabcd": rest}
