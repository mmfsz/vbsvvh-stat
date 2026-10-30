"""vbsvvh-stat: datacards, expected limits, plots and impacts for the VBS VVH channels.

    vbsvvh-stat validate CHANNEL.yaml
    vbsvvh-stat run RESULT.yaml [--stage cards|limits] [--plots-only] [--jobs N]
    vbsvvh-stat impacts RESULT.yaml [--fit F] [--signal S] [--plots-only] [--jobs N]
"""

import argparse
import json
import shutil
import sys

from . import cards, combine, plots, provenance
from .config import Analysis, ConfigError, load_channel, load_result
from .naming import Dictionary


def _validate(args):
    analysis = Analysis()
    channel = load_channel(args.channel)
    from .inputs import validate
    errors, warnings, summary = validate(channel, analysis)
    print(f"channel {channel['tag']}: {channel['shapes']}")
    print("\n".join(summary))
    for w in warnings:
        print(f"warning: {w}")
    for e in errors:
        print(f"error: {e}")
    print(f"{len(errors)} error(s), {len(warnings)} warning(s)")
    return 1 if errors else 0


def _fit_eras(index):
    return {fit: sorted({e for r in runs for p in index["members"][r] for e in index["planes"][p]["eras"]})
            for fit, runs in index["fits"].items()}


def _run(args):
    analysis = Analysis()
    result = load_result(args.result, analysis)
    out = result["results_dir"]
    meta = out / "fits.json"
    if args.plots_only:
        index = json.loads(meta.read_text())
        lims = json.loads((out / "limits.json").read_text())
    else:
        out.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(result["path"], out / "result.yaml")
        (out / "provenance.json").write_text(json.dumps(provenance.capture(result), indent=2))
        index = cards.write_cards(result, analysis, Dictionary(analysis.master))
        meta.write_text(json.dumps({k: index[k] for k in ("planes", "members", "fits")}, indent=2))
        if args.stage == "cards":
            return 0
        lims = combine.limits(result, index, jobs=args.jobs)
        if args.stage == "limits":
            return 0
    plots.plot_limits(lims, index["fits"], _fit_eras(index), analysis, out / "plots")
    print(f"results in {out}")
    return 0


def _impacts(args):
    analysis = Analysis()
    result = load_result(args.result, analysis)
    out = result["results_dir"]
    lims = json.loads((out / "limits.json").read_text())
    todo = {args.fit: [args.signal]} if args.fit and args.signal else result["impacts"]
    if not todo:
        sys.exit("no impacts requested: give --fit and --signal, or 'impacts:' in the result")
    for fit, signals in todo.items():
        for signal in signals:
            path = combine.impacts(result, fit, signal, lims, jobs=args.jobs, plots_only=args.plots_only)
            print(f"impacts {fit} {signal}: {path}")
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(prog="vbsvvh-stat", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="command", required=True)
    v = sub.add_parser("validate", help="check a channel's shapes.root + channel.yaml")
    v.add_argument("channel")
    v.set_defaults(func=_validate)
    r = sub.add_parser("run", help="cards, expected limits and plots of a result")
    r.add_argument("result")
    r.add_argument("--stage", choices=["cards", "limits"], help="stop after this stage")
    r.add_argument("--plots-only", action="store_true", help="redraw the plots from limits.json")
    r.add_argument("--jobs", type=int, default=4)
    r.set_defaults(func=_run)
    i = sub.add_parser("impacts", help="nuisance impacts (needs `run` first)")
    i.add_argument("result")
    i.add_argument("--fit")
    i.add_argument("--signal")
    i.add_argument("--plots-only", action="store_true", help="redraw the plots from impacts.json")
    i.add_argument("--jobs", type=int, default=4)
    i.set_defaults(func=_impacts)
    args = ap.parse_args(argv)
    try:
        return args.func(args)
    except ConfigError as exc:
        sys.exit(f"error: {exc}")


if __name__ == "__main__":
    sys.exit(main())
