"""What a result was made from, written to ``<results>/provenance.json``.

No local path is recorded: the inputs are identified by their sha256, the
environment by its CMSSW release and Combine/CombineHarvester tags, so the file
can be handed out with the results.
"""

import hashlib
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

from .config import REPO


def _git(path, *args):
    try:
        return subprocess.run(["git", *args], cwd=path, capture_output=True, text=True,
                              timeout=20).stdout.strip() or None
    except (OSError, subprocess.SubprocessError):
        return None


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _arg(arg):
    """A command-line argument; a path becomes relative to the repository, or its file name
    if it lies outside."""
    path = Path(arg).resolve()
    if not path.exists():
        return arg
    return str(path.relative_to(REPO)) if path.is_relative_to(REPO) else path.name


def capture(result):
    cmssw = os.environ.get("CMSSW_BASE")
    return {
        "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "argv": ["vbsvvh-stat"] + [_arg(a) for a in sys.argv[1:]],
        "vbsvvh_stat": {"commit": _git(REPO, "rev-parse", "HEAD"),
                        "dirty": bool(_git(REPO, "status", "--porcelain"))},
        "systematics_dictionary": _git(REPO / "systematics", "rev-parse", "HEAD"),
        "cmssw": os.environ.get("CMSSW_VERSION"),
        "combine": _git(f"{cmssw}/src/HiggsAnalysis/CombinedLimit", "describe", "--tags", "--always") if cmssw else None,
        "combineharvester": _git(f"{cmssw}/src/CombineHarvester", "describe", "--tags", "--always") if cmssw else None,
        "result": dict(result["raw"], channels=sorted(result["raw"]["channels"])),
        "channels": {tag: [{"shapes_sha256": _sha256(i["shapes"]), "provenance": i["provenance"]}
                           for i in ch["inputs"]]
                     for tag, ch in result["channels"].items()},
    }
