"""Nuisance names against the CMS dictionary (systematics/systematics_master.yml).

The CMS datacard CI (``check_names.py``) needs every nuisance of a card in an
analysis-local ``systematics.yml`` with a class. A name the master dictionary
knows takes the master's class; any other is ``custom`` with a description.
Matching follows ``check_names.py``: a master entry is a regex, ``fullmatch``ed.
"""

import re

import yaml

from .config import ANALYSIS_PREFIX, load_yaml

#: Descriptions of the nuisances the tool itself adds, by name pattern.
_TOOL_MADE = (
    (re.compile(rf"{ANALYSIS_PREFIX}_.+_mcstat_.+"), "MC statistical uncertainty of one process in one bin"),
    (re.compile(rf"{ANALYSIS_PREFIX}_.+_abcd_closure"), "ABCD non-closure of the data-driven background in region A"),
)


class Dictionary:
    def __init__(self, master_path):
        master = load_yaml(master_path)
        self.patterns = [(re.compile(p), cls) for cls, entries in master.items() for p in entries]

    def cls(self, name):
        for pattern, cls in self.patterns:
            if pattern.fullmatch(name):
                return cls
        return None


def declaration(names, analysis, dictionary):
    """The CI's ``systematics.yml`` body for the nuisances ``names``."""
    out = {}
    for name in sorted(names):
        cls = dictionary.cls(name)
        if cls:
            out[name] = {"class": cls}
            continue
        key = analysis.nuisance_entry(name)
        desc = (analysis.nuisances.get(key) or {}).get("description") if key else None
        desc = desc or (analysis.lumi.get(name) or {}).get("description")
        desc = desc or next((d for p, d in _TOOL_MADE if p.fullmatch(name)), None)
        if not desc:
            raise ValueError(f"nuisance {name!r} is neither in the CMS dictionary nor described "
                             "in nuisances.yaml / lumi.yaml")
        out[name] = {"class": "custom", "description": desc}
    return out


def write_declaration(path, names, analysis, dictionary):
    with open(path, "w") as f:
        yaml.safe_dump(declaration(names, analysis, dictionary), f, sort_keys=True, width=100)
