"""provenance.json can be handed out: no local path in it."""

import json

import yaml

from conftest import ABCD, abcd_hists, write_channel
from vbsvvh_stat.config import REPO, load_result
from vbsvvh_stat.provenance import capture


def _strings(obj):
    if isinstance(obj, dict):
        for v in obj.values():
            yield from _strings(v)
    elif isinstance(obj, list):
        for v in obj:
            yield from _strings(v)
    elif isinstance(obj, str):
        yield obj


def test_no_absolute_path(tmp_path, monkeypatch, analysis):
    ch = write_channel(tmp_path / "ch", "T", abcd_hists(["2018"]), ABCD)
    cfg = tmp_path / "result.yaml"
    cfg.write_text(yaml.safe_dump({"name": "test", "signals": ["Signal_C2V_1p5_C3_1p0"],
                                   "channels": {"T": {"input": str(ch)}}}))
    monkeypatch.setattr("sys.argv", [str(REPO / "vbsvvh_stat" / "__main__.py"), "run", str(cfg)])
    prov = capture(load_result(cfg, analysis))
    assert not [s for s in _strings(prov) if s.startswith("/")], json.dumps(prov, indent=1)
    assert prov["argv"][:2] == ["vbsvvh-stat", "run"]
    assert prov["result"]["channels"] == ["T"]
    assert [list(i) for i in prov["channels"]["T"]] == [["shapes_sha256", "provenance"]]
