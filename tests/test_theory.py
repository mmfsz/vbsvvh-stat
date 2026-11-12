"""Signal theory variations: one normalisation for every channel, chosen by the result."""

import json

import numpy as np
import pytest
import yaml

from conftest import write_channel
from vbsvvh_stat.cards import write_cards
from vbsvvh_stat.config import ConfigError, load_channel, load_result
from vbsvvh_stat.inputs import validate
from vbsvvh_stat.naming import Dictionary
from vbsvvh_stat.theory import build_table, factors, parse_sample

SIG = "Signal_C2V_1p5_C3_1p0"
MUF = "QCDscale_fac_vbsvvh"
#: Generator-level ratios [Up, Down] per era, as theory_inclusive.yaml holds them.
TABLE = {MUF: {SIG: {"2017": [0.88, 1.16], "2018": [0.90, 1.12]}}}


def lnn(card_text, name):
    """{(bin, process): (down, up)} of an lnN row, read from the card text."""
    last_bin = bins = procs = None
    for line in card_text.splitlines():
        f = line.split()
        if f and f[0] == "bin":
            last_bin = f[1:]
        if f and f[0] == "process" and procs is None:
            bins, procs = last_bin, f[1:]
        if f and f[0] == name:
            out = {}
            for b, p, c in zip(bins, procs, f[2:]):
                if c != "-":
                    dn, up = (float(x) for x in c.split("/"))
                    out[(b, p)] = (dn, up)
            return out
    return {}


def channel(tmp_path, theory, eras=(("2018", 2.0),), acc=1.02):
    """One signal region; the varied yield = nominal x inclusive ratio x acceptance change."""
    h = {}
    for era, y in eras:
        r_up, r_dn = TABLE[MUF][SIG][era]
        h[f"{era}/SR/{SIG}"] = y
        h[f"{era}/SR/ttbar"] = 5.0
        h[f"{era}/SR/data_obs"] = 7.0
        if theory == "raw":
            up, dn = y * r_up * acc, y * r_dn / acc
        else:
            up, dn = y * acc, y / acc
        h[f"{era}/SR/{SIG}__{MUF}Up"] = up
        h[f"{era}/SR/{SIG}__{MUF}Down"] = dn
    return write_channel(tmp_path / f"ch_{theory}", "T", h, theory=theory)


def cards(tmp_path, monkeypatch, analysis, ch, **extra):
    monkeypatch.setenv("VBSVVH_STAT_RESULTS", str(tmp_path / "results"))
    analysis.theory = TABLE
    cfg = {"name": "t", "signals": [SIG], "channels": {"T": {"input": str(ch)}}, "unblind": True, **extra}
    path = tmp_path / "result.yaml"
    path.write_text(yaml.safe_dump(cfg))
    result = load_result(path, analysis)
    index = write_cards(result, analysis, Dictionary(analysis.master), log=lambda *_: None)
    return index["cards"][SIG]["Run2"][0].read_text(), index


def test_raw_delivery_is_written_acceptance_only(tmp_path, monkeypatch, analysis):
    text, _ = cards(tmp_path, monkeypatch, analysis, channel(tmp_path, "raw"))
    dn, up = lnn(text, MUF)[("T_SR_Run2", SIG)]
    assert up == pytest.approx(1.02) and dn == pytest.approx(1 / 1.02)
    assert "delivered raw, written acceptance" in text


def test_flat_acceptance_leaves_no_nuisance(tmp_path, monkeypatch, analysis):
    text, _ = cards(tmp_path, monkeypatch, analysis, channel(tmp_path, "raw", acc=1.0))
    assert lnn(text, MUF) == {}


def test_raw_delivery_written_total(tmp_path, monkeypatch, analysis):
    text, _ = cards(tmp_path, monkeypatch, analysis, channel(tmp_path, "raw"), theory="total")
    dn, up = lnn(text, MUF)[("T_SR_Run2", SIG)]
    assert up == pytest.approx(0.90 * 1.02) and dn == pytest.approx(1.12 / 1.02)


def test_acceptance_delivery_written_total(tmp_path, monkeypatch, analysis):
    text, _ = cards(tmp_path, monkeypatch, analysis, channel(tmp_path, "acceptance"), theory="total")
    dn, up = lnn(text, MUF)[("T_SR_Run2", SIG)]
    assert up == pytest.approx(0.90 * 1.02) and dn == pytest.approx(1.12 / 1.02)


def test_both_deliveries_give_the_same_card(tmp_path, monkeypatch, analysis):
    raw, _ = cards(tmp_path / "a", monkeypatch, analysis, channel(tmp_path / "a", "raw"))
    acc, _ = cards(tmp_path / "b", monkeypatch, analysis, channel(tmp_path / "b", "acceptance"))
    assert lnn(raw, MUF) == pytest.approx(lnn(acc, MUF))


def test_each_era_divided_by_its_own_ratio_before_the_sum(tmp_path, monkeypatch, analysis):
    eras = (("2017", 1.0), ("2018", 3.0))
    text, index = cards(tmp_path, monkeypatch, analysis, channel(tmp_path, "raw", eras=eras))
    dn, up = lnn(text, MUF)[("T_SR_Run2", SIG)]
    assert up == pytest.approx(1.02) and dn == pytest.approx(1 / 1.02)
    assert index["summary"]["T"]["Run2"]["theory"] == {"delivered": "raw", "written": "acceptance"}


def test_missing_ratio_stops_the_cards(tmp_path, monkeypatch, analysis):
    ch = channel(tmp_path, "raw", eras=(("2018", 2.0),))
    TABLE_SAVED = TABLE[MUF][SIG].pop("2018")
    try:
        with pytest.raises(ConfigError, match="no generator-level ratio"):
            cards(tmp_path, monkeypatch, analysis, ch)
    finally:
        TABLE[MUF][SIG]["2018"] = TABLE_SAVED


def test_undeclared_theory_is_an_input_error(tmp_path, analysis):
    path = channel(tmp_path, "raw")
    cfg = yaml.safe_load(path.read_text())
    cfg.pop("theory")
    path.write_text(yaml.safe_dump(cfg))
    errors = validate(load_channel(path), analysis)[0]
    assert any("declare in channel.yaml" in e for e in errors)


def test_background_variations_are_untouched(tmp_path, monkeypatch, analysis):
    h = {"2018/SR/" + SIG: 2.0, "2018/SR/ttbar": 5.0, "2018/SR/data_obs": 7.0,
         "2018/SR/ttbar__ps_fsrUp": 5.5, "2018/SR/ttbar__ps_fsrDown": 4.5}
    analysis.theory = {"ps_fsr": {SIG: {"2018": [0.5, 2.0]}}}
    ch = write_channel(tmp_path / "ch", "T", h, theory="raw")
    monkeypatch.setenv("VBSVVH_STAT_RESULTS", str(tmp_path / "results"))
    path = tmp_path / "result.yaml"
    path.write_text(yaml.safe_dump({"name": "t", "signals": [SIG], "channels": {"T": {"input": str(ch)}},
                                    "unblind": True}))
    result = load_result(path, analysis)
    index = write_cards(result, analysis, Dictionary(analysis.master), log=lambda *_: None)
    dn, up = lnn(index["cards"][SIG]["Run2"][0].read_text(), "ps_fsr")[("T_SR_Run2", "T_ttbar")]
    assert up == pytest.approx(1.1) and dn == pytest.approx(0.9)


def test_factors():
    assert factors(TABLE, MUF, SIG, "2018", "acceptance", "acceptance") is None
    assert factors(TABLE, MUF, SIG, "2018", "raw", "total") is None
    assert factors(TABLE, MUF, SIG, "2018", "raw", "acceptance") == pytest.approx((1 / 0.90, 1 / 1.12))
    assert factors(TABLE, MUF, SIG, "2018", "acceptance", "total") == pytest.approx((0.90, 1.12))


@pytest.mark.parametrize("name,expected", [
    ("VBSWZH_c2v1p5_c3_1p0_UL18", ("VBSWZH", SIG, "2018")),
    ("VBSWWH_OS_c2v1p0_c3_10p0_UL16APV", ("VBSWWH_OS", "Signal_C2V_1p0_C3_10p0", "2016preVFP")),
    ("VBSWWH_SS_c2v1p0_c3_1p0_UL16", ("VBSWWH_SS", "Signal_C2V_1p0_C3_1p0", "2016postVFP")),
    ("VBSZZH_c2v1p5_c3_1p0_Run3Summer24_ext1", ("VBSZZH", SIG, "2024")),
    ("VBSZZH_c2v1p5_c3_1p0_Run3Summer24Summer24for2025", ("VBSZZH", SIG, "2025")),
    ("TTToHadronic_UL18", None),
])
def test_parse_sample(name, expected):
    assert parse_sample(name) == expected


def _runs_file(path, scale_up, scale_dn, sumw):
    import uproot
    lhe = np.ones((1, 9))
    lhe[0, 5], lhe[0, 3] = scale_up, scale_dn
    with uproot.recreate(path) as f:
        f["Runs"] = {"genEventSumw": np.array([sumw]), "LHEScaleSumw": lhe,
                     "PSSumw": np.array([[1.01, 0.99, 0.98, 1.02]])}


def test_build_table_weights_processes_by_cross_section(tmp_path):
    procs = {"VBSWWH_OS": (0.80, 1.20), "VBSWWH_SS": (0.90, 1.10), "VBSWZH": (0.90, 1.10), "VBSZZH": (0.90, 1.10)}
    for p, (u, d) in procs.items():
        d_ = tmp_path / f"{p}_c2v1p5_c3_1p0_UL18"
        d_.mkdir()
        _runs_file(d_ / "a.root", u, d, 3.0)
        _runs_file(d_ / "b.root", u, d, 1.0)
    xsec = {"run2": {f"{p}_c2v1p5_c3_1p0": (3.0 if p == "VBSWWH_OS" else 1.0) for p in procs}}
    table = build_table([tmp_path], xsec, log=lambda *_: None)
    assert table[MUF][SIG]["2018"] == pytest.approx([(3 * 0.8 + 3 * 0.9) / 6, (3 * 1.2 + 3 * 1.1) / 6])
    assert table["ps_fsr"][SIG]["2018"] == pytest.approx([0.99, 1.02])
    assert table["ps_isr"][SIG]["2018"] == pytest.approx([1.01, 0.98])
    assert table["QCDscale_ren_vbsvvh"][SIG]["2018"] == pytest.approx([1.0, 1.0])
