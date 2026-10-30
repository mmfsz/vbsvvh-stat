"""Cards: the numbers Combine reads, parsed back with Combine's own DatacardParser."""

import json
import math
import subprocess
from optparse import OptionParser

import pytest
import yaml

from conftest import ABCD, abcd_hists, th1, write_channel
from vbsvvh_stat.cards import write_cards
from vbsvvh_stat.config import ConfigError, load_result
from vbsvvh_stat.naming import Dictionary

SIG = "Signal_C2V_1p5_C3_1p0"


def parse(path):
    from HiggsAnalysis.CombinedLimit.DatacardParser import addDatacardParserOptions, parseCard
    parser = OptionParser()
    addDatacardParserOptions(parser)
    options, _ = parser.parse_args([])
    with open(path) as f:
        return parseCard(f, options)


def kappa(dc, name, bin_, proc):
    """(down, up) of an lnN entry."""
    for syst in dc.systs:
        if syst[0] == name:
            k = syst[4][bin_][proc]
            return tuple(k) if isinstance(k, list) else (1 / k, k)
    raise KeyError(name)


def make_result(tmp_path, monkeypatch, analysis, channels, **extra):
    monkeypatch.setenv("VBSVVH_STAT_RESULTS", str(tmp_path / "results"))
    cfg = {"name": "test", "signals": [SIG],
           "channels": {tag: {"input": str(path)} for tag, path in channels.items()}, **extra}
    path = tmp_path / "result.yaml"
    path.write_text(yaml.safe_dump(cfg))
    result = load_result(path, analysis)
    return result, write_cards(result, analysis, Dictionary(analysis.master), log=lambda *_: None)


def test_era_sum_dilutes_era_tagged_and_sums_correlated_nuisances(tmp_path, monkeypatch, analysis):
    # README example: 2017 = 1.1, 2018 = 2.1, a 9% 2018-only nuisance -> 1.0590625 on the sum.
    h = {}
    for era, y, fsr in (("2017", 1.1, 1.05), ("2018", 2.1, 1.10)):
        h[f"{era}/SR/{SIG}"] = y
        h[f"{era}/SR/ttbar"] = 5.0
        h[f"{era}/SR/data_obs"] = 7.0
        h[f"{era}/SR/{SIG}__ps_fsrUp"] = y * fsr
        h[f"{era}/SR/{SIG}__ps_fsrDown"] = y / fsr
    h[f"2018/SR/{SIG}__CMS_scale_met_unclustered_energy_2018Up"] = 2.1 * 1.09
    h[f"2018/SR/{SIG}__CMS_scale_met_unclustered_energy_2018Down"] = 2.1 / 1.09
    ch = write_channel(tmp_path / "ch", "T", h)
    result, index = make_result(tmp_path, monkeypatch, analysis, {"T": ch}, unblind=True)
    dc = parse(index["cards"][SIG]["Run2"][0])
    b = "T_SR_Run2"
    assert dc.exp[b][SIG] == pytest.approx(3.2)
    assert dc.exp[b]["T_ttbar"] == pytest.approx(10.0)
    assert dc.obs[b] == pytest.approx(14.0)
    assert kappa(dc, "CMS_scale_met_unclustered_energy_2018", b, SIG)[1] == pytest.approx(1.0590625, rel=1e-6)
    assert kappa(dc, "CMS_scale_met_unclustered_energy_2018", b, SIG)[0] == pytest.approx((2.1 / 1.09 + 1.1) / 3.2, rel=1e-6)
    assert kappa(dc, "ps_fsr", b, SIG)[1] == pytest.approx((1.1 * 1.05 + 2.1 * 1.10) / 3.2, rel=1e-6)
    k17, k18 = analysis.lumi_kappas("2017"), analysis.lumi_kappas("2018")
    lumi = "lumi_13TeV_15161718_l"
    assert kappa(dc, lumi, b, SIG)[1] == pytest.approx((1.1 * k17[lumi] + 2.1 * k18[lumi]) / 3.2, rel=1e-6)
    assert kappa(dc, lumi, b, "T_ttbar")[1] == pytest.approx((k17[lumi] + k18[lumi]) / 2, rel=1e-6)
    assert kappa(dc, f"CMS_HIG26013_T_mcstat_SR_Run2_{SIG}", b, SIG)[1] == pytest.approx(1 + math.sqrt(0.032) / 3.2, rel=1e-5)


def test_abcd_block(tmp_path, monkeypatch, analysis):
    ch = write_channel(tmp_path / "ch", "T", abcd_hists(["2018"], bkg=(2.0, 20.0, 30.0, 330.0)), ABCD)
    result, index = make_result(tmp_path, monkeypatch, analysis, {"T": ch})
    path = index["cards"][SIG]["Run2"][0]
    dc = parse(path)
    assert dc.exp["T_A_Run2"]["T_bkg"] == 1.0
    names = {s[0] for s in dc.systs}
    # c = nA nD / (nB nC) = 2 * 330 / 600 = 1.1, MC variances 0.01 n per region
    err = 1.1 * math.sqrt(sum(0.01 / n for n in (2.0, 20.0, 30.0, 330.0)))
    assert kappa(dc, "CMS_HIG26013_T_Run2_abcd_closure", "T_A_Run2", "T_bkg")[1] == pytest.approx(
        1 + math.sqrt(0.1 ** 2 + err ** 2))
    assert not any("_mcstat_" in n and n.endswith("_bkg") for n in names)
    assert all(s[4]["T_A_Run2"].get("T_bkg", 0) == 0 for s in dc.systs if s[0].startswith("lumi"))
    text = path.read_text()
    n = 21.0
    assert f"CMS_HIG26013_T_Run2_abcd_B rateParam T_B_Run2 T_bkg 21 [{n - 3 * math.sqrt(n):.8g},{n + 3 * math.sqrt(n):.8g}]" in text
    assert ("CMS_HIG26013_T_Run2_abcd_A rateParam T_A_Run2 T_bkg (@0*@1/@2) "
            "CMS_HIG26013_T_Run2_abcd_B,CMS_HIG26013_T_Run2_abcd_C,CMS_HIG26013_T_Run2_abcd_D") in text
    assert index["summary"]["T"]["Run2"]["closure"]["delta"] == pytest.approx(0.1)


def test_closure_override(tmp_path, monkeypatch, analysis):
    ch = write_channel(tmp_path / "ch", "T", abcd_hists(["2018"]), ABCD)
    _, index = make_result(tmp_path, monkeypatch, analysis, {"T": ch}, closure={"T": {"Run2": 0.25}})
    assert index["summary"]["T"]["Run2"]["closure"] == {"delta": 0.25, "error": None, "method": "config",
                                                         "kappa": 1.25}
    with pytest.raises(ConfigError, match="the non-closure must be a number"):
        make_result(tmp_path, monkeypatch, analysis, {"T": ch}, closure={"T": {"Run2": "data_sideband"}})


def test_substitute_builds_a_scaled_relabelled_era(tmp_path, monkeypatch, analysis):
    hists = abcd_hists(["2024"], nuis=[("CMS_scale_met_unclustered_energy_{era}", 1.02)])
    ch = write_channel(tmp_path / "ch", "T", hists, ABCD)
    _, index = make_result(tmp_path, monkeypatch, analysis, {"T": ch}, runs={"Run3": {"substitute": {2025: 2024}}})
    dc = parse(index["cards"][SIG]["Run3"][0])
    r = analysis.luminosity["2025"] / analysis.luminosity["2024"]
    assert dc.exp["T_A_Run3"][SIG] == pytest.approx(1.0 + r)
    share = r / (1 + r)
    assert kappa(dc, "CMS_scale_met_unclustered_energy_2025", "T_A_Run3", SIG)[1] == pytest.approx(1 + share * 0.02, rel=1e-6)
    assert kappa(dc, "CMS_scale_met_unclustered_energy_2024", "T_A_Run3", SIG)[1] == pytest.approx(1 + (1 - share) * 0.02, rel=1e-6)
    assert kappa(dc, "lumi_13p6TeV_2025", "T_A_Run3", SIG)[1] == pytest.approx(1 + share * 0.05, rel=1e-6)
    assert index["planes"]["Run3"]["eras"] == ["2024", "2025"]
    assert "21 [" not in (index["cards"][SIG]["Run3"][0]).read_text()  # data in B scaled too
    assert f"T_B_Run3 T_bkg {21 * (1 + r):.8g}" in index["cards"][SIG]["Run3"][0].read_text()


def test_eras_of_an_explicit_run_leave_the_default_run(tmp_path, monkeypatch, analysis):
    ch = write_channel(tmp_path / "ch", "T", abcd_hists(["2024", "2025"]), ABCD)
    _, index = make_result(tmp_path, monkeypatch, analysis, {"T": ch},
                           runs={"Run3_2024": {"eras": [2024]}, "Run3_2025": {"eras": [2025]}})
    assert set(index["planes"]) == {"Run3_2024", "Run3_2025"}


def test_only_planes_in_a_fit_get_cards(tmp_path, monkeypatch, analysis):
    ch = write_channel(tmp_path / "ch", "T", abcd_hists(["2018", "2024"]), ABCD)
    result, index = make_result(tmp_path, monkeypatch, analysis, {"T": ch}, fits={"Run3": ["Run3"]})
    assert set(index["planes"]) == {"Run3"} and set(index["members"]) == {"Run3"}
    assert [c.name for c in (result["results_dir"] / "cards" / SIG).iterdir()] == ["T_Run3.txt"]


def test_substitute_into_an_unimplemented_era(tmp_path, monkeypatch, analysis):
    ch = write_channel(tmp_path / "ch", "T", abcd_hists(["2024"]), ABCD)
    with pytest.raises(ConfigError, match="not implemented yet"):
        make_result(tmp_path, monkeypatch, analysis, {"T": ch}, runs={"Run3": {"substitute": {2023: 2024}}})


def test_names_carry_the_channel_tag_but_signals_do_not(tmp_path, monkeypatch, analysis):
    ch = write_channel(tmp_path / "ch", "MET1FJ", abcd_hists(["2018"]), ABCD)
    _, index = make_result(tmp_path, monkeypatch, analysis, {"MET1FJ": ch})
    dc = parse(index["cards"][SIG]["Run2"][0])
    assert set(dc.bins) == {f"MET1FJ_{r}_Run2" for r in "ABCD"}
    assert set(dc.processes) == {SIG, "MET1FJ_bkg"}


def test_shape_channel_combined_with_an_abcd_channel(tmp_path, monkeypatch, analysis):
    shape = {}
    for era in ("2017", "2018"):
        shape[f"{era}/SR/{SIG}"] = th1([0.5, 1.0, 1.5], [0.01, 0.01, 0.01])
        shape[f"{era}/SR/ttbar"] = th1([10.0, 6.0, 3.0], [0.5, 0.4, 0.3])
        shape[f"{era}/SR/data_obs"] = th1([10.0, 7.0, 4.0], [10.0, 7.0, 4.0])
        for d, k in (("Up", 1.03), ("Down", 0.97)):
            shape[f"{era}/SR/ttbar__ps_fsr{d}"] = th1([10.0 * k, 6.0, 3.0 / k], [0.5, 0.4, 0.3])
            shape[f"{era}/SR/{SIG}__ps_fsr{d}"] = th1([0.5 * k, 1.0, 1.5], [0.01, 0.01, 0.01])
    ch_shape = write_channel(tmp_path / "shape", "S", shape)
    ch_abcd = write_channel(tmp_path / "abcd", "T", abcd_hists(["2017", "2018"], nuis=[("ps_fsr", 1.02)]), ABCD)
    result, index = make_result(tmp_path, monkeypatch, analysis, {"S": ch_shape, "T": ch_abcd})
    cards = sorted(index["cards"][SIG]["Run2"])
    text = {c.name: c.read_text() for c in cards}
    assert "shapes * * " in text["S_Run2.txt"] and "S_SR_Run2 autoMCStats 10 0 1" in text["S_Run2.txt"]
    assert (result["results_dir"] / "cards" / "S_Run2.root").exists()
    work = tmp_path / "ws"
    work.mkdir()
    with open(work / "combined.txt", "w") as out:
        subprocess.run(["combineCards.py", "-S"] + [f".={c}" for c in cards], cwd=work, stdout=out, check=True)
    subprocess.run(["text2workspace.py", "combined.txt", "-o", "ws.root", "-m", "125"], cwd=work, check=True,
                   capture_output=True)
    assert (work / "ws.root").exists()


def test_a_core_nuisance_must_be_in_every_channel(tmp_path, monkeypatch, analysis):
    with_fsr = write_channel(tmp_path / "a", "T", abcd_hists(["2018"], nuis=[("ps_fsr", 1.02)]), ABCD)
    without = write_channel(tmp_path / "b", "U", abcd_hists(["2018"]), ABCD)
    with pytest.raises(ConfigError, match="U does not carry ps_fsr, a core nuisance"):
        make_result(tmp_path, monkeypatch, analysis, {"T": with_fsr, "U": without})


def test_a_nuisance_outside_the_core_is_in_the_channels_it_applies_to(tmp_path, monkeypatch, analysis):
    analysis.nuisances["CMS_eff_m_id_syst"] = {}
    with_muons = write_channel(tmp_path / "a", "T", abcd_hists(["2018"], nuis=[("CMS_eff_m_id_syst", 1.02)]), ABCD)
    without = write_channel(tmp_path / "b", "U", abcd_hists(["2018"]), ABCD)
    _, index = make_result(tmp_path, monkeypatch, analysis, {"T": with_muons, "U": without})
    cards = {c.name: c.read_text() for c in index["cards"][SIG]["Run2"]}
    assert "CMS_eff_m_id_syst" in cards["T_Run2.txt"] and "CMS_eff_m_id_syst" not in cards["U_Run2.txt"]


def _two_inputs(tmp_path, monkeypatch, analysis, eras_a, eras_b, abcd_b=ABCD):
    a = write_channel(tmp_path / "a", "T", abcd_hists(eras_a), ABCD)
    b = write_channel(tmp_path / "b", "T", abcd_hists(eras_b), abcd_b)
    monkeypatch.setenv("VBSVVH_STAT_RESULTS", str(tmp_path / "results"))
    path = tmp_path / "result.yaml"
    path.write_text(yaml.safe_dump({"name": "test", "signals": [SIG],
                                    "channels": {"T": {"input": [str(a), str(b)]}}}))
    result = load_result(path, analysis)
    return result, write_cards(result, analysis, Dictionary(analysis.master), log=lambda *_: None)


def test_one_channel_from_two_inputs_one_per_run(tmp_path, monkeypatch, analysis):
    result, index = _two_inputs(tmp_path, monkeypatch, analysis, ["2017", "2018"], ["2024"])
    assert sorted(index["fits"]) == ["Run2", "Run2Run3", "Run3"]
    summary = json.loads((result["results_dir"] / "cards" / "summary.json").read_text())
    assert summary["T"]["Run2"]["eras"] == ["2017", "2018"] and summary["T"]["Run3"]["eras"] == ["2024"]


def test_an_era_in_two_inputs(tmp_path, monkeypatch, analysis):
    with pytest.raises(ConfigError, match=r"era\(s\) \['2018'\] of .* also in another input"):
        _two_inputs(tmp_path, monkeypatch, analysis, ["2017", "2018"], ["2018", "2024"])


def test_inputs_of_a_channel_share_the_abcd_block(tmp_path, monkeypatch, analysis):
    with pytest.raises(ConfigError, match="differ in their abcd block"):
        _two_inputs(tmp_path, monkeypatch, analysis, ["2018"], ["2024"], dict(ABCD, background="qcd"))


def _fit(work, args):
    return subprocess.run(["combine", "-d", "ws.root", "-m", "125"] + args, cwd=work,
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True).stdout


def test_counting_cards_are_blinded_unless_unblind(tmp_path, monkeypatch, analysis):
    ch = write_channel(tmp_path / "ch", "T", abcd_hists(["2018"]), ABCD)
    result, index = make_result(tmp_path, monkeypatch, analysis, {"T": ch})
    card = index["cards"][SIG]["Run2"][0]
    text = card.read_text()
    assert "# BLINDED" in text and not any(l.startswith("observation") for l in text.splitlines())
    assert json.loads((result["results_dir"] / "cards" / "summary.json").read_text())["T"]["Run2"]["blinded"]
    work = tmp_path / "ws"
    work.mkdir()
    subprocess.run(["text2workspace.py", str(card), "-o", "ws.root", "-m", "125"], cwd=work, check=True,
                   capture_output=True)
    assert "Expected 50.0%: r <" in _fit(work, ["-M", "AsymptoticLimits", "-t", "-1", "--run", "blind"])
    observed = _fit(work, ["-M", "AsymptoticLimits", "-n", ".obs"])
    assert "No observed data 'data_obs' in the workspace" in observed and "Observed Limit" not in observed

    (tmp_path / "results").rename(tmp_path / "blinded")
    _, index = make_result(tmp_path, monkeypatch, analysis, {"T": ch}, unblind=True)
    dc = parse(index["cards"][SIG]["Run2"][0])
    assert dc.obs["T_A_Run2"] == 0.0 and dc.obs["T_B_Run2"] == pytest.approx(21.0)


def test_unblind_must_be_a_boolean(tmp_path, monkeypatch, analysis):
    ch = write_channel(tmp_path / "ch", "T", abcd_hists(["2018"]), ABCD)
    with pytest.raises(ConfigError, match="'unblind' must be true or false"):
        make_result(tmp_path, monkeypatch, analysis, {"T": ch}, unblind="yes")
