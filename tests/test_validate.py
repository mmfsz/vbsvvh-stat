"""validate(): every contract violation is reported with its own message."""

import pytest

from conftest import ABCD, abcd_hists, th1, write_channel
from vbsvvh_stat.config import load_channel
from vbsvvh_stat.inputs import validate


def _errors(tmp_path, analysis, hists, abcd=ABCD, raw_keys=()):
    ch = load_channel(write_channel(tmp_path, "T", hists, abcd, raw_keys))
    return validate(ch, analysis)[0]


def test_good_channel_passes(tmp_path, analysis):
    hists = abcd_hists(["2017", "2018"], nuis=[("CMS_scale_met_unclustered_energy_{era}", 1.02), ("ps_fsr", 1.05)])
    assert _errors(tmp_path, analysis, hists) == []


def test_summary_never_shows_data(tmp_path, analysis):
    ch = load_channel(write_channel(tmp_path, "T", abcd_hists(["2018"]), ABCD))
    summary = validate(ch, analysis)[2]
    assert summary and not any("data_obs" in line for line in summary)


@pytest.mark.parametrize("mutate,message", [
    (lambda h: h.pop("2018/A/Signal_C2V_1p5_C3_1p0__ps_fsrDown"), "only one of Up/Down"),
    (lambda h: h.pop("2018/B/data_obs"), "no data_obs"),
    (lambda h: h.update({"2018/A/Signal_C2V_1p5_C3_1p0__CMS_madeupUp": 1.1,
                         "2018/A/Signal_C2V_1p5_C3_1p0__CMS_madeupDown": 0.9}), "not in nuisances.yaml"),
    (lambda h: h.update({"2018/A/Signal_C2V_9p9_C3_1p0": 1.0}), "not in signals.yaml"),
    (lambda h: h.update({"2018/A/CMS_HIG26013_x": 1.0}), "reserved"),
    (lambda h: h.update({"2018/D/data_obs": 0.0}), "no data in ABCD region D"),
    (lambda h: h.update({"2018/A/bkg": th1([1.0, 2.0])}), "binning differs"),
    (lambda h: h.update({"2022/A/Signal_C2V_1p5_C3_1p0": 1.0}), "era 2022 is not implemented yet"),
    (lambda h: h.update({"2018/A/bkg__QCDscale_fac_vbsvvhUp": 2.2,
                         "2018/A/bkg__QCDscale_fac_vbsvvhDown": 1.8}), "act on the signal only"),
])
def test_bad_inputs_are_caught(tmp_path, analysis, mutate, message):
    hists = abcd_hists(["2017", "2018"], nuis=[("ps_fsr", 1.05)])
    mutate(hists)
    errors = _errors(tmp_path, analysis, hists)
    assert any(message in e for e in errors), errors


def test_missing_sumw2(tmp_path, analysis):
    hists = abcd_hists(["2018"])
    raw = hists.pop("2018/A/bkg")
    errors = _errors(tmp_path, analysis, hists, raw_keys=[("2018/A/bkg", th1(raw, sumw2=False))])
    assert any("no Sumw2" in e for e in errors), errors


def test_correlated_nuisance_missing_in_one_era(tmp_path, analysis):
    hists = abcd_hists(["2017", "2018"], nuis=[("ps_fsr", 1.05)])
    for r in "ABCD":
        for d in ("Up", "Down"):
            hists.pop(f"2017/{r}/Signal_C2V_1p5_C3_1p0__ps_fsr{d}")
    errors = _errors(tmp_path, analysis, hists)
    assert any("correlated nuisance ps_fsr missing in era(s) ['2017']" in e for e in errors), errors


def test_energy_tagged_nuisance_is_correlated_within_its_run_only(tmp_path, analysis):
    hists = {**abcd_hists(["2017", "2018"], nuis=[("CMS_res_j_13TeV", 1.03)]),
             **abcd_hists(["2024", "2025"], nuis=[("CMS_res_j_13p6TeV", 1.03)])}
    assert _errors(tmp_path, analysis, hists) == []
    for r in "ABCD":
        for d in ("Up", "Down"):
            hists.pop(f"2017/{r}/Signal_C2V_1p5_C3_1p0__CMS_res_j_13TeV{d}")
    errors = _errors(tmp_path, analysis, hists)
    assert any("correlated nuisance CMS_res_j_13TeV missing in era(s) ['2017']" in e
               for e in errors), errors


def test_abcd_on_a_shape_region(tmp_path, analysis):
    hists = abcd_hists(["2018"])
    for p in ("Signal_C2V_1p5_C3_1p0", "bkg", "data_obs"):
        hists[f"2018/A/{p}"] = th1([1.0, 1.0], [1.0, 1.0])
    errors = _errors(tmp_path, analysis, hists)
    assert any("must be a 1-bin (counting) region" in e for e in errors), errors


def test_unknown_key_in_nuisances_yaml(tmp_path):
    import shutil
    from vbsvvh_stat.config import REPO, Analysis, ConfigError
    for f in ("signals.yaml", "lumi.yaml"):
        shutil.copy(REPO / f, tmp_path / f)
    (tmp_path / "nuisances.yaml").write_text("ps_fsr: {cor: true}\n")
    with pytest.raises(ConfigError, match=r"ps_fsr: unknown key\(s\) \['cor'\]"):
        Analysis(tmp_path)
