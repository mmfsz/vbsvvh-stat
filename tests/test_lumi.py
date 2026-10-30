"""lumi.yaml against the LUM POG numbers it was copied from."""

import math

import numpy as np
import pytest

#: LUM POG Run 3 covariance (Dec 2025), in %^2, order 2022, 2023, 2024.
RUN3_COV = np.array([[1.9023, 0.2336, 0.2742],
                     [0.2336, 1.6403, 0.8959],
                     [0.2742, 0.8959, 2.5791]])
RUN3_YEARS = ("2022", "2023", "2024")


def _components(analysis, years):
    """(n_nuisances, n_years) matrix of the relative uncertainties in %."""
    rows = [[100 * (e["kappa"].get(y, 1.0) - 1) for y in years] for e in analysis.lumi.values()]
    return np.array([r for r in rows if any(r)])


def _total(analysis, year):
    return math.sqrt(sum((e["kappa"][year] - 1) ** 2 for e in analysis.lumi.values() if year in e["kappa"]))


@pytest.mark.parametrize("year,total", [("2016", 0.0123), ("2017", 0.0082), ("2018", 0.0084)])
def test_run2_year_totals(analysis, year, total):
    assert _total(analysis, year) == pytest.approx(total, abs=5e-5)


def test_run3_reproduces_the_covariance(analysis):
    m = _components(analysis, RUN3_YEARS)
    assert m.T @ m == pytest.approx(RUN3_COV, abs=0.02)


def test_run3_subsets_reproduce_the_published_subset_tables(analysis):
    # The POG's minimal tables for two years; our table restricted to the same
    # years must give the same covariance (it may use one more nuisance).
    published = {
        ("2022", "2023"): np.array([[1.38, 0.17], [0.0, 1.27]]),
        ("2023", "2024"): np.array([[1.15, 0.0], [0.56, 1.61]]),
    }
    for years, table in published.items():
        m = _components(analysis, years)
        assert m.T @ m == pytest.approx(table.T @ table, abs=0.03)


def test_run2_dilution_reproduces_the_integrated_table(analysis):
    # POG "three years integrated" column for 2016+2017+2018, from the per-year
    # table diluted by the Golden Legacy luminosity shares.
    lumi = {"2016": 36.31, "2017": 42.07, "2018": 59.56}
    total = sum(lumi.values())
    integrated = {"lumi_13TeV_1516_l": 1.0031, "lumi_13TeV_151617_l": 1.0018, "lumi_13TeV_15161718_l": 1.0064}
    for name, k in integrated.items():
        kappa = analysis.lumi[name]["kappa"]
        diluted = 1 + sum(lumi[y] / total * (kappa.get(y, 1.0) - 1) for y in lumi)
        assert diluted == pytest.approx(k, abs=1e-4)


def test_2025_alone(analysis):
    assert analysis.lumi_kappas("2025") == {"lumi_13p6TeV_2025": 1.05}


def test_both_2016_halves_share_the_2016_nuisances(analysis):
    assert analysis.lumi_kappas("2016preVFP") == analysis.lumi_kappas("2016postVFP")
