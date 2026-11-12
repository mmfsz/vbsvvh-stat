# vbsvvh-stat

Datacards, expected limits, limit plots and nuisance impacts for the VBS VVH
channels, from one `shapes.root` + `channel.yaml` per channel (see `INPUTS.md`).

Everything runs in one environment, a CMSSW area with Combine and
CombineHarvester at the versions recommended by the Combine developers.

## Installation

Once, on an EL9 machine with CVMFS:

```bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=el9_amd64_gcc13
cmsrel CMSSW_16_0_0
cd CMSSW_16_0_0/src
cmsenv
git -c advice.detachedHead=false clone --depth 1 --branch v11.1.0 \
    https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
git -c advice.detachedHead=false clone --branch v4.0.0 \
    https://github.com/cms-analysis/CombineHarvester.git CombineHarvester
scram b clean; scram b -j4
cd -
git clone --recurse-submodules https://github.com/mmfsz/vbsvvh-stat.git vbsvvh-stat
```

In every new shell:

```bash
export VBSVVH_STAT_CMSSW=/path/to/CMSSW_16_0_0
source vbsvvh-stat/setup.sh      # cmsenv there, puts vbsvvh-stat on PATH
```

## A result

A result is one YAML in `configs/`; git ignores every config there except
`configs/example.yaml`, so each analyzer keeps their own. For example:

```yaml
name: MET1FJ_Run2Run3
channels:
  MET1FJ: {input: /path/to/channel.yaml}
  # MET1FJ: {input: [/path/to/run2/channel.yaml, /path/to/run3/channel.yaml]}  # eras from several deliveries
signals: [Signal_C2V_1p5_C3_1p0, Signal_C2V_1p0_C3_10p0, Signal_C2V_1p0_C3_1p0]
impacts: {Run2: [Signal_C2V_1p5_C3_1p0], Run2Run3: [Signal_C2V_1p5_C3_1p0]}
# runs:    {Run3: {substitute: {2025: 2024}}}   # projection: 2025 built from 2024
# closure: {MET1FJ: {Run2: 0.12}}              # ABCD non-closure instead of the MC value
# fits:    {Run2Run3: [Run2, Run3]}             # instead of the default fits
# unblind: true                                 # observed data in the cards (default: blinded)
# theory: total                                 # signal theory rows with the cross section (default: acceptance)
```

By default the eras of each run are summed into one plane, one card per (channel,
run, signal), and the fits are each run alone plus all runs together. `runs:`
changes that per run: `substitute` (projections), `split: true` (one plane per
era), or a run of its own with explicit `eras:`. With `fits:`, only the planes of
those fits get cards.

A channel can take its eras from several deliveries (for example Run 2 and Run 3
exported separately): `input` is then a list. Each era must come from one of them
only, and they must have the same tag and `abcd` block. Each era goes to its run as
usual.

The results are blinded by default (see "Blinding" below).

```bash
vbsvvh-stat validate path/to/channel.yaml       # check a channel's inputs
vbsvvh-stat run configs/<result>.yaml           # cards, limits, plots
vbsvvh-stat impacts configs/<result>.yaml       # impacts listed under impacts:
```

`run --stage cards` stops after the cards; `run --plots-only` and `impacts --plots-only`
redraw the plots.

## Outputs

```
results/<name>/                 result.yaml, provenance.json, fits.json, limits.json
  cards/<signal>/<tag>_<run>.txt  datacards; cards/systematics.yml for the CMS name check;
  cards/summary.json              eras, substitutions and ABCD closure per channel and run
  plots/limits_<signal>.*         one row per fit
  plots/limits_<fit>.*            one row per signal hypothesis
  impacts/impacts_<fit>_<signal>.*         every nuisance
  impacts/impacts_<fit>_<signal>_abcd.*    ABCD rateParams and non-closure only
  impacts/impacts_<fit>_<signal>_noabcd.*  every other nuisance, on its own axis
  work/<fit>/<signal>/            combined card, workspace, Combine output and logs
```

`VBSVVH_STAT_RESULTS` overrides `results/`.

## How the cards are built

### Eras summed into one plane

By default the eras of a run are summed into one set of bins (one ABCD plane per
run). The yield of a summed bin is the sum of the eras' yields, and so is its
varied yield for each nuisance: every era contributes its varied yield if it
carries the nuisance, and its nominal yield if it does not. The kappa on the
summed bin is the summed varied yield divided by the summed nominal yield, for
Up and Down separately.

- A nuisance correlated between eras (no era in its name, e.g. `ps_fsr`) is
  carried by every era, so its kappa is the sum of all varied yields over the sum
  of all nominal yields.
- A nuisance of one era (e.g. `CMS_scale_met_unclustered_energy_2018`) varies only
  that era's yield, so its kappa is diluted by that era's share of the bin. For a
  9% effect in 2018 on a bin with 2017 = 1.1 and 2018 = 2.1 events:

      kappa = (1.09 * 2.1 + 1.1) / (2.1 + 1.1) = 1.0590625 = 1 + (2.1 / 3.2) * 0.09

Shape regions are summed bin by bin in the same way.

### Signal theory variations

The QCD scale and parton-shower variations of the signal (`acceptance: true` in
`nuisances.yaml`) move both its cross section and its acceptance. One result writes
them one way in every channel: `theory: acceptance` (the default; the signal strength
is relative to the generator cross section, which is factored out) or `theory: total`
(the cross-section part kept). A channel declares how it delivered them
(`channel.yaml`, `theory: raw | acceptance`), and the cards convert each era before
the eras are summed:

| delivered | result | each era's Up/Down |
|---|---|---|
| raw | acceptance | divided by the generator-level ratio |
| acceptance | total | multiplied by it |
| raw | total | unchanged |
| acceptance | acceptance | unchanged |

The ratios are in `theory_inclusive.yaml`, per nuisance, signal hypothesis and era:
the varied over the nominal generator-level cross section (`Runs/LHEScaleSumw`,
`Runs/PSSumw` of the NanoAOD skims), summed over the four VBS VVH processes with
their cross sections (`theory_xsec.yaml`). Rebuild it when the signal samples change:

```bash
vbsvvh-stat theory-table /path/to/Run2_Sig_skims /path/to/Run3_Sig_skims
```

A conversion that needs a ratio the table lacks stops the cards. For the VBS VVH
signal at LO the factorisation-scale ratio is about 0.90/1.12 at the SM point and
0.88/1.16 at kappa_2V = 1.5; the renormalisation scale is exactly 1 (no alpha_S in a
pure-electroweak LO matrix element).

### ABCD

Counting regions only. The data-driven background has rate 1 and absolute-yield
rateParams for B, C, D, initialised at the observed count N and allowed in
[N - 3 sqrt N, N + 3 sqrt N]; A is B·C/D. A non-closure lnN acts on the background
in A. From the summed MC background of the plane, delta = |1 - nA nD / (nB nC)|
with MC statistical error sigma, and the lnN is 1 + sqrt(delta^2 + sigma^2), so the
uncertainty is never smaller than the MC statistical precision of the closure test.
A value given in the result is used as 1 + delta.

### Blinding

Every fit is an Asimov fit (`-t -1`): Combine builds the dataset from the model alone,
so the blind results need no data. A blinded counting card (the default) therefore has
no `observation` line; its workspace holds no `data_obs`, and an observed fit stops with
"No observed data 'data_obs' in the workspace". The ABCD B, C, D counts stay in the
`rateParam` lines, which set the background prediction. Shape cards keep
`observation -1` and their `data_obs` histograms. `unblind: true` in the result writes
the observed data in every card.

### Luminosity

`lumi.yaml` holds the LUM POG tables as published. Each nuisance gives a kappa per
year; an era takes the column of its year (both 2016 halves are 2016, 2022EE is
2022, 2023BPix is 2023), and the per-era variations are summed like any other.

- Run 2: the "new recommendation scheme" (Cholesky), 3 parameters for 2016-2018.
  Diluted with the luminosity shares it reproduces the POG's integrated table.
- Run 3, 2022-2024: the POG's three-year table, with its generic `lumi_1/2/3` named
  by the years each acts on (`lumi_13p6TeV_222324`, `_2324`, `_2024`, CMS dictionary
  names). Keeping only the columns of the years present gives exactly the
  covariance of those years, so the same table serves any subset.
- 2025: 5%, uncorrelated with everything (preliminary).

The Run 2 Cholesky names and `lumi_13p6TeV_2025` are not in the CMS dictionary yet;
`systematics.yml` declares them `custom` with a description.

### Projections

`runs: {Run3: {substitute: {2025: 2024}}}` builds 2025 from the 2024 inputs: every
histogram (MC, variations, data) scaled by L2025/L2024, nuisances tagged `_2024`
renamed `_2025`, and the 2025 luminosity uncertainty attached. The stand-in then
enters the Run 3 plane like a real era. Its MC is the 2024 MC again, so the MC
statistical error of the closure is underestimated in a projection.

## Analysis-wide tables

- `signals.yaml`: the signal hypotheses and their plot labels.
- `nuisances.yaml`: every nuisance a channel may deliver.
- `lumi.yaml`: luminosities and the LUM POG uncertainty tables (see "Luminosity" above).
- `theory_inclusive.yaml`: generator-level ratios of the signal theory variations, built
  by `vbsvvh-stat theory-table` from the skims and `theory_xsec.yaml` (see "Signal
  theory variations" above).
- `systematics/`: the CMS nuisance dictionary (submodule; update with
  `git submodule update --remote systematics`).

## Adding a nuisance

A channel can only deliver nuisances listed in `nuisances.yaml`. To add one, for
example the muon identification efficiency:

1. **Name.** Look for the source in the CMS dictionary,
   `systematics/systematics_master.yml`, and use its name. Only a source that the
   dictionary does not cover gets a name of its own, with a `description`. A custom
   uncertainty derived specifically for this analysis starts with `CMS_HIG26013_`,
   like the ones the tool adds (`CMS_HIG26013_<channel>_mcstat_*`,
   `CMS_HIG26013_<channel>_<run>_abcd_closure`).
2. **Correlation.** The name decides it: the same name is one nuisance, fully
   correlated wherever it appears, in every era and channel; different names are
   independent. A source uncorrelated between eras carries the era in its name. The
   muon identification scale factors have a systematic part correlated between eras
   (`CMS_eff_m_id_syst`) and a statistical part per era (`CMS_eff_m_id_stat_2018`, ...).
   The working point is not part of the name: as in the CMS dictionary, channels using
   different working points carry the same nuisances, fully correlated.
3. **Add the entries** to `nuisances.yaml`, `{era}` standing for every era:

   ```yaml
   CMS_eff_m_id_syst: {}
   CMS_eff_m_id_stat_{era}: {}
   ```

   Do not mark it `core`: the core (the nuisances every channel must carry, such as
   JES, JER, pileup, b tagging) is decided by the analysis contacts. A nuisance
   outside the core is carried only by the channels it applies to.
4. **Deliver** its Up/Down histograms in `shapes.root` (see `INPUTS.md`), then
   `vbsvvh-stat validate` the channel and `run` a result combining it with the
   other channels.

## Tests

```bash
python3 -m pytest -q tests
```
