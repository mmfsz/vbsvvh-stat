# Inputs for the VBS VVH statistical analysis

Each channel delivers two files in one directory: `shapes.root` and `channel.yaml`.
From them `vbsvvh-stat` builds the Combine datacards of every channel in the same way,
so that the channels can be fitted together.

Check a delivery before sending it:

```bash
source setup.sh                      # in the vbsvvh-stat repository
vbsvvh-stat validate path/to/channel.yaml
```

It exits with an error message for every problem found, and prints the MC yields
per era and region (data are never printed).

## `shapes.root`

One TH1 per era, region and process, with the sum of squared weights stored
(`TH1::Sumw2()`, or `hist` with `Weight` storage when writing with uproot):

| key | content |
|---|---|
| `<era>/<region>/<process>` | nominal |
| `<era>/<region>/<process>__<nuisance>Up` | systematic variation up |
| `<era>/<region>/<process>__<nuisance>Down` | systematic variation down |
| `<era>/<region>/data_obs` | data |

- **Eras**: `2016preVFP`, `2016postVFP`, `2017`, `2018`, `2024`, `2025`. All the eras
  of the channel go in the same file. `2022`, `2022EE`, `2023`, `2023BPix` are not
  implemented yet and are rejected.
- **Regions**: named freely. A region with one bin is a counting region; a region
  with several bins is fitted as a shape and needs the same binning in every era
  and for every process.
- **Signals**: one process per signal hypothesis, named exactly as in the
  analysis-wide list (`signals.yaml`): `Signal_C2V_1p0_C3_1p0`, `Signal_C2V_1p5_C3_1p0`,
  `Signal_C2V_1p0_C3_10p0`. The same name in every channel is what makes one signal
  strength per hypothesis across all channels.
- **Backgrounds**: every other process is an MC background (or a data-driven ABCD
  background, see below). Name them freely; the tool prefixes them with the channel
  tag in the cards, so a background is never shared between channels by accident.
- **Names**: letters, digits and single underscores, starting with a letter; `__`
  only as the separator before a nuisance name.
- **Nuisances**: only names from the analysis-wide list `nuisances.yaml` (CMS
  dictionary names where one exists). A source that is not correlated between eras
  carries the era in its name (`CMS_scale_met_unclustered_energy_2018`); a source
  correlated between eras has no era in its name and must be present in every era
  of the process.
  The nuisances marked `core: true` are carried by every channel; the others by the
  channels they apply to. To add a nuisance, see "Adding a nuisance" in the README.
  A nuisance named after the signal (`*_vbsvvh`) goes on the signal processes only.
- **Signal theory variations** (`QCDscale_fac_vbsvvh`, `QCDscale_ren_vbsvvh`, `ps_isr`,
  `ps_fsr`: `acceptance: true` in `nuisances.yaml`) are best delivered **raw**: the
  varied yield from the generator weight, cross section included, with nothing divided
  out. Say so in `channel.yaml` (`theory: raw`). The tool then writes them as the result
  asks, acceptance-only by default, dividing each era by the same variation's
  generator-level ratio (`theory_inclusive.yaml`), so every channel uses the same
  denominator. A channel that already divided by the generator-level ratio declares
  `theory: acceptance`. A channel that divided by its own post-skim yield (varied
  over nominal events after the skim, not the generator cross section) declares
  `theory: skim` and delivers that yield as an extra region named in `theory_norm`
  (signal nominal and Up/Down, no data); the tool multiplies it back out and divides
  by the generator-level ratio instead. The skim keeps high-x events, whose muF weight
  is further from 1 than the generator-level mean, so a post-skim denominator removes
  part of the acceptance change too (AllHad2FJ region A: ±3.5% instead of ±6–7.5%).
  The PDF row is not covered: its
  acceptance-only form needs each eigenvector divided before the quadrature sum, which
  only the channel can do.
- **Normalisation uncertainties** (cross sections, trigger, ...) are delivered like
  the others, as Up/Down histograms equal to nominal × (1 ± x).
- **Not to be delivered** (the tool adds them, and the names are reserved):
  luminosity (`lumi_*`), MC statistics (`CMS_HIG26013_<channel>_mcstat_*`), ABCD
  parameters and non-closure (`CMS_HIG26013_<channel>_<run>_abcd_*`).

## `channel.yaml`

```yaml
channel: MET1FJ                                     # channel tag
abcd: {background: bkg, A: A, B: B, C: C, D: D}     # only for an ABCD channel
theory: raw                                          # raw | acceptance | skim, see below
# theory_norm: all_events                           # only for theory: skim
provenance: {...}                                    # optional, free-form
```

- `channel`: the tag used in every name of this channel in the cards.
- `abcd`: for a channel whose background is estimated with the ABCD method from
  data. `background` is the process estimated from data; `A` (the signal region),
  `B`, `C`, `D` are the four counting regions of the plane, with A = B·C/D. The
  background histograms are optional: if present they are its MC prediction, used
  for the non-closure uncertainty; if absent, the non-closure is given in the
  result configuration.
- `theory`: how the signal theory variations were delivered: `raw` (varied yield over
  nominal yield, cross section included), `acceptance` (each era divided by the
  generator-level ratio of the same variation) or `skim` (each era divided by the
  varied over nominal signal yield of the region named in `theory_norm`). Required
  when the signal carries one.
- `theory_norm`: with `theory: skim` only, the region holding the post-skim signal
  yield the variations were divided by. It is read for the conversion and not
  written to the cards.
- `provenance`: anything that identifies how the inputs were made; it is copied into the results.
