Title: Characterization Observing Conditions

# Characterization Observing Conditions

This family shows the *observing conditions* at every attempted
characterization of an Earth-like planet -- working angle, delta magnitude,
phase, and SNR -- and which of those attempts failed. The question it answers
is not "how many planets did the mission characterize" but "what did the hard
ones look like, and is there a pattern to the failures."

The data source is `reduce-earth-char-list.csv.gz`, produced by the yield
reduction in `reduce_drms.py`. Plots are made by `plot_drm_earth_chars.py`,
are named `det-earth-char-wa-dmag-*.png` and `det-earth-char-hist-phi*.png`,
and appear in the *Characterization Observing Plots* section of the ensemble
page.

The planet class is the configurable "earthlike" class, so a scenario that
redefines it in `config-reduce.json` -- to a Sub-Neptune population, say --
gets that name in the plot titles instead.

## The data

**One row per planet per characterization observation.** A row is an
*attempt*, not a planet: an Earth-like planet whose star is characterized
three times contributes three rows, and a star with two Earth-like planets
contributes two rows per observation. Counts on these plots are therefore
counts of attempts, and do not match the yield counts elsewhere on the page.

The columns used are the observing conditions `WA` (mas), `dMag`, `phi`,
`char_SNR`, and `MV` (the *star's* V magnitude), plus four flags:
`is_success`, `n_success`, `is_deep`, and `is_promo`.

**Phi** is the Lambertian phase function, between 0 and 1. It is already
folded into `dMag`; a planet at an unfavorable phase is dimmer by
2.5 log<sub>10</sub> Phi magnitudes, which is why Phi gets two plots of its
own.

## Outcomes

Every row falls in exactly one of three outcomes, and the distinction between
the last two matters for reading the scatter plots:

| Condition | Outcome | Meaning |
|-----------|---------|---------|
| `is_success > 0` | **success** | this planet was characterized on this observation, in some band |
| `is_success == 0`, `n_success > 0` | **split** | this planet failed, but another Earth-like planet at the same star succeeded on the same observation |
| otherwise | **failure** | nothing was characterized |

Splits are not a curiosity: a multi-planet system is characterized as a
system, so whenever an observation succeeds for one planet and not another,
every planet in it lands in this middle category.

## Target populations

Two flags describe the *star*, and they are made exclusive before plotting:

* **deep-dive** -- `is_deep`: the star was on the deep-dive list.
* **promoted** -- `is_promo` and not `is_deep`: promoted to the
  characterization list, with deep-dive taking precedence.

Only the promoted population gets scatter plots -- the deep-dive panels were
retired in 2024 -- which is why every WA/dMag filename carries `-promo-`. The
histograms stack both populations.

**If the scenario has no deep-dive targets, the blue series is empty** and its
legend entry reads "(N/A)". That is the common case, not a failure.

## The "easy" box

Every scatter plot draws a red dashed rectangle at WA > 70 mas and dMag < 26,
and an inset counts the successes and failures inside it. These bounds are
fixed constants, not read from the scenario's instrument parameters, so on a
mission whose inner working angle or contrast floor differs markedly the box
is a fixed reference mark rather than a statement about feasibility. The
"failed" count in the inset covers only the plotted (promoted) population; the
"successful" count covers all of them.

## The scatter plots

All four share a layout: working angle on x, delta magnitude on y, with the
y-axis clipped at zero and the x-axis left to autoscale.

* **Grey `+` markers** are every *successful* characterization in the
  ensemble, across all populations. They are context -- the cloud the mission
  actually managed.
* **Colored dots** are the unsuccessful promoted attempts, shaded by the
  quantity the plot is named for.

One asymmetry is worth knowing. The SNR plot separates the two bad outcomes,
drawing failures as dots and splits as `x` markers, and its legend counts them
apart. The phi and vmag plots merge them: their colored dots are every
attempt that was not a success, splits included. The same plot region can
therefore look denser on the phi plot than on the SNR plot.

The shadings:

* **snr** -- `char_SNR`, the SNR actually achieved. Low values at small WA and
  large dMag are the expected shape.
* **phi** -- log<sub>10</sub> Phi, on a color scale pinned to [-2, 0], so
  everything at Phi < 0.01 saturates at the bottom. Rows with Phi exactly 0
  land there too.
* **vmag** -- the planet's apparent visual magnitude, `MV + dMag`. This is the
  one plot that depends on how bright the *host star* is.
* **zzz-dmag** (extra) -- the y coordinate is a *phase-adjusted* dMag, showing
  where each planet would sit if it were observed at a reference phase of
  Phi = 0.7. Planets already above that phase are not moved, and the
  correction is capped at 5 magnitudes. Comparing it against the plain phi
  plot separates "too dim to characterize" from "caught at a bad phase." The
  `zzz` in the name only sorts it last on the page.

## The phase histograms

Both histogram only the **failed** characterizations, stacked by population
(blue deep-dive below, orange promoted above).

* **hist-phi** bins Phi linearly, 25 bins across [0, 1]. A second y-axis on
  the right carries a line per population giving the *failure probability* in
  each bin -- failures over attempts at that phase. The bars say where the
  failures are; the lines say whether that phase is actually hostile or merely
  well-populated. Read the lines where the bars are tall enough to support
  them.
* **hist-phi-log** bins the same failures by 2.5 log<sub>10</sub> Phi, the
  magnitude penalty the phase imposes, in 0.4-magnitude bins from -6 to 0. The
  leftmost bar is a catch-all for everything below -6, including the Phi = 0
  rows. There are no rate lines on this one.

## The plots

Five are always produced; one more with `mode_op` set to `+`.

| Plot | Content | Extra? |
|------|---------|:------:|
| `wa-dmag-promo-snr`      | WA/dMag, shaded by achieved char SNR; splits drawn separately |   |
| `wa-dmag-promo-phi`      | WA/dMag, shaded by log Phi, scale pinned to [-2, 0]          |   |
| `wa-dmag-promo-vmag`     | WA/dMag, shaded by planet apparent V magnitude (`MV + dMag`) |   |
| `hist-phi`               | failed chars vs. Phi, with per-bin failure probability       |   |
| `hist-phi-log`           | failed chars vs. 2.5 log Phi, the phase magnitude penalty    |   |
| `wa-dmag-promo-zzz-dmag` | WA/phase-adjusted-dMag, shaded by Phi                        | x |

To get the extra plot for every family, `make S=... graphics-extra`. To ask
for it for this family alone -- and have it stick, rather than depending on
which target was last run -- add this to the `config-reduce.json` of the
scenario, or of the family or experiment above it:

```

  "graphics": {
    "mode_op": {
      "earth_chars": "+"
    }
  }
```

`earth_chars` is this family's name in the plot registry; `plot_drm_driver.py
--list` shows them all.

## Caveats

* **These are attempts, not planets.** Nothing here is a yield. A frequently
  revisited star dominates the point count.
* **Successes and failures are drawn from different populations.** The grey
  successes span the whole ensemble; the colored failures are promoted
  targets only. Do not read a success:failure ratio off the plot.
* **`char_SNR` is the achieved SNR, not the threshold.** The plots do not draw
  the scenario's SNR requirement, so where the success boundary should fall is
  not marked.
* **A planet can appear as both.** Success is per observation, so the same
  planet can be a grey `+` from one visit and a colored failure from another.
