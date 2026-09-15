Title: Radius/Luminosity Plots

# Radius/Luminosity Plots

These are bar charts of yield, broken out by planet type. Each bar is one of
the 15 radius/insolation bins -- five radius classes, each split into hot,
warm, and cold -- averaged over the ensemble, with an error bar giving the
spread across simulations. They answer "what kind of planets did this mission
get, and how many."

The data comes from `reduce-radlum.csv` (the binned histograms) and
`reduce-earth.csv` (the single exo-Earth tally), both written by the yield
reduction in `reduce_drms.py`. Plots are named `det-radlum-*.png` and appear
in the *Radius/Luminosity Plots* section of the ensemble page.

## Reading the axis

The x-axis runs through the 15 bins in a fixed order: five groups of three,
one group per radius class, each labeled Hot / Warm / Cold. The pale gray
stripes behind the bars mark the radius groups, since only the insolation
label is written out.

**The leftmost bar is not one of the 15.** It is set off in green, labeled
with the short name of the earthlike class, and counts that class separately.
The class is a region of the plane rather than a bin, so it cuts across the
grid and its count is not the sum of anything to its right --
see [Tailoring the Earthlike Planet Class](earthlike-tailoring.html) for where
its bounds come from and how to move them.

Bars are stacked wherever a count has two components, and the inset text block
on each plot says what the segments mean.

## Unique, All, and Strict

Three of these plots count characterizations, and they differ along *two*
independent axes. Getting them confused is the easiest mistake to make here.

**Repeats.** A planet characterized on three separate visits can count once or
three times.

* *Unique* counts each planet once, on its first successful characterization.
  The reduction keeps a running set (`set_chars_uniq`) and only tallies
  planets not already in it.
* *All* counts every characterization event, repeats included. These are the
  `xchar` quantities in the reduction, parallel to `xdet` for detections.

**Bands.** An observation may characterize in several spectral bands at once,
and EXOSIMS reports a status per band: `1` for a full characterization, `-1`
for a partial one.

* *Full* means status `1` in **any** band. A planet fully characterized in
  blue but only partially in red still counts as full.
* *Partial* means status `-1`.
* *Strict* means status `1` in **every** band -- the reduction accumulates
  `charizations_strict` across the bands of one observation and keeps only
  planets that scored in all of them.

So the three plots are:

| Plot | Repeats | Bands |
|------|---------|-------|
| `radlum-char`        | unique | full (lower segment) + partial (upper segment) |
| `radlum-char-all`    | **with repeats** | full + partial, stacked the same way |
| `radlum-char-strict` | unique | **all bands full**; no partials at all |

Note that `-all` and `-strict` differ on *both* counts, so they are not a
nested pair. The nesting that does hold is the one the reduction comments
state, within the unique counts for the earthlike class:

```

    exoE_char_strict  <=  exoE_char_full  <=  chars_earth_unique
```

where the rightmost term is `exoE_char_full + exoE_char_part`, the total the
plain `radlum-char` plot shows.

Two consequences worth keeping in mind. A single-band scenario has nothing for
strictness to bite on, so `radlum-char-strict` collapses toward the full part
of `radlum-char` and the plot carries little information. And strict
characterizations have **no SNR** attached -- the SNR differs between bands, so
the reduction records none, which is why the SNR plot has no strict variant.

## The plots

Seven, always produced; this family has no extras.

| Plot | Shows |
|------|-------|
| `radlum-population`  | occurrence rate: planets of each type per star, over the whole target list |
| `radlum-det`         | unique detections, split by detection mode (blue / other) |
| `radlum-det-all`     | detections with repeats, same split |
| `radlum-char`        | unique characterizations, full + partial stacked |
| `radlum-char-all`    | characterizations with repeats, full + partial stacked |
| `radlum-char-strict` | unique characterizations, full in every band |
| `radlum-char-snr`    | mean characterization SNR per bin, full chars only |

`radlum-population` is the sanity check of the set: it is a property of the
simulated universe and the target list, not of the schedule, so it is where
you confirm planets are being generated at the intended rate. Its leftmost
green bar is the observed eta for the earthlike class.

## Caveats

* **Error bars are +/- 1 sigma across the ensemble**, not standard errors of
  the mean. They describe how much one simulation differs from another.
* **The population plot is normalized per star**; the yield plots are counts
  per simulation. Do not read a throughput by dividing one plot by another by
  eye -- the binned throughputs live on the radius/SMA bin plots instead.
* **A missing column skips a plot rather than failing.** Older reductions
  lack some of these quantities; the driver says which plot it skipped and
  why.
* **Bins with no planets are empty, not zero-yield.** A class absent from the
  simulated population produces nothing to characterize, which is a different
  statement from a class the mission failed to reach.
