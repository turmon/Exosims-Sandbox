Title: Radius/Luminosity Plots

# Radius/Luminosity Plots

These are bar charts of yield, broken out by planet type. Each bar is one of
the 15 radius/insolation bins -- five radius classes, each split into hot,
warm, and cold -- averaged over the ensemble, with an error bar giving the
spread across simulations. They answer "what kind of planets did this mission
get, and how many."

The data comes from `reduce-radlum.csv` (the binned histograms) and
`reduce-earth.csv` (the focused tally for exo-Earths), both written by the 
reduction code `reduce_drms.py`. Plots are named `det-radlum-*.png` and appear
in the *Radius/Luminosity Plots* section of the ensemble webpage.

## Reading the axis

The x-axis runs through the 15 bins in a fixed order: five groups of three,
one group per radius class, each labeled Hot / Warm / Cold. The pale gray
stripes behind the bars mark the radius groups, since only the insolation
label is written out.

The leftmost bar is set off in green, labeled
with the short name of the earthlike class, and counts that class separately.
The class is a region of the plane rather than a bin, so it cuts across the
grid and its count is not the sum of anything to its right --
see [Tailoring the Earthlike Planet Class](earthlike-tailoring.html) for where
its bounds come from and how to move them.

Bars are stacked wherever a count has two components, and the inset text block
on each plot says what the segments mean.

## Unique, All, and Strict

Three of these plots count characterizations, and they differ in two
separate ways: Repeats and Bands.

**Repeats.** A planet characterized on three separate visits can count once or
three times.

* *Unique* counts each planet once, on its first successful characterization.
  The reduction keeps a running set (`set_chars_uniq`) and only tallies
  planets not already there.
* *All* counts every successful characterization, repeats included. These are the
  `xchar` quantities in the reduction, parallel to `xdet` for detections.

**Bands.** An observation may characterize in several spectral bands at once,
and EXOSIMS reports band-by-band status: `1` for a full characterization, `-1`
for a partial char within that band. Here's how we manage that detail:

* *Full* means status `1` in **any** band. A planet fully characterized in
  blue but only partially in red still counts as full.
* *Partial* means all status was `-1`. This is mutually exclusive with Full.
* *Strict* means status `1` in **every** band. The reduction accumulates
  `charizations_strict` across bands of that observation and tallies only
  planets that scored in all of them.

So the three plots are:

| Plot | Repeats | Bands |
|------|---------|-------|
| `radlum-char`        | unique | full (lower segment) + partial (upper segment) |
| `radlum-char-all`    | **with repeats** | full + partial, stacked the same way |
| `radlum-char-strict` | unique | **all bands full**; no partials at all |

Note that `-all` and `-strict` differ in *both* ways, so they are not a
nested pair. The nesting that does hold (see the reduction code, 
within the unique counts for the earthlike class):

```

    exoE_char_strict  <=  exoE_char_full  <=  chars_earth_unique
                                           =  exoE_char_full + exoE_char_part
```

Note, the rightmost term (`exoE_char_full + exoE_char_part`) is
the total that the plain `radlum-char` plot shows.

Two consequences worth noting:

* "Strict" offers nothing special for a single-band scenario,
  so `radlum-char-strict` equals the "full" segments of 
  of `radlum-char`, and the plot carries no extra information. 
* Strict characterizations have no SNR attached -- the SNR differs between bands, so
  the reduction records none. Thus, the SNR plot has no strict variant.

## Plot Summary

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
Simulated Universe and the target list, not of the schedule, so it is where
you confirm planets are being generated at the intended rate. Its leftmost
green bar is the observed eta for the earthlike class.

## Clarifications

* **Error bars are +/- 1 sigma across the ensemble**, not standard errors of
  the mean. They describe how much one simulation differs from another.
  Divide by the square root of the ensemble size to get the standard error
  of the mean.
* **The population plot is normalized per star**; the yield plots are counts
  per simulation. Do not read a throughput by dividing one plot by another by
  eye -- the binned throughputs live on the radius/SMA bin plots instead.
* **Missing reduction data skips a plot rather than failing.**  Summary CSVs
  from years-old reductions may lack some quantities. The driver says which 
  plot it skipped and why.
