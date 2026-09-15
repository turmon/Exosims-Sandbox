Title: Radius/SMA Bin Plots

# Radius/SMA Bin Plots

These plots put the same 5x3 planet classes on the radius vs.
luminosity-scaled-SMA plane, drawn to scale. Where the radius/luminosity bar
charts show *how many*, these show *where* -- each bin is a rectangle in its
true position, carrying its number inside it.

Two different programs contribute to this section, which is worth knowing
because they behave differently:

* `util/rad-sma-rectangle-bin-plot.py`, driven by
  `util/rad-sma-rectangle-plot-driver.sh`, draws the **binned tables** from
  `reduce-radlum.csv` and `reduce-earth.csv`. These are written as both `.png`
  and `.pdf`.
* `plot_drm_rad_sma_chars.py`, one of the gallery modules, draws the
  **scatter and density** of individual characterizations from
  `reduce-earth-char-list.csv.gz`. PNG only.

All are named `det-rad-sma-*.png` and appear in the *Radius/SMA Bin Plots*
section of the ensemble page.

## The plane

The x coordinate is the luminosity-scaled semi-major axis, `a / sqrt(L)`, not
the raw orbit. That is the coordinate the bins are defined in, and it is what
makes planets around different stars comparable. Both axes are logarithmic.

Drawn behind the data are the 5x3 bin rectangles and the outline of the
earthlike class. That class is a region rather than a bin, so it straddles bin
boundaries -- its placement, and its name in the titles, come from
`config-reduce.json`, described in
[Tailoring the Earthlike Planet Class](earthlike-tailoring.html).

## Counts and throughputs

The binned tables come in two kinds, and the distinction is the same one the
radius/luminosity plots make between full and strict:

* **Counts** (`char_full`, `char_strict`) -- the mean number of unique
  characterizations per simulation falling in each bin. *Full* means the
  planet scored a full characterization in at least one spectral band;
  *strict* means it scored one in every band. Strict is always a subset, so
  its numbers never exceed the full ones.
* **Throughputs** (`char_tput_full`, `char_tput_strict`) -- the proportion of
  the planets *present* in that bin that were characterized. The reduction
  forms this as

```

    h_RpL_char_tput_full = h_RpL_char_full / (Nstar * h_RpL_population)
```

  where `h_RpL_population` is the per-star occurrence rate, so the denominator
  is the expected number of planets of that type across the target list.

A throughput is a rate on an absolute scale, which makes it the right plot for
comparing scenarios; a count depends on how many planets the universe happened
to generate. Where a bin holds no planets the throughput is 0/0 and shows as
**NaN** rather than zero -- an empty cell means "that planet type was not
present", which is a different statement from "the mission reached none of
them".

For the deeper version of this question -- throughput as a smooth function of
position rather than per bin, and the separation of scheduler targeting from
response -- see the Radius/SMA Density plots, which estimate the same
conditional probabilities without binning.

## The plots

Seven, always produced; this family has no extras. The first five are the
binned tables, and are written as PNG and PDF both.

| Plot | Shows |
|------|-------|
| `rad-sma-population`       | occurrence rate per bin: planets of that type per star |
| `rad-sma-char-full`        | mean unique characterizations per bin, full in any band |
| `rad-sma-char-strict`      | the same, requiring full in every band |
| `rad-sma-char-tput-full`   | characterized / present, per bin |
| `rad-sma-char-tput-strict` | the same, strict |
| `rad-sma-scatter-point`    | one point per successful characterization, pooled over the ensemble |
| `rad-sma-scatter-density`  | the same points as a kernel density estimate |

The two scatter plots are drawn from the characterization list rather than the
binned histograms, so they show individual events and are not averaged per
simulation. They pool every simulation in the ensemble, which is why their
point counts are much larger than the per-bin means beside them.

## Caveats

* **The binned tables are per-simulation means**; the scatter plots are
  ensemble totals. The two are not on the same footing and their numbers
  should not be compared directly.
* **NaN is not zero.** See above -- it marks an absent planet population, and
  the distinction matters most in the corner bins.
* **Strict characterizations carry no SNR**, since SNR differs across bands.
  Nothing in this section is SNR-weighted.
* **The bin grid is fixed at 5x3.** It can be moved via the `RpL_bins` group
  of `config-reduce.json`, but not reshaped; the reduction and both plotting
  programs assume that shape.
