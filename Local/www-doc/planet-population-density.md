Title: Planet-Population Densities and Throughputs

# Planet-Population Densities and Throughputs

This series of plots shows *where in the radius/SMA plane* a mission's
planets are, and what fraction of them it reaches. Every plot in the series
lives on the same axes as the Kopparapu bin plots: planet radius against
luminosity-scaled semi-major axis, both logarithmic. The 5x3 bins and the
outline of the earthlike region are drawn behind the data as context.

The data source is `reduce-planet-population.csv.gz`, one row per planet per
simulation, produced by `per_planet_yield()` in `reduce_drms.py`. The columns
used here are `sma_scaled`, `radius`, and four flags: `det_ok`, `char_ok`,
`star_det_obs`, `star_char_obs`.

**Luminosity-scaled SMA.** The x coordinate is `a / sqrt(L_star)`, not the raw
orbit. That is the coordinate the bins and the earthlike class are defined in,
and it is what makes planets around different stars comparable. Raw SMA does
not belong on these axes.

## The populations

Four planet populations, and one *star* variable that sits behind them:

* **A** -- all planets in the table. Note this is not the whole simulated
  universe: the reduction keeps only planets around stars the mission actually
  visited, because the remainder would all read `det_ok = char_ok = 0`.
* **D** -- `det_ok`: the planet was successfully detected at least once.
* **C** -- `char_ok`: the planet was successfully characterized at least once.
* **B** -- *bycatch*, C \ D: characterized but never detected. The name is a
  fishing analogy: these planets are swept up as a side effect of working a
  star for other reasons.
* **S** -- the star. `star_det_obs` and `star_char_obs` say whether the star
  was *observed* that way, whatever came of it. S is not a planet property,
  and it is the scheduler's choice.

**The funnel is not a chain.** A -> D -> C would be two filters in series, but
characterization also draws directly from planets that were never detected, so
C has two feeder streams:

```

                 +--------------------------------------------+
                 |  S : the star                              |
                 |  star_det_obs / star_char_obs              |
                 |  (scheduler's choice; not a planet         |
                 |   property, but it gates both paths)       |
                 +---------------------+----------------------+
                        |                            |
                        v                            v
     A  ------------->  D  ---------------------->  C n D
    all    P(det|A)   det_ok      P(char|D)       "targeted catch"
  planets    |                                        ^
             |  1 - P(det|A)                          |  both are in C
             v                                        v
            A \ D  ---------------------------->  B = C \ D
         not detected      P(char|not det)        "bycatch"
```

Because S gates both arrows out of the star box, a plot of P(char | planet
present) mixes two unrelated processes: whether the scheduler chose to work
this star at all, and whether a planet of this radius and SMA shows up once it
did. The two are separable, and separating them is the point of the extra
plots below.

## Coverage

Every planet falls in exactly one of four cells. The table gives the cell
names, which population each belongs to, and -- for orientation -- typical
shares from one 100-simulation ensemble:

| `det_ok` | `char_ok` | cell            | in A | in D | in C | example share |
|:--------:|:---------:|-----------------|:----:|:----:|:----:|--------------:|
|    0     |     0     | missed          |  x   |      |      |         0.713 |
|    1     |     0     | detected only   |  x   |  x   |      |         0.194 |
|    0     |     1     | **bycatch (B)** |  x   |      |  x   |         0.058 |
|    1     |     1     | targeted catch  |  x   |  x   |  x   |         0.035 |

Read across: the *density* plots show where each column's population lies. The
*throughput* plots show the ratio between two columns, as a function of
position. The bycatch row is the one that A -> D -> C reasoning misses; in the
example ensemble it is 62% of all characterizations.

## Densities and throughputs

Two kinds of plot, with different estimators.

**Density.** A kernel density estimate of one population, in units of
probability per dex^2. It answers "where are these planets?" and integrates to
one, so densities of different populations are *not* directly comparable in
height -- each is normalized to its own count.

**Throughput.** A conditional probability P(numerator | position), also called
a throughput here: of the planets at this radius and SMA in the denominator
population, what fraction are also in the numerator population? These are
bounded in [0, 1], and all such plots use the same absolute 0-1 color scale, so
their colors mean one thing across the whole series.

A throughput is *not* computed as a ratio of two densities. Two kernel fits
choose their bandwidth and orientation from their own samples, so their ratio
is unbounded, blows up where the denominator thins out, and is not a
probability. Instead the numerator and denominator are kernel sums over the
same sample with the same kernel,

        P(num | x)  =  sum_i w_i(x) num_i  /  sum_i w_i(x)

which is the Nadaraya-Watson estimator of the indicator and lies in [0, 1] by
construction. Where the denominator holds less than a set fraction of its peak
kernel weight, the cell is left blank rather than shown as a ratio of two
nearly-zero numbers.

## The plots

Six are always produced; four more with `mode_op` set to `+`
(`make S=... graphics-plus`). All are named `det-planet-pop-*.png` and appear
in the *Radius/SMA Densities* section of the ensemble page.

| Plot | Quantity | Extra? |
|------|----------|:------:|
| `density-all`        | density of A, all planets in the table               |     |
| `density-det`        | density of D, detected planets                       |     |
| `density-char`       | density of C, characterized planets                  |     |
| `tput-all2det`       | P(detected \| planet present)                        |     |
| `tput-all2char`      | P(characterized \| planet present)                   |     |
| `tput-det2char`      | P(characterized \| detected) -- the targeted path only |   |
| `tput-all2starchar`  | P(star observed for char. \| planet present) -- *targeting* | x |
| `tput-starchar2char` | P(characterized \| star observed for char.) -- *response*   | x |
| `tput-nodet2char`    | P(characterized \| not detected) -- the *bycatch rate*      | x |
| `tput-char2det`      | P(detected \| characterized) -- *provenance* of the chars    | x |

Two identities tie the series together, and are worth checking on any new
scenario:

* `all2starchar` x `starchar2char` = `all2char`, in the aggregate. The first
  factor is scheduling, the second is response to a planet at that radius and
  SMA. If the targeting map is flat, target selection is planet-agnostic; if
  it is not, some planet property is driving the scheduler.
* `all2char` = `det2char` x `all2det` + `nodet2char` x (1 - `all2det`), in the
  aggregate: the targeted and bycatch channels, weighted by how much of the
  population each draws from.

## Caveats

* **The population is restricted to observed stars.** "All planets" is what the
  mission had the chance to observe, not the simulated universe. Totals from
  this file are not population totals.
* **The flags are cumulative over the mission.** `det_ok` means *ever*
  detected, so P(char | det) does not imply the detection came first. A
  question about ordering needs first-detection and first-char times, which
  this table does not carry.
* **Large populations are sampled.** The estimate is capped at a fixed number
  of points, drawn with a fixed seed so plots reproduce; the title says when
  this happened. Broad structure is stable across seeds, fine detail in thin
  regions is not.
* **The kernel smooths across boundaries.** Density appears just outside the
  earthlike region and outside the populated part of the plane. The outline is
  drawn over the density so the boundary stays visible.
