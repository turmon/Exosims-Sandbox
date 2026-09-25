Title: Tailoring the Earthlike Planet Class

# Tailoring the Earthlike Planet Class

Several plot families single out the "Earthlike" planet class
and give it a dedicated count, label, or title.
Both *the numerical definition of the class * and *what it is called* 
are set per scenario, in `config-reduce.json`, under the `earthlike` key.

The class is used for more than labeling. `RpLBins.is_earthlike()` decides
which planets are counted in the exo-Earth tallies during reduction, so
changing the geometry changes the counts, not just the captions.

## Where the file goes

`config-reduce.json` is read from the scenario directory. If it is not there,
exactly **one** level up is tried, and only if that parent is a `.fam` or
`.exp` directory. There is no further walk up the directory tree.

That one level of inheritance lets a family set the class once for all
its members. A scenario's own file does not merge with its parent's -- the
first file found wins entirely, so a scenario that overrides a family must
restate everything it wants.

The plots announce which file they used, once per run:

```

    rad_sma_common.py: Loaded reduction config: sims/MyFam.fam/config-reduce.json
```

## Planet-property Bounds

Five attributes place the class in the radius / luminosity-scaled-SMA plane.
The defaults describe the historical Earthlike class:

| Attribute | Default | Meaning |
|-----------|--------:|---------|
| `Earth_Rp_lo`     | 0.80  | lower planet-radius bound, Earth radii |
| `Earth_Rp_hi`     | 1.40  | upper planet-radius bound, Earth radii |
| `Earth_SMA_lo`    | 0.95  | inner bound on luminosity-scaled SMA, AU |
| `Earth_SMA_hi`    | 1.67  | outer bound on luminosity-scaled SMA, AU |
| `Earth_Rp_scaled` | true  | slope the lower radius bound (see below) |

**`Earth_Rp_scaled` is the one with a surprise.** Left true, the lower
radius bound is not a constant: it is `Earth_Rp_lo / sqrt(a)`, evaluated at
each planet's scaled SMA. This is the SAG13 sense of "earthlike", and it
matches the EXOSIMS `SurveySimulation` prototype. The region it cuts out is
therefore not a rectangle -- its floor slopes down as SMA grows, which is why
the outline drawn on the radius/SMA plots has the "Nevada" shape.

Set it false and the bounds become axis-parallel, giving a plain rectangle.
That is the setting to use when repurposing the machinery for a class
where the SAG13 scaling has no meaning.

## Names

The same group carries the display names, which are *not* bin geometry and
never affect a count:

| Key | Default | Used for |
|-----|---------|----------|
| `name`        | Earth     | singular noun, in titles |
| `name_plural` | Earths    | count noun |
| `name_adj`    | Earthlike | class-membership adjective |
| `name_short`  | Earth     | axis labels and other tight spots |
| `symbol`      | `\oplus`  | LaTeX subscript, as in eta |

You can usually simply set `name` by hand. *The rest are derived from it when
omitted:* the plural gains an `s`, the adjective gains `-like` unless the name
already contains a space, the short name copies the name, and the symbol is
built from the short name. Override a derived one only where the derivation
reads badly.

## An example

Redefining the class as an approximate Sub-Neptune population, with
axis-parallel bounds:

```

{
  "earthlike": {
    "_comment": "Default _scaled is true. Change to false for axis-parallel boundaries.",
    "Earth_Rp_scaled": false,
    "Earth_Rp_lo": 1.80,
    "Earth_Rp_hi": 4.00,
    "Earth_SMA_lo": 0.10,
    "Earth_SMA_hi": 18.257,
    "name": "Sub-Neptune",
    "name_short": "Sub-Nep"
  }
}
```

Every plot that names the class now says "Sub-Neptune", the leftmost column of
the radius/luminosity plots is labeled "Sub-Nep", and the class outline on the
radius/SMA plots moves to the new bounds.

Keys beginning with `_` are ignored, so `_comment` is the place to say why a
scenario differs. In some files, `_original_Earth_Rp_scaled` is left in the file
as a convention for recording the value you replaced.

## The 5x3 bins

A second group, `RpL_bins`, sets the bin edges (the number of edges is
one more than the number of bins):

| Attribute | Meaning |
|-----------|---------|
| `Rp_bins` | 6 planet-radius edges, giving 5 radius bins |
| `L_bins`  | 5x4 luminosity edges, giving 3 insolation bins per radius |

These are checked on load and must keep that 5x3 shape; the plots and the
reduction both assume it. Most scenarios leave them alone -- retuning the
class bounds is common, moving the bin grid is not.

## After changing it

**Geometry changes require re-reducing.** `is_earthlike()` runs during
reduction, so the exo-Earth counts in `reduce-earth.csv` and the binned
histograms in `reduce-radlum.csv` are already committed to the old class by
the time any plot is drawn. Re-run `make S=... reduce`, then the graphics.

**Name-only changes need just the graphics**, since the names ride through
`reduce-info.csv` to the plot titles.

A misspelled attribute in the config file is not silently dropped. 
Names that no group recognizes are collected and reported:

```

    rad_sma_common.py: Warning: Unused attribute(s) in .../config-reduce.json: Earth_Rp_low
```

## Other reduction-configuration information in this file

`config-reduce.json` can carry other settings. 

`reduce_info_extras`: 
Extracts named cells from the reduction code and places them into
the ensemble-summary `reduce-info.csv`. Once there, the rest of
the reduction code folds these extra values into the overall
`reduce-yield-plus.csv` sunnary file. 

`graphics.mode_op`: 
Requests a plot family's "extra" plots to be made by default when
`make graphics` is run. Example:
```
    {
      "graphics": {
        "_comment": "Make extra (mode_op '+') promotion plots",
        "mode_op": {
          "promote": "+"
        }
      }
    }
```
This is similar to running `make graphics-extra`, but more granular,
because we can request extra plots from only one plot family.

Both of these settings are separate from the "earthlike" class definition.

