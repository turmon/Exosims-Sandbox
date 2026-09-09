# Reduce DRM Tools -- Components for data reduction

A registry-driven, modular system for reducing EXOSIMS DRM ensembles.

These notes are incomplete -- adapted from another project.

## Architecture

### Directory Layout

```
util/reduce_drm_tools/
  __init__.py               Package init
  utils.py                  Config-file loading, unit stripping
  PlanetBins.py             RpLBins: numeric planet binning, is_earthlike()
  PlanetNames.py            PlanetNames: display names for the earthlike class
```

The driver lives one level up:

```
util/reduce_drms.py    Registry, drm loading, dispatch loop
```

### Data Flow

1.`EXOSIMS` produces "drm" files: `sims/SCENARIO/drm/NNN.pkl`
2. Pickles are loaded
3. Histograms are made for each pickle, iterating over its contents
4. Histograms are Averaged across the ensemble.
5. `reduce_drms.py` writes CSV files: `sims/SCENARIO/reduce-TYPE.csv`

Steps 2 and 3 are done in separate streams of control using 
`multiprocessing` because file loads are much faster when done 
in parallel.

The file-naming convention uses printf-style templates with two `%s` slots:

- **Destination template**: `src_tmpl % (TYPE, "csv")` -- e.g.,
  `"sims/scenario/reduce-%s.%s" % ("times", "csv")`
  yields `sims/scenario/reduce-times.csv`


### Reduction Registry

Dispatch is table-driven via `PLOT_REGISTRY` in `plot_drm_driver.py`. Each
entry is a dict with some keys:

| Key          | Type       | Description                                       |
|--------------|------------|---------------------------------------------------|
| `name`       | str        | Short identifier (used with `--only`/`--skip`)    |
| `module`     | str        | Module name in `drm-reduce-tools``                 |
| `enabled`    | bool       | Whether the plot runs by default                  |

Current registry entries:

| name           | csv_files                       |
|----------------|---------------------------------|
| yield_times    | yield-time                      |
| fuel_used      | times                           |
| events         | events                          |
| event_counts   | event-counts, earth-char-count  |
| visit_times    | visit-time                      |
| time_used      | times                           |
| promote        | promote, promote-hist           |
| star_targets   | star-target                     |
| earth_chars    | earth-char-list                 |
| radlum         | radlum, earth                   |


### The `config-reduce.json` File

Per-scenario reduction customization, an optional file living in the scenario
directory (`sims/SCENARIO/config-reduce.json`).  It is loaded by
`utils.load_reduce_config(dirname)`, which looks in `dirname` and -- if
`dirname`'s parent is a `.fam` or `.exp` directory -- one level up, so a whole
family or experiment can share one file.  Absence is normal, not an error.

Top-level keys:

| Key                  | Consumed by            | Meaning                                     |
|----------------------|------------------------|---------------------------------------------|
| `earthlike`          | `PlanetBins.py`, `PlanetNames.py` | Definition and naming of the earthlike planet class |
| `RpL_bins`           | `PlanetBins.py`        | `Rp_bins` / `L_bins` radius-luminosity grid  |
| `reduce_info_extras` | `reduce_drms.py`       | Extra columns to append to `reduce-info.csv` |
| `graphics`           | `plot_drm_driver.py`   | Which plot families make their extra plots  |

Any key starting with `_` is a comment and is ignored.

#### The `earthlike` group

Two kinds of entry.  The **numeric** ones (`Earth_Rp_scaled`, `Earth_Rp_lo`,
`Earth_Rp_hi`, `Earth_SMA_lo`, `Earth_SMA_hi`) re-define which planets
`RpLBins.is_earthlike()` counts, and thus the numbers in the `reduce-*.csv`
files.  The **display-name** ones (`name`, `name_plural`, `name_adj`,
`name_short`, `symbol`) say what to call that class in plot titles, axis
labels, table headers, and HTML captions -- see `PlanetNames.py`.  They affect
no numbers, so changing a name needs only `make graphics`/`make html`, not a
re-reduction.

```json
{
  "earthlike": {
    "_comment": "a Sub-Neptune population: axis-parallel Rp/SMA box",
    "Earth_Rp_scaled": false,
    "Earth_Rp_lo": 1.00, "Earth_Rp_hi": 3.50,
    "Earth_SMA_lo": 1.796, "Earth_SMA_hi": 18.257,

    "name":        "Sub-Neptune",
    "name_plural": "Sub-Neptunes",
    "name_adj":    "Sub-Neptune-like",
    "name_short":  "SubNep",
    "symbol":      "\\mathrm{SN}"
    }
}
```

All five names are optional.  Given only `name`, the rest are derived
(`name + "s"`, `name + "-like"`, `name`, and an upright `\mathrm{...}` of the
short name).  Given none, the historical `Earth` / `Earths` / `Earthlike` /
`Earth` / `\oplus` are used, so un-customized scenarios are unaffected.

To check what a scenario resolves to:

```
$ util/reduce_drm_tools/PlanetBins.py  sims/SCENARIO     # the numbers
$ util/reduce_drm_tools/PlanetNames.py sims/SCENARIO     # the names
$ diff <(util/reduce_drm_tools/PlanetNames.py sims/SCENARIO) \
       <(util/reduce_drm_tools/PlanetNames.py)           # just the differences
```

Note that the CSV column names, output filenames, and dict keys keep their
historical `earth`/`exoE` spellings regardless of the configured name: they are
data plumbing, not labels.

#### The `graphics` group

Several plot families have *extra* plots, made only when `mode.op` contains
`+`: `make S=... graphics-extra` asks for all of them, plain `make S=...
graphics` for none.  A scenario that wants one family's extras, and only that
family's, says so here:

```json
{
  "graphics": {
    "_comment": "the planet-population throughput maps, and nothing else",
    "mode_op": {
      "planet_pop": "+"
      }
    }
}
```

The families that have extras to ask for are `star_targets`, `promote`,
`yield_times`, `event_counts`, `visit_times`, `time_used`, `earth_chars`, and
`planet_pop`; naming any other plot is harmless but does nothing.

Keys are `fnmatch` patterns over the driver's plot names -- `plot_drm_driver.py
--list` prints them -- and the **first match in file order wins**, so put the
specific ones first and use `"*"` as a scenario-wide default:

```json
  "mode_op": {"planet_pop": "+", "*": ""}
```

An entry applies to the plots it names, whatever the command line said, so the
config can turn a family's extras on during a plain `make graphics` and equally
turn them off during `make graphics-extra`.  Precedence, highest first: a
plot's own `mode` in `PLOT_REGISTRY` (code), this file, then `--mode_op`.

A pattern matching no plot draws a warning: it is a typo, and would otherwise
do nothing quietly.  This group is read by the graphics driver alone; it
changes no reduced numbers, so `make graphics` is enough to see its effect.


## Usage

### Driver Usage

```
python plot_drm_driver.py SRC_TMPL DEST_TMPL [options]
```

Options:
- `--only NAME` -- run only the named plot
- `--skip NAME` -- skip the named plot (repeatable)
- `--list` -- list all registered plots and exit
- `--mode_op OP` -- set the global `mode['op']` string
- `--pdf` -- also write PDF output
- `-v` / `--verbose` -- increase verbosity (repeatable)
- `-q` / `--quiet` -- minimal output

Typical Makefile invocation:

```bash
util/plot_drm_driver.py \
    "sims/MyScript/reduce-%s.%s" \
    "sims/MyScript/gfx/det-%s.%s"
```

Run a single plot:

```bash
util/plot_drm_driver.py \
    "sims/MyScript/reduce-%s.%s" \
    "sims/MyScript/gfx/det-%s.%s" \
    --only fuel_used -v
```

## Adding a New Plot Module


## Adding a New Reduction Module

### Step-by-step

1. **Create the file**: `util/plot_drm_gallery/plot_drm_NEWNAME.py`.
   Copy the skeleton from an existing simple module like `plot_drm_fuel_used.py`.

2. **Implement the plot function**: `plot_drm_NEWNAME(reduce_info, plot_data,
   dest_tmpl, mode)` following the contract below.

3. **Add `main()` and `__main__` block** for standalone use. Use
   `cs.load_csv_files()` to load CSVs and `pd.read_csv()` for
   `reduce-info.csv`.

4. **Register in PLOT_REGISTRY** in `plot_drm_driver.py`:
   ```python
   {
       'name': 'NEWNAME',
       'module': 'plot_drm_NEWNAME',
       'function': 'plot_drm_NEWNAME',
       'csv_files': ['your-csv-id'],
       'enabled': True,
       'mode': {},
   },
   ```

5. **Test standalone**, then test via the driver:
   ```bash
   # standalone
   util/plot_drm_gallery/plot_drm_NEWNAME.py SRC_TMPL DEST_TMPL
   # via driver
   util/plot_drm_driver.py SRC_TMPL DEST_TMPL --only NEWNAME
   ```

### Reduction Subclass Contract

**Signature:**
```python
def plot_drm_NEWNAME(reduce_info, plot_data, dest_tmpl, mode):
```

**Parameters:**
- `reduce_info` (dict) -- metadata from `reduce-info.csv`
- `plot_data` (list[DataFrame]) -- pre-loaded CSVs, in the order listed in
  `csv_files` in the registry entry
- `dest_tmpl` (str) -- output path template with two `%s` placeholders
- `mode` (dict) -- operation settings (`op`, `verbose`, `ext_list`)

**Return value:**
- `list[str]` -- list of basenames written (from `tracker.get_files()`)
- `[]` (empty list) -- plot was skipped (e.g., missing data columns)
- `None` -- hard error occurred



