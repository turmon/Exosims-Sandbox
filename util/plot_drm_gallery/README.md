# Plot Gallery -- Modular Plots for Ensembles

A registry-driven, modular plotting system for reduced EXOSIMS DRM ensembles.
Each module is a self-contained plot generator that can run standalone or be
dispatched by the central driver (`plot_drm_driver.py`).

## Architecture

### Directory Layout

```
util/plot_drm_gallery/
  __init__.py               Package init
  common_style.py           Shared PlotTracker class, title helper, CSV loader
  plot_drm_earth_chars.py   Earth characterization scatter plots
  plot_drm_event_counts.py  Detection/characterization event count bar charts
  plot_drm_events.py        Event timeline plots
  plot_drm_fuel_used.py     Fuel use and delta-v vs. time
  plot_drm_promote.py       Promotion funnel plots
  plot_drm_radlum.py        Radius-luminosity scatter/histogram
  plot_drm_star_targets.py  Star target observation plots
  plot_drm_time_used.py     Time allocation breakdown
  plot_drm_visit_times.py   Visit time distribution plots
  plot_drm_yield_times.py   Yield vs. time plots
  apply_planet_overlay.py   Planet overlay utility (not a plot module)
```

The driver lives one level up:

```
util/plot_drm_driver.py    Registry, CSV loading, dispatch loop
```

### Data Flow

1. `reduce_drms.py` produces CSV files: `sims/SCENARIO/reduce-TYPE.csv`
2. The driver (or standalone module) loads CSVs into pandas DataFrames
3. Plot functions receive DataFrames, produce matplotlib figures
4. Figures are saved to disk as PNGs (and optionally PDFs) via `PlotTracker`

The file-naming convention uses printf-style templates with two `%s` slots:

- **Source template**: `src_tmpl % (TYPE, "csv")` -- e.g.,
  `"sims/scenario/reduce-%s.%s" % ("times", "csv")`
  yields `sims/scenario/reduce-times.csv`

- **Dest template**: `dest_tmpl % (plot_stem, ext)` -- e.g.,
  `"sims/scenario/gfx/det-%s.%s" % ("fuel", "png")`
  yields `sims/scenario/gfx/det-fuel.png`

### Plot Registry

Dispatch is table-driven via `PLOT_REGISTRY` in `plot_drm_driver.py`. Each
entry is a dict with these keys:

| Key          | Type       | Description                                       |
|--------------|------------|---------------------------------------------------|
| `name`       | str        | Short identifier (used with `--only`/`--skip`)    |
| `module`     | str        | Module name in `plot_drm_gallery`                 |
| `function`   | str        | Function name to call within that module          |
| `csv_files`  | list[str]  | CSV IDs to load (passed as `plot_data`)           |
| `enabled`    | bool       | Whether the plot runs by default                  |
| `mode`       | dict       | Per-plot mode overrides (merged with global mode) |

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

### The `mode` Dict

Every plot function receives a `mode` dict -- a package to flow information
downward to the plot code. Standard keys:

- `op` (str) -- operation mode string. Empty string is normal; `"+"` requests
  extra/optional plots.
- `verbose` (int) -- verbosity level. 0 = quiet, 1 = normal, 2+ = debug.
- `ext_list` (list[str]) -- file extensions to write, e.g. `['png']` or
  `['png', 'pdf']`.

The driver constructs an overall mode dict from CLI flags, then merges any
per-plot `mode` overrides from the registry entry before calling the plot
function.

### The `reduce_info` Dict

Metadata from the top-level summary file `reduce-info.csv`, converted to a dict. 

Keys we use include `experiment` (scenario name) and `ensemble_size`
(number of runs/DRMs). Passed to every plot function for use in
titles via `cs.plot_make_title(reduce_info)`.


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

### Standalone Module Usage

Each plot module has its own `main()` and can run independently:

```bash
util/plot_drm_gallery/plot_drm_fuel_used.py \
    "sims/MyScript/reduce-%s.%s" \
    "sims/MyScript/gfx/det-%s.%s"
```

Options are the same subset: `--mode_op`, `-v`/`--verbose`, `-q`/`--quiet`.

Standalone modules use `cs.load_csv_files()` to load their own CSVs (the
driver does this centrally when dispatching).


## Adding a New Plot Module

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

### Plot Function Contract

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

**Must use** `PlotTracker` for all file output so files are tracked correctly.


### Styling Conventions

All plot modules should follow these conventions for visual consistency:

- **Figure size**: `(8.5, 5)`
- **Title**: two lines -- `cs.plot_make_title(reduce_info)` on top (scenario
  metadata), descriptive title below. Bold, ~12pt (`fontsize=11*1.1`).
- **Axis labels**: bold (`fontweight='bold'`)
- **Tick labels**: size 13
- **Grid**: on (`ax.grid(True)`)


## Associating Plots with Functions and Data

`PlotTracker` in `common_style.py` is the centralized funnel for
writing figures.  This allows tracking all output filenames, and
the function they came from.

`PlotTracker` allows us to generate a metadata table (in the standard `gfx`
directory with the plots) called `det_plot_list.csv`. This maps graphics
filenames to (1) the plot function that produced them; (2) the data they 
used. Here are some lines for example:

```
graphics_file,routine,csv_files
det-time-det-allplan-month.png,plot_drm_yield_times,yield-time
det-event-count-slew.png,plot_drm_event_counts,event-counts;earth-char-count
```

This table allows you to solve problems like:

- If a problem is seen in `det-time-det-allplan-month.png`, it can be
  found in the `plot_drm_yield_times` function.

- If you want to reproduce the plot in `det-event-count-slew.png`, you can
  find the data in `reduce-event-counts.csv` and `reduce-earth-char-count.csv`.

