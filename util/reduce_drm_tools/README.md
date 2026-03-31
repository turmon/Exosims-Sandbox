# Reduce DRM Tools -- Components for data reduction

A registry-driven, modular system for reducing EXOSIMS DRM ensembles.

These notes are incomplete -- adapted from another project.

## Architecture

### Directory Layout

```
util/reduce_drm_tools//
  __init__.py               Package init
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


### The `config-reduce` Dict

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



