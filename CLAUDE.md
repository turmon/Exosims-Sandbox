# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

EXOSIMS Sandbox is a post-processing and analysis framework for EXOSIMS (Exoplanet Observation Simulation). It manages the full pipeline: running simulations, reducing DRM (Design Reference Mission) output to CSV summaries, generating matplotlib plots, and serving results via a web interface.

The `EXOSIMS` symlink (gitignored) must point to a checkout of the EXOSIMS source code. Simulation scripts (JSON parameter files) go in `Scripts/`. All generated output goes under `sims/SCENARIO/`.

## Key Commands

All data reduction and plotting is orchestrated through `make` with `S=SCENARIO` (the basename of a JSON script file, without `.json`):

```bash
# Full pipeline: reduce + plot + generate HTML
make S=MyScript html

# Individual steps
make S=MyScript reduce          # Reduce DRMs to CSV files
make S=MyScript graphics        # Generate matplotlib plots from CSVs
make S=MyScript html-only       # Regenerate HTML without re-reducing

# Single-DRM visualizations (N = 1, 2, 5, 10, 20, 50, 100, or T for all)
make S=MyScript path-movie-5    # Generate 5 path movies
make S=MyScript obs-timeline-10 # Generate 10 observation timelines
make S=MyScript keepout-5       # Generate 5 keepout plots

# Experiment targets (S = experiment name, not single script)
make S=MyExperiment exp-reduce
make S=MyExperiment exp-html

# Web server
make html-start                 # Start Apache httpd
make html-status                # Show running servers
make html-stop                  # Stop server

# Force re-run with -B flag
make -B S=MyScript graphics
```

Running simulations is separate from make: use `add-sims.sh` to submit simulation runs, which invokes `Local/sandbox_driver.py`.

## Architecture

### Data Flow

1. **Run**: `add-sims.sh` → `Local/sandbox_driver.py` → EXOSIMS → `sims/SCENARIO/drm/*.pkl` + `spc/*.spc`
2. **Reduce**: `util/reduce_drms.py` loads pickled DRMs, computes statistics → `sims/SCENARIO/reduce-*.csv`
3. **Plot**: `util/plot_drm_driver.py` reads CSVs, dispatches to modular plot functions → `sims/SCENARIO/gfx/*.png`
4. **Present**: `util/html-summary.py` generates browsable HTML → `sims/SCENARIO/html/`

### Directory Layout

- **`Local/`** — Mission-specific code: simulation driver (`sandbox_driver.py`), custom EXOSIMS modules (`EXOSIMS_local/`), web server config and assets (`www-resources/`, `www-service/`)
- **`util/`** — All analysis/reduction/plotting utilities, called by the Makefile
- **`util/plot_drm_gallery/`** — Modular matplotlib plotting modules, one per plot type
- **`Scripts/`** — JSON parameter files for EXOSIMS (gitignored except `sims.json`)
- **`Experiments/`** — Parameter sweep infrastructure with JSON transform configs and seed files
- **`sims/`** — All generated output (gitignored): DRMs, CSVs, plots, HTML per scenario
- **`docs/`** — Static documentation website

### Plot Gallery System

Plots are table-driven via `PLOT_REGISTRY` in `util/plot_drm_driver.py`. Each entry references a module in `util/plot_drm_gallery/` (e.g., `plot_drm_yield_times.py`). To add a new plot type:

1. Create `util/plot_drm_gallery/plot_drm_NEWNAME.py` with a function `plot_drm_NEWNAME(args, mode)`
2. Register it in `PLOT_REGISTRY` in `util/plot_drm_driver.py`
3. Shared styling helpers are in `util/plot_drm_gallery/common_style.py`

### Data Reduction

`util/reduce_drms.py` implements a map/reduce pattern over DRM ensembles with ~20 analysis modes (yield, funnel, timeline, fuel, delta_v, etc.). It uses multiprocessing for parallelization. `util/reduce_drm_sets.py` handles multi-ensemble reduction.

### Experiment System

Experiments in `Experiments/` define parameter sweeps across multiple ensembles. They use JSON transform files to vary parameters and seed lists (`seed100.txt`, `seed1000.txt`) for reproducibility. The Makefile has `exp-*` targets to process all ensembles in an experiment.

## Conventions

- Python 3 throughout; key dependencies: numpy, scipy, pandas, matplotlib, astropy
- Most utilities are standalone scripts with argparse and `-h` help
- Commit messages use prefixes: `Enhance:`, `Fix:`, `Add:`
- The `S=` variable to make accepts a script basename, a `.json` path, or a `sims/` directory path (the Makefile normalizes it)
- Set `OPENBLAS_NUM_THREADS` to control thread usage during parallel reduction
