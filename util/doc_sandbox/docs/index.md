# EXOSIMS Sandbox

Our environment for running batches of Exosims simulations
is called the EXOSIMS Sandbox.
It has four main functional pieces:
simulation execution, data reduction, 
graphics generation, and webpage generation.
Also, there is some support for
ensembles of related simulations for parameter-tuning, 
utilities for data extraction, 
and a webserver for inspection of results.

In general, the workflow is to create an Exosims script, 
run an ensemble of simulations that dumps output files to disk,
reduce data from those simulations to
a set of CSV files, and generate plots.

Most of the workflow execution is controlled by a `Makefile`, and
thus invoked by running `make` from the shell with an appropriate target.
The components making up this tower of abstraction are described [elsewhere](execution.md).
The `make` mechanism both feeds off of and imposes a file structure,
which we describe next.

## Sandbox File Layout

We run simulations tailored to many mission scenarios.
We show HabEx here for example, other simulations have
the same layout.

    HabEx/
        README 
        Makefile             -- controls execution
        Scripts/             -- JSON scripts for Exosims
		  .../HabEx.fam/     -- family of related scripts
		  .../HabEx.fam/H4m.json -- typical input script 
        sims/                -- Ensembles of DRM outputs
		  .../HabEx.fam/H4m/ -- DRMs and data from above script
        add-sims.sh          -- adds more Exosims runs to sims/
        Experiments/         -- groups of ensembles
		  seed100.txt        -- standard set of 100 PRNG seeds
        util/                -- mission-generic utilities
        Local/               -- Sandbox internal configuration
	Python-venvs/
		README.txt           -- how to make a Python virtual environment
		exosims-3.2.0        -- venv containing EXOSIMS 3.2.0


The main conventions here are:

* EXOSIMS input scripts are all under the `Scripts` directory.
All runs using that script (an *ensemble*)
have output in the corresponding directory
underneath `sims`.
This convention allows all utilities to place reduced data, graphics,
and html index pages into standard locations.
* Python virtual environments are used for EXOSIMS model runs and
data reduction.
* Families of related scripts can be placed under a directory ending
in `.fam`; these families can be nested.
* For use cases like parameter sweeps, several related ensembles will
be generated for different settings.  We call this an *experiment*,
and Experiment scripts are contained in a directory ending in `.exp`, with
an index file `s_index.json` describing the parameter settings.
For more, see `HOWTO-run-experiment.txt`.
* The `Makefile` has targets (beginning in `exp-`)
applicable to *both* Families and
Experiments, such that they process all the members at once.
* Most reduction and graphics code lives in `util` and `Local`;
the latter is intended to be mission-tailored, but the tailoring
has turned out to be minimal.

### Simulation Output Files

The subdirectories of `sims` are of particular interest, because they store
Exosims outputs like DRMs, reduced data like CSV files, and graphical
outputs like plots and movies:

	drm/                   -- DRM outputs (time-ordered observation lists)
		.../NNN.pkl        -- a specific run as a python pickle
		...                -- (and a bunch more runs)
    spc/                   -- SPC (star-planet configuration) files
		.../NNN.spc        -- a specific SimulatedUniverse as a pickle
		...
	log/                       -- run logs
		console/NNN.out   -- console log from a run
		environ/NNN.txt   -- runtime environment from a run
		outspec/NNN.json  -- outspec from a run
    reduce-info.csv       -- CSV files summarizing the ensemble
    reduce-radlum.csv
    reduce-times.csv
    reduce-visits.csv, etc.
    reduce-pool.json
    html/                  -- html page(s) indexing graphics below
    gfx/                   -- various plots as images 
    tbl/                   -- tabulated data from reduction
    path/                  -- single-DRM movies/graphics
		.../NNN.mp4        -- observing sequence movie
		.../NNN-*.png      -- various graphics
	path-ens/              -- ensemble-wide observing sequences

Each Exosims run produces a single DRM, stored in `drm/` according to the
random-number seed, and a corresponding star-planet configuration
(the `SimulatedUniverse`), stored in `spc/`.
Depending on how it was run, it may produce a text runlog in `log/console`, and
a JSON outspec in `log/outspec`.

This set of outputs (DRM, SPC, logs) is made by
a tailored `SurveyEnsemble` class
(`Local/EXOSIMS_local/IPClusterEnsembleJPL2.py`),
and a shell-level driver and "run_one" (`Local/sandbox_driver.py`).
These have been generalized and updated from the Exosims standard files,
which are documented in EXOSIMS as subclasses following the
`SurveyEnsemble` prototype.

As mentioned, after an ensemble is run, a reduction script
is typically invoked that generates the CSV output files,
and subsequent graphics commands use the ensemble summaries,
or in some cases just single DRMs.

## Data Reduction

We separated data reduction from plotting so that reduction of
the ensemble could be done once, and then various plots could
be remade and tweaked quickly from the reduced data.
The python driver is `reduce_drms.py`.
It contains some re-usable code for loading the ensemble of DRMs
and for computing summaries (means, standard errors,
medians, quantiles, histograms).
With only one exception, `reduce_drms` produces averages over
the ensemble.  Note that this does not limit us to mean values!
For example, a histogram of yields is also an average over the ensemble:
bin number "N" of a yield histogram
contains the average, across the ensemble,
of this 0/1 variable: "N" exo-Earths were detected in the simulation.
 
Data reduction happens at the ensemble level, by reading in
the set of typically ~100 DRMs and corresponding SPCs, and putting
out CSV files such as the following non-exclusive list:

    reduce-info.csv    -- metadata
    reduce-times.csv   -- temporal summaries like fuel use 
    reduce-radlum.csv  -- histograms segmented by radius/luminosity
    reduce-visits.csv  -- histograms of revisits
	reduce-earth.csv   -- exo-Earth detection counts
 
Currently, there are 14 such files.

The CSV format is good at storing vectors, so the above CSV files are
all vectorized along some index set: time, radius/luminosity, or 
revisit-count.
This format enforces a certain discipline in how a summary should be
done: within a bin defined by the index set
("month 15 of the simulation", "radius in range X, luminosity in range Y"),
counts that fall in each bin are accumulated across DRMs.
Then, means, standard errors, and quantiles can be computed from
the counts.
So, all CSVs have two columns indicating the lower and upper bin
boundaries, and for each quantity
of interest, there is one column each for the mean, standard error, 
each quantile needed, and the number of ensemble members comprising
that mean.

As an example, the fields now stored in `reduce_times.csv` include:

    h_det_time_lo,         -- bin boundaries
    h_det_time_hi,
    h_det_time_all_mean,   -- cumulative detections
    h_det_time_all_std, 
    h_det_time_all_nEns, 
    h_det_time_unq_mean,   -- unique detections only
    h_det_time_unq_std, 
    h_det_time_unq_nEns, 
    h_det_time_rev_mean,   -- revisits only
    h_det_time_rev_std, 
    h_det_time_rev_nEns, 
    h_time_fuel_all_mean,  -- fuel use
    h_time_fuel_all_std, 
    h_time_fuel_all_nEns, 
    h_time_fuel_slew_mean, -- fuel used for slews
    h_time_fuel_slew_std, 
    h_time_fuel_slew_nEns, 
    h_time_fuel_keep_mean, -- fuel for station-keeping
    h_time_fuel_keep_std, 
    h_time_fuel_keep_nEns

In all, there are 89 fields in that particular file, and about 1100 summarized quantities
across all the files.

One exception to this format convention is the file `reduce-earth-char-list.csv`, which
simply contains a list of all attempted characterizations across the entire
ensemble (all observations for all simulations).

## Graphical Output

Graphical outputs are generated by Matlab and by Python.
Python is preferable, being un-encumbered by licensing, 
more powerful, and in the same language as EXOSIMS itself.
However, we are comfortable with Matlab plots, and the reduced
data in CSV is easy to load into Matlab, so we adopted a combined approach.

### Graphical Output Files: Ensemble From Matlab
   
To make the whole-ensemble
plots from the CSV files above,
you run `make S=(SCRIPT) graphics`. 
This runs a driver (in the shell) that invokes Matlab,
which reads the CSVs and invokes plotting m-files
(all in `Local/Matlab/mfile`) for each
plot flavor.  The Matlab functions include, for example:

	plot_drms_script.m       -- driver script
	plot_drm_det_times.m     -- detections-vs-time
	plot_drm_fuel_use.m      -- fuel-vs-time
	plot_drm_radlum.m        -- radius/luminosity
	plot_drm_signal_end.m    -- signals success by writing a file

The resulting output graphics are put into files like the following in 
`sims/<script>/gfx/`:

	gfx/det-detects.png
	gfx/det-cume-detects.png
	gfx/det-fuel.png
	gfx/det-radlum-det.png
	gfx/det-radlum-char.png
	gfx/det-radlum-det-all.png

Currently, 90 files are generated, almost all PNG image files.

See the [binned ensemble plot gallery](binned-ensemble-gallery.md).

### Graphical Output Files: Ensemble From Python
  
Other plots show the tour of characterizations
made by starshade missions.
These are made by the python script
 `util/ens-path-graphics.py`
 which is invoked by `util/ens-path-summary.sh`.
 The `make` target is `path-ensemble`. 

The resulting output files are put into files in 
`sims/<script>/path-ens`:

	path-map.png                  -- map-format plot of activity
	path-adjacency-lon.png    -- slews, targets ordered by lon 
	path-adjacency-lat.png    -- slews, targets ordered by lat
	path-visits.csv                 -- mean number of visits per star 
	path-slews.csv                 -- total slews between pairs of stars

See the [ensemble path plot gallery](path-ensemble-gallery.md).

The latter two CSV files are used by a javascript widget in the generated HTML
page to show a zoomable plot of slews across the ensemble.


### Graphical Output Files: Single DRM 

A final set of plots shows the observations, keepout, or slews for a single
DRM.
They are useful in analyzing observation-scheduling behavior, 
verifying that keepout constraints are honored, and related questions.
They are made by `plot-keepout-and-obs.py`, `plot-timeline.py`,
and `keepout_path_graphics.py`.

For a given DRM with a certain (integer) SEED,
these files are placed in the directory
`sims/<script>/path/`, and named like:

	SEED-obs-timeline.png  -- time-format plot of all observations
	SEED-obs-keepout-all.png -- time-format plot of keepout (all obs.)
	SEED-obs-keepout-char.png -- time-format plot of keepout (char only)
	SEED-obs-timeline-info.csv
	SEED.mp4 -- movie of all detection and char observations
	SEED-final.png -- final frame of above movie

Some other plots and summaries are placed in the directory
`sims/<script>/path/SEED-cume`, including:

	path-visits.csv
	path-slews.csv

which are used by the same javascript widget described above, to
show a zoomable plot of slews for just that DRM.

See the [single-drm gallery](single-drm-gallery.md). 


## Web Output

These diverse plots are stored as files in given locations, but are easier
to grasp when displayed in context.  The web output helps with this: it is
a web page devoted to that particular ensemble of runs.  To make the web page
plots, you run `make S=(SCRIPT) html`.  This became our default 
diagnostic path in understanding results from an ensemble.

The script that performs this is `util/html-summary.py`.   When pointed at
a particular ensemble directory, it locates *all*  graphical files below
that point.  For each located file, it pattern-matches the filename against
a table of patterns that gives the skeleton of a web page.  It then
iterates over the filled-in skeleton to generate a page divided into
sections such as "Resource Use" and "Planet Yield".

Additionally, it generates a top-level index into all the ensembles.  This
index page serves as a dashboard for the sandbox.
