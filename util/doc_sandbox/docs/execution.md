# Pipeline Execution Scheme

At the top level, you need to perform simulations, and perform
data reduction and plotting. 
For the first, you use `add-sims.sh`, and execution of all
data reduction, plotting, and webpage generation is
is controlled by standard unix `make`.

This page describes how `add-sims.sh` and `make` 
are invoked in standard Sandbox usage patterns.

## Quick start

The top-level workflow to run 100 simulations from a canned
set of seeds (using N-fold parallelism on the local computer),
reduce their output, and
produce a summarizing webpage is simply:

	# Create the input script, Scripts/example.json
	# Enter a venv with: source ../Python-venvs/ENV/bin/activate
	add-sims.sh Scripts/example.json Experiments/seed100.txt
    make S=example reduce
	make S=example html

If you go away and forget whether the graphics are up to date because
you might have added more simulations, you can just:

    make S=example html 
 
The graphical plots and the html page containing them
will be refreshed if there were new simulations added
since the last time the webpage was made, otherwise, the `make` returns
immediately.

## Core Functions

### Fundamental operating concept

The short-circuiting of unneeded processing just described
is possible because the `Makefile` includes the dependency network
of the data and graphics.  The `html` target depends on the 
the `reduce` target, which finally depends on the `sims/example/drm`
directory, where the simulation outputs are placed.
Dependency management is a powerful feature of `make`.

This works automatically whenever downstream data (webpage) needs to change in
response to new upstream data (simulations).
However, if the underlying graphics code changes ("put the titles in 14 point") you
need to force the webpage refresh by supplying `-B` to `make`,
which forces the rebuild (mnemonic: B => force build):

    make -B S=example html 

The chain of data artifacts (DRM's, reduced data, graphics, html pages)
is connected by driver scripts throughout.
Thus, `make` calls a shell-script driver to do data reductions and produce graphics.
Typically the shell-script will enforce the naming conventions
on inputs and outputs, and then call a Matlab or Python script to do
the actual processing.
So, there are roughly three levels of abstraction: the `make` target,
the shell driver script invoked by `make`,
and the 'doing' routine, in matlab or python.

The `Makefile` lists all the targets at the top of the file, with an
explanation.  This document adds some context to that helpful summary.

### Adding simulations to the ensemble

Adding simulations is done outside of make with the `add-sims.sh`
shell script driver.
Its basic usage is simply:

    add-sims.sh SCRIPT N

where `SCRIPT` is the JSON script for Exosims,
and `N` is the number of sims
to add to the ensemble tied to `SCRIPT`.
More commonly, `N` is actually the filename of a 
standardized list of seeds, e.g., `Experiments/seed100.txt`.

This pushes down to call a driver script for
Exosims, `Local/sandbox_driver.py`.
The main function of this driver is to put result files in the proper place,
and to perform logging.
In particular, the observing sequence ("DRM") for each sim is named for
the seed and goes in one directory, and the star-planet configuration
for that sim ("SPC"), also named for the seed, goes in a separate directory.

In the above usage case, a number of independent jobs 
(maximum of about 40) is run in parallel 
at the shell level on the local host only using GNU `parallel`.

Two parallelization options to add-sims are noteworthy:

    -Z => run in parallel across mustang2/3/4, aftac1/2/3 (100-way parallel)
    -S => run in parallel across mustang2/3/4 (64-way parallel)

These are available as long as you have passwordless ssh set up
between the various machines above.
(Test this with `ssh mustang4 ls`
while logged in to mustang2, for example.)
In the above cases, the independent jobs are farmed out to the 
indicated machines by GNU `parallel`.

Scripts typically need to generate cache files when run the first
time.
We support this by allowing the `N` option above to refer not to a 
seed filename, but to give a seed itself:

    =SEED => run a single simulation with the given integer SEED
    =0 => initialize caches, but exit before running the simulation

For the single-seed cases above, the usual simulation chatter is
sent to standard output rather than being logged, to make it easier
to diagnose first-run problems.

So typically, when generating results from a new script, you run:

    add-sims.sh SCRIPT =0

first, to initialize the cache files, and only then run:

    add-sims.sh -Z SCRIPT Experiments/seed100.txt

to run 100 Monte Carlo simulations.

More options and further details in the output of
`add-sims.sh -h`, or in the `add-sims.sh` header.


### Generating basic result summaries

As described in the top-level documentation, data from simulation
runs is reduced by averaging across the ensemble and dumping the results
to CSV files.  The reduction itself is triggered by invoking

    make S=example reduce 

Typically the reduction does not have to be directly invoked, because the
generation of other products invokes reduction as a side-effect, if it is
needed.  The reduction script is multi-threaded internally.

Graphics are generated by
Matlab and Python scripts that generally only
access the CSV files.
The top-level invocation for the basic graphical suite is:

    make S=example graphics 
 
This is not multi-threaded and takes a few minutes.
As with reduction, this is not typically directly invoked.

The graphics can be summarized and organized into a single-script HTML
file (with linked graphics) using:

    make S=example html 
 
As noted, this is the usual top-level command for reduction/graphics as
well.  There is an overall summary or "Table of Contents" page,
for all ensembles, that is also
updated when this is run.


## Extra Functions

### Additional result summaries

Additional graphical products can be made for a sample of single DRMs,
showing things like instrument selection (versus time) and slews (in a lat/lon format).
Python code is used for these graphics, and DRMs are directly accessed
rather than whole-ensemble reduced data.

The `Makefile` targets for these are:

   `path-movie-N`: make "N" path-movies and final frames 
   `path-final-N`: make "N" final frames, only
   `obs-timeline-N`: make "N" observing-target timelines
   `keepout-N`: make "N" keepout-vs-time plots

where `N` is 1, 2, 5, 10, 20, etc.  If the `-N` part is omitted, 10 is used.
These run independently across `N` DRMs, so they can be parallelized by using the
`-j P` option to `make` for `P`-way parallelization.  Typical usage is

    make -j 10 S=example path-movie  
 
which makes 10 path movies (the default number) in parallel.

The dependencies in the `Makefile` are incomplete for these graphics.
That is, `make` knows how to build them, but the `html` target doesn't
automatically know to re-make the html index file when they are changed.
So, when any of these additional graphics are made or changed, you should
update just the html index with

    make S=example html-only
 
This re-makes the html file (which is quick) but does not refresh the
basic graphic suite (which takes minutes).  The html-maker will notice the
additional graphics generated as above and include them in the new index.



### Webserver support

Results are produced as graphic files and as a webpage that summarizes the
graphics.  You could view the html in several ways, but one easy way is to
start a server within the simulations directory.  A simple Python server
(`python -m http.server 8000`) is too slow and does not support some
video functions, so we use httpd (Apache).

To control this server, the `Makefile` has these targets:

*  `html-start`: start httpd server

*  `html-stop`: stop httpd server

*  `html-status`: show running servers, if any, by inspecting their log
	files.

Using `make html-status` tells you where to point your browser.
The server typically can stay up for months without trouble.


### Tuning experiment support

We had occasional need to run multiple linked ensembles to understand
performance of schedulers in response to different tuning parameters.
For instance, a scheduler might need to trade off slew time against
Brown completeness ("should we spend time slewing to a far-off
high-completeness target, or
integrate longer on a nearby target with smaller completeness").
So we seek to compute and maximize yield
across a selected set of tuning parameters.

This was handled by creating a series of Exosims input scripts, one for
each parameter setting, running
`add-sims.sh` for each, and performing data reduction as noted above.
This amounts to an outer loop around the above process.

Regarding execution, we would generate a list of the input script
filenames, and farm this list out to the machines available (aftac1,
aftac2, aftac3).  The detailed process is described in the file 
`Experiments/run-experiment-howto.txt`.

One important detail has to do with Monte Carlo noise.  When tuning the
scheduler, we want to randomize over multiple simulated universes to
provide robust parameter choices.  But, when comparing two parameter
settings, we want to use the same *set* of simulated universes for both, to
suppress the additional noise that would arise if the two settings used
independent simulated universes.  (This is a simple instance of stratified
sampling.)  So the sets of random number seeds used for each ensemble
across the tuning experiment should be the same.  This is one purpose of the
standardized seed files in the `Experiments` directory.

