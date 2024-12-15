# EXOSIMS Sandbox

This is the code making up an EXOSIMS "Sandbox" -- a directory containing
scripts (JSON-encoded parameter files) that drive EXOSIMS simulation runs,
source code that performs runs and analysis, and support code for display
of results over a browser interface.

All the simulation scripts used by the Standards Team for the STDT final
report reside here.  In principle, multiple such Sandbox directories could
exist for separate families of simulations.

## Process

To use the sandbox for runs and analysis:

1. The EXOSIMS directory (a symbolic link) should point to Exosims source code (see also
under File Contents, below).

2. Put your Exosims input script, a json file, under Scripts.

3. Use `add-sims.sh` to add simulations to the ensemble run (see also File Contents, below).

4. Invoke `make S=<script> html`, where `<script>`
is the script name above, to reduce data, make plots, and update the web
page for that simulation.

5.  Reduced data and plots are placed in the `sims/<script>` directory.  The plots can be
viewed directly.

6. Plots can also be viewed on the generated webpage.  Start a server with
`make html-serve`, then invoke `make html-status` and point your
browser to the indicated location.

## Documentation

Further documentation is available on how execution works and how 
products are generated is available [on github](https://turmon.github.io/Exosims-Sandbox/),
or the [JPL github](https://github.jpl.nasa.gov/pages/turmon/EXOSIMS-sandbox/).

## File Contents

The contents here include:

* Files/links/directories you might want to alter
  + `EXOSIMS`     

    Symlink to the source code of Exosims used.  The file `EXOSIMS/EXOSIMS/__init__.py` should exist.

  + `Scripts`     

    .json scripts for Exosims input, placed here by convention.

* Files/links/directories you might want to inspect or run
  + `Makefile`    

    Controls data reduction and startup of ipython parallel engines.  See the Makefile header for actions.

  + `add-sims.sh`

    Driver script to add more ensemble members.  It contains usage instructions.

  + `sims/*`

    Dumped Exosims results, categorized by script-file root
    name.  These generated files are not part of the source code distribution.

* Files/links/directories that are mostly infrastructure

  + `Local`       

    Local modules to be imported by Exosims, including the run_one()
    method, and some graphical summarization scripts.

  + `util`       

    Utility scripts, mostly called at the direction of the `Makefile`.

  + `Experiments`

    Infrastructure for running coupled sets of Exosims runs for parameter tuning.

  + `ipyparallel`

    iPython-parallel per-user, per-machine configuration and lock files

## Extras

Optionally, use `make ipp-create` to create an iPython-parallel profile, 
and use `make ipp-start` to start an iPython-parallel cluster. 
Our own usage does not typically use this facility, so this is generally
not necessary.
 
