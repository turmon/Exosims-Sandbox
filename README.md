# EXOSIMS Sandbox

This is an Exosims ensemble directory containing scripts (JSON-encoded
parameter files) that drive Exosims
simulation runs, source code supporting runs and analysis, and support code
for display of results over a browser interface.

All simulation scripts for the Standards Team for the STDT final report
reside here.
In principle, multiple such directories could exist for separate families of simulations.

## Process

To use the sandbox for runs and analysis:

1. The EXOSIMS directory should point to Exosims source code (see also
under File Contents, below).

2. Put your Exosims input script, a json file, under Scripts.

3. Use `add-sims.sh` to add simulations to the ensemble run (see also File Contents, below).

4. Invoke `make S=<script> html`, where `<script>`
is the script name above, to reduce data, make plots, and update the web
page for that simulation.

5.  Reduced data and plots are placed in the `sims/<script>` directory.  The plots can be
viewed directly.

6. Plots can also be viewed on the generated webpage.  Start a server with `make html-serve`.

Further documentation is available on how execution works and how 
products are generated is available [on github](https://turmon.github.io/Exosims-Sandbox/),
or the [JPL github](https://github.jpl.nasa.gov/pages/turmon/EXOSIMS-sandbox/).

## File Contents

The contents here include:

* Files/links/directories you might want to alter
  + EXOSIMS     

    Symlink to the source code of Exosims used.  The file EXOSIMS/EXOSIMS/__init__.py should exist.

  + Scripts     

    .json scripts for Exosims input, placed here by convention.

* Files/links/directories you might want to inspect or run
  + Makefile    

    Controls data reduction and startup of ipython parallel engines.  See the Makefile header for actions.

  + add-sims.sh 

    Driver script to add more ensemble members.  It contains usage instructions.

  + sims/*      -- Dumped Exosims results, categorized by script-file root name.

* Files/links/directories that are mostly infrastructure
  + Local       

    Local modules for Exosims, including the run_one() method.

  + util       

    utility scripts, mostly called from Makefile.

  + ipyparallel 

    ipython parallel per-user, per-machine configuration and lock files

## Extras

Optionally, use `make ipp-create` to create an iPython-parallel profile, 
and use `make ipp-start` to start an iPython-parallel cluster. 
Our own usage does not typically use this facility, so this is generally
not necessary.
 
