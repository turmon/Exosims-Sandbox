# EXOSIMS Sandbox

This is an Exosims ensemble directory containing scripts, source code, and dumped results.

In principle, multiple such directories would exist for separate families of simulations.  
We chose to have all simulations for the Standards Team for the STDT final report reside here.

## Process

Before use:
  (1) The EXOSIMS directory should point to Exosims source code (see below).
  (optional) Use "make ipp-create" to create an iPython-parallel profile.
  (optional) Use "make ipp-start" to start an iPython-parallel cluster.

To use:
  (2) Put your Exosims input script, a json file, under Scripts.
  (3) Use "add-sims.sh" to add simulations to the ensemble run (see below).
  (4) Use "make S=<script> html", where <script> is the script name in (2), to reduce data, make plots, and update the web page.

Reduced data and plots are placed in the sims/<script> directory.  The plots can be
viewed directly, or on the generated webpage:
  (5) View summary plots by starting a server: "make html-serve"

Further documentation is available on how execution works and how 
products are generated is available [on github](https://turmon.github.io/Exosims-Sandbox/),
or the [JPL github](https://github.jpl.nasa.gov/pages/turmon/EXOSIMS-sandbox/).

## Contents

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

