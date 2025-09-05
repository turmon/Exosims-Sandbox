# Detection parameter sweep 2025 August

Example parameter sweep setup.


## Introduction

The design here is implemented in "explode.py" which
iterates over various detection parameters

The cachedir is set in "explode.py" 
to use a dedicated area for this Experiment.

## Usage

See the Makefile for how this basic usage translates to commands.

The steps are:
* Find an appropriate template script and place it in "exosims_template.json"

* Edit "explode.py" to set the range of n_det_remove and n_det_min, and the 
cache directory.

* Edit this file if desired.

* Template out the scripts:
>> make scripts

* Inspect the resulting scripts in the Staging/ directory

* Iterate if needed

* When satisfied, move the scripts to the working directory (up one level):
>> make finalize


