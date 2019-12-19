#!/bin/sh -f
# copy current code snapshot

# options:
#   -h follow symlinks
#   ignore-failed-read passes over read-prohibited files in silence
#   exclude sims, Exosims source, ipyparallel

d=$(date +%Y-%m-%d_%H-%M)
fn=$HOME/Sandbox-Clones/Sandbox_$d.tgz
tar czf "$fn" -h --ignore-failed-read --exclude sims --exclude EXOSIMS --exclude ipyparallel --exclude spc --exclude exps .
ls -l "$fn"
