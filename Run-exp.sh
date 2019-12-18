#!/bin/bash

# EXOSIMS path
#OPTS="-p /proj/exep/rhonda/exosims_1.4Walker/EXOSIMS"
OPTS=

host=$(hostname | sed 's/.jpl.nasa.gov//')
# number of parallel jobs
if [[ $host == aftac1 ]]; then
  procs=26
else
  procs=34
fi

for s in $(cat Run-exp.$host); do
  echo "*** Processing $s"
  date
  ./add-sims.sh $OPTS -S Experiments/seed100.txt -P $procs $s 1
done
