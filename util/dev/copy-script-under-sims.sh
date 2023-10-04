#!/bin/bash
#
# update sims/ dir to have a copy of the script.json under
# each entry, so that the index.html page can link to it
# this is now done automatically by reduction, but doing
# it manually like this allows any new "make html-only"
# to pick it up as well.
#
# turmon aug 2023

# sample matches
#  sims/ExampleExp.exp/s_c1=0.99_c3=0.01/html
#  sims/HabEx_2.4m_Coro_ljsDDPC_20180722/html

# (this is the command I used to make /var/tmp/htmls,
# it takes a long time to run)
# for hdir in $(find sims -name html) ; do

for hdir in $(cat /var/tmp/htmls) ; do
  if [ ! -r $hdir/index.html ]; then
    echo "# No index.html: $hdir"
    continue
  fi
  # sims/ExampleExp.exp/s_c1=0.99_c3=0.01
  sim_root=${hdir%/html}
  # ExampleExp.exp/s_c1=0.99_c3=0.01
  scenario=${sim_root#sims/}
  # Scripts/ExampleExp.exp/s_c1=0.99_c3=0.01.json
  script=Scripts/$scenario.json
  # sims/ExampleExp.exp/s_c1=0.99_c3=0.01/reduce-script.json
  jf=sims/$scenario/reduce-script.json
  #
  canary=sims/$scenario/reduce-info.csv
  # 
  if [ ! -r $script ]; then
    echo "# No script: $script"
    continue
  fi
  if [ ! -r $canary ]; then
    echo "# No reduction: $script"
    continue
  fi
  if [ -r $jf ]; then
    echo "# Already: $script"
    continue
  fi
  echo cp --preserve=timestamps $script $jf
done
