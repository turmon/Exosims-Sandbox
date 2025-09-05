
The enclosed is a way to generate scripts systematically, using a
simply python for-loop approach.
 
Directory structure like:
 
$ ls -lR Sweep_det_202508.exp/ 
Sweep_det_202508.exp/:
total 0
drwxrwsr-x 2 turmon exosims 106 Aug  7 09:18 Design/
 
Sweep_det_202508.exp/Design:
total 24K
-rw-rw-r-- 1 turmon exosims  475 Aug  7 09:15 Makefile
-rw-rw-r-- 1 turmon exosims  779 Aug  7 09:18 README.md
-rw-rw-r-- 1 turmon exosims  11K Aug  7 09:00 exosims_template.json
-rwxrwxr-x 1 turmon exosims 2.7K Aug  7 09:05 explode.py*
 
Place this under your Scripts/ directory in your laptop Sandbox 
(so you would have a path .../Scripts/Sweep_det_202508.exp/...).
 
Take a look at the file README.md in the Design/ directory.
 
The approach is: use the Design/ directory as a work area.
There is a template script file (exosims_template.json). I used a
random script I had laying around. You would want to use one of the
same scripts you have been using that has correct pathnames, etc.,
your template. Just replace the one here with yours.

The explode.py script loads this template and iterates over a bunch
of n_det_remove and n_det_min values, writing one script for each
tuple.

By following the directions in README.md, you will end up
with a bunch of scripts. You can then do:

$ add-sims.sh Scripts/Sweep_det_202508.exp/s_det_min1_rm1.json =0
$ add-sims.sh -j 4 Scripts/Sweep_det_202508.exp/s_det_min1_rm1.json Experiments/seed46.txt

for each one, to get ensembles. And then you can reduce them all
in one command with:

$ make S=Scripts/Sweep_det_202508.exp exp-reduce 
$ make S=Scripts/Sweep_det_202508.exp exp-html-only

to reduce the data and make a summary webpage.

Note the exp-html-only rather than exp-html. The former does not
make all the graphics for each ensemble (doing so takes a long time with
exp-html).  If you want to do the reduction on a partway-done set
of ensembles (e.g., you've made 6 ensembles and are curious about
how they look), you can do the two commands above and it will just
process what is there.

If you want a couple of sample graphics, you can make them in the
individual run:

$ make S=Scripts/Sweep_det_202508.exp/s_det_min1_rm1.json html

and it will just make that *one* ensembles graphics.


