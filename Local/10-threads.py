# thread limits for remote engines -- ipython-parallel setup
# 
# this file is copied from Local/ upon profile creation
#
# Rationale:
# This limit on OPENBLAS_NUM_THREADS is an experiment (2/2018) to stop a bug on
# aftac3 (only) instances.  I'm not sure if this is the appropriate place
# to inject the environment variable limiting # threads.
#  | make ipp-status
#  | util/ipcluster-check.py --file ipyparallel/rhonda/aftac3/security/ipcontroller-client.json
#  | OpenBLAS blas_thread_init: pthread_create: Resource temporarily unavailable
#  | OpenBLAS blas_thread_init: RLIMIT_NPROC 4096 current, 1032001 max
# See: https://github.com/HazyResearch/deepdive/issues/595
# Verification on aftac3 with:
#  | for i in $(seq 80); do python -c 'import numpy' &  done
# Yet, setting:
#  | export OPENBLAS_NUM_THREADS=8
# stops it.  
# If this equals 64, then you can see trouble again (80*64 
# possible threads created).
# Note that (see ulimit) the per-user process limit is 4096, and
# it seems that if (# numpy's) x (# threads) > 4096, then this
# limit can be reached.  But, it depends on timing of process creation.
#
# This setting must be done *before* numpy is imported.
#
import os
os.environ["OPENBLAS_NUM_THREADS"] = "8"

