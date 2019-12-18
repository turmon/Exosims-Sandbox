# exosims path for remote ipython-parallel engines
# 
# this file is copied from Local/ upon profile creation

import sys
# ./EXOSIMS points to the main EXOSIMS package directory
#   in general, ./EXOSIMS/EXOSIMS/MissionSim.py should exist, so 
#   that import EXOSIMS.MissionSim works.
sys.path.append("./EXOSIMS")
# ./Local points to site-local modules
#    in general, ./Local/EXOSIMS_local/__init__.py should exist,
#    so that import EXOSIMS_local.LocalModule works, and in particular,
#    so that module dependencies (in a json script) like: 
#      "SurveyEnsemble": "EXOSIMS_local.IPClusterEnsembleJPL"
#    will work.
#    Note that the working directory is also added to sys.path,
#    so this may duplicate that, depending on how python is called.
sys.path.append("./Local")

