#!/usr/bin/env python
#
# Check EXOSIMS input script "specs" for keywords that are unused.
#
# Usage:
#   input_script_check.py script.json
#
# Notes:
#   * JPL scripts use a localized runner, so you will need to put the sim-runner
#     in your PYTHONPATH with:
#
#     PYTHONPATH=Local util/input_script_check.py ...
#    
#   * You can also give a list of scripts, and they will be cycled through
#   and checked one-by-one.
#
# Use with -h for more usage help.

# turmon 2023-09 adapted from an original by Dmitry Savransky
#
# notes: an alternative to requiring PYTHONPATH=Local is to just
# use the prototype SurveyEnsemble

import json
import sys
import argparse
import copy
from typing import List, Dict, Union, Any, Tuple

import EXOSIMS.MissionSim
from EXOSIMS.util.get_module import get_module
from EXOSIMS.util.keyword_fun import get_all_mod_kws, check_opticalsystem_kws


def parse_mods(specs: Dict[str, Any]) -> Dict[str, type]:
    """Check for presence of all required modules in input specs and return list of
    module class types.

    Args:
        specs (str or dict):
            Either full path to JSON script or an :ref:`sec:inputspec` dict

    Returns:
        dict:
            dict of all module classes along with MissionSim
    """

    req_mods = [
        "StarCatalog",
        "PlanetPopulation",
        "PlanetPhysicalModel",
        "OpticalSystem",
        "ZodiacalLight",
        "BackgroundSources",
        "PostProcessing",
        "Completeness",
        "TargetList",
        "SimulatedUniverse",
        "Observatory",
        "TimeKeeping",
        "SurveySimulation",
        "SurveyEnsemble",
    ]

    if "modules" not in specs:
        print(f'Fatal: "modules" not present in specs', file=sys.stderr)
        sys.exit(2)
    mods = {}
    for k in req_mods:
        if k not in specs["modules"]:
            print(f'Fatal: Required module {k} not present', file=sys.stderr)
            sys.exit(2)
        else:
            mods[k] = get_module(specs["modules"][k], k, silent=True)

    mods["MissionSim"] = EXOSIMS.MissionSim.MissionSim

    return mods


def check_for_unused_kws(
    specs: Union[Dict[str, Any], str]
) -> Tuple[List[str], Dict[str, Any]]:
    """Check input specification for consistency with module inputs

    Args:
        specs (str or dict):
            Either full path to JSON script or an :ref:`sec:inputspec` dict

    Returns:
        tuple:
            unused (list):
                List of unused keywords
            specs (dict):
                Original input (useful if read from disk)
    """

    if isinstance(specs, str):
        with open(specs, "r") as f:
            try:
                specs = json.loads(f.read())
            except:
                print(f'Fatal: Error opening specs file "{specs}"', file=sys.stderr)
                raise

    mods = parse_mods(specs)
    allkws, allkwmods, ukws, ukwcounts = get_all_mod_kws(mods)

    unused = list(set(specs.keys()) - set(ukws))
    if "modules" in unused:
        unused.remove("modules")
    # unused keywords starting with _ are allowed (comments)
    unused = [k for k in unused if not k.startswith('_')]

    return unused, specs


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Check an input script for spurious entries."
    )
    parser.add_argument("path", nargs="*", help="Pathname(s) of input script(s)")
    parser.add_argument("-o", "--optical", action="store_false", help="Omit OpticalSystem load and check")
    parser.add_argument("-l", "--list", action="store_true", help="list spurious keys without text")
    args = parser.parse_args()

    exit_status = 0

    for f in args.path:
        # check top-level keys
        unused, specs = check_for_unused_kws(f)
        if len(unused) > 0:
            exit_status = max(1, exit_status)
            if args.list:
                print("\n".join(unused))
            else:
                if len(args.path) > 1:
                    print(f'For {f}:')
                print(
                    "These input keywords were not used in any "
                    "module init:\n\t{}".format("\n\t".join(unused))
                    )

        # check the optical system keys
        if args.optical:
            try:
                OS = get_module(
                    specs["modules"]["OpticalSystem"], "OpticalSystem", silent=True
                )(**copy.deepcopy(specs))
                out = check_opticalsystem_kws(specs, OS)
                if out != "":
                    exit_status = max(1, exit_status)
                    print(f"\n{out}")
            except:  # noqa: E722
                exit_status = max(2, exit_status)
                print(
                    "Could not instantiate OpticalSystem with this script, "
                    "likely due to missing files.",
                    file=sys.stderr
                )
    # exit
    sys.exit(exit_status)
