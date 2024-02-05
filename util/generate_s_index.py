#!/usr/bin/env python
r'''
generate_s_index.py: generate the s_index.json file for an Experiment

usage:
  `generate_s_index.py [ -o outfile ] SCRIPT [...]`

where:

*  SCRIPT ... is a list of JSON scripts, or 
*  SCRIPT is a directory containing JSON scripts
    (if a directory, the directory is scanned for all its EXOSIMS scripts)
*  -o allows naming a specific output .json file

Only `.json` files that look like EXOSIMS input scripts are indexed.

The output file (-o) typically would not be given. If not, its path is 
inferred from the directory of the first given SCRIPT (or from the 
script-directory, if that style of argument is used). 

The full output filename will be the script directory name as determined
above, plus `s_index.json`.

Typical usage:

*  Usually best:
    `util/generate_s_index.py Scripts/ExampleExp.exp`
*  Also OK:
    `util/generate_s_index.py Scripts/ExampleExp.exp/*.json`

NOTE:
This is a "primitive" indexer that does not do anything clever about the actual
Experiment parameters. The `s_index.json` file lists all viable EXOSIMS scripts, 
but *the index does not break out any specific parameter values*.

That is, if iterating over telescope diameter and contrast, or scheduler coefficients,
symbolic names of these parameters will be entered in the s_index.json file generated
by other tools, like json-xform.py. Such index files will look like this:

``` json
{
  "diam": 5,
  "contrast": 1e-10,
  "metric": "A",
  "IWA": 2.0,
  "script_name": "H5_C1e-10_baseA_IWA2.0.json",
  "run_name": "H5_C1e-10_baseA_IWA2.0"
},

```

An entry produced by the present routine just has an arbitrary script number:

``` json
{
  "index": 11,
  "script_name": "H5_C1e-10_baseA_IWA2.0.json",
  "run_name": "H5_C1e-10_baseA_IWA2.0"
},
```

Thus, this is a backup mechanism to allow Sandbox tools to work without complaint. 
The downstream effect (not catastrophic) is that Experiment-wide tables like 
`reduce-yield-plus.csv` will not have row-by-row parameter information other than 
the textual Ensemble name (e.g., `H5_C1e-10_baseA_IWA2.0`).

'''
# turmon aug 2023

import argparse
import sys
import os
import glob
import json
from datetime import datetime

# global verbosity mode, also usable for debugging print's
VERBOSITY = 0

########################################
###
###  Utility Functions
###
########################################

def ensure_permissions(fn):
    r'''Ensure correct permissions on the named data file.  
    We use rw-rw-r-- = 664 (octal), to allow group-write.'''
    try:
        os.chmod(fn, 0o664)
    except OSError:
        pass # e.g., don't own the file


########################################
###
###  Main routines
###
########################################

def exosims_sniff_test(script):
    r'''True if a script looks like EXOSIMS. No errors thrown.'''
    with open(script, 'r') as fp:
        try:
            specs = json.load(fp)
            return ('modules' in specs) and ('Observatory' in specs['modules'])
        except:
            return False
            

def dump(args, index):
    if not args.scripts and args.outfile == '.':
        print(f'{args.progname}: No scripts and assumed outfile: have no filename for output')
        return
    # translate -o . to a script filename if needed
    if args.outfile == '.':
        root = os.path.dirname(args.scripts[0])
        outfile = os.path.join(root, 's_index.json')
        args.outinfo = os.path.join(root, 's_index.txt')
    else:
        outfile = args.outfile
    # ensure enclosing dir exists
    directory = os.path.dirname(outfile)
    if directory and not os.path.exists(directory):
        os.makedirs(directory)
    # finally, write it
    with open(outfile, 'w') as fp:
        json.dump(index, fp, indent=2)
    ensure_permissions(outfile)
    print(f'{args.progname}: "{outfile}" now indexes {len(index)} scripts.')
    # put in a note claiming the file
    if args.outinfo:
        with open(args.outinfo, 'w') as fp:
            print(f'{args.progname} by {os.environ["USER"]} indexed {len(index)} scripts on {datetime.now()}',
                      file=fp)
        ensure_permissions(args.outinfo)


def load(args):
    # iterate over given scripts
    index = []
    skipped = 0
    for i, script in enumerate(args.scripts):
        # skip the special index file, always
        if script == 's_index.json':
            continue
        # give a warning if the named file isn't there, but don't die
        if not os.path.exists(script):
            print(f'{args.progname}: Named file {script} not found. Skipping.', file=sys.stderr)
            skipped = skipped + 1
            continue
        # only include if it looks like Exosims
        if args.auto and not exosims_sniff_test(script):
            skipped = skipped + 1
            continue
        # OK, it's in the index
        script_name = os.path.basename(script)
        run = os.path.splitext(script_name)[0]
        index.append({
            "script_name": script_name,
            "run_name": run,
            "index": i})
    # don't spook them
    #if skipped > 0:
    #    print(f'{args.progname}: Skipped {skipped} files.')
    return index
    

def main(args):
    index = load(args)
    dump(args, index)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generate basic s_index.json file for Experiment.",
                                     epilog='Typically you should use directory-name input. For more on usage, see the top of the file.')
    parser.add_argument('scripts', metavar='SCRIPT', nargs='*', default=[], help='JSON script files, or one script-dir')
    parser.add_argument('-o', '--outfile', metavar='FILE', help='output file name (default = SCRIPT dir)',
                            default='.')
    #parser.add_argument('-a', action='store_false', dest='auto',
    #                        default=True, help='index the named file only if it looks like an EXOSIMS script')
    parser.add_argument('-v', help='verbosity', action='count', dest='verbose',
                            default=0)
    args = parser.parse_args()
    
    # set umask in hopes that files will be group-writable
    os.umask(0o002)

    VERBOSITY = args.verbose

    # program name, for convenience
    args.progname = os.path.basename(sys.argv[0])
    
    # just leave this on always
    args.auto = True

    # check for dir
    args.outinfo = ''
    if len(args.scripts) > 0 and os.path.isdir(args.scripts[0]):
        dir_mode = True
        # outfile is in the dir itself
        args.outfile = os.path.join(args.scripts[0], 's_index.json')
        args.outinfo = os.path.join(args.scripts[0], 's_index.txt')
        # skip non-Exosims .json's
        args.auto = True
        # locate the other .json's
        #   the filter() is un-needed (it would be filtered later) but
        #   doing it here allows the correct number to be reported
        #   the sorted() ensures the ordering will at least be consistent
        #   from run-to-run
        args.scripts = sorted(
            filter(
                lambda s: not s.endswith('/s_index.json'),
                glob.glob(os.path.join(args.scripts[0], '*.json'))))
        print(f'{args.progname}: Full-directory mode ({len(args.scripts)} JSONs found)')

    # call the main
    main(args)
    sys.exit(0)
