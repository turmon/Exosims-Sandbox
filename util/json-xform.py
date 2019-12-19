#!/usr/bin/env python
#
# json-xform: Transform an EXOSIMS JSON script into a set of such scripts.
#
# Usage:
#   json-xform [-v] [-o FILE] SCRIPT XSCRIPT ...
#
# where:
#   SCRIPT is the JSON script used by EXOSIMS
#   XSCRIPT is a JSON script giving the transform
# and, optionally:
#   -o FILE to specify where files should go; FILE must
#      have one %s to receive a file index or name.
#   -v for some verbose output
# The input SCRIPT (first required argument) is transformed according
# to the given XSCRIPT(s), and written to the given file template.
# Also, an index file summarizing the output files and the changing
# variables is written to an index file.
#
# The XSCRIPT is in JSON, and encodes one of 3 types of transforms.
#   update =>
#     Replace values in selected fields in SCRIPT with whichever
#     ones were given in XSCRIPT.  Example:
#       XSCRIPT = {"dMagLim": 23}
#     means to replace dMagLim within SCRIPT with the above value.
#   sequence =>
#     Fields in XSCRIPT are given by lists, and for each list
#     entry, a separate version of SCRIPT is generated using that
#     particular entry.  Example:
#       XSCRIPT = {"dMagLim": [23, 24, 25]}
#     means to generate 3 versions of SCRIPT, one for each value above.
#   cross =>
#     Similar to sequence, but the Cartesian product of the lists
#     for each field is used.  Example:
#       XSCRIPT = {"dMagLim": [23, 24, 25], "minComp": [0, 0.1]}
#     means to generate 6 versions of SCRIPT, one for each combination
#     of values.
# By default, a transform is taken to be "update".  To specify
# otherwise, insert the property "xform_type" into the top-level of
# the XSCRIPT, e.g. XSCRIPT = {"xform_type": "sequence", ...}.  (This
# property should have been in the sequence and cross examples above.)
#
# These transform generalize in the following ways:
#   For all:
#     By default, the specified values just take the place of the old value in
#     the controlling SCRIPT.
#     - You can make the specified values be interpreted as giving proportions,
#     not absolute quantities, by including the directive:
#       XSCRIPT = {"xform_action": {"ppFact": "relative"} ...}
#     within XSCRIPT.  This is particularly helpful for sequence and cross.  In
#     this case, a subsequent {..."ppFact": 2, ...} is taken to be "use 2x
#     the base script's ppFact." 
#     - For vector values, you can "mask" some substitutions by including:
#       XSCRIPT = {"xform_action": {"coeffs": "masked"} ...}
#     and in this case a NaN in the given substitution will pass the original
#     value through, with no substitution in that entry.  So if
#       coeffs = [[1, NaN], [2, NaN]],
#     then 1 and 2 will be (in turn) substituted into the first entry of
#     the vector "coeffs", but the second entry will be left as it was.
#   update:
#     If multiple fields are given in XSCRIPT, of course all given
#     values are replaced.
#     If the fieldname in XSCRIPT begins with "+", the corresponding
#     value is added in to the field's value.  Example:
#       XSCRIPT = {"+modules": {"SurveyEnsemble": "IPClusterEnsemble"}}
#     replaces just the SurveyEnsemble module entry, leaving the other
#     modules intact. Similarly:
#       XSCRIPT = {"+scienceInstruments": [{"Name": "new-imager", "QE": 0.6}]}
#     adds a new science instrument to the existing list of instruments.
#   sequence:
#     If multiple fields are specified, the sequence is drawn from each, wrapping
#     around if lengths are unequal.  Example:
#       XSCRIPT = {"dMagLim": [23, 24, 25], "minComp": [0, 0.1]}
#     generates three SCRIPTs, with (23, 0), (24, 0.1), (25, 0).  In particular,
#     a length-1 list will act like a scalar constant of that value.
#
# Transforms can be chained.  For instance, the first transform ("update") adds a
# field, and a subsequent transform ("sequence") then generates a sequence of scripts,
# where each starts from the updated script given by applying the first transform.
# As another example, a sequence on one variable followed by a sequence on another
# variable is equivalent to a cross product of the two variables.
# An XSCRIPT file be just one bare transform, or it can contain a list of transforms,
# which are applied in the sequence given, and of course multiple XSCRIPT files can
# be supplied as arguments.
#
# Many output files can be created, so file naming is important.  The template
# in the -o option gives the general naming pattern (e.g., -o script%s.json).
# For each output script, its particular script-name is plugged in to the %s 
# slot in the template supplied by -o.
# By default, the script-name is just the index number of the generated script.
# For a single sequence transform of length 6, the script-names would then be
# 0, 1, ..., 5, and the output filenames would be script0.json, etc.
# The (optional) xform_name parameter in the XSCRIPT gives great control over
# naming.  It is a string given to the Python format() method, with a reference
# dictionary composed of the field->value mapping in the script.  So for example,
#   "xform_name": "_ppLarge"
# would give a script-name of _ppLarge, and a filename of script_ppLarge.json.
# Alternatively, using the "magic" {} specifier of the format mini-language,
#   "xform_name": "_pp{ppFact:4.2f}"
# would pick up the ppFact value of the script, and format it as a real number with
# two digits to the right of the decimal point, resulting in files named like
# "script_pp0.10.json", "script_pp0.80.json", etc.  If the XSCRIPT is:
#  { "xform_type": "cross",
#    "xform_name": "life{missionLife:.0f}Xmag{dMagLim:.1f}",
#    "missionLife": [1, 2, 3, 5.0, 6],
#    "dMagLim": [22, 23] }
# then the file names are:
#   script_life1Xmag22.0.json, script_life1Xmag23.0.json, script_life2Xmag22.0.json, ...
#   script_life5Xmag23.0.json, script_life6Xmag22.0.json, script_life6Xmag23.0.json
#
# When transforms are chained, names are also chained, with an intervening _
# character.  This allows a series of several sequence transforms to emerge
# with correct names.  If, in the end, more than one file is written with
# the same name (a name collision), a warning is issued.

# turmon
#   oct 2018, created

from __future__ import print_function
import argparse
import sys
import os 
import json
import math
import copy
import itertools
from collections import OrderedDict, Counter

def load_scripts(args):
    r'''Load the root script, and the transform scripts, from given files.'''
    d_root = json.load(open(args.script), object_pairs_hook=OrderedDict)
    xforms = []
    for xname in args.xforms:
        xform_next = json.load(open(xname), object_pairs_hook=OrderedDict)
        if type(xform_next) == OrderedDict:
            xform_next = [xform_next] # bare dict is shorthand for one-item list
        if type(xform_next) == list and all(type(item) == OrderedDict for item in xform_next):
            xforms.extend(xform_next)
        else:
            raise RuntimeError('xform script must be a dict or list-of-dicts: %s' % xname)
    # set expected keys in each xform if missing
    for xform in xforms:
        xform.setdefault('xform_type', 'update')
        xform.setdefault('xform_action', {})
        name = '' if xform['xform_type'] == 'update' else '{index:d}'
        xform.setdefault('xform_name', name)
    return d_root, xforms

def update_script(d, xform):
    r'''Return a fresh script starting from d, but with (key,values) updated as in xform.'''
    # create a copy so starting script is unchanged
    d_new = copy.deepcopy(d)
    # replace values as needed
    for f in xform:
        if f.startswith('xform_') and f != 'xform_name':
            # for internal use, do not propagate into the script
            continue
        elif f.startswith('+'):
            # special case: augment d_new[f_base] with xform[f]
            f_base = f[1:] # strip the + to get the destination fieldname
            if f_base not in d_new:
                raise RuntimeError('Cannot augment a non-existent field: %s' % f_base)
            if type(d_new[f_base]) == type(xform[f]) == OrderedDict:
                d_new[f_base] = update_script(d_new[f_base], xform[f]) # recurse
            elif type(d_new[f_base]) == type(xform[f]) == list:
                d_new[f_base].extend(xform[f])
            else:
                # e.g., if d_new[f_base] is an int or a string
                raise RuntimeError('Do not know how to augment the field %s' % f_base)
        else:
            # typical case is a plain update
            d_new[f] = xform[f]
    return d_new

def unmask(val0, vals):
    try:
        len(val0)
    except:
        raise RuntimeError('Masked transformantions only work on vector quantities')
    val1 = [] # return value is a list, one for each entry in "vals"
    for v_subs in vals:
        if len(val0) != len(v_subs):
            raise ValueError('Masked substitutions only work on equal-length vectors')
        val_new = v_subs
        for i in range(len(val_new)):
            if math.isnan(val_new[i]):
                val_new[i] = val0[i] # retain old value
        val1.append(val_new)
    return val1

def apply_one_xform(d, xform):
    r'''Apply one transform to one dict, returning the list of resulting dicts.'''
    # shorthand names
    xform = copy.deepcopy(xform) # we need to modify xform within here
    d_name = d.get('xform_name', '')
    if d_name: d_name += '_'
    x_type = xform['xform_type']
    x_name = xform['xform_name']
    x_action = xform['xform_action']
    x_fields = [key for key in xform.keys() if not key.startswith('xform_')]
    # implement relative (proportional) specification of values
    for f in x_fields:
        if f in x_action and x_action[f] == 'relative':
            try: 
                xform[f] = [d[f]*scale for scale in xform[f]]
            except TypeError:
                print('Error: Could not use relative values on the field "%s".' % f)
                raise
        if f in x_action and x_action[f] == 'masked':
            xform[f] = unmask(d[f], xform[f])
    # apply the selected transform type to d
    #   one branch must be taken, and it must define:
    #   d_new = list of new dicts
    #   f_new = set of fields with varying values
    if x_type == 'update':
        d_new_1 = update_script(d, xform)
        d_new_1['xform_name'] = d_name + x_name.format(**d_new_1)
        d_new = [d_new_1]
        f_new = set() # NB: values do not vary, so do not include fields
    elif x_type == 'sequence':
        n_seq = max(len(xform[f]) for f in x_fields)
        d_new = []
        # generate a list of n_seq new dictionaries, with keys in fields
        for i in range(n_seq):
            # this is the set of keys in this dict
            d_up = OrderedDict()
            for f in x_fields:
                d_up[f] = xform[f][i % len(xform[f])]
            d_new_1 = copy.deepcopy(d)
            d_new_1.update(d_up)
            d_new_1['xform_name'] = d_name + x_name.format(index=i, **d_new_1)
            d_new.append(d_new_1)
        f_new = set(x_fields)
    elif x_type == 'cross':
        values = [xform[key] for key in x_fields]
        d_new = []
        for i, value in enumerate(itertools.product(*values)):
            d_up = OrderedDict()
            for f, v in zip(x_fields, value):
                d_up[f] = v
            d_new_1 = copy.deepcopy(d)
            d_new_1.update(d_up)
            d_new_1['xform_name'] = d_name + x_name.format(index=i, **d_new_1)
            d_new.append(d_new_1)
        f_new = set(x_fields)
    else:
        raise ValueError('Do not know the transform type "%s"' % xform_type)
    return d_new, f_new

def apply_xforms(d_list, xforms):
    r'''Apply the list of transforms to each script in d_list.
    Return the list of transformed scripts, and the set of fields within
    those scripts that vary.'''
    d_now = d_list
    f_now = set()
    for xform in xforms:
        # each single script d in d_now results in a list of new scripts
        d_new = []
        for d in d_now:
            d_new_1, f_new_1 = apply_one_xform(d, xform)
            d_new.extend(d_new_1)
            f_now.update(f_new_1)
        d_now = d_new
    return d_now, f_now

def dump_scripts(args, scripts, fields):
    r'''Dump the list of scripts to the designated outfile template.'''
    files_written = Counter()
    if not args.outfile:
        return files_written
    # one dict per script, to summarize the script's info in an index file
    script_summary = []
    for index, s in enumerate(scripts):
        # remove this name key before writing
        name = s.pop('xform_name', str(index))
        # substitute space -> _ in case "name" had spaces
        outname = (args.outfile % name).replace(' ', '_')
        if files_written[outname] > 0:
            print('Warning: Name clash on writing "%s"' % outname)
        files_written[outname] += 1
        # record the final file basename, and the fields
        summ1 = OrderedDict()
        summ1['script_name'] = os.path.basename(outname)
        summ1['run_name'] = os.path.splitext(os.path.basename(outname))[0]
        for f in fields:
            summ1[f] = s[f]
        script_summary.append(summ1)
        # write the file
        with open(outname, 'w') as f:
            if args.verbose:
                print('Dumping to: %s' % outname)
            json.dump(s, f, indent=2)
    # write the summary file
    outname = args.outfile % 'index'
    with open(outname, 'w') as f:
        if args.verbose:
            print('Dumping index to: %s' % outname)
        json.dump(script_summary, f, indent=2)
    return files_written
        
def main(args):
    r'''Load script and transforms, apply transforms, write results.'''
    d_root, xforms = load_scripts(args)
    scripts, fields = apply_xforms([d_root], xforms)
    files_written = dump_scripts(args, scripts, fields)
    return files_written

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Transform SCRIPT according to XFORMs.", epilog='')
    parser.add_argument('script', metavar='SCRIPT', default='',
                            help='json script file for EXOSIMS')
    parser.add_argument('xforms', metavar='XFORM', nargs='*', default='',
                            help='transform script(s) in json')
    parser.add_argument('--outfile', '-o', type=str, default='script_xform%s.json',
                            help='Output file template (with %%s).')
    parser.add_argument('-v', default=False, action='store_true', 
                            dest='verbose', help='Verbosity.')
    args = parser.parse_args()
    args.progname = os.path.basename(sys.argv[0])
    
    # set umask in hopes that files/dirs will be group-writable
    os.umask(0002)

    if not os.path.exists(args.script):
        raise IOError('Given script file (%s) does not exist.' % args.script)
    if args.outfile:
        if '%s' not in args.outfile:
            raise ValueError('Need a %%s somewhere in the outfile template "%s"' % args.outfile)
        # ensure the directory
        outdir = os.path.dirname(args.outfile % 'test')
        if len(outdir) > 0 and not os.path.isdir(outdir):
            os.makedirs(outdir, 0775)
    # do it
    files_written = main(args)
    print('%s: Wrote %d new scripts into %d different files.' %
              (args.progname, sum(files_written.values()), len(files_written)))
    sys.exit(0)
