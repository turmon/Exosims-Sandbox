#!/usr/bin/env python
"""drm-ls.py: Summarize DRMs to the terminal similar to "ls".

## Usage
```
  drm-ls.py [-lscDi] DIR_OR_DRM ...
```

Simplest usage:
```
  drm-ls.py sims/ENSEMBLE
```

## Options
```
 -l gives long-format output (extra columns)
 -s gives summary output: rollups only (no per-DRM output)
 -c gives CSV output instead of tabular output
 -i gives run info, from the environment logs (alone: instead of the listing)
 -D gives diagnostic output: a traceback for DRMs that fail to load
 -h gives help
```

Each argument is a Sandbox directory, or a DRM file.  Sandbox directories
follow these conventions, and the DRMs within them are found accordingly:

+ `X/drm/` exists: X is an ensemble, and its DRMs are `X/drm/*.pkl`
+ `X` is named `*.fam` or `*.exp`: a family or experiment, whose
  subdirectories are examined in turn, by these same rules
+ `X` is `sims` itself: treated as a family, so all DRMs are listed
+ otherwise, X is ignored (e.g., `Analysis/` within an experiment)

For convenience, an ensemble's `drm/` directory may also be given.  DRM files
given explicitly (e.g., `sims/ENS/drm/17*.pkl`) are listed as-is, grouped
under their ensemble.

## Listing

Output is a tree, indented like the `tree` utility, following the Sandbox
hierarchy: sims, families, experiments, ensembles, and (unless `-s`) the DRMs
within each ensemble, labeled by seed.  The line for each sims, family,
experiment, or ensemble is a rollup over all the DRMs below it: their count
(Ndrm), and the per-DRM mean of each column.  The columns are `Nobs`, the
number of observations (DRM entries) of any kind, and, with `-l`:

+ `Ndet_ok`: successful detections, summed over observations (so a planet
  re-detected on a revisit counts again)
+ `Nchar_ok`: successful (full, not partial) characterizations, summed likewise
+ `Nstar_det`: number of distinct stars with at least one successful detection

With `-c`, the same rows are given as CSV, with the full path and the kind of
each row (sims, family, experiment, ensemble, drm) in place of the tree.

A DRM that cannot be loaded -- typically, a pickle made with packages (or
versions) not present in this python installation -- is skipped with a
warning, and omitted from the rollups.  Warnings are collapsed to one line
per ensemble for each distinct error.

## Run info

The run info (`-i`) is a similar tree.  For each node, it gives the number
of ensembles below it (Nscen, for sims, families, and experiments), the
number of DRMs below it (Nens, noting how many have logs, if not all), and
summarizes the run start time, the EXOSIMS version, and the EXOSIMS path,
across those DRMs, as recorded in the environment log for each DRM (for
`ENS/drm/SEED.pkl`, this is `ENS/log/environ/SEED.txt`, or for older
ensembles, the header of `ENS/run/outseed_SEED.txt`).  Values that differ
among the runs are listed along with their counts (for the EXOSIMS path,
one per line).  Given alone, `-i` does not examine the DRMs themselves.
With `-l` or `-s`, the run-info tree comes before the listing.
"""


import argparse
import sys
import os
import os.path
import glob
import time
import traceback
import pickle
from collections import defaultdict
import numpy as np
import astropy.units as u
#from astropy.time import Time


# unpickling python2/numpy pickles within python3 requires this
PICKLE_ARGS = {'encoding': 'latin1'}

# global modes
CSV_OUTPUT = False

# column order for DRM summaries
KEY_ORDER = ['Nobs', 'Ndet_ok', 'Nchar_ok', 'Nstar_det']


############################################################
#
# Sandbox tree
#
############################################################

class Node(object):
    r'''One Sandbox directory: sims, family, experiment, or ensemble.

    Ensembles hold DRM filenames (drms); the others hold child Nodes (children).
    After loading, stats holds one summary dict per DRM, keyed by filename.'''
    def __init__(self, name, path, kind):
        self.name = name  # label in the tree
        self.path = path  # directory path
        self.kind = kind  # ensemble, family, experiment, or sims
        self.children = []
        self.drms = []
        self.stats = {}

    def all_drms(self):
        r'''All DRM filenames at or below this node.'''
        return self.drms + [fn for c in self.children for fn in c.all_drms()]

    def n_ensembles(self):
        r'''Number of ensembles at or below this node.'''
        return (self.kind == 'ensemble') + sum(c.n_ensembles() for c in self.children)

    def all_stats(self):
        r'''All per-DRM summary dicts at or below this node.'''
        return list(self.stats.values()) + [s for c in self.children for s in c.all_stats()]


def container_kind(d):
    r'''Return the kind of Sandbox container d is, or None if not one.'''
    if d.endswith('.fam'):
        return 'family'
    if d.endswith('.exp'):
        return 'experiment'
    if os.path.basename(d) == 'sims':
        return 'sims'
    return None


def ensemble_node(ens_dir, name):
    r'''Return the Node for ensemble directory ens_dir, with all its DRMs.'''
    node = Node(name, ens_dir, 'ensemble')
    node.drms = sorted(glob.glob(os.path.join(ens_dir, 'drm', '*.pkl')))
    return node


def sandbox_node(d, name):
    r'''Return the Node for Sandbox directory d, or None if d is not one.'''
    if os.path.isdir(os.path.join(d, 'drm')):
        return ensemble_node(d, name)
    kind = container_kind(d)
    if kind and os.path.isdir(d):
        node = Node(name, d, kind)
        # descend, ignoring non-Sandbox subdirectories
        for sub in sorted(os.scandir(d), key=lambda e: e.name):
            if sub.is_dir():
                child = sandbox_node(sub.path, sub.name)
                if child is not None:
                    node.children.append(child)
        return node
    return None


def args_to_trees(arglist):
    r'''Return a list of Nodes, one tree per argument (or per group of DRM files).'''
    trees = []
    for arg in arglist:
        arg = os.path.normpath(arg)
        if os.path.isfile(arg):
            # explicit DRM file: group with the previous one, if same ensemble
            ens_dir = os.path.dirname(os.path.dirname(arg))
            if not (trees and trees[-1].kind == 'ensemble' and trees[-1].path == ens_dir
                    and trees[-1].explicit):
                node = Node(ens_dir, ens_dir, 'ensemble')
                node.explicit = True
                trees.append(node)
            trees[-1].drms.append(arg)
            continue
        if os.path.basename(arg) == 'drm' and os.path.isdir(arg):
            # an ensemble's drm/ directory, given directly
            node = ensemble_node(os.path.dirname(arg), os.path.dirname(arg))
        else:
            node = sandbox_node(arg, arg)
        if node is None:
            print('%s: Warning: skipping %s: not a DRM, ensemble, family, or experiment'
                  % (os.path.basename(sys.argv[0]), arg), file=sys.stderr)
            continue
        node.explicit = False
        trees.append(node)
    return trees


def tree_prefix(depth):
    r'''Prefix for the line naming a node at the given depth.'''
    return '' if depth == 0 else '|  ' * (depth - 1) + '|- '


############################################################
#
# DRM loading and summaries
#
############################################################

def load_drm(fn):
    r"""Return (drm, None), with drm the list of observations in a DRM pickle.

    If the DRM cannot be loaded -- e.g., a pickle made with packages (or
    versions) not present in this python installation -- instead return
    (None, (reason, traceback-text)).  The caller reports it."""
    try:
        with open(fn, 'rb') as f:
            drm = pickle.load(f, **PICKLE_ARGS)
    except Exception as e:
        # anything can go wrong within an unpickle: return it, and go on
        return None, ('could not load (%s: %s)' % (type(e).__name__, e), traceback.format_exc())
    # we read something: examine it to verify it is a DRM.
    # we take a valid DRM to be either:
    #   (1) an empty list (could be a non-DRM, but tough luck)
    #   (2) a non-empty list containing dictionaries with
    #       detection statuses
    # non-DRMs can include any un-pickled thing, including
    # dictionaries, etc.
    if isinstance(drm, list):
        if len(drm) == 0:
            return drm, None # case 1 above
        if len(drm) > 0 and isinstance(drm[0], dict) and ('star_ind' in drm[0]):
            return drm, None # case 2 above
    return None, ('loaded, but not a DRM', '')


def detail_summarize(drm):
    r'''More detailed DRM summary.'''
    # successful detections - an array entry for each observation
    dets = np.array([np.sum(np.maximum(0, obs['det_status'])) if 'det_status' in obs else 0
                     for obs in drm])
    n_det = np.sum(dets)
    # successful characterizations - an array entry for each observation
    chars = np.array([np.sum(np.maximum(0, obs['char_status'])) if 'char_status' in obs else 0
                      for obs in drm ])
    n_char = np.sum(chars)
    # stars we visited, in order
    all_stars = np.array([obs['star_ind'] for obs in drm])
    all_star_det = set(all_stars[dets > 0]) # de-duplicate
    n_star_det = len(all_star_det)
    # package and return
    result = dict(Nobs=len(drm), Ndet_ok=n_det, Nchar_ok=n_char, Nstar_det=n_star_det)
    return result


def skip_report(node, failures):
    r'''Report the DRMs within node that failed to load, one line per distinct reason.

    failures is a list of (fn, (reason, traceback-text)).'''
    by_reason = defaultdict(list)
    for fn, why in failures:
        by_reason[why[0]].append((fn, why[1]))
    for reason, fails in by_reason.items():
        if len(fails) == 1:
            what = fails[0][0]
        else:
            what = '%d DRMs in %s' % (len(fails), node.path)
        print('%s: Warning: skipping %s: %s' % (os.path.basename(sys.argv[0]), what, reason),
              file=sys.stderr)
        if DIAGNOSE and fails[0][1]:
            # one traceback per group, from its first DRM
            print(fails[0][1], end='', file=sys.stderr)


def load_tree_stats(node, verbosity):
    r'''Load and summarize the DRMs in the tree at node.  Return count of unreadable DRMs.'''
    failures = []
    for fn in node.drms:
        drm, why = load_drm(fn)
        if drm is None:
            failures.append((fn, why))
        elif verbosity > 1:
            node.stats[fn] = detail_summarize(drm) # inspect the drm
        else:
            node.stats[fn] = dict(Nobs=len(drm)) # don't look within the drm
    skip_report(node, failures)
    n_bad = len(failures)
    for child in node.children:
        n_bad += load_tree_stats(child, verbosity)
    return n_bad


def tree_rows(node, keys, summary, depth=0):
    r'''Return table rows (path, kind, label, Ndrm, values) for the tree at node.'''
    stats = node.all_stats()
    if stats:
        means = ['%.2f' % np.mean([s[k] for s in stats]) for k in keys]
    else:
        means = [''] * len(keys)
    rows = [(node.path, node.kind, tree_prefix(depth) + node.name, str(len(stats)), means)]
    for child in node.children:
        rows.extend(tree_rows(child, keys, summary, depth + 1))
    if not summary:
        for fn in node.drms:
            if fn in node.stats:
                seed = os.path.splitext(os.path.basename(fn))[0]
                values = ['%.0f' % node.stats[fn][k] for k in keys]
                rows.append((fn, 'drm', tree_prefix(depth + 1) + seed, '', values))
    return rows


def table_print(trees, verbosity, summary):
    r'''Print the DRM summary table, one tree per Node in trees.'''
    keys = KEY_ORDER if verbosity > 1 else KEY_ORDER[:1]
    rows = [r for node in trees for r in tree_rows(node, keys, summary)]
    if CSV_OUTPUT:
        print(','.join(['path', 'kind', 'Ndrm'] + keys))
        for path, kind, _, ndrm, values in rows:
            print(','.join([path, kind, ndrm] + values))
        return
    # plain text: label column left-justified, numeric columns right-justified
    header = ['DRM', 'Ndrm'] + keys
    table = [header] + [[label, ndrm] + values for _, _, label, ndrm, values in rows]
    widths = [max(len(r[j]) for r in table) for j in range(len(header))]
    for r in table:
        line = r[0].ljust(widths[0])
        line += ''.join('  ' + r[j].rjust(widths[j]) for j in range(1, len(r)))
        print(line)


############################################################
#
# Run info (environment logs)
#
############################################################

def environ_filenames(fn):
    r'''Return the candidate run-environment logs for DRM filename fn, newest format first.'''
    ens_dir = os.path.dirname(os.path.dirname(os.path.abspath(fn)))
    seed = os.path.basename(fn).split('.')[0]
    return [os.path.join(ens_dir, 'log', 'environ', seed + '.txt'),
            os.path.join(ens_dir, 'run', 'outseed_%s.txt' % seed)]


def load_environ(fn):
    r'''Return the run environment for DRM filename fn as a dict, or None if no log.

    The log is either log/environ/SEED.txt, with "key: value" lines, or (older)
    run/outseed_SEED.txt, with the same lines prefixed by "# ", followed by the seed.'''
    for log_fn in environ_filenames(fn):
        try:
            with open(log_fn) as f:
                lines = f.readlines()
            break
        except OSError:
            pass
    else:
        return None
    env = {}
    for line in lines:
        if line.startswith('#'):
            line = line[1:]
        key, sep, value = line.partition(':')
        if sep:
            env[key.strip()] = value.strip()
    # older logs give the time as, e.g., "Sun May 10 19:24:13 2026": standardize it
    if 'time' in env:
        try:
            env['time'] = time.strftime('%Y-%m-%d %H:%M:%S',
                                        time.strptime(env['time'], '%a %b %d %H:%M:%S %Y'))
        except ValueError:
            pass
    return env


def info_lines(envs):
    r'''Return (key, value) lines summarizing the run environments envs.

    A value is a string, or a list of strings, to be printed one per line.'''
    lines = []
    if not envs:
        return lines
    # run start time: the range, since it differs for every run
    times = sorted(env['time'] for env in envs if 'time' in env)
    if times:
        span = times[0] if times[0] == times[-1] else '%s to %s' % (times[0], times[-1])
        lines.append(('time', span))
    # others: typically constant, but list each distinct value, with counts
    # (paths are long, so several of them go one per line)
    for key in ('EXOSIMS_version', 'EXOSIMS_path'):
        counts = defaultdict(int)
        for env in envs:
            counts[env.get(key, '(missing)')] += 1
        if len(counts) == 1:
            values = list(counts)[0]
        else:
            values = ['%s (%d)' % (v, n) for v, n in sorted(counts.items())]
            if key != 'EXOSIMS_path':
                values = ', '.join(values)
        lines.append((key, values))
    return lines


def info_print(node, envs_of, depth=0):
    r'''Print the run-info tree at node.  envs_of maps DRM filename to its environment.'''
    prefix = '# ' if CSV_OUTPUT else ''
    fns = node.all_drms()
    envs = [envs_of[fn] for fn in fns if envs_of[fn] is not None]
    indent = prefix + '|  ' * (depth + 1)
    print('%s%s%s' % (prefix, tree_prefix(depth), node.name))
    # Nscen = number of ensembles (drm/ directories) below a family/experiment/sims
    if node.kind != 'ensemble':
        print('%sNscen = %d' % (indent, node.n_ensembles()))
    # Nens = number of DRMs, noting any that lack logs
    logs_note = '' if len(envs) == len(fns) else '  (%d with logs)' % len(envs)
    print('%sNens = %d%s' % (indent, len(fns), logs_note))
    for key, value in info_lines(envs):
        if isinstance(value, list):
            print('%s%s:' % (indent, key))
            for v in value:
                print('%s    %s' % (indent, v))
        else:
            print('%s%s: %s' % (indent, key, value))
    for child in node.children:
        info_print(child, envs_of, depth + 1)


############################################################
#
# Main routine
#
############################################################

def main(args):
    r'''Main routine: load and process DRMs.'''

    global CSV_OUTPUT
    CSV_OUTPUT = args.csv_output
    # this DIAGNOSE is intended to be for error-reporting to diagnose
    # issues with the code or the files, but not for ordinary use
    global DIAGNOSE
    DIAGNOSE = args.DIAGNOSE

    # Find the DRMs, as trees following the Sandbox hierarchy
    trees = args_to_trees(args.drm)
    if not trees:
        return

    # Run info from the run-environment logs, if requested
    if args.info:
        envs_of = {fn: load_environ(fn)
                   for node in trees for fn in node.all_drms()}
        for node in trees:
            info_print(node, envs_of)
        # -i alone (no -l or -s) is info only
        if args.verbose == 1 and not args.summary:
            return

    # Summarize each DRM, and roll up
    n_bad = sum(load_tree_stats(node, args.verbose) for node in trees)
    table_print(trees, args.verbose, args.summary)
    if n_bad:
        print('%s: Warning: %d file(s) could not be loaded as DRMs, and are omitted'
              % (os.path.basename(sys.argv[0]), n_bad), file=sys.stderr)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Summarize EXOSIMS DRM(s).",
                                     epilog='')
    parser.add_argument('drm', metavar='DIR_OR_DRM', nargs='+',
                            help='Sandbox directory (ensemble, family, experiment), or DRM file')
    parser.add_argument('-l', '--long', help='long-format listing',
                      dest='verbose', action='count', default=1)
    parser.add_argument('-D', '--diagnose', help='diagnostic output: traceback for DRMs that fail to load (one per group)',
                      dest='DIAGNOSE', action='count', default=0)
    parser.add_argument('-s', '--summary', help='rollups only (no per-DRM output)',
                      dest='summary', action='store_true', default=False)
    parser.add_argument('-c', '--csv', help='CSV output', default=False,
                      dest='csv_output', action='store_true')
    parser.add_argument('-i', '--info', help='tree of run info, from environment logs', default=False,
                      dest='info', action='store_true')
    args = parser.parse_args()

    main(args)
    sys.exit(0)
