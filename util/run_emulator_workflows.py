#!/usr/bin/env -S uv run
# /// script
# requires-python = ">=3.13"
# dependencies = [
#   "ExoEmulator[inverse]",
# ]
#
# [tool.uv.sources]
# ExoEmulator = { path = "../../Tools/ExoEmulator", editable = true }
# ///
'''Train, validate, and plot yield emulators and inverses from the CLI.

Actions:
    design    Scatter/cornerplots of training and validation data only.
    emulator  Full emulator workflow: scatter plots + emulator fit + PDP plots.
    inverse   Emulator and inverse-yield parameter posterior plots.
    all       Design, emulator, and inverse (everything).

Run with:
    # Inline single-job form:
    demo/run_emulator_workflows.py -a ACTION [-a ACTION ...] \\
           -d DATA_DIR [-d DATA_DIR ...] [-v VALID_DIR] \\
           [-b BAND[,...]] [-o OUTPUT_DIR] [-C ROOT_DIR]

    # Workflow file form (one or more files):
    demo/run_emulator_workflows.py WORKFLOW.json [WORKFLOW2.json ...] \\
           [-b BAND[,...]] [-C ROOT_DIR]

Reads config-workflow.json from the current working directory to determine the
emulator class, band, training method, constructor arguments, yield targets,
and MCMC parameters.  In file mode, each job dict can override any key from
config-workflow.json; config-workflow.json is optional if all required keys are
supplied in the job dict.
'''
# Design note: workflows, jobs, actions
#
# A *workflow file* is a JSON file containing a list of *job dicts*.  Each job
# dict describes one complete run: which data directories to use, which spectral
# band, which output directory, and which actions to perform.  Multiple workflow
# files can be given on the command line; they are executed in order.
#
# A *job* is a single dict passed to run_workflows().  It is assembled from two
# sources: (1) the job dict loaded from the workflow file (or built inline from
# CLI flags), and (2) a per-invocation job_template derived from CLI flags
# (-d, -v, -b, -C, -a).  CLI flags in the template supersede matching keys in
# the file-supplied job dict, so they act as global overrides.
#
# An *action* (design | emulator | inverse | all) selects which processing
# steps run_workflows() will execute for a given job.  In file mode the actions
# come from the job dict; in inline mode they come from -a flags (defaulting to
# 'all' when -a is omitted).
import matplotlib
matplotlib.use('Agg')

import os
import sys
import json
import time
import datetime
import importlib
import argparse
import contextlib
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# ---------------------------------------------------------------------------
# Class registry: maps config class_name -> (module, DataClass, EmulatorClass)
# ---------------------------------------------------------------------------

PROG = os.path.basename(sys.argv[0])

CONFIG_FILE = 'config-workflow.json'

CLASS_REGISTRY = {
    'YokohamaUltimateV2YieldEmulator': (
        'ExoEmulator.YokohamaUltimateV2',
        'YokohamaUltimateV2EmulatorData',
        'YokohamaUltimateV2YieldEmulator'),
    'YokohamaExtendedV2YieldEmulator': (
        'ExoEmulator.YokohamaExtendedV2',
        'YokohamaExtendedV2EmulatorData',
        'YokohamaExtendedV2YieldEmulator'),
}


# ---------------------------------------------------------------------------
# Working directory management
# ---------------------------------------------------------------------------

# temporarily change working directory
@contextlib.contextmanager
def push_dir(d):
    old_dir = os.getcwd()
    if d is not None:
        os.chdir(d)
    try:
        yield
    finally:
        os.chdir(old_dir)

# ---------------------------------------------------------------------------
# Capture logs of called functions
# ---------------------------------------------------------------------------

class Prefixer:
    """A file-like object that prefixes all written text."""
    def __init__(self, prefix, original_stdout):
        self.prefix = prefix
        self.original_stdout = original_stdout

    def write(self, text):
        """Writes the prefixed text to the original stdout."""
        # the \n that follows print() is a separate write()
        if text and text != '\n':
            self.original_stdout.write(self.prefix + text)
        else:
            self.original_stdout.write(text)
        self.flush()

    def flush(self):
        """Flushes the original stdout."""
        self.original_stdout.flush()
    
    def __getattr__(self, attr):
        """Delegates all other attribute access to the original stdout."""
        return getattr(self.original_stdout, attr)

@contextlib.contextmanager
def prefix_stdout(prefix):
    """
    A context manager to temporarily prefix strings written to sys.stdout.
    """
    # Save the current stdout
    old_stdout = sys.stdout
    try:
        # Replace sys.stdout with the custom Prefixer object
        sys.stdout = Prefixer(prefix, old_stdout)
        yield
    finally:
        # Ensure the original sys.stdout is restored when exiting the 'with' block
        sys.stdout = old_stdout
        

# ---------------------------------------------------------------------------
# Save-figure helper
# ---------------------------------------------------------------------------

_written_files = []

def save_fig(fig, path):
    '''Save fig to path, close it, and record the path.'''
    fig.savefig(path, bbox_inches='tight', dpi=150)
    plt.close(fig)
    _written_files.append(path)

def save_table(df, path):
    # Set index=False to exclude the pandas index column
    markdown_table = df.to_markdown(index=False)
    with open(path, 'w') as fp:
        print(markdown_table, file=fp)
    _written_files.append(path)


# ---------------------------------------------------------------------------
# Workflow functions
# ---------------------------------------------------------------------------

YIELD_BASENAME = Path('reduce-yield-plus.csv')


def run_design(bands, data_files, valid_files, out_dir, DataClass):
    '''Scatter/cornerplots of training (and optionally validation) data.'''
    if not data_files and not valid_files:
        print(f'{PROG}: No data: skipping design workflow.')
        return

    for band in bands:
        print(f'{PROG}: === Mode: {band} [design] ===')

        with prefix_stdout('||  '):
            # we have at least one of data or validation
            e_data = DataClass(band)

            if data_files:
                e_data.load_training_files(data_files)
                fig, _ = e_data.univariate_training()
                save_fig(fig, out_dir / f'scatter_univariate_train_{band}.png')
                fig, _ = e_data.cornerplot_training()
                save_fig(fig, out_dir / f'scatter_cornerplot_train_{band}.png')

            if valid_files:
                e_data.load_validation_files(valid_files)
                fig, _ = e_data.univariate_validation()
                save_fig(fig, out_dir / f'scatter_univariate_valid_{band}.png')
                fig, _ = e_data.cornerplot_validation()
                save_fig(fig, out_dir / f'scatter_cornerplot_valid_{band}.png')

            # if both: combo plot
            if data_files and valid_files:
                fig, _ = e_data.cornerplot_combined()
                save_fig(fig, out_dir / f'scatter_cornerplot_combo_{band}.png')


def run_emulation(bands, data_files, valid_files, out_dir,
                  DataClass, EmulatorClass, class_args,
                  emulator_method, emulator_method_args,
                  emulator_workflow_options=None):
    '''Scatter plots + emulator fit + section/PDP plots.'''
    if not data_files:
        print(f'{PROG}: No training data: skipping emulator workflow.')
        return

    opts = emulator_workflow_options or {}
    for band in bands:
        print(f'{PROG}: === Mode: {band} [emulator {emulator_method if emulator_method else "(default)"}] ===')

        with prefix_stdout('||  '):
            e_data = DataClass(band)
            e_data.load_training_files(data_files)
            if valid_files:
                e_data.load_validation_files(valid_files)

            ye = EmulatorClass(band, **class_args)
            ye.attach_data(e_data)
            ye.train(method=emulator_method, method_args=emulator_method_args)

            # The class does not maintain the method used for training,
            # just the (fixed) default, sigh
            if emulator_method:
                ye.label_tag = f'{band} ({emulator_method})'
            else:
                ye.label_tag = f'{band} ({ye.train_method})'

            fig = ye.training_plots()
            save_fig(fig, out_dir / f'emulator_fit_train_{band}.png')

            if valid_files:
                fig = ye.validation_plots()
                save_fig(fig, out_dir / f'emulator_fit_valid_{band}.png')

            df = ye.qoi_fit_metrics({'Method':emulator_method, 'Band':band})
            print(df)
            save_table(df, out_dir / f'emulator_fit_metrics_{band}.md')

            p0 = np.array(
                tuple((ye.l_bound[cn] + ye.u_bound[cn]) / 2.0 for cn in ye.col_names),
                dtype=[(cn, 'float64') for cn in ye.col_names])

            fig, *_ = ye.plot_2d_around(p0)
            save_fig(fig, out_dir / f'section_cornerplot_around_{band}.png')

            fig, *_ = ye.partial_dependence_univariate_plot()
            save_fig(fig, out_dir / f'pdepend_univariate_{band}.png')

            if not opts.get('skip_pdp', False):
                fig, _ = ye.partial_dependence_plot(common_diag_scale=True)
                save_fig(fig, out_dir / f'pdepend_cornerplot_{band}.png')


def run_grid_search(bands, data_files, valid_files, out_dir,
                    DataClass, EmulatorClass, class_args,
                    emulator_method, param_grid):
    '''Parameter grid search for a given emulator method.'''
    if not data_files:
        print(f'{PROG}: No training data: skipping gridsearch workflow.')
        return

    for band in bands:
        print(f'{PROG}: === Mode: {band} [grid_search {emulator_method if emulator_method else "(default)"}] ===')

        with prefix_stdout('||  '):
            e_data = DataClass(band)
            e_data.load_training_files(data_files)
            if valid_files:
                e_data.load_validation_files(valid_files)

            ye = EmulatorClass(band, **class_args)
            ye.attach_data(e_data)

            ye.grid_search(method=emulator_method, param_grid=param_grid)

            # The class does not maintain the method used for training,
            # just the (fixed) default, sigh
            if emulator_method:
                ye.label_tag = f'{band} ({emulator_method})'
            else:
                ye.label_tag = f'{band} ({ye.train_method})'

            df_gs = ye.grid_search_summary
            save_table(df_gs, out_dir / f'grid_search_all_metrics_{band}.md')

            fig = ye.training_plots()
            save_fig(fig, out_dir / f'grid_search_fit_train_{band}.png')

            if valid_files:
                fig = ye.validation_plots()
                save_fig(fig, out_dir / f'grid_search_fit_valid_{band}.png')

            df = ye.qoi_fit_metrics({'Method':emulator_method, 'Band':band})
            print(df)
            save_table(df, out_dir / f'grid_search_fit_metrics_{band}.md')



def run_inverse(bands, data_files, valid_files, out_dir,
                DataClass, EmulatorClass, class_args,
                emulator_method, emulator_method_args, 
                inverse_args, run_pymc_args, qoi_targets):
    '''Emulator fit (inv variant) + PyMC MCMC plots per QOI target.'''
    if not data_files:
        print(f'{PROG}: No training data: skipping inverse workflow.')
        return
    from ExoEmulator.YieldInverse import YieldInverse

    for band in bands:
        print(f'{PROG}: === Band: {band} [inverse {emulator_method if emulator_method else "(default)"}] ===')

        with prefix_stdout('||  '):
            e_data = DataClass(band)
            e_data.load_training_files(data_files)
            if valid_files:
                e_data.load_validation_files(valid_files)

            ye = EmulatorClass(band, **class_args)
            ye.attach_data(e_data)
            ye.train(method=emulator_method, method_args=emulator_method_args)

            # The class does not maintain the method used for training,
            # just the (fixed) default, sigh
            if emulator_method:
                ye.label_tag = f'{band} ({emulator_method})'
            else:
                ye.label_tag = f'{band} ({ye.train_method})'

            fig = ye.training_plots()
            save_fig(fig, out_dir / f'emulator_inv_fit_train_{band}.png')

            if valid_files:
                fig = ye.validation_plots()
                save_fig(fig, out_dir / f'emulator_inv_fit_valid_{band}.png')

            df = ye.qoi_fit_metrics({'Method':emulator_method, 'Band':band})
            print(df)
            save_table(df, out_dir / f'emulator_inv_metrics_{band}.md')

            yi = YieldInverse(ye, out_dir=out_dir, **inverse_args)
            run_label = yi.run_label
            file_label = yi.file_label

            grand_summary = []
            for yval in qoi_targets:
                rv = yi.run_pymc(yval, **run_pymc_args)
                save_fig(rv['fig'], out_dir / f'inverse_posterior_cornerplot_{file_label}_QOI{yval:g}.pdf')
                summ1 = rv['summary']
                #summ1.insert(0, 'col_names', ye.col_names)
                summ1.insert(0, 'QOI_target', yval)
                summ1.insert(0, 'RunLabel', run_label)
                grand_summary.append(summ1)
            # MCMC diagnostics to markdown table
            grand_summary_df = pd.concat(grand_summary, ignore_index=True)
            save_table(grand_summary_df,
                       out_dir / f'inverse_posterior_summary_{file_label}.md')

            fig = yi.plot_combined_posteriors(qoi_targets)
            if fig is not None:
                save_fig(fig, out_dir / f'inverse_posterior_cornerplot_{file_label}_combo.pdf')


# ---------------------------------------------------------------------------
# Single-job workflow runner
# ---------------------------------------------------------------------------

# Keys that control paths/dirs rather than emulator config
#_PATH_KEYS = {'data', 'valid', 'cd', 'output'}
_PATH_KEYS = {'cd'}


def run_workflow_job(pool):
    '''Execute one job in one workflow using a fully-resolved configuration pool.

    Parameters
    ----------
    pool : dict
        Fully merged dict (config-workflow.json base + workflow job overrides + CLI
        flags). Path keys (cd) are excluded. Contains:
        Control keys: actions (list[str]), band (str), override (bool),
                      run_tag (str), job_parent (str), job_number (int)
        Data keys: data (list[str], required), valid (list[str], optional),
                   output (str, optional; default 'plots')
        Emulator keys: emulator_workflow_options (dict, optional) — per-job
                       emulator controls; supported keys: skip_pdp (bool).
        Inverse keys: qoi_targets (list[float], required for inverse/all).
        Config keys: any key from config-workflow.json.
    '''

    class_name = pool.get('class_name', '')
    class_args = pool.get('class_args', {})
    emulator_method = pool.get('emulator_method', '')
    emulator_method_args = pool.get('emulator_method_args', {})
    param_grid = pool.get('param_grid', {})

    # --- Resolve bands ---
    config_band = pool.get('band', 'VIS')
    bands = [b.strip() for b in config_band.split(',') if b.strip()]

    # --- Resolve emulator/data class ---
    if class_name not in CLASS_REGISTRY:
        knowns = ', '.join(CLASS_REGISTRY.keys())
        print(f'{PROG}: Notice: unknown class_name "{class_name}". Valid names: {knowns}. Skipping job.')
        return 0
    module_name, ed_attr, ye_attr = CLASS_REGISTRY[class_name]
    module = importlib.import_module(module_name)
    DataClass = getattr(module, ed_attr)
    EmulatorClass = getattr(module, ye_attr)

    # --- Prepare paths, which can be empty ---
    data_dirs = pool.get('data', [])
    data_files = [Path(d) / YIELD_BASENAME for d in data_dirs]
    valid_dirs = pool.get('valid', [])
    valid_files = [Path(d) / YIELD_BASENAME for d in valid_dirs]

    # --- Output directory ---
    if pool['override']:
        # do not enforce any directory-naming conventions
        out_dir = Path(pool['output'] or 'gfx')
        sentinel_dir = None # "make"-driven updates irrelevant
    else:
        # enforce naming convention within sims/SCENARIO/Analysis
        # output to: TAG/OUTDIR/gfx/*.png
        out_dir = Path(pool['workflow_tag'] + '.wfl', pool['output'], 'gfx')
        sentinel_dir = Path(pool['workflow_tag'] + '.wfl')
    out_dir.mkdir(parents=True, exist_ok=True)

    actions = pool.get('actions', [])

    # --- Dispatch ---
    if any(v in ('design', 'all') for v in actions):
        run_design(bands, data_files, valid_files, out_dir, DataClass)

    if any(v in ('emulator', 'all') for v in actions):
        emulator_workflow_options = pool.get('emulator_workflow_options', {})
        run_emulation(bands, data_files, valid_files, out_dir,
                      DataClass, EmulatorClass, class_args,
                      emulator_method, emulator_method_args,
                      emulator_workflow_options=emulator_workflow_options)

    if any(v in ('grid_search', ) for v in actions):
        emulator_workflow_options = pool.get('emulator_workflow_options', {})
        run_grid_search(bands, data_files, valid_files, out_dir,
                        DataClass, EmulatorClass, class_args,
                        emulator_method, param_grid)
        
    if any(v in ('inverse', 'all') for v in actions):
        # --- Inverse config (only required for inverse/all) ---
        # TODO: normalize this so we set up appropriate defaults
        inverse_args = pool.get('inverse_args', {})
        run_pymc_args = pool.get('run_pymc_args', {})
        if 'qoi_targets' not in pool:
            print(f'{PROG}: Notice: config must contain a "qoi_targets" key for inverse/all. Skipping job.')
            return 0
        qoi_targets = pool['qoi_targets']
        run_inverse(bands, data_files, valid_files, out_dir,
                    DataClass, EmulatorClass, class_args,
                    emulator_method, emulator_method_args,
                    inverse_args, run_pymc_args, qoi_targets)
    # --- Record result ---
    out_spec = out_dir.parent / Path('run-spec.json')
    with open(out_spec, 'w') as spec_file:
        json.dump(pool, spec_file, indent=4)

    # Sentinel file for re-"make" (workflow.wfl/gfx/graphics-info.txt)
    sentinel = out_dir / 'graphics-info.txt'
    with open(sentinel, 'w') as fh:
        fh.write(f'user: {os.environ["USER"]}\n')
        fh.write(f'program: {PROG}\n')
        fh.write(f'workflow: {pool.get("workflow_tag", "unknown")}\n')
        fh.write(f'file: {pool.get("job_parent", "unknown")}\n')
        fh.write(f'date: {datetime.datetime.now().isoformat(timespec="seconds")}\n')
    # all-job sentinel (may be re-written throughout execution)
    # known filesystem location: FOO.{exp,fam}/Analysis/graphics-info.txt
    # so that Makefile knows where it is
    if sentinel_dir is not None:
        sentinel_link = sentinel_dir / 'graphics-info.txt'
        if os.path.lexists(sentinel_link):
            os.unlink(sentinel_link)
        os.symlink(sentinel.relative_to(sentinel_dir), sentinel_link)

    return 1



# ---------------------------------------------------------------------------
# Main helpers
# ---------------------------------------------------------------------------

def load_workflows(args):
    '''Return list of (wf_file, workflow_tag, jobs) from CLI args.

    In file mode, reads each workflow JSON file. In inline mode, builds a
    single synthetic workflow from the CLI flags directly.
    '''
    def extract_workflow_tag(filename):
        '''If possible, extract TAG from <dir/dir/workflow-TAG.json>, otherwise ""'''
        if not filename:
            return ''
        file_base = os.path.basename(filename)
        if file_base.endswith('.json') and file_base.startswith('workflow-'):
            return file_base[9:-5]
        return ''

    if args.workflows:
        workflows_info = []
        for wf_path in args.workflows:
            with open(wf_path) as fh:
                jobs = json.load(fh)
            if not isinstance(jobs, list):
                sys.exit(f'{PROG}: Error: {wf_path} must contain a list of job dicts.')
            workflows_info.append((wf_path, extract_workflow_tag(wf_path), jobs))
    elif args.data or args.valid:
        # --- Inline mode: args ==> a single job, single workflow ---
        actions = args.actions or ['all']
        jobs = [{
            'data':     args.data,
            'valid':    args.valid,
            'band':     args.band,
            'output':   args.output,
            'override': args.override,
            'actions':  actions,
        }]
        workflows_info = [(None, '', jobs)]
    else:
        # no data or workflow: no work
        workflows_info = []
    return workflows_info


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Train, validate, and plot designs, yield emulators and inverse yield.')
    parser.add_argument('workflows', nargs='*', metavar='WORKFLOW.json',
                        help='One or more JSON workflow files (list of job dicts each)')
    parser.add_argument('-a', '--action', dest='actions', action='append',
                        default=[],
                        choices=['design', 'emulator', 'grid_search', 'inverse', 'all'],
                        metavar='ACTION',
                        help='Action to run (design|emulator|inverse|all); '
                             'repeat for multiple')
    parser.add_argument('-d', '--data', action='append', default=[],
                        metavar='DIR', dest='data',
                        help='Directory with reduce-yield-plus.csv; '
                             'repeat to stack multiples')
    parser.add_argument('-v', '--valid', action='append', default=[],
                        metavar='DIR', dest='valid',
                        help='Directory with validation reduce-yield-plus.csv; '
                             'repeat to stack multiples')
    parser.add_argument('-b', '--band', default=None, metavar='BAND[,BAND,...]',
                        help='Comma-separated spectral bands (default: from config, else VIS)')
    parser.add_argument('-o', '--output', default='', metavar='OUTDIR',
                        help='Output directory for plots (default: Analysis/TAG/OUTDIR/gfx)')
    parser.add_argument('-O', '--override', action="store_true", help='Override -o subdirectory conventions')
    parser.add_argument('-C', '--cd', '--chdir', default=None, metavar='DIR',
                        help='Change to DIR before processing (other paths then relative to DIR)')
    args = parser.parse_args()

    # ensure group-write of generated files
    os.umask(0o002)

    # Build job_template from CLI overrides that apply across all jobs/files
    job_template = {}
    if args.cd is not None:
        job_template['cd'] = args.cd
    if args.data:
        job_template['data'] = args.data
    if args.valid:
        job_template['valid'] = args.valid
    if args.actions:
        job_template['actions'] = args.actions
    if args.band is not None:
        job_template['band'] = args.band

    workflows_info = load_workflows(args)

    t0 = time.monotonic()
    n = 0
    # Loop over workflows
    for wf_file, workflow_tag, jobs in workflows_info:
        per_wf_template = {**job_template, 'workflow_tag': workflow_tag}
        if wf_file is not None:
            # must define the override key in File Mode
            per_wf_template['override'] = False
            if args.cd == '':
                # -C '' => working dir set to dirname(workflow file)
                per_wf_template['cd'] = os.path.dirname(wf_file)
        for job_num, job in enumerate(jobs):
            t00 = time.monotonic()
            if 'title' not in job:
                job['title'] = f'Job %d' % (job_num+1) if wf_file else 'Manual'
            print(f'{PROG}: {"="*25}')
            print(f'{PROG}: > Workflow "{workflow_tag}": {len(jobs)} jobs')
            print(f'{PROG}: > Beginning job "{job["title"]}" -- {job_num+1}/{len(jobs)}')
            cd = per_wf_template.get('cd', job.get('cd', None))
            with push_dir(cd):
                config_path = Path(CONFIG_FILE)
                if config_path.exists():
                    with config_path.open() as fh:
                        cfg = json.load(fh)
                else:
                    print(f'{PROG}: No base config loaded.')
                    cfg = {}
                pool = cfg.copy()
                pool.update({k: v for k, v in job.items() if k not in _PATH_KEYS})
                pool.update({k: v for k, v in per_wf_template.items() if k not in _PATH_KEYS})
                pool['job_parent'] = wf_file or 'manual'
                pool['job_number'] = job_num
                # this will be a directory name and a (brief) section head
                if 'output' not in pool:
                    pool['output'] = f'job_{job_num+1}'
                n += run_workflow_job(pool)
            elapsed = time.monotonic() - t00
            print(f'{PROG}: > Finished job {job_num+1}/{len(jobs)}: {elapsed:.1f} sec')

    # --- Summary ---
    if _written_files:
        print(f'{PROG}: Wrote files:')
        for p in _written_files:
            print(f'    {p}')
    else:
        print(f'{PROG}: Wrote no files.')

    elapsed = time.monotonic() - t0
    print(f'{PROG}: Ran {n} workflow job(s) in {elapsed:.1f} sec')
    print(f'{PROG}: Done.')

if __name__ == '__main__':
    main()
