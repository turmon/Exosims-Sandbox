#!/usr/bin/env -S uv run
# /// script
# requires-python = ">=3.13"
# dependencies = [
#   "markdown>=3.7",
# ]
# ///
'''Generate two-level HTML summaries for emulator workflow directories.

Run with:
    demo/emulator_html_summary.py [options] [DIR ...]

Example:
    demo/emulator_html_summary.py -R . -S inverse.css demo/DemoData/Yokohama_Extended.fam/v2-v3-v23-nf.exp

Each DIR must be a .fam or .exp directory. With -r/--recursive, matching
subdirectories are also processed. Directories without an Analysis/ subdir
are skipped silently.

For each processed DIR, writes index.{md,html} in DIR/Analysis/ and at
each DIR/Analysis/<TAG>.wfl/. Only subdirectories ending in .wfl are
treated as workflow output directories; the semantic tag used in titles
and links is the directory name with the .wfl suffix stripped.
'''
import json
import argparse
import os
import re
import sys
from pathlib import Path

import time

import markdown

PROG = os.path.basename(sys.argv[0])

# ---------------------------------------------------------------------------
# CSS hook
# ---------------------------------------------------------------------------

# Hook: literal HTML injected into <head> of the generated HTML file.
# Edit this string to point to your CSS file.
HEAD_EXTRA = '<link rel="stylesheet" href="emulator.css">'

CONFIG_FILE = 'config-workflow.json'

# ---------------------------------------------------------------------------
# File-pattern → section config
# ---------------------------------------------------------------------------

SECTIONS = [
    {
        'heading': 'Experimental Sampling',
        'intro': 'Scatter plots of training and validation data across all parameters.',
        'img_globs': ['scatter_*.png'],
        'md_globs': [],
    },
    {
        'heading': 'Emulator Fit',
        'intro': 'Emulator fit diagnostics (true vs. predicted, residuals) on training and validation data.',
        'img_globs': ['emulator_fit_*.png', 'grid_search_fit_*.png'],
        'md_globs': ['emulator_fit_metrics_*.md', 'grid_search_fit_metrics_*.md'],
    },
    {
        'heading': 'Emulator Cross-Section',
        'intro': '2D cross-section corner plots of the emulator surface around a reference point.',
        'img_globs': ['section_cornerplot_around_*.png'],
        'md_globs': [],
    },
    {
        'heading': 'Emulator Diagnostics',
        'intro': 'Partial dependence plots showing univariate and pairwise parameter effects.',
        'img_globs': ['pdepend_*.png'],
        'md_globs': [],
    },
    {
        'heading': 'Emulator Grid Search',
        'intro': 'Emulator parametric grid search on cross-validated training data. (Top-ranked model is shown above.)\n\nClick table column heads to sort.',
        'img_globs': [],
        'md_globs': ['grid_search_all_metrics_*.md'],
    },
    {
        'heading': 'Inverse Emulator Fit',
        'intro': 'Emulator fit diagnostics on training and validation data (for inverse).',
        'img_globs': ['emulator_inv_fit_*.png'],
        'md_globs': ['emulator_inv_metrics_*.md'],
    },
    {
        'heading': 'Inverse Posterior',
        'intro': 'Parameter posterior distributions, one per QOI target value.',
        'img_globs': ['inverse_posterior_cornerplot_*_QOI*.pdf'],
        'md_globs': ['inverse_posterior_summary_*.md'],
    },
    {
        'heading': 'Inverse Posterior: Combined',
        'intro': 'Combined posterior distributions summarizing all QOI target values.',
        'img_globs': ['inverse_posterior_*_combo.pdf'],
        'md_globs': [],
    },
]

# ---------------------------------------------------------------------------
# Per-file description
# ---------------------------------------------------------------------------

# Module-level table: (compiled regex, format string). Groups are positional.
_FILE_DESC_PATTERNS = [
    (r'^scatter_univariate_train_(.+)$',
     'Univariate scatter plots of training data for {0}.'),
    (r'^scatter_cornerplot_train_(.+)$',
     'Corner scatter plot of training data for {0}.'),
    (r'^scatter_univariate_valid_(.+)$',
     'Univariate scatter plots of validation data for {0}.'),
    (r'^scatter_cornerplot_valid_(.+)$',
     'Corner scatter plot of validation data for {0}.'),
    (r'^scatter_cornerplot_combo_(.+)$',
     'Combined corner scatter plot of training and validation data for {0}.'),
    (r'^emulator_fit_train_(.+)$',
     'Emulator fit diagnostics on training data for {0}.'),
    (r'^emulator_fit_valid_(.+)$',
     'Emulator fit diagnostics on validation data for {0}.'),
    # grid search
    (r'^grid_search_fit_train_(.+)$',
     'Emulator search, best-fit diagnostics on training data for {0}.'),
    (r'^grid_search_fit_valid_(.+)$',
     'Emulator search, best-fit diagnostics on validation data for {0}.'),
    #
    (r'^section_cornerplot_around_(.+)$',
     '2D cross-section corner plots around a reference point for {0}.'),
    (r'^pdepend_univariate_(.+)$',
     'Univariate partial dependence plots for {0}.'),
    (r'^pdepend_cornerplot_(.+)$',
     'Corner plot of partial dependence for {0}.'),
    (r'^emulator_inv_fit_train_(.+)$',
     'Emulator fit diagnostics on training data for {0} (inverse mode).'),
    (r'^emulator_inv_fit_valid_(.+)$',
     'Emulator fit diagnostics on validation data for {0} (inverse mode).'),
    (r'^inverse_posterior_(.+)_(QOI.+)$',
     'Parameter posterior for QOI target value {1} in {0}.'),
    (r'^inverse_posterior_(.+)_combo$',
     'Combined posterior distributions across all QOI targets in {0}.'),
]


def _file_description(fname):
    '''Return a short description string derived from filename stem.'''
    stem = Path(fname).stem
    for pattern, template in _FILE_DESC_PATTERNS:
        m = re.match(pattern, stem)
        if m:
            return template.format(*m.groups())
    return ''


# ---------------------------------------------------------------------------
# Config loading
# ---------------------------------------------------------------------------

def _load_config(base_dir):
    '''Return (class_name, band, qoi_targets) from config-workflow.json, or None values.'''
    config_path = Path(base_dir) / 'Analysis' / CONFIG_FILE
    if not config_path.exists():
        return None, None, None
    try:
        with config_path.open() as fh:
            cfg = json.load(fh)
        return (
            cfg.get('class_name', ''),
            cfg.get('band', ''),
            cfg.get('qoi_targets', []),
        )
    except Exception:
        return None, None, None


# ---------------------------------------------------------------------------
# Markdown generation helpers
# ---------------------------------------------------------------------------

def slugify(text):
    return re.sub(r'[^a-z0-9]+', '-', text.lower()).strip('-')


def detect_present_sections(gfx_dir):
    '''Return list of heading strings for sections that have matching files in gfx_dir.'''
    result = []
    for section in SECTIONS:
        found = any(any(gfx_dir.glob(g)) for g in section['img_globs'])
        found = found or any(any(gfx_dir.glob(g)) for g in section['md_globs'])
        if found:
            result.append(section['heading'])
    return result


def build_nav_html(h2_items):
    '''Build a fixed <nav> bar HTML string.

    Parameters
    ----------
    h2_items : list of (label, h2_slug, [section_heading, ...])
    '''
    items = []
    for label, h2_slug, sections in h2_items:
        sub = '\n'.join(
            f'      <li><a href="#{h2_slug}-{slugify(s)}">{s}</a></li>'
            for s in sections
        )
        items.append(
            f'<li><a href="#{h2_slug}">{label}</a>\n'
            f'  <ul class="level-2">\n{sub}\n  </ul>\n</li>'
        )
    inner = '\n'.join(f'  {it}' for it in items)
    return f'<nav>\n<ul class="level-1">\n{inner}\n</ul>\n</nav>'


def _md_image_entries(gfx_dir, link_prefix, h2_slug=None):
    '''Return Markdown lines for SECTIONS that match files in gfx_dir.

    Parameters
    ----------
    gfx_dir : Path
        Directory to scan for image files.
    link_prefix : str
        Prepended to each filename to form the relative link, e.g. 'gfx/' or
        'HGB/gfx/'.
    h2_slug : str or None
        When provided, H3 headings get an explicit anchor id of
        ``{h2_slug}-{section_slug}``.
    '''
    lines = []
    for section in SECTIONS:
        matched_imgs = []
        seen = set()
        for glob in section['img_globs']:
            for p in sorted(gfx_dir.glob(glob)):
                if p.name not in seen:
                    seen.add(p.name)
                    matched_imgs.append(p.name)

        matched_mds = []
        for glob in section['md_globs']:
            for p in sorted(gfx_dir.glob(glob)):
                if p not in matched_mds:
                    matched_mds.append(p)

        if not matched_imgs and not matched_mds:
            continue
        heading_id = f' {{#{h2_slug}-{slugify(section["heading"])}}}' if h2_slug else ''
        lines.append(f'### {section["heading"]}{heading_id}')
        lines.append('')
        lines.append(section['intro'])
        lines.append('')

        for md_path in matched_mds:
            lines.append(f'<!-- begin {md_path.name} -->')
            lines.append(md_path.read_text())
            lines.append(f'<!-- end {md_path.name} -->')
            lines.append('')

        for fname in matched_imgs:
            link = link_prefix + fname
            desc = _file_description(fname)
            if desc:
                lines.append(f'**{fname}:** {desc}')
            else:
                lines.append(f'**{fname}:**')
            lines.append('')
            lines.append(f'<a class="plot-thumb" href="{link}"><img src="{link}" alt="{fname}"></a>')
            lines.append('')
    return lines


def identify_jobs(tag_dir):
    '''Return job_dirs for a TAG directory.

    Returns a dict mapping JOB -> job spec's, provided
    appropriate contents are found.
    '''
    job_specs = dict()
    for d in tag_dir.iterdir():
        if not d.is_dir():
            continue
        if not (d / 'gfx').is_dir():
            continue
        # the job's run-spec file
        spec_file = d / 'run-spec.json'
        if not spec_file.exists():
            continue
        with open(spec_file, 'r') as fp:
            job_specs[d] = json.load(fp)
    # return in name-sorted order
    return dict(sorted(job_specs.items(), key=lambda x: x[0].name))


# ---------------------------------------------------------------------------
# Page builders
# ---------------------------------------------------------------------------

def titlify(job_dir, spec):
    '''Printable title from a job directory Path and its spec dict.'''
    if 'title' in spec:
        return f"{spec['title']} ({job_dir.name})"
    else:
        return f"{job_dir.name}"

def build_job_markdown(tag_dir):
    '''Build and return Markdown text for one TAG directory.

    tag_dir/gfx/ present): single-level sections under H2.
    '''
    dname = tag_dir.name
    tag = dname[:-4] if dname.endswith('.wfl') else dname

    job_specs = identify_jobs(tag_dir)

    h2_items = [(
        titlify(job_path, spec),
        slugify(job_path.name),
        detect_present_sections(job_path / 'gfx')
        )
        for job_path, spec in job_specs.items()
        ]

    lines = [f'# Workflow: {tag}', '']
    lines.append(build_nav_html(h2_items))
    lines.append('')
    lines.append('Upward to all <a href="../index.html">Workflows</a>')
    lines.append('')
    lines.append(f'Workflow consists of {len(job_specs)} ' +
                 f'job{"s" if len(job_specs) > 1 else ""}' +
                 f' [(workflow spec)](../workflow-{tag}.json)')

    # Table of contents at top as basic list
    lines.append('')
    for job_path, job_spec in job_specs.items():
        lines.append(f'  - {titlify(job_path, job_spec)} ([link]({tag}/index.html#{slugify(job_path.name)}))')
        lines.append('')


    for job_path, spec in job_specs.items():
        job = job_path.name
        job_slug = slugify(job)
        # needs a helper function to handle edge cases
        job_title = titlify(job_path, spec)
        lines.append(f'## Job: {job_title} {{#{job_slug}}}')
        lines.append(f'Full job specification: [run-spec.json]({job}/run-spec.json)')
        lines.append('')
        lines.extend(_md_image_entries(job_path / 'gfx', f'{job}/gfx/', h2_slug=job_slug))

    return '\n'.join(lines), job_specs


def build_workflow_markdown(base_dir, tags, workflow_job_specs):
    '''Build and return Markdown text for the workflow index page.

    Writes to <base_dir>/Analysis/index.md. Links point upward to
    config-workflow.json and workflow-<TAG>.json files in base_dir.
    '''
    lines = [f'# Workflow: {base_dir.name}', '']
    lines.append('Upward to <a href="../index.html">Experiment/Family</a>')
    lines.append('')

    config_path = base_dir / 'Analysis' / CONFIG_FILE
    if config_path.exists():
        lines.append(f'Workflow Template: [{CONFIG_FILE}]({CONFIG_FILE})')
        lines.append('')

    for tag in tags:
        lines.append(f'## Workflow: {tag}')
        lines.append('')
        lines.append(f'Down to [all workflow results]({tag}.wfl/index.html), or by job:')
        lines.append('')
        for job_path, job_spec in workflow_job_specs[tag].items():
            lines.append(f'  - {titlify(job_path, job_spec)} ([direct link]({tag}.wfl/index.html#{slugify(job_path.name)}))')
            lines.append('')
        workflow_spec = base_dir / 'Analysis' / f'workflow-{tag}.json'
        if workflow_spec.exists():
            lines.append(f'Workflow Specification: [workflow-{tag}.json](workflow-{tag}.json)')
            lines.append('')

    return '\n'.join(lines)


# ---------------------------------------------------------------------------
# HTML generation
# ---------------------------------------------------------------------------

def build_html(md_text, title, head_extra=None):
    '''Convert Markdown text to a full HTML document string.'''
    body = markdown.markdown(md_text, extensions=['tables', 'attr_list'])
    if head_extra is None:
        head_extra = HEAD_EXTRA
    head_extra = f'\n  {head_extra}' if head_extra else ''
    # attr_list to set class is unsupported for python-markdown tables
    # See: "attr-list" docs -> Limitations -> Implied Elements
    # "There is no way to use an attribute list to define attributes
    # on implied elements, including [...] table [...]"
    # Surrounding in <div> does not work either. 2026/03.
    new_body = body.replace('<table>', '<table class="sortable">')

    return (
        f'<!DOCTYPE html>\n'
        f'<html>\n'
        f'<head>\n'
        f'  <meta charset="utf-8">\n'
        f'  <title>{title}</title>{head_extra}\n'
        f'</head>\n'
        f'<body>\n'
        f'{new_body}\n'
        f'</body>\n'
        f'</html>\n'
    )


# ---------------------------------------------------------------------------
# Output helper
# ---------------------------------------------------------------------------

def save_md_and_html(md_path, md_text, stylesheets=(), javascripts=()):
    '''Write md_path and the corresponding .html file.

    Parameters
    ----------
    stylesheets : sequence of Path
        Absolute paths to CSS files. When provided, overrides HEAD_EXTRA;
        links are computed relative to the HTML file's directory.
    javascripts : sequence of Path
        Absolute paths to JS files. When provided, overrides HEAD_EXTRA;
        links are computed relative to the HTML file's directory.
    '''
    md_path.write_text(md_text)
    # print(f'{PROG}: Wrote: {md_path}')

    title_line = next(
        (l for l in md_text.splitlines() if l.startswith('# ')), '# Report'
    )
    title = title_line.lstrip('# ').strip()


    html_dir = md_path.parent
    links = [
        f'  <link rel="stylesheet" href="{os.path.relpath(css, html_dir)}">'
            for css in stylesheets
        ] + [
        f'  <script src="{os.path.relpath(js, html_dir)}"></script>'
            for js in javascripts
        ]
    if links:
        head_extra = '\n'.join(links)
    else:
        head_extra = None  # falls back to HEAD_EXTRA in build_html

    html_path = md_path.with_suffix('.html')
    html_path.write_text(build_html(md_text, title, head_extra=head_extra))
    print(f'{PROG}: Wrote: {html_path}')
    return 1


# ---------------------------------------------------------------------------
# Multi-dir helpers
# ---------------------------------------------------------------------------

def iter_targets(base_dir, recursive):
    """Yield .fam/.exp dirs rooted at base_dir to index.

    base_dir itself is yielded iff its name ends with '.fam' or '.exp'.
    With recursive=True, matching subdirectories are also yielded (depth-first).
    Non-matching directories are silently skipped.
    """
    if not (base_dir.name.endswith('.fam') or base_dir.name.endswith('.exp')):
        return
    yield base_dir
    if recursive:
        for child in sorted(base_dir.iterdir()):
            if child.is_dir():
                yield from iter_targets(child, recursive)


def process_dir(base_dir, stylesheets, javascripts):
    """Generate index files for base_dir if it has an Analysis/ subdir.

    Returns the number of HTML files written.
    """
    analysis_dir = base_dir / 'Analysis'
    if not analysis_dir.is_dir():
        print(f'{PROG}: Skipping {base_dir}')
        return 0
    print(f'{PROG}: Indexing {base_dir}')
    n = 0
    workflow_job_specs = dict()
    tags = sorted(d.name[:-4] for d in analysis_dir.iterdir()
                  if d.is_dir() and d.name.endswith('.wfl'))
    for tag in tags:
        wfl_dir = analysis_dir / (tag + '.wfl')
        md_text, job_specs = build_job_markdown(wfl_dir)
        workflow_job_specs[tag] = job_specs
        n += save_md_and_html(wfl_dir / 'index.md', md_text, stylesheets, javascripts)
    md_text = build_workflow_markdown(base_dir, tags, workflow_job_specs)
    n += save_md_and_html(analysis_dir / 'index.md', md_text, stylesheets, javascripts)
    return n


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Generate multi-page HTML summary for available emulator workflows.')
    parser.add_argument('dir', metavar='DIR', nargs='*',
                        help='Workflow root directory (0 or more)')
    parser.add_argument('-r', '--recursive', action='store_true',
                        help='Descend recursively into .fam/.exp subdirectories')
    parser.add_argument('-R', '--resources', default=None, metavar='WWW_DIR',
                        help='Directory containing CSS stylesheet files')
    parser.add_argument('-S', '--css', action='append', default=[], metavar='STYLE',
                        help='CSS stylesheet filename relative to --resources; may be repeated')
    parser.add_argument('-J', '--javascript', action='append', default=[], metavar='JS',
                        help='JS script filename relative to --resources; may be repeated')
    args = parser.parse_args()

    stylesheets = []
    if args.css:
        resource_dir = Path(args.resources) if args.resources else Path('.')
        for style in args.css:
            css = (resource_dir / style).resolve()
            if not css.exists():
                print(f'Warning: stylesheet not found: {css}', file=sys.stderr)
            stylesheets.append(css)
    javascripts = []
    if args.javascript:
        resource_dir = Path(args.resources) if args.resources else Path('.')
        for script in args.javascript:
            js = (resource_dir / script).resolve()
            if not js.exists():
                print(f'Warning: script not found: {js}', file=sys.stderr)
            javascripts.append(js)

    t0 = time.monotonic()
    n = 0
    for d in args.dir:
        for target in iter_targets(Path(d), args.recursive):
            n += process_dir(target, stylesheets, javascripts)
    elapsed = time.monotonic() - t0
    print(f'{PROG}: Wrote {n} indexes in {elapsed:.1f} sec')
    print(f'{PROG}: Done.')


if __name__ == '__main__':
    main()
