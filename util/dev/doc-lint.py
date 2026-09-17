#!/usr/bin/env python3
r'''
doc-lint.py: check that documentation blocks will render as markdown

The Sandbox utility scripts carry their user documentation in a module
docstring (.py) or a leading #-comment block (.sh). `mkdocs` renders those
blocks as Markdown, by way of util/doc_sandbox/gen_pages.py, so a block that
reads fine in a terminal can still render badly on the documentation page.

This checks for the ways that goes wrong:

+ `fence`: an unbalanced ``` fence, or a closing fence carrying an info
  string (```shell), which is not a close at all -- the rest of the block
  is swallowed into the code block.
+ `indent`: a run of option descriptions ("-o FILE gives ...") indented
  1-3 spaces and not inside a ``` fence. Markdown needs 4 spaces, or a
  fence, to keep the line breaks; at 1-3 spaces the whole option list
  reflows into one run-on paragraph.
+ `nodoc`: no documentation block at all.
+ `usage`: a usage line carrying bracketed flags (`prog [-v VERB] ARG`)
  that is neither inside a ``` fence nor inline backticks. mkdocs reads
  `[-v VERB]` as a cross-reference, fails to resolve it, and warns on
  every build.
+ `nosummary`: the generated script index would show an empty description
  for this script. gen_pages.py builds that description by finding the
  script's own name inside the block ("drm-ls.py: summarize DRMs ..."), so
  a block that never names the script, or names it differently, yields a
  bare entry on the Script Index page.
+ `filename`: the block's leading "name.py:" disagrees with the real
  filename, usually left behind by a rename.
+ `topics`: the two plot-topic lists disagree -- TOPICS in
  Local/www-doc/Makefile, the *.md files there, and the "Plot Reference"
  nav in util/doc_sandbox/mkdocs.yml all have to name the same pages.

Synopsis:
```
   util/dev/doc-lint.py [-q] [-c CHECK] [FILE ...]
```

With no FILE, checks the standard set: util/*.py, util/*.sh,
util/plot_drm_gallery/*.py, and the Sandbox-root driver scripts.

The `nosummary` and `filename` checks apply only to the scripts that
gen_pages.py actually puts in the Script Index (util/*.py, util/*.sh and
the root driver scripts); the gallery modules are not indexed, so there is
no index entry for them to get wrong.

+ `-q` reports only a count, no detail
+ `-c CHECK` limits to one check, by the names above

Exits nonzero if anything is reported, so it can gate a build.

Typical usage:
  $ make doc-lint
  $ util/dev/doc-lint.py -c fence
  $ util/dev/doc-lint.py util/drm-ls.py

turmon 2026
'''

import sys
import os
import ast
import re
import argparse
from pathlib import Path

# repo root, found relative to this script (util/dev/doc-lint.py)
ROOT = Path(__file__).resolve().parent.parent.parent
# driver scripts at the Sandbox root that are documented like util/ scripts
SANDBOX_SCRIPTS = ['add-sims.sh', 'exp-add-sims.sh']
# a line that looks like an option/flag description, e.g. "-o FILE gives ..."
OPTION_LINE = re.compile(r'^\s{1,3}\[?-{1,2}\w')
# a bracketed flag, e.g. "[-v VERB]" -- mkdocs reads this as a cross-reference
BRACKET_FLAG = re.compile(r'\[\s*-{1,2}\w')
# the plot doc pages live here; see the "topics" check
WWW_DOC = 'Local/www-doc'
MKDOCS_YML = 'util/doc_sandbox/mkdocs.yml'

CHECKS = ['fence', 'indent', 'usage', 'nodoc', 'nosummary', 'filename', 'topics']


def doc_block(path):
    r'''Return (block, offset) for a .py or .sh file, or (None, 0).

    Mirrors what util/doc_sandbox/gen_pages.py extracts, so that what is
    checked here is what actually gets rendered. `offset` converts a
    line number within the block to a line number within the file.
    '''
    text = path.read_text(errors='replace')
    if path.suffix == '.py':
        try:
            tree = ast.parse(text)
            doc = ast.get_docstring(tree, clean=False)
        except SyntaxError:
            return None, 0
        if doc is not None:
            # line 1 of the block is the line the opening quote is on
            return doc, tree.body[0].lineno - 1
    # .sh files, and .py files documented with a leading #-block
    lines, offset = [], 0
    for n, l in enumerate(text.splitlines(), 1):
        if l.startswith('#!'):
            continue
        # "##" marks end-of-block, as does any non-# line
        if l.startswith('##') or not l.startswith('#'):
            break
        if not lines:
            offset = n - 1
        lines.append(l[2:] if l.startswith('# ') else l[1:])
    return ('\n'.join(lines), offset) if lines else (None, 0)


def check_fence(block, offset=0):
    r'''Unbalanced fences, and closing fences that carry an info string.'''
    # line numbers here are block-relative; main() shifts them to file lines
    out, open_info, open_line = [], None, 0
    for n, l in enumerate(block.splitlines(), offset + 1):
        m = re.match(r'\s*(```+)\s*(\S*)', l)
        if not m:
            continue
        info = m.group(2)
        if open_info is None:
            open_info, open_line = info, n
        elif info:
            # a fence with an info string opens a block; it never closes one
            out.append((n, f'closing fence has an info string "```{info}" '
                           f'(opened line {open_line}); rest of block is swallowed'))
            open_info, open_line = info, n
        else:
            open_info = None
    if open_info is not None:
        out.append((open_line, 'fence opened here is never closed'))
    return out


def check_indent(block, offset=0):
    r'''Runs of shallowly-indented option lines, which markdown reflows.

    Keyed on the option lines themselves rather than on a "where:"-style
    heading, because plenty of these blocks are introduced by ordinary
    prose ("...and optionally:") that no heading pattern would catch.
    '''
    out, in_fence, run = [], False, []

    def flush():
        # a run is only a problem if it actually describes options, and if
        # the options themselves sit at 1-3 spaces (deeper continuation
        # lines are part of the same list and do not break it)
        opts = [(n, l) for n, l in run if OPTION_LINE.match(l)]
        if len(run) >= 2 and opts:
            indent = min(len(l) - len(l.lstrip()) for _, l in opts)
            if 1 <= indent <= 3:
                out.append((opts[0][0],
                            f'option list is indented {indent} spaces '
                            f'({len(run)} lines); markdown needs 4, '
                            f'or a ``` fence'))
        run.clear()

    for n, l in enumerate(block.splitlines(), offset + 1):
        if re.match(r'\s*```', l):
            in_fence = not in_fence
            flush()
            continue
        if in_fence:
            continue
        if not l.strip():
            # a blank line does not break a run: option lists are often spaced
            continue
        if l.startswith((' ', '\t')):
            run.append((n, l))
        else:
            flush()
    flush()
    return out


def check_usage(block, offset=0):
    r'''Bracketed flags outside a fence, which mkdocs treats as autorefs.'''
    out, in_fence = [], False
    for n, l in enumerate(block.splitlines(), offset + 1):
        if re.match(r'\s*```', l):
            in_fence = not in_fence
            continue
        if in_fence:
            continue
        # inline code spans are safe: mkdocs does not resolve refs inside them
        if BRACKET_FLAG.search(re.sub(r'`[^`]*`', '', l)):
            out.append((n, 'usage line has bracketed flags outside a ``` fence; '
                           'mkdocs will read them as cross-references and warn'))
    return out


def check_nosummary(block, path):
    r'''The script index entry would have an empty description.

    This mirrors find_oneliner() in util/doc_sandbox/gen_pages.py: the
    description is whatever follows the script's own name in the block.
    A blank first line is fine -- the whole block is searched.
    '''
    stem = path.stem
    for l in block.split('\n'):
        if stem not in l:
            continue
        line1 = l[l.find(stem) + len(stem):].strip(':\n')
        for ext in ('.sh', '.py'):
            if line1.startswith(ext):
                line1 = line1[3:]
        if line1.startswith(':'):
            line1 = line1[1:].strip()
        if line1:
            return []
        break
    return [(1, f'script index entry will be empty: no "{stem}: <summary>" '
                f'line in the documentation block')]


def check_filename(block, path):
    m = re.match(r'\s*([\w.-]+\.(?:py|sh))\s*:', block)
    if m and m.group(1) != path.name:
        return [(1, f'block names "{m.group(1)}" but the file is "{path.name}"')]
    return []


def check_topics():
    r'''The plot-topic lists must agree across three places.'''
    out = []
    mk = ROOT / WWW_DOC / 'Makefile'
    yml = ROOT / MKDOCS_YML
    if not mk.is_file() or not yml.is_file():
        return out
    m = re.search(r'^TOPICS\s*=\s*(.*)$', mk.read_text(), re.MULTILINE)
    if not m:
        return [(0, f'{WWW_DOC}/Makefile: no TOPICS variable found')]
    topics = set(m.group(1).split())
    on_disk = {p.name for p in (ROOT / WWW_DOC).glob('*.md')
               if not p.name.startswith('index')}
    in_nav = {f'{s}.md' for s in
              re.findall(r'plots/([\w.-]+)\.md', yml.read_text())}
    for missing in sorted(on_disk - topics):
        out.append((0, f'{missing} exists in {WWW_DOC} but is not in TOPICS'))
    for missing in sorted(topics - on_disk):
        out.append((0, f'{missing} is in TOPICS but not in {WWW_DOC}'))
    for missing in sorted(topics - in_nav):
        out.append((0, f'{missing} is in TOPICS but not in the mkdocs "Plot Reference" nav'))
    for missing in sorted(in_nav - topics):
        out.append((0, f'{missing} is in the mkdocs nav but not in TOPICS'))
    return out


def indexed_files():
    r'''Scripts that gen_pages.py puts in the Script Index.'''
    files = sorted((ROOT / 'util').glob('*.py'))
    files += sorted((ROOT / 'util').glob('*.sh'))
    files += [ROOT / n for n in SANDBOX_SCRIPTS]
    return [f for f in files if f.is_file() and f.name != '__init__.py']


def default_files():
    files = indexed_files()
    files += sorted((ROOT / 'util' / 'plot_drm_gallery').glob('*.py'))
    return [f for f in files if f.is_file() and f.name != '__init__.py']


def main():
    parser = argparse.ArgumentParser(
        description='Check that documentation blocks render as markdown.')
    parser.add_argument('files', metavar='FILE', nargs='*',
                        help='files to check (default: the standard set)')
    parser.add_argument('-q', action='store_true', dest='quiet',
                        help='report only a count')
    parser.add_argument('-c', metavar='CHECK', dest='check', choices=CHECKS,
                        help='run only this check (%s)' % ', '.join(CHECKS))
    args = parser.parse_args()

    want = (lambda c: args.check in (None, c))
    files = [Path(f) for f in args.files] if args.files else default_files()
    # index-entry checks only make sense for scripts that get an index entry
    indexed = {f.resolve() for f in indexed_files()}

    findings = []
    for path in files:
        rel = os.path.relpath(path, ROOT)
        if rel.startswith('..'):
            # a file from outside the Sandbox: show it as given
            rel = str(path)
        block, offset = doc_block(path)
        if block is None:
            if want('nodoc'):
                findings.append((rel, 0, 'nodoc', 'no documentation block'))
            continue
        for check, fn in (('fence', check_fence), ('indent', check_indent),
                          ('usage', check_usage)):
            if want(check):
                findings += [(rel, n, check, msg) for n, msg in fn(block, offset)]
        if path.resolve() in indexed:
            for check, fn in (('nosummary', check_nosummary),
                              ('filename', check_filename)):
                if want(check):
                    findings += [(rel, n + offset, check, msg)
                                 for n, msg in fn(block, path)]
    if want('topics'):
        findings += [('(plot topics)', n, 'topics', msg) for n, msg in check_topics()]

    if not args.quiet:
        for rel, n, check, msg in findings:
            where = f'{rel}:{n}' if n else rel
            print(f'{where}: {check}: {msg}')
    print(f'doc-lint.py: {len(findings)} finding(s) over {len(files)} file(s).')
    return 1 if findings else 0


if __name__ == '__main__':
    sys.exit(main())
