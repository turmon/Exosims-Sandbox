#!/usr/bin/env python
r''' 
pre_mkdocs.py: generate files _before_ mkdocs run

At present, this generates simple template-based
placeholder .md files for python utilities that 
are present in the util_detail/ directory

This is NOT an mkdocs plugin!

It is run before mkdocs to make .md files that mkdocs
consumes when it is eventually run.

There is probably a native mkdocs way to generate these
files, e.g. from an __init.py__ containing them

'''

import glob
from pathlib import Path

# module to index
MODULE = 'util_detail'
# filename of overall script index
INDEX_FILE = 'docs/impl_index.md'
# directory name for generated files
INDEX_DIR = 'implementation'

# directory to look for scripts in
ROOT_DIR = Path(f'{MODULE}')


def gen_markup(script):
    global MODULE

    md_lines = []
    md_lines.append(f'# {script.stem}')
    md_lines.append(f'Implementation details of {script}.')
    md_lines.append(f'::: {MODULE}.{script.stem}')
    return '\n'.join(md_lines)

###
### Main routine
###

def main():
    global INDEX_DIR

    # index file line-by-line (no newlines)
    index = []
    index.append('---')
    index.append(f'title: Implementation Details - Index')
    index.append('---')

    # Python script portion
    index.append('')
    index.append('## Python Scripts')
    index.append('')

    index.append('''Detailed pages giving deeper views
    of a few functions in the `Sandbox` project code are here.''')
    index.append('')
    index.append('')

    # look through "featured scripts" directory only
    for script in sorted(ROOT_DIR.glob('*.py')):
        # a stop-list
        if script.stem in ['__init__']:
            continue
        out = str(Path(f'docs/{INDEX_DIR}') / script.stem) + '.md'
        md = gen_markup(script)
        with open(out, "w") as f:
            f.write(md)
        # hyperlink to documentation
        index.append(f'- [{script.stem}]({INDEX_DIR}/{script.stem}.md)')

    # write the index file out
    with open(INDEX_FILE, "w") as fp:
        fp.write('\n'.join(index))

#
# "do it"
#
main()
print('pre_mkdocs: Done.')

