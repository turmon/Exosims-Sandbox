r''' 
gen_pages.py: extract documentation blocks into markdown

plug-in to extract documentation blocks from shell and
python scripts.

Designed to run as an mkdocs plugin. It runs under "mkdocs build"
and programmatically looks over all such scripts and extracts 
their top-level documention block.

One feels there must be a built-in capability to
do some of this automatically, but I could not find it.
'''

import os
import glob
import ast
from pathlib import Path
import docstring_to_markdown
import mkdocs_gen_files

# Informationally: the most basic usage of mkdocs_gen_files
# with mkdocs_gen_files.open("gen_example.md", "w") as f:
#     print("Hello, *world*!", file=f)
# mkdocs_gen_files.set_edit_path("gen_example.md", "gen_pages.py")

# verbosity
VERBOSE = False
# directory to look for scripts in
ROOT_DIR = Path('..')
# directory for detailed views of script implementations
IMPL_DIR = Path('util_detail')
# filename of overall script index
INDEX_FILE = 'script_index.md'

def find_oneliner(block, title):
    r'''Find one-line description of form title: Description

    Both the block-comment .sh and .py files have a one-line
    description at the top, like:
       drm_process: process drms in a good way
    This routine looks for such descriptions in text.''' 
    line1 = ''
    for l in block.split('\n'):
        if title in l:
            line1 = l[l.find(title)+len(title):].strip(':\n')
            if line1.startswith('.sh'):
                line1 = line1[3:]
            if line1.startswith('.py'):
                line1 = line1[3:]
            if line1.startswith(':'):
                line1 = line1[1:].strip()
            break
    return line1
        
def gen_metadata(title):
    '''Generate metadata for a clean-looking title

    Function names like drm_process would be manipulated
    by the usual title maker - this forces the filename
    to be used as the title.'''
    chunk = []
    chunk.append('---\n')
    chunk.append(f'title: {title}\n')
    chunk.append('---\n')
    return chunk

def get_doc_block_sh(script, title):
    # line1 is the 
    line1 = ''
    with open(script, 'r') as orig:
        chunk = gen_metadata(title)
        for l in orig.readlines():
            # skip shebang lines
            if l.startswith('#!'):
                continue
            # the ## marks end-of-block
            if l.startswith('##'):
                break
            # non-# (including blank) is end-of-block
            if not l.startswith('#'):
                break
            if not line1 and title in l:
                line1 = l[l.find(title)+len(title):].strip(':\n')
                if line1.startswith('.sh'):
                    line1 = line1[3:]
                if line1.startswith('.py'):
                    line1 = line1[3:]
                if line1.startswith(':'):
                    line1 = line1[1:].strip()
            if l.startswith('# '):
                chunk.append(l[2:])
            elif l.startswith('#'):
                chunk.append(l[1:])
    block = ''.join(chunk)
    # utility to convert google-style doc-blocks to markdown
    # (is mostly a no-op, but can recognize headings like Args:
    # and Note:)
    block_md = docstring_to_markdown.google_to_markdown(block)
    return block_md, line1

def get_doc_block_py(fn, stem):
    r'''Pulls out a #-delimited doc block from .py files'''
    # we can use the same approach as for shell scripts
    doc, line1 = get_doc_block_sh(fn, stem)
    return doc, line1

def get_doc_py(script, stem):
    # ast approach only parses out the __doc__ string,
    # it does not execute the code like importlib
    # approaches would
    with open(script) as f:
        doc = ast.get_docstring(ast.parse(f.read()))
    if doc:
        line1 = find_oneliner(doc, stem)
        metadata = ''.join(gen_metadata(stem))
        return metadata + doc, line1
    doc, line1 = get_doc_block_py(script, stem)
    if doc:
        return doc, line1
    return 'No documentation found.', stem
    
###
### Main routine
###

def main():
    # index file line-by-line (no newlines)
    index = []
    index.append('---')
    index.append(f'title: Script Index')
    index.append('---')

    # Shell script portion
    index.append('')
    index.append('## Shell Scripts')
    index.append('')

    for script in sorted(ROOT_DIR.glob('*.sh')):
        out = script.stem + '.md'
        if VERBOSE: 
            print(f'{script} -> {out}')
        md, line1 = get_doc_block_sh(script, script.stem)
        with mkdocs_gen_files.open(out, "w") as f:
            f.write(md)
        # hyperlink to documentation
        index.append(f'- [{script.stem}]({script.stem}.md): {line1}')

    # Python script portion
    index.append('')
    index.append('## Python Scripts')
    index.append('')

    for script in sorted(ROOT_DIR.glob('*.py')):
        # if the implementation dir has an entry, find it
        # and generate a link in the index
        details_file = IMPL_DIR / script.name
        if VERBOSE:
            print(f'looking for {details_file}...')
        has_detail = os.path.isfile(details_file)
        if has_detail:
            detail_info = f' Implementation [details](implementation/{details_file.stem}.md)'
        else:
            detail_info = ''
        # make the basic link to the documentation block
        out = script.stem + '.md'
        if VERBOSE: 
            print(f'{script} -> {out}')
        md, line1 = get_doc_py(script, script.stem)
        with mkdocs_gen_files.open(out, "w") as f:
            f.write(md)
        # hyperlink to documentation
        index.append(f'- [{script.stem}]({script.stem}.md): {line1}{detail_info}')

    # write the index file out
    with mkdocs_gen_files.open(INDEX_FILE, "w") as fp:
        fp.write('\n'.join(index))

#
# "do it"
#
main()
