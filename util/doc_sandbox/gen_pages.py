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
import re
from pathlib import Path
import docstring_to_markdown
import docstring_to_markdown.google
import mkdocs_gen_files

# Informationally: the most basic usage of mkdocs_gen_files
# with mkdocs_gen_files.open("gen_example.md", "w") as f:
#     print("Hello, *world*!", file=f)
# mkdocs_gen_files.set_edit_path("gen_example.md", "gen_pages.py")

# verbosity
VERBOSE = False
# directory to look for scripts in (the util/ dir)
ROOT_DIR = Path('..')
# the Sandbox root, for the top-level driver scripts
SANDBOX_DIR = Path('../..')
# driver scripts at the Sandbox root that belong in the script index
SANDBOX_SCRIPTS = ['add-sims.sh', 'exp-add-sims.sh']
# directory holding the hand-written detail pages (":::" mkdocstrings stubs)
IMPL_DIR = Path('docs/implementation')
# filename of overall script index
INDEX_FILE = 'script_index.md'
# plot documentation: source (rendered separately by Local/www-doc/Makefile
# for the web UI; rendered here so it also reaches the exported/github site)
PLOT_DOC_DIR = Path('../../Local/www-doc')
# subdirectory of the mkdocs tree that the plot docs are written to
PLOT_OUT_DIR = 'plots'
# image types copied alongside the plot docs (.pptx source art is skipped)
PLOT_MEDIA_SUFFIXES = ('.png', '.jpg', '.jpeg', '.gif')

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

# google_to_markdown() always emits one final section, even for a block with
# no Google-style headings at all -- it lands as an empty "#### " heading plus
# an empty list item at the foot of the page.  Matches only the empty-titled
# section, so real trailing sections (#### Note, etc.) are left alone.
EMPTY_SECTION_RE = re.compile(r'\n+#### *\n+- *\n?$')

def google_to_markdown(block):
    r'''Convert a doc block to markdown, minus the parser's empty last section.'''
    return EMPTY_SECTION_RE.sub('\n', docstring_to_markdown.google.google_to_markdown(block))

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
    # original pre-2026, fixed namespace issue 04/2026:
    #  utility to convert google-style doc-blocks to markdown
    #  (is mostly a no-op, but can recognize headings like Args:
    #  and Note:)
    block_md = google_to_markdown(block)
    # turmon 04/2026: above may be obsolete? Consider:
    # block_md = docstring_to_markdown.convert(block)
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
    
def gen_plot_docs():
    r'''Render the Local/www-doc plot pages into the mkdocs tree.

    These pages describe individual plots and tabulations. They are
    *also* rendered, independently, by Local/www-doc/Makefile (markdown
    + m4) for the sandbox web UI, where they are deep-linked next to the
    plots themselves. Rendering them here too is what puts them into the
    exportable (github pages) documentation, which otherwise would not
    carry any plot documentation at all.

    Source-of-truth stays in Local/www-doc; nothing is written back there.
    '''
    if not PLOT_DOC_DIR.is_dir():
        # a sandbox without the www-doc tree still builds
        return
    for page in sorted(PLOT_DOC_DIR.glob('*.md')):
        # by convention "index*" files are the generated index, which the
        # mkdocs nav replaces
        if page.name.startswith('index'):
            continue
        title, body = split_plot_metadata(page.read_text())
        # the www-doc pages cross-link each other by rendered .html name;
        # mkdocs resolves links by source (.md) name instead
        body = re.sub(r'\]\(([^)/:#]+)\.html([)#])', r'](\1.md\2', body)
        out = f'{PLOT_OUT_DIR}/{page.stem}.md'
        if VERBOSE:
            print(f'{page} -> {out}')
        with mkdocs_gen_files.open(out, "w") as f:
            f.write(''.join(gen_metadata(title or page.stem)) + body)
    # supporting imagery, referenced as Media/foo.png from the pages above
    media = PLOT_DOC_DIR / 'Media'
    if media.is_dir():
        for art in sorted(media.iterdir()):
            if art.suffix.lower() not in PLOT_MEDIA_SUFFIXES:
                continue
            with mkdocs_gen_files.open(f'{PLOT_OUT_DIR}/Media/{art.name}', "wb") as f:
                f.write(art.read_bytes())


def split_plot_metadata(text):
    r'''Split a "Title: ..." metadata header off a www-doc page.

    The www-doc pages carry a python-markdown Meta-Data block --
    "<tag>: <value>" lines terminated by a blank line (see
    Local/www-doc/util/meta-extract.sh). Returns (title, body).
    '''
    lines = text.splitlines(keepends=True)
    title, n = '', 0
    for n, l in enumerate(lines):
        if not l.strip():
            n += 1
            break
        m = re.match(r'([A-Za-z][\w-]*):\s*(.*)', l)
        if not m:
            # no metadata block at all; keep the whole file
            return '', text
        if m.group(1).lower() == 'title':
            title = m.group(2).strip()
    return title, ''.join(lines[n:])


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

    # util/*.sh plus the named driver scripts at the Sandbox root
    # (add-sims.sh is step 3 of the README workflow, so it belongs here)
    shell_scripts = sorted(
        list(ROOT_DIR.glob('*.sh')) +
        [SANDBOX_DIR / n for n in SANDBOX_SCRIPTS if (SANDBOX_DIR / n).is_file()],
        key=lambda s: s.stem)
    for script in shell_scripts:
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
        # stop-list
        if script.name in ['__init__.py']:
            continue
        # if a hand-written detail page exists, link to it from the index
        details_file = IMPL_DIR / (script.stem + '.md')
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

    # plot documentation (shared source with the sandbox web UI)
    gen_plot_docs()

#
# "do it"
#
main()
