#!/usr/bin/env python
#
# The given directory (or directories) within sims/, containing an ensemble of DRMs 
# and corresponding graphical outputs, is summarized to HTML.  Note: invoke without
# the sim/ prefix.
#
# Usage:
#   html-summary.py [-r] [-i] [SIM ...]
#
# Optionally (-i), a top-level index.html will be updated as a table of contents
# for all the simulations.  You generally want this.
#
# Options are:
#  -i => generate the top-level index.html after generating any named SIM's.
#  -r => performs recursive descent (otherwise, only named files/dirs are examined)
#  -h => help
#
# Simplest usage:
#  $ html-summary.py -i SIM
#   where SIM is one of the sim ensemble directories, like HabEx_4m_TSDD_top200DD.
#   Re-generate the HTML and (with -i) the top-level index pointing to the HTML.
# Other usage:
#  $ html-summary.py -i SIM.exp
#   vs.
#  $ html-summary.py -r -i SIM.exp
#   where SIM.exp is a multi-ensemble experiment.  The call with -r summarizes
#   all ensembles under SIM.exp; without, it generates only the top-level summary
#   file, which lists a one-line overview of each ensemble.
#  $ html-summary.py -r
#   regenerates *all* indexes for everything in sims/*.
#
# Notes and Specifics
# * The output is written to the file:
#     sims/SIM/index.html
#   This file can then be viewed by starting an HTTP server on an unused port,
#   as explained below under "server".
# * The generated HTML file links directly to whatever graphical outputs are present
#   below the sims/SIM/... directory, so "make S=SIM graphics" should be run before
#   running this script.  On the other hand, if the graphics are re-made *after*
#   running this script, then a re-load of the HTML will show the new graphic files,
#   because the referenced files are new.
#
# HTTP Server: To see the HTML files and the images they point to, an HTTP server
# must be started to send them to a browser.  The method we use is to start a server 
# that listens on a non-standard port, and the browser can be directed to that port
# by giving a special URL qualifier.
# The server invocation wrapper is html-serve.sh, and it can be invoked by:
#   util/html-serve.sh status   -- to see if there is one already running
#   util/html-serve.sh start    -- to start one if not
# See that file for more.
# 
# For more on usage, use the -h option.
# Some options may be described there but not documented here.

# author:
#  Michael Turmon, JPL, 2018

# Apologies: the "business logic" is entwined with the "presentation logic."

import argparse
import sys
import os
import glob
import re
import csv
#import pickle
import cPickle as pickle
import StringIO
from collections import defaultdict, OrderedDict
import numpy as np
import astropy.units as u
#from astropy.time import Time

# image to use in case something we expect is not found
DUMMY_IMAGE = '/Local/www-resources/image-not-found.png'

# section heads, one for each category of graphic
SECTION_HEADS = {
        'radlum': '',
        'rad-sma': '',
        'duration': '''X-axis shows event duration.  Note that x-axis range varies between plots 
        to accomodate large variations in duration.
        Frequency values between plots when x-axis units are the same are comparable,
        but if units change, the frequencies are not directly comparable.
        Off-scale durations are not shown.''',
        'event-count': '',
        'yield': '''Time axis is mission clock time.''',
        'cume': '''Time axis is mission clock time.''',
        'perstar-det':  '',
        'perstar-char': '',
        'promote': '''Funnel analysis from detection to characterization, 
        separated into deep dive targets and promoted targets.  
        <p>Promotions in the
        <em>table</em> are determined
        using the promoted_stars variable output at the end of the simulation.
        Those in the <em>plots</em> are determined using detection results
        aggregated by the reduction code as the simulation proceeds.
        <p>
        Some x-axis units in this section are in terms of cumulative
        detection integration time rather than mission clock time.''',
        'earth-char': '''Promotions here are all determined 
        using the promoted_stars variable output at the end of the simulation.''',
        'path': '',
         }

def ensure_permissions(fn):
    r'''Ensure correct permissions on the named data file.  
    We use rw-rw-r-- = 664 (octal), to allow group-write.'''
    try:
        os.chmod(fn, 0o664)
    except OSError:
        pass # e.g., don't own the file

class ChangeDir:
    """Context manager for changing the current working directory."""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)
    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)
    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)
        

# This class is used for the overall index page, the single-sim 
# summary page, the per-path summary pages, and (recursively) 
# for TOC's within these pages.

class HTML_helper(object):
    r'''Primitive HTML-writer.

    Could have used a templater, but wanted to keep it self-contained.
    The HTML is put into a stringIO buffer so that it can be retrieved
    later and manipulated.  The use for this is to plunk a navbar and a
    table of contents in to the file, and these elements are accumulated 
    as <hN> elements are inserted into the main body.
    To do this, subsidiary HTML_helper's are created.'''
    ##
    ## CONSTANT DATA FOR THE CLASS
    ##
    # fixed string: HTML5 header
    opener = '''<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <!-- autogenerated by html-summary.py -->
    <title>%s</title>
    <link rel="stylesheet" href="/Local/www-resources/ensemble-reports.css">
    <link rel="shortcut icon" type="image/png" sizes="16x16" href="/Local/www-resources/favicon.png">
    <!-- Plotly.js -->
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <script src="/Local/www-resources/sorttable.js"></script>
</head>
<body>
'''
    trailer = '''
</body>
</html>
'''
    # table of contents (listing <hN> elements) is surrounded by below strings.
    # TOC placeholder inserted by calling self.toc_here(), a toc sub-object is
    # filled in during the run, and then at dump-time, the string for the sub-object 
    # is inserted at the placeholder point.  Same for nav-bar, below.
    toc_opener   = '<!-- toc goes below -->\n'
    toc_sentinel = '<!-- toc goes above -->\n'
    # nav-bar (listing <hN> elements) is surrounded by below strings.
    nav_opener   = '<!-- nav list goes below -->\n'
    nav_sentinel = '<!-- nav list goes above -->\n'
    ##
    ## CLASS METHODS
    ##
    def __init__(self, outfile, title, level=2):
        r'''Open a HTML file "outfile" with the given title.
        If outfile is None, a buffer is accumulated for later use, but no external 
        file will be created - obtain it with getvalue().  In this case, no HTML
        header or footer are sent out.'''
        global DIAGNOSE
        if DIAGNOSE:
            print 'Writing HTML to', outfile
        self.title = title
        self.outfile = outfile
        self.f_external = open(outfile, 'w') if outfile else None
        # a temporary buffer to write all HTML elements into
        self.f = StringIO.StringIO()
        # counter for HTML element IDs (shared by sub-objects)
        self.id_count = 0
        # indent level
        self.level = level
        # if we have a TOC (resp. nav) object, this is non-None
        self.toc = None
        self.nav = None
        # are we in a list?
        # (NB: lists used only for nav-bar at the moment)
        self.list_level = 0
        self.list_type = 'UnknownList' # ul, etc.
        self.list_attrs = {}
        # are we in the middle of a table?
        self.in_table = False
    def __enter__(self):
        r'''Support with ... as usage pattern.'''
        return self
    def __exit__(self, exc_type, exc_value, traceback):
        # close anything in progress
        self.close()
        # obtain the rendered HTML
        buf = self.f.getvalue()
        # insert TOC if we had one
        if self.toc:
            buf_orig = buf
            pos = buf_orig.find(self.toc_sentinel)
            buf = buf_orig[:pos] + self.toc.getvalue() + buf_orig[pos:]
        # insert nav-bar if we had one
        if self.nav:
            buf_orig = buf
            pos = buf_orig.find(self.nav_sentinel)
            buf = buf_orig[:pos] + self.nav.getvalue() + buf_orig[pos:]
        # we are done with the intermediate StringIO buffer
        self.f.close()
        if self.f_external:
            # write the external file: opener, body, trailer
            self.f_external.write(self.opener % self.title)
            self.f_external.write(buf)
            self.f_external.write(self.trailer)
            self.f_external.close()
            ensure_permissions(self.outfile)
    def indent(self, bump=0):
        self.level += bump
        return ' '*(self.level*2)
    def getvalue(self):
        r'''Return the string that is the rendered HTML for this object.'''
        return self.f.getvalue()
    def close(self):
        r'''Close the HTML object including rendering the sub-objects (TOC and nav).'''
        if self.nav:
            self.nav.close()
        if self.toc:
            self.toc.close()
        # for navbar: close out navbar's <ul> list, which come from sections,
        # which are terminated implicitly by the end of all sections
        while self.list_level > 0:
            self.f.write(self.indent(  ) + '</li>\n')
            self.f.write(self.indent(-1) + '</%s>\n' % self.list_type)
            self.list_level -= 1
        assert not self.in_table, 'Terminated HTML construction within a table'
    def toc_here(self, txt):
        r'''Insert TOC marker here, later to be replaced with a header-based TOC.'''
        # mark where the toc will be inserted later, once finished
        self.f.write(self.indent() + self.toc_opener)
        self.f.write(self.indent() + self.toc_sentinel)
        # establish the inner html-holder for the toc
        # it has no file, and thus no header or footer
        self.toc = HTML_helper(None, 'NoName', level=self.level)
        self.toc.header(txt, level=3)
    def nav_here(self):
        r'''Insert NAV marker here, later to be replaced with a header-based navbar.'''
        # mark where the nav-bar will be inserted later, once finished
        self.f.write(self.indent() + '<nav>\n')
        self.f.write(self.indent() + self.nav_opener)
        self.f.write(self.indent() + self.nav_sentinel)
        self.f.write(self.indent() + '</nav>\n')
        # establish the inner html-holder for the nav
        # it has no file, and thus no header or footer
        self.nav = HTML_helper(None, 'NoName', level=self.level)
    ##
    ## CLASS METHODS FOR HTML ELEMENTS
    ##
    def table_top(self, cols, elem_class=None):
        # allows multiple calls for multi-line table header (<th> elements)
        # by maintaining "in_table" property
        attrs = ['']*len(cols) # unused at present
        class_text = ' class="%s"' % elem_class if elem_class else ''
        if not self.in_table:
            self.f.write(self.indent() + ('<table%s>\n' % class_text) )
            self.indent(1)
            self.in_table = True
        self.f.write(self.indent() + '<tr>\n')
        self.f.write(''.join(['%s<th%s>%s</th>\n' %
                                  (self.indent(), attrs[i], cols[i]) for i in range(len(cols))]))
        self.f.write(self.indent() + '</tr>\n')
    def table_end(self):
        assert self.in_table, 'Mismatched table_top/table_end sequence'
        self.f.write(self.indent(-1) + '</table>\n')
        self.in_table = False
    def table_row(self, cols):
        self.f.write(self.indent() + '<tr>\n')
        self.f.write(''.join(['%s<td>%s</td>\n' % (self.indent(), col) for col in cols]))
        self.f.write(self.indent() + '  </tr>\n')
    def header(self, txt, level=2):
        self.id_count += 1
        marker = ('sec_%d' % self.id_count) if (self.toc or self.nav) else '' # toc id's would clash w/ outer id's
        attr_text = (' id="%s"' % marker) if marker else ''
        self.f.write(self.indent() + '<h%d%s>%s</h%d>\n' % (level, attr_text, txt, level))
        # enter the header in the TOC
        if self.toc:
            self.toc.paragraph('&nbsp;' * ((level-1)*4)) # primitive indentation
            self.toc.link('#' + marker, txt)
        # enter the header in the nav-bar
        if self.nav and level <= 3:
            attrs = {'class':('level-%d' % level)}
            self.nav.ul(attrs)
            attrs = {'class':'dir'} if level == 2 else {}
            self.nav.li('', level, attrs)
            # special-case the top-level header name (hacky - requires a <H1> in the doc)
            self.nav.link('#' + marker, txt if level > 1 else 'Sections')
    def ul(self, attrs={}):
        self.list_type = 'ul'
        if attrs:
            self.list_attrs = ' ' + ' '.join(['%s="%s"' % (key, val) for key, val in attrs.iteritems()])
        else:
            self.list_attrs = ''
    def li(self, txt, level, attrs):
        # establish the LI attributes, if any
        attr_txt = ''
        if attrs:
            attr_txt = ' ' + ' '.join(['%s="%s"' % (key, val) for key, val in attrs.iteritems()])
        # close prior LI, etc.
        if level == self.list_level:
            self.f.write(self.indent() + '</li>\n')
        elif level > self.list_level:
            self.f.write(self.indent(1) + '<%s%s>\n' % (self.list_type, self.list_attrs)) # new list
        elif level < self.list_level:
            self.f.write(self.indent(  ) + '</li>\n') # close existing list item
            self.f.write(self.indent(-1) + '</%s>\n' % self.list_type) # close existing list
        # write the LI
        self.f.write(self.indent() + '<li%s>\n' % attr_txt)
        if txt:
            self.f.write(self.indent() + '%s\n' % txt)
        self.list_level = level
    def paragraph(self, txt=''):
        self.f.write(self.indent() + '<p>' + txt + '\n')
    def text(self, txt, br=False):
        self.f.write(self.indent() + '%s\n' % txt)
        if br:
            self.f.write(self.indent() + '<br>\n')
    def div(self, txt, **attr):
        attr_txt = ''
        if attr:
            attr_txt = ' ' + ' '.join(['%s="%s"' % (key, val) for key, val in attr.iteritems()])
        self.f.write(self.indent(  ) + '<div%s>\n' % attr_txt)
        self.f.write(self.indent( 1) + '%s\n' % txt)
        self.f.write(self.indent(-1) + '</div>\n')
    def link(self, href, content, inner=False):
        line = '<A href="%s">%s</A>' % (href, content)
        if not inner:
            self.f.write(self.indent() + line + '\n')
        else:
            return line
    def script(self, path, literal=False):
        r'''Javascript literal or external file.'''
        if literal:
            if isinstance(path, basestring):
                path = [path]
            self.f.write(self.indent() + '<script type="text/javascript">\n')
            self.indent(1)
            self.f.write(self.indent() + ('\n'+self.indent()).join(path) + '\n')
            self.f.write(self.indent(-1) + '</script>\n')
        else:
            self.f.write(self.indent() + '<script src="%s"></script>\n' % path)
    def image(self, img, width=None, inner=False):
        src = '.' + img if img else DUMMY_IMAGE
        attr = ('width=%d%% ' % width) if width else ''
        # we set the image (.png, typically) background here
        line = '<A href="%s"><img %s src="%s" style="background:white"></A>\n' % (src, attr, src)
        if not inner:
            self.f.write(self.indent() + line)
        else:
            return line
    def video(self, src):
        self.f.write(self.indent()   + '<video controls>\n')
        self.f.write(self.indent(1)  + '<source src="%s" type="video/mp4">\n' % ('.'+src))
        self.f.write(self.indent()   + 'Your browser does not support the video tag.\n')
        self.f.write(self.indent(-1) + '</video>\n')



class SimSummary(object):
    # NOTE: The following two lists are the main hook for adding new plot families to the
    # generated HTML:
    #  -- The "graphics_map" associates a file pattern (e.g., /det-radlum*) to
    #  a symbolic tag ("radlum").
    #  -- The "graphics_show" gives a list of symbolic tags to put into HTML sections (in the
    #  order given), and the displayed section names.
    # list-of-pairs: (str,tag) means if filename contains 'str', it is a graphic of type 'tag'
    graphics_map = [('/det-perstar-det', 'perstar-det'),
                        ('/det-perstar-char', 'perstar-char'),
                        ('/det-radlum', 'radlum'),
                        ('/det-rad-sma', 'rad-sma'),
                        ('/det-duration', 'duration'),
                        ('/det-event-count', 'event-count'),
                        ('/det-time', 'yield'),
                        ('/det-cume', 'cume'), ('/det-detects', 'cume'), # these plots are now in yield
                        ('/det-fuel', 'cume'),
                        ('/det-delta-v', 'cume'),
                        ('/det-obstime', 'cume'),
                        ('/det-promote', 'promote'),
                        ('/det-phist', 'promote'),
                        ('/det-earth-char-', 'earth-char'),
                        ('/det-earth-char-count', 'earth-char'),
                        ('/path-ens/path-', 'path')]
    # (tag, str, target) means the as-shown-in-html name for graphics of type 'tag' is 'str',
    # and it can be remade with "make target"
    # this is a list because it is in order of presentation
    graphics_show = [
        ('radlum', 'Radius/Luminosity', 'reduce'),
        ('rad-sma', 'Radius/SMA', 'reduce'),
        ('duration', 'Event Duration', 'reduce'),
        ('event-count', 'Event Count', 'reduce'),
        ('yield', 'Mission Yield vs. Time', 'reduce'),
        ('cume', 'Mission Resources vs. Time', 'reduce'),
        ('perstar-det',  'Per-Star Detection', 'reduce'),
        ('perstar-char', 'Per-Star Characterization', 'reduce'),
        ('promote', 'Target Promotion', 'reduce'),
        ('earth-char', 'Earth Characterizations', 'reduce'),
        ('path', 'Full-Ensemble Path', 'path-ensemble')]
    # list-of-pairs: (str,tag) means if filename contains 'str', it is a table of type 'tag'
    tables_map = [('/table-funnel', 'promote'),
                      ]
    # (tag: str) means the as-shown-in-html name for table of type 'tag' is 'str'
    tables_show = {'promote': 'Target Promotion',
                       }
    def __init__(self, sim_dir, name=''):
        self.name = name if name else sim_dir
        self.Ndrm = 0
        self.graphics = defaultdict(list) # dict-of-list
        self.tables = defaultdict(list) # dict-of-list
        self.path_graphics = defaultdict(dict) # dict-of-dict
        self.csvs = []

    def add(self, filename):
        if '/path/' in filename:
            self.add_path_graphic(filename)
        elif filename.endswith('.png') or filename.endswith('.pdf'):
            self.add_graphic(filename)
        elif filename.endswith('.md'):
            self.add_table(filename)
        elif '/drm/' in filename:
            self.add_drm(filename)
        elif filename.endswith('.csv'):
            self.add_csv(filename)
        # skip some explicitly
        elif '/log/' in filename:
            pass
        elif '/spc/' in filename:
            pass
        elif '/sys/' in filename:
            pass
        elif '/run/' in filename:
            pass
        
    def add_path_graphic(self, filename):
        r'''Separate path graphics into types.  Unlike add_graphic, this is not table-driven.

        TODO: Getting a little cumbersome at this point.
        '''
        if '-frames/' in filename:
            pass # skip frame-by-frame graphics
        elif filename.endswith('-final.png'):
            seed = re.sub(r'.*/([0-9]+)-final\.png', r'\1', filename)
            self.path_graphics[seed]['final'] = filename
        elif filename.endswith('.mp4'):
            seed = re.sub(r'.*/([0-9]+)\.mp4', r'\1', filename)
            self.path_graphics[seed]['movie'] = filename
        elif filename.endswith('obs-timeline.png'):
            seed = re.sub(r'.*/([0-9]+)-obs-timeline\.png', r'\1', filename)
            self.path_graphics[seed]['obs-timeline'] = filename
        elif filename.endswith('obs-keepout-all.png'):
            seed = re.sub(r'.*/([0-9]+)-obs-keepout-all\.png', r'\1', filename)
            self.path_graphics[seed]['obs-keepout-all'] = filename
        elif filename.endswith('obs-keepout-char.png'):
            seed = re.sub(r'.*/([0-9]+)-obs-keepout-char\.png', r'\1', filename)
            self.path_graphics[seed]['obs-keepout-char'] = filename
        elif filename.endswith('.png'):
            seed = re.sub(r'.*/path/([0-9]+)-cume/.*\.png', r'\1', filename)
            # split out corona/starshade
            key = 'shade' if filename.endswith('-shade.png') else 'corona'
            try:
                self.path_graphics[seed][key].append(filename)
            except KeyError:
                self.path_graphics[seed][key] = [filename]

    def add_graphic(self, filename):
        if filename.endswith('.pdf'):
            return
        for fn_str, tag in self.graphics_map:
            if fn_str in filename:
                self.graphics[tag].append(filename)
                break

    def add_table(self, filename):
        for fn_str, tag in self.tables_map:
            if fn_str in filename:
                self.tables[tag].append(filename)
                break

    def add_drm(self, filename):
        if filename.endswith(('.pkl', '.drm')):
            self.Ndrm += 1
        return
        
    def add_csv(self, filename):
        self.csvs.append(filename)
        return
        
    def render_table(self, hh, t_file):
        r'''Renders a table, in a markdown-format file, as html.

        Table is expected to be in the form:

        | colhead1 | colhead2 | ... | colheadN |
        | ---      | ---      | ...            |
        | data 1_1 | data 1_2 | ...            |
        ...
        | data 2_1 | data 2_2 | ...            |
        <newline>
        
        with an optional header (# Title) in front of the column heads.
        Markdown puts a |---|---| style line between the header rows 
        and the data rows (allowing for multi-line heads), and that this 
        code requires such a separation.
        '''
        with open(t_file, 'r') as fp:
            # 0 = not in a table; 1 = in table header; 2 = in table body
            in_table = 0
            for line_orig in fp:
                line = line_orig.strip()
                if len(line) == 0:
                    # blank line -> end table
                    if in_table > 0: hh.table_end()
                    in_table = 0
                elif line.startswith('#'):
                    # section header line -> end table (if any) + write section header
                    if in_table > 0: hh.table_end()            
                    in_table = 0
                    # level = 4 is styled as plain/bold in the HTML
                    hh.header(line.lstrip('#'), level=4)
                elif line.startswith('|') and '---' in line:
                    # table head ends and table body begins -> nothing to write
                    in_table = 2
                elif line.startswith('|') and in_table <= 1:
                    # header row of table (first line, or continuing)
                    in_table = 1 # in case it was first line
                    cols = line.split('|')
                    hh.table_top(cols[1:-1])
                elif line.startswith('|') and in_table > 1:
                    # table body
                    cols = line.split('|')
                    hh.table_row(cols[1:-1])
                else:
                    print 'Unexpected table line <%s>' % line
            # be sure to close the table out
            if in_table > 0:
                hh.table_end()            

    def render(self, filename, uplink='Sandbox'):
        #print self.name
        with HTML_helper(filename, self.name) as hh:
            # nav-bar
            hh.nav_here()
            # title - h1 tag
            hh.header(os.path.basename(self.name), level=1)
            # navigation link
            hh.link('../../', 'Up to %s Root' % uplink)
            # summary
            hh.header('Ensemble Summary')
            hh.paragraph('Ensemble: %s' % self.name)
            hh.paragraph('Size: %d sims\n' % self.Ndrm)
            # table of contents
            hh.toc_here('Contents')
            # overall images
            for tag, shown, target in self.graphics_show:
                # overall section header
                hh.header(shown + ' Plots')
                # optional caption
                if SECTION_HEADS[tag]:
                    hh.paragraph(SECTION_HEADS[tag])
                # tables, if any, are at top (ad hoc at the moment)
                if tag in self.tables_show and len(self.tables[tag]) > 0:
                    # print 'Rendering', self.tables[tag]
                    hh.header('Tables for %s' % shown, level=3)
                    for table in self.tables[tag]:
                        self.render_table(hh, table)
                # plot widgets for star-target info
                if tag == 'perstar-det':
                    hh.header('Interactive Detection Plot Widget', level=3)
                    hh.div('<!-- det plot goes here -->', id='detPlotDiv', style='width: 900px; height: 700px;')
                    hh.div('Detection QOI for Plot Shading: <select class="det_qoi_select"> </select>')
                elif tag == 'perstar-char':
                    hh.header('Interactive Characterization Plot Widget', level=3)
                    hh.div('<!-- char plot goes here -->', id='charPlotDiv', style='width: 900px; height: 700px;')
                    hh.div('Characterization QOI for Plot Shading: <select class="char_qoi_select"> </select>')
                    # assume dets always accompany chars, and insert the script at this point
                    hh.script('/Local/www-resources/star-target-plots.js')
                elif tag == 'path' and self.graphics[tag]:
                    hh.header('Interactive Ensemble Path Widget', level=3)
                    hh.div('<!-- ensemble path plot goes here -->',
                               id='slewPlotDiv', style='width: 1000px; height: 600px;')
                    # set up path to the data CSVs, and label mode, for the JS viewer
                    hh.script('var ens_path_root = "../path-ens";', literal=True)
                    hh.script([
                        'var ens_path_root = "../path-ens";',
                        'var ens_path_mode = "ensemble";'],
                              literal=True)
                    hh.script('/Local/www-resources/ens-path-plots.js')
                plots = self.graphics[tag]
                if not plots:
                    hh.paragraph('No such plots.  Generate with: make S=... %s' % target)
                else:
                    hh.table_top(('Filename', 'Plot Preview'), elem_class='gfx')
                    for p in sorted(plots):
                        img = hh.image(p, width=100, inner=True)
                        hh.table_row((os.path.basename(p), img))
                    hh.table_end()
            # specific path movies
            hh.header('Single-Run Paths')
            path_movie_seeds = sorted(self.path_graphics.keys())
            if not path_movie_seeds:
                hh.paragraph('No path movies yet.  Fix with: make S=... path-movie-N\n')
            else:
                hh.table_top(('Seed', 'Final Frame'), elem_class='gfx')
                for seed in path_movie_seeds:
                    info = self.path_graphics[seed]
                    img_target = info['final'] if 'final' in info else None
                    alink = hh.link('path-%s.html' % seed, seed, inner=True)
                    img = hh.image(img_target, width=100, inner=True)
                    hh.table_row((alink, img))
                hh.table_end()

    def render_path(self, seed, filename):
        r'''Write one simulation's path information to html.'''
        #print 'Rendering path to', filename
        info = self.path_graphics[seed]
        with HTML_helper(filename, 'Path from %s' % self.name) as hh:
            # title (as h1)
            hh.header('Observing Path Summary', level=1)
            # navigation link
            hh.link('index.html', 'Return to Ensemble Root')
            # summary
            hh.header('Overview')
            hh.paragraph('Ensemble: %s\n' % self.name)
            hh.paragraph('Seed: %s\n' % seed)
            hh.toc_here('Contents')
            # path timeline
            hh.header('Observing Timeline')
            if 'obs-timeline' in info:
                img_target = info['obs-timeline']
                hh.image(img_target, width=70)
            else:
                hh.paragraph('Timeline not available.  Fix with: make ... obs-timeline-N')
            # path keepout
            hh.header('Keepout and Observations')
            if 'obs-keepout-all' in info:
                img_target = info['obs-keepout-all']
                hh.image(img_target, width=70)
                img_target = info['obs-keepout-char']
                hh.image(img_target, width=70)
            else:
                hh.paragraph('Keepout timeline not available.  Fix with: make ... keepout-N')
            # path summary
            hh.header('Path Overview')
            img_target = info['final'] if 'final' in info else None
            hh.image(img_target, width=70)
            ### insert the single-drm path widget
            hh.header('Interactive Observational Tour Widget', level=3)
            if 'movie' in info:
                # (movie is a proxy for the path CSV data)
                hh.div('<!-- tour path plot goes here -->',
                           id='slewPlotDiv', style='width: 1000px; height: 600px;')
                # point the path viewer, below, to the right directory
                hh.script([
                    'var ens_path_root = "../path/%s-cume";' % seed,
                    'var ens_path_mode = "single";'],
                              literal=True)
                hh.script('/Local/www-resources/ens-path-plots.js')
            else:
                hh.paragraph('Data for widget not available.  Fix with: make ... path-movie-N')
            # path movie
            hh.header('Path Movie of This Tour')
            if 'movie' in info:
                hh.video(info['movie'])
            else:
                hh.paragraph('No path movie available.  Fix with: make ... path-movie-N')
            # cumulative images
            if 'corona' in info:
                hh.header('Cumulative Coronagraph Keepout for This Tour')
                hh.table_top(('Filename', 'Plot Preview'), elem_class='gfx')
                for img in sorted(info['corona']):
                    line = hh.image(img, width=100, inner=True)
                    hh.table_row((os.path.basename(img), line))
                hh.table_end()
            if 'shade' in info:
                hh.header('Cumulative Occulter Keepout for This Tour')
                hh.table_top(('Filename', 'Plot Preview'), elem_class='gfx')
                for img in sorted(info['shade']):
                    line = hh.image(img, width=100, inner=True)
                    hh.table_row((os.path.basename(img), line))
                hh.table_end()

    def render_paths(self, filename_tmpl):
        r'''Write path information, typically from multiple seeds, to separate html.'''
        for seed in self.path_graphics.keys():
            self.render_path(seed, filename_tmpl % seed)
                    
def index_ensemble(args, path_sim):
    r'''Make an index for one ensemble.'''
    if not os.path.isdir(os.path.join(path_sim, 'drm')):
        return # not an ensemble -- e.g., a recursive call into an arb. dir
    outdir = 'html'
    print 'Indexing:', path_sim
    with ChangeDir(path_sim):
        # ensure the output directory
        if not os.path.isdir(outdir):
            os.makedirs(outdir, 0o775)
        # gather information
        sim_info = SimSummary('.', name=path_sim)
        for root, dirs, files in os.walk('.'):
            level = root.replace(path_sim, '').count(os.sep)
            # (delete certain "dirs" to avoid walking them)
            for fn in files:
                fn_full = os.path.join(root, fn)
                if fn.startswith('.'):
                    continue # metadata files often of form ._FOO
                #print 'Processing', fn_full
                sim_info.add(fn_full)
        # render information
        fn = os.path.join(outdir, 'index.html')
        sim_info.render(fn)
        for seed in sim_info.path_graphics.keys():
            sim_info.render_paths(os.path.join(outdir, 'path-%s.html'))

def exp_summary(d):
    r'''Summarize one experiment directory as an OrderedDict of strings.'''
    # For now, this summary (for multi-ensembles) just uses the same
    # summary file as the single-ensemble case.
    rv = sim_summary(d)
    # could glob for these, but not worth the time
    unknown_value = '(N/A)' # ()'s help sort order
    rv['path_count'] = unknown_value
    rv['gfx_count'] = unknown_value
    return rv

def sim_summary(d):
    r'''Summarize one simulation directory as an OrderedDict of strings.

    If the argument is None, instead, return the header corresponding to the list.'''
    # do not display floats to 12 digits of precision
    def fmt_float(x):
        return ('%.3f' % float(x));
    def fmt_int(x):
        return ('%d' % int(x));
    def fmt_str(x):
        return x
    # properties we want to record: key, label, formatter
    props = [
        ('ensemble_size', 'Ens. Size', fmt_str),
        ('runtime', 'Reduction Date', fmt_str),
        ('user', 'User', fmt_str),
        ('detections_earth_unique', 'Earths (Det.)', fmt_float),
        ('chars_earth_unique', 'Earths (Char.)', fmt_float),
        ('chars_earth_strict', 'Earths (Strict)', fmt_float),
        ('gfx_count', 'Ens. Graphs', fmt_int),
        ('path_count', "Path Summ's", fmt_int),
        ('path_gfx', "Path Graphs", fmt_int),
        ]
    # Return the header columns if None was input
    if d is None:
        rv = OrderedDict()
        for dtag, dname, _ in props: rv[dtag] = dname
        return rv
    # number of unique paths that were indexed in any way (fn = SEED-foo-bar.png)
    path_count = len(set(fn.split('-')[0] for fn in glob.glob(os.path.join(d, 'path/[0-9]*-*.*'))))
    # number of path graphics
    path_gfx = (len(glob.glob(os.path.join(d, 'path/[0-9]*.png'))) +
                len(glob.glob(os.path.join(d, 'path/[0-9]*.mp4'))))
    # total graphic file counts (non-path)
    gfx_count  = (len(glob.glob(os.path.join(d, 'gfx/*.png'))) +
                  len(glob.glob(os.path.join(d, 'path-ens/*.png'))))
    # grab the metadata for d from the CSV written by reduce-drms.py
    info_fn = os.path.join(d, 'reduce-info.csv')
    unknown_value = '(N/A)' # ()'s help sort order
    try:
        with open(info_fn) as f:
            info_items = csv.DictReader(f);
            info = info_items.next() # it is a 1-line csv
        # insert the extra data, just retrieved above
        info['path_count'] = path_count
        info['path_gfx'] = path_gfx
        info['gfx_count'] = gfx_count
        # make the return value
        rv = OrderedDict()
        for dtag, dname, dfmt in props:
            if dtag in info:
                rv[dtag] = dfmt(info[dtag])
            else:
                rv[dtag] = unknown_value
    except IOError:
        # order does not matter
        return {dtag: unknown_value for dtag, _, _ in props}
    return rv

def index_all(args, startpath, title, uplink=None):
    r'''Make a summary table of all ensembles below startpath.

    If args.recursive, descend recursively and summarize child directories.'''
    print 'Indexing:', startpath
    # output HTML
    filename = os.path.join(startpath, 'index.html')
    with HTML_helper(filename, title) as hh:
        # title as H1
        hh.header(title, level=1)
        # navigation link
        if uplink:
            hh.link('../', 'Up to %s Root' % uplink)
        # table of individual sims
        hh.header('Ensembles')
        # make the table be sortable so that the JS sorter knows about it
        hh.table_top(['Name'] + sim_summary(None).values(), elem_class='sortable')
        for root, dirs, files in os.walk(startpath):
            # break below => loops once => "dirs" is just top-level subdirs of startpath
            for d in sorted(dirs):
                # link to the sim below, summarize the sim in a list of a few properties
                if d.endswith('.exp'):
                    alink = hh.link('%s/index.html' % d, d, inner=True)
                    properties = exp_summary(os.path.join(root, d))
                    hh.table_row([alink] + properties.values())
                elif os.path.isdir(os.path.join(root, d, 'drm')):
                    alink = hh.link('%s/html/index.html' % d, d, inner=True)
                    properties = sim_summary(os.path.join(root, d))
                    # format the properties as a row
                    hh.table_row([alink] + properties.values())
                if args.recurse:
                    # recurse (1 level max) down into sims/d -- no-op if no drm/ there
                    index_ensemble(args, os.path.join(startpath, d))
            break # important: goes only one level deep!
        hh.table_end()

############################################################
#
# Main routine
#
############################################################

def main(args):
    r'''Main routine: Summarize files below a given location to HTML.'''

    # This DIAGNOSE is intended to be for error-reporting to diagnose
    # issues deep within this code, but not for ordinary use
    # We want to be able to get to it anywhere, so it's global.
    global DIAGNOSE
    DIAGNOSE = args.DIAGNOSE

    # this supports running with: CMD -i SIM, which will generate the
    # index for SIM, and then regenerate the global index
    for sim in args.sim:
        sim_dir = os.path.join('sims', sim)
        if sim.endswith('.exp'):
            index_all(args, sim_dir, 'Experiment Root: ' + sim.split('/')[-1], uplink='Sandbox')
        else:
            index_ensemble(args, sim_dir)
    # re-do top-level index to sync with results of loop above
    if args.index:
        # don't let index_all() recurse into sims/* if args were named; if args given,
        # the call here is only to generate the overall root-level index.
        if len(args.sim) > 0: args.recurse = False
        index_all(args, 'sims', 'Simulation Ensemble Root')
    
    
# not all options are currently used
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Summarize EXOSIMS sims to HTML.",
                                     epilog='')
    parser.add_argument('sim', metavar='SIM', nargs='*', default=[],
                            help='sim ensemble directory or directories')
    parser.add_argument('-i', '--index', help='generate index (non-recursive by default)',
                      dest='index', action='store_true', default=False)
    parser.add_argument('-r', '--recurse', help='recursive over all sims',
                      dest='recurse', action='store_true', default=False)
    parser.add_argument('-D', '--diagnose', help='diagnostic (debugging) output',
                      dest='DIAGNOSE', action='count', default=0)
    args = parser.parse_args()
    
    # set umask in hopes that files/dirs will be group-writable
    os.umask(0o002)

    main(args)
    sys.exit(0)



