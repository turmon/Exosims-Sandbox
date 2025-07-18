#!/usr/bin/env python
#
# html-summary.py: Generate HTML-format summary of graphics and tables
# 
# The given directory (or directories) within `sims/`, containing an ensemble of DRMs 
# and corresponding graphical outputs, is summarized to HTML.
# Note that arguments should not contain the `sims/` prefix.
#
# Usage:
#   `html-summary.py [-r] [-i] [SIM ...]`
#
# Optionally (using `-i`), a top-level index.html is updated as a table of contents
# for all the simulations. You typically want this.
# Also for `-i`, any intermediate index.html's between `sims/index.html` and
# `sims/SIM/index.html` will also be made; this is needed when SIM is nested ("families").
#
# Options are:
# 
# +  -i => generate intermediate index.html's after generating any named SIM indexes.
# +  -r => performs recursive descent (otherwise, only named SIMS are indexed)
# +  -h => help
#
# Simplest usage:
# ```
#  $ html-summary.py -i SIM
# ```
#   where SIM is one of the sim ensemble directories, like HabEx_4m_TSDD.
#   Re-generates the HTML for this ensemble and the top-level index.
# Other usage:
# ```
#  $ html-summary.py -i SIM.exp
#  $ html-summary.py -r -i SIM.exp
# ```
#   where SIM.exp is a multi-ensemble experiment.  The usage with "-r" summarizes
#   all ensembles under SIM.exp; without "-r", it generates only the top-level summary
#   file, which lists a one-line overview of each ensemble.
# ```
#  $ html-summary.py -r
# ```
#   regenerates *all* indexes for everything in sims/*.
#
# Notes and Specifics
# 
# * The output is written to the file:
#     sims/SIM/index.html
#   This file can be viewed by starting an HTTP server on an unused port, see below.
# * The generated HTML file links directly to whatever graphical outputs are present
#   below the sims/SIM/... directory.  Use "make S=SIM html" to make the graphics files
#   and then call this script to make the graphics viewable through the webpage.
#
# HTTP Server: To see the HTML files and the images they point to, an HTTP server
# must be started to send them to a browser.  We start a server listening on a
# non-standard port, and the browser is directed to that port with a special URL qualifier.
# The server invocation wrapper is html-serve.sh, and it can be invoked by:
#   make html-status   -- to see if there is one already running
#   make html-start    -- to start one if not
# See the Makefile, or util/html-serve.sh, for more.
# 
# For more on usage, use the -h option.
# Some options may be described there but not documented here.

# author:
#  Michael Turmon, JPL, 2018, 2020, 2022

# Apologies:
# * The "business logic" is entwined with the "presentation logic."
# * A superior approach would have been to generate this page in
# markdown, and then convert the markdown to html with python's markdown library.
# We wouldn't have had as fine control over the HTML (actually a minor plus),
# and we wouldn't have to worry about matching tags, for a huge win.
# Notes: profiling (2023-10) reveals the glob.glob's in sim_summary dominate 
# runtime, see below.

import argparse
import sys
import os
import glob
import json
import time
import re
import csv
import html
#import pickle
import six.moves.cPickle as pickle
from collections import defaultdict, OrderedDict
from io import StringIO
from pathlib import Path

# path to resources
WWW_RES = Path('/Local/www-resources')
# doc directory within this
WWW_DOC = WWW_RES / 'doc'

# image to use in case something we expect is not found
DUMMY_IMAGE = WWW_RES / 'image-not-found.png'

# section heads, one for each category of graphic
SECTION_HEADS = {
        'radlum': '',
        'rad-sma': '''"Throughput" is the proportion of planets present, that were characterized
            in the way indicated. Not-a-number entries in this plot correspond to 0/0 conditions
            arising due to planet types not present in the selected population.
            The "Popuation" plot allows verification that planets are being generated at the 
            correct rate, e.g., verification of eta<sub>Earth</sub>.''',
        'duration': '''X-axis shows event duration.  Note that x-axis range varies between plots 
            to accomodate large variations in duration.
            Frequency values between plots when x-axis units are the same are comparable,
            but if units change, the frequencies are not directly comparable.
            Off-scale durations are not shown.''',
        'event-count': '''All detection and characterization counts are attempts, 
             without regard to success.''',
        'visit-time': '''Counts in these plots represent the number of target stars visited or re-visited. 
             <p>Time axis is mission clock time. Time-axis for the various lines is offset slightly to 
             reduce overplotting of multiple time series.''',
        'yield': '''<b>Plots:</b> Time axis is mission clock time. 
             The time axis for the various quantities (All/Unique/Revisit)
             is offset slightly to 
             reduce overplotting of multiple time series.
             <p>
             <b>Characterization plots:</b> If the observing scenario
             contains multiple observing modes, the "revisit" line will
             count the second and following modes as revisits. For example,
             if there are characterization observations in VIS and 3 NUV wavelengths,
             the 3 subsequent NUV observations will be counted as revisits.
             <p>
             <b>Plot/Table relationship:</b> Table entries for
             Detected (cumulative) give the final
             achieved value shown on the Detections plots.''',
        'cume': '''Time axis is mission clock time. 
             Time-axis for the various lines is offset slightly to 
             reduce overplotting of multiple time series.''',
        'perstar-det':  '',
        'perstar-char': '',
        'det-funnel': '''Progression of detection observations, separated into
            all observed stars, and only promoted stars.
            <p>
            Attempts means all detection attempts on that target. 
            Each attempt ends in either failure (further subdivided into 
            one of three modes), or success.
            Repeat detection attempts result in a count of successful detections.
            ''',
        'promote': '''Funnel analysis from detection to characterization, 
            separated into deep dive targets and promoted targets.  
            <p>Promotions in the
            <em>table</em> are determined
            using the promoted_stars variable output at the end of the simulation.
            Those in the <em>plots</em> are determined using detection results
            aggregated by the reduction code as it processes the DRM.
            <p>
            <em>The reduction code is separate from the EXOSIMS scheduler. 
            The reduction code's post-hoc processing, seen in the plots, 
            may not accurately reproduce scheduler promotion behaviors.</em>
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
        

############################################################
#
# HTML Generation
#
# This class is used for the overall index page, the single-sim 
# summary page, the per-path summary pages, and (recursively) 
# for TOC's within these pages.
############################################################

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
    <!-- tabulator.js -->
    <link href="https://unpkg.com/tabulator-tables@6.3.1/dist/css/tabulator.min.css" rel="stylesheet">
    <script type="text/javascript" src="https://unpkg.com/tabulator-tables@6.3.1/dist/js/tabulator.min.js"></script>
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
            print('Writing HTML to', outfile)
        self.title = title
        self.outfile = outfile
        self.f_external = open(outfile, 'w') if outfile else None
        # a temporary buffer to write all HTML elements into
        self.f = StringIO()
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
        # where are we in a table? (0=outside; 1=thead; 2=tbody; 3=tfoot)
        self.in_table = 0
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
        assert self.in_table == 0, 'Terminated HTML construction within a table'
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
        # by maintaining "in_table" => 0=outside; 1=thead; 2=tbody; 3=tfoot
        attrs = ['']*len(cols) # unused at present
        class_text = ' class="%s"' % elem_class if elem_class else ''
        if self.in_table == 0:
            self.f.write(self.indent() + ('<table%s>\n' % class_text) )
            self.f.write(self.indent(1) + '<thead>\n')
            self.in_table = 1
        elif self.in_table != 1:
            assert False, 'Got table_top within a table body or foot'
        self.f.write(self.indent() + '<tr>\n')
        self.f.write(''.join(['%s<th%s>%s</th>\n' %
                                  (self.indent(), attrs[i], cols[i]) for i in range(len(cols))]))
        self.f.write(self.indent() + '</tr>\n')
    def table_end(self):
        if self.in_table == 0:
            assert False, 'Mismatched table_top/table_end sequence'
        elif self.in_table == 1:
            self.f.write(self.indent() + '</thead>\n') # empty is OK
        elif self.in_table == 2:
            self.f.write(self.indent() + '</tbody>\n')
        elif self.in_table == 3:
            self.f.write(self.indent() + '</tfoot>\n')
        self.f.write(self.indent(-1) + '</table>\n')
        self.in_table = 0
    def table_row(self, cols):
        if self.in_table == 0:
            assert False, 'Got table_row outside table_top/table_end'
        elif self.in_table == 1:
            self.f.write(self.indent() + '</thead>\n')
            self.f.write(self.indent() + '<tbody>\n')
            self.in_table = 2 # thead -> tbody
        elif self.in_table == 2:
            pass # <tr> within tbody
        elif self.in_table == 3:
            pass # <tr> within tfoot is OK
        self.f.write(self.indent() + '<tr>\n')
        self.f.write(''.join(['%s<td>%s</td>\n' % (self.indent(), col) for col in cols]))
        self.f.write(self.indent() + '  </tr>\n')
    def table_foot(self):
        if self.in_table == 1:
            self.f.write(self.indent() + '</thead>\n') # head->foot OK
        elif self.in_table == 2:
            self.f.write(self.indent() + '</tbody>\n') # body->foot
        else:
            assert False, 'Unexpected table_foot'
        self.f.write(self.indent() + '<tfoot>\n')
        self.in_table = 3
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
            self.list_attrs = ' ' + ' '.join(['%s="%s"' % (key, val) for key, val in attrs.items()])
        else:
            self.list_attrs = ''
    def li(self, txt, level, attrs):
        # establish the LI attributes, if any
        attr_txt = ''
        if attrs:
            attr_txt = ' ' + ' '.join(['%s="%s"' % (key, val) for key, val in attrs.items()])
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
    def paragraph(self, txt='', br=False):
        self.f.write(self.indent() + '<p>' + txt + '\n')
        if br:
            self.f.write(self.indent() + '<br>\n')
    def text(self, txt, br=False):
        self.f.write(self.indent() + '%s\n' % txt)
        if br:
            self.f.write(self.indent() + '<br>\n')
    def div(self, txt, **attr):
        attr_txt = ''
        if attr:
            attr_txt = ' ' + ' '.join(['%s="%s"' % (key, val) for key, val in attr.items()])
        self.f.write(self.indent(  ) + '<div%s>\n' % attr_txt)
        self.f.write(self.indent( 1) + '%s\n' % txt)
        self.f.write(self.indent(-1) + '</div>\n')
    def span(self, txt, **attr):
        '''Returns a string but does not render to the HTML doc.'''
        attr_txt = ''
        if attr:
            attr_txt = ' ' + ' '.join(['%s="%s"' % (key, val) for key, val in attr.items()])
        rv = ('<span%s>' % attr_txt) + txt + '</span>'
        return rv
    def link(self, href, content, inner=False, **attr):
        attr_txt = ''
        if attr:
            attr_txt = ' ' + ' '.join(['%s="%s"' % (key, val) for key, val in attr.items()])
        line = f'<A href="{href}"{attr_txt}>{content}</A>'
        if not inner:
            self.f.write(self.indent() + line + '\n')
        else:
            return line
    def script(self, path, literal=False):
        r'''Javascript literal or external file.'''
        if literal:
            # bare str -> singleton list for join below
            if isinstance(path, str):
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


############################################################
#
# Summarize one Ensemble Simulation to HTML
#
############################################################

# container class for plot information
class GraphicsDescription(object):
    # title: the as-shown-in-html name for graphics of this type
    # target: graphics can be remade with "make S=... target"
    # infotype: the type of info (sometimes it is both tables and plots)
    # filename: for descriptive webpage (parallel with TableDescription below)
    title = 'Plot' # always replaced
    target = ''
    infotype = 'Plots' # sometimes Tables + Plots
    filename = None
    def __init__(self, title, target, infotype=None, filename=None):
        self.title = title
        self.target = target
        if infotype:
            self.infotype = infotype
        if filename:
            self.filename = filename
    
# container class for table information
class TableDescription(object):
    title = 'Table'
    text = 'This is a table.'
    filename = ''
    def __init__(self, title, text, filename):
        self.title = title
        self.text = text
        self.filename = filename
    

class SimSummary(object):
    r'''Summarize one ensemble to HTML.

    Algorithm: First step, load all "interesting" files into various
    dictionaries depending on their filename (the "add" method), following
    the patterns in the graphics_map variable.
    Second step, the "render" method, iterates over a list of HTML-page sections, 
    encoded in the graphics_show variable, extracting any graphics files that 
    were added above and inserting those files into an HTML template 
    for the ensemble.'''
    # do not index content of dirs having these names
    # ('drm' not skipped b/c we use it to count ensemble size)
    sim_dir_no_index = set(('export', 'log', 'spc', 'sys', 'log_sim', 'run', 'html'))
    # NOTE: The following two lists are the main hook for adding new plot families to the
    # generated HTML:
    #  -- "graphics_map" associates a file pattern (e.g., /det-radlum*) to
    #  a symbolic tag ("radlum").
    #  -- "graphics_show" gives a list of symbolic tags to put into HTML sections (in the
    #  order given), and the displayed section names.
    #
    # graphics_map: maps filename -> tag
    # list-of-pairs: (str,tag) means if filename contains 'str', it is a graphic of type 'tag'
    graphics_map = [('/det-perstar-det', 'perstar-det'),
                        ('/det-perstar-char', 'perstar-char'),
                        ('/det-radlum', 'radlum'),
                        ('/det-rad-sma', 'rad-sma'),
                        ('/det-duration', 'duration'),
                        ('/det-event-count', 'event-count'),
                        ('/det-visit-time', 'visit-time'),
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
    # tag: GD() means that graphics of type "tag" have the given description
    # Note: the order here gives the order of presentation!
    graphics_show = {
        'radlum':      GraphicsDescription('Radius/Luminosity',          'graphics'),
        'rad-sma':     GraphicsDescription('Radius/SMA',                 'graphics'),
        'duration':    GraphicsDescription('Event Duration',             'graphics'),
        'event-count': GraphicsDescription('Event Count',                'graphics'),
        'visit-time':  GraphicsDescription('Visits vs. Time',            'graphics'),
        'yield':       GraphicsDescription('Yield vs. Time',             'graphics', infotype='Information'),
        'cume':        GraphicsDescription('Mission Resources vs. Time', 'graphics'),
        'perstar-det': GraphicsDescription('Per-Star Detection',         'graphics', filename=WWW_DOC/'per-star-metrics.html'),
        'perstar-char':GraphicsDescription('Per-Star Characterization',  'graphics', filename=WWW_DOC/'per-star-metrics.html'),
        'promote':     GraphicsDescription('Target Promotion',           'graphics', infotype='Information'),
        'earth-char':  GraphicsDescription('Earth Characterizations',    'graphics'),
        'path':        GraphicsDescription('Full-Ensemble Path',         'path-ensemble')
        }
    # tables_map: maps filename -> tag
    # (str,tag) means if filename contains 'str', it is a table of type 'tag'
    # and will be shown in-line with the associated graphics in "graphics_show".
    # descriptive text is associated with the tag in "tables_show" below
    tables_map = [
        ('/table-funnel', 'promote'),
        ('/table-det-funnel', 'yield'), # with yield-vs-time
        ]
    # tag: TD() means the as-shown-in-html information for table of type 'tag' is TD
    tables_show = {
        'promote': TableDescription(
            'Target Promotion',
            '''These tables show the promotion-to-characterization process for
            the given star populations.''',
            WWW_DOC / 'promotion-tabulation.html'),
        'yield': TableDescription(
            'Detection Progress',
            '''These tables show detection attempts and their outcome (failure/success)
            across the given star populations.''',
            WWW_DOC / 'detection-tabulation.html')
        }

    def __init__(self, sim_dir, name=''):
        self.name = name if name else sim_dir
        self.Ndrm = 0
        self.graphics = defaultdict(list) # dict-of-list
        self.tables = defaultdict(list) # dict-of-list
        self.path_graphics = defaultdict(dict) # dict-of-dict
        self.csvs = []
        self.readme_info = '' # one-line summary
        self.readme_file = '' # filename

    def add(self, filename):
        r'''Dispatcher: adds filename to the correct category of item.'''
        # see also: self.sim_dir_no_index for skipped dirs that
        # do not reach this routine
        if '/path/' in filename:
            self.add_path_graphic(filename)
        elif filename.endswith('.png') or filename.endswith('.pdf'):
            self.add_graphic(filename)
        elif filename.endswith(('/README.md', '/README.html')):
            self.add_readme(filename)
        elif filename.endswith('.md'):
            self.add_table(filename)
        elif '/drm/' in filename:
            self.add_drm(filename)
        elif filename.endswith('.csv'):
            self.add_csv(filename)
        
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
        elif filename.endswith('obs-timeline-part2.png'):
            # part2 covers up to 10 year mission duration
            seed = re.sub(r'.*/([0-9]+)-obs-timeline-part2\.png', r'\1', filename)
            self.path_graphics[seed]['obs-timeline-part2'] = filename
        elif filename.endswith('star-obs-trace.png'):
            seed = re.sub(r'.*/([0-9]+)-star-obs-trace\.png', r'\1', filename)
            self.path_graphics[seed]['star-obs-trace'] = filename
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
        
    def add_readme(self, filename):
        # in the function call, need to zero-out the "root" argument
        # d is the Ensemble base directory, and d is the path down
        meta, meta_fn = obtain_metadata_info('', os.path.dirname(filename))
        self.readme_info = meta
        self.readme_file = meta_fn

    def render_table(self, hh, t_file):
        r'''Renders a table, in a markdown-format file, as html.

        Table is expected to be in the form:

        | colhead1 | colhead2 | ... | colheadN |
        | ---      | ---      | ...            |
        | data 1_1 | data 1_2 | ...            |
        ...
        | data 2_1 | data 2_2 | ...            |
        <newline>
        coda
        
        with an optional header (# Title) in front of the column heads.
        Markdown puts a |---|---| style line between the header rows 
        and the data rows (allowing for multi-line heads), and this 
        code requires such a separation. The ending blank line is required.
        The optional coda at the end is a text block entered as one
        HTML paragraph.
        Another approach: use the python markdown library.
        '''
        # slots for post-table markup (e.g. comments after the table)
        coda = []
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
                elif line.startswith('[//]'):
                    # markdown-format comment
                    pass
                elif in_table == 0:
                    # default case is text to appear after the table:
                    # it will appear as a single paragraph
                    coda.append(line)
                else:
                    # we are in-table, but line does not match table content
                    print('Unexpected table line <%s>' % line)
            # be sure to close the table out
            if in_table > 0:
                hh.table_end()            
            # insert the coda, if any, as one paragraph
            if coda:
                hh.paragraph('\n'.join(coda))

    def render(self, filename, uplink):
        r'''Main driver for rendering a directory to HTML files.'''
        # print('Rendering', self.name)
        # second argument to the class is the HTML doc title
        with HTML_helper(filename, 'Sim: ' + os.path.basename(self.name)) as hh:
            # nav-bar
            hh.nav_here()
            # title - h1 tag
            hh.header(os.path.basename(self.name), level=1)
            # navigation link
            hh.link('../../', 'Up to %s' % uplink)
            # fixed link to documentation
            hh.paragraph('Sandbox ' + hh.link(WWW_RES/'doc_sandbox/', 'documentation', inner=True))
            # summary
            hh.header('Ensemble Summary')
            hh.paragraph(f'Ensemble: {self.name}')
            hh.paragraph(f'Size: {self.Ndrm} sims')
            hh.paragraph('Source JSON ' + hh.link('../reduce-script.json', 'script', inner=True))
            # table of contents
            hh.toc_here('Contents')
            # overall images
            for tag, g_desc in self.graphics_show.items():
                target = g_desc.target
                # overall section header
                hh.header(f'{g_desc.title} {g_desc.infotype}')
                # optional caption
                if SECTION_HEADS[tag]:
                    hh.paragraph(SECTION_HEADS[tag])
                # tables, if any, are at section top (ad hoc at the moment)
                if tag in self.tables_show and len(self.tables[tag]) > 0:
                    # print(f'Rendering {self.tables[tag]}')
                    # "t_desc" abstracts the table explanation
                    t_desc = self.tables_show[tag]
                    hh.header(f'Tables for {t_desc.title}', level=3)
                    if t_desc.text:
                        hh.paragraph(t_desc.text, br=True)
                    if t_desc.filename:
                        hh.paragraph(hh.link(t_desc.filename,
                                    'Rules for counting', inner=True) +
                            ' in these tables',
                            br=True)
                    for table in self.tables[tag]:
                        self.render_table(hh, table)
                # plot widgets for star-target info
                if tag == 'perstar-det':
                    hh.header('Interactive Detection Plot Widget', level=3)
                    hh.div('<!-- det plot goes here -->', id='detPlotDiv', style='width: 900px; height: 700px;')
                    hh.div('Detection QOI for Plot Shading: <select class="det_qoi_select"> </select>')
                    if g_desc.filename:
                        hh.paragraph(hh.link(g_desc.filename, 'Plot description and specifics', inner=True), br=True)
                elif tag == 'perstar-char':
                    hh.header('Interactive Characterization Plot Widget', level=3)
                    hh.div('<!-- char plot goes here -->', id='charPlotDiv', style='width: 900px; height: 700px;')
                    hh.div('Characterization QOI for Plot Shading: <select class="char_qoi_select"> </select>')
                    if g_desc.filename:
                        hh.paragraph(hh.link(g_desc.filename, 'Plot description and specifics', inner=True), br=True)
                    # assume dets always accompany chars, and insert the script at this point
                    hh.script(WWW_RES/'star-target-plots.js')
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
                    hh.script(WWW_RES/'ens-path-plots.js')
                plots = self.graphics[tag]
                if not plots:
                    if target:
                        hh.paragraph(f'No such plots.  Generate with: make S=... {target}')
                else:
                    hh.table_top(('Filename', 'Plot Preview'), elem_class='gfx')
                    for p in sorted(plots):
                        img = hh.image(p, width=100, inner=True)
                        hh.table_row((os.path.basename(p), img))
                    hh.table_end()
            # single-run plots: path movies, keepout, obs-timelines
            hh.header('Single-Run Paths and Timelines')
            path_movie_seeds = sorted(self.path_graphics.keys())
            if not path_movie_seeds:
                hh.paragraph('No path movies or observation timelines yet.')
                hh.paragraph('Generate with: make S=... path-movie-N or make S=... obs-timeline-N\n')
            else:
                # (plot-count incorrect if there are coronagraph/starshade keepout plots.
                # they are dicts-of-dict underneath path_graphics[seed] - not worth it to fix)
                num_plot = sum(len(p) for p in self.path_graphics.values())
                hh.paragraph(f'Found {num_plot} plots across {len(path_movie_seeds)} runs.')
                hh.table_top(('Seed', 'Representative Figure'), elem_class='gfx')
                for seed in path_movie_seeds:
                    info = self.path_graphics[seed]
                    # find a thumbnail image for the seed, or None
                    thumb_tries = {'final': 100, 'star-obs-trace': 50, 'obs-timeline': 50}
                    img_target, img_width = None, 100
                    for thumb, width in thumb_tries.items():
                        if thumb in info:
                            img_target = info[thumb]
                            img_width = width
                            break
                    # make up and insert table row
                    alink = hh.link('path-%s.html' % seed, seed, inner=True)
                    img = hh.image(img_target, width=img_width, inner=True)
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
                hh.paragraph('Observing timeline plot ' +
                    hh.link(WWW_DOC/'obs-timeline.html', 'format description', inner=True),
                    br=True)
                img_target = info['obs-timeline']
                hh.image(img_target, width=70)
                if 'obs-timeline-part2' in info:
                    img_target = info['obs-timeline-part2']
                    hh.image(img_target, width=70)
            else:
                hh.paragraph('Timeline not available.  Generate with: make ... obs-timeline-N')
            # star-obs-trace
            hh.header('Star-Observation Trace')
            if 'star-obs-trace' in info:
                hh.paragraph('Star-observation trace plot ' +
                    hh.link(WWW_DOC/'star-obs-trace.html', 'format description', inner=True),
                    br=True)
                img_target = info['star-obs-trace']
                hh.image(img_target, width=70)
            else:
                hh.paragraph('Star-Observation trace not available.  Generate with: make ... obs-timeline-N')
            # path keepout
            hh.header('Keepout and Observations')
            if 'obs-keepout-all' in info:
                img_target = info['obs-keepout-all']
                hh.image(img_target, width=70)
                img_target = info['obs-keepout-char']
                hh.image(img_target, width=70)
            else:
                hh.paragraph('Keepout timeline not available.  Generate with: make ... keepout-N')
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
                hh.script(WWW_RES/'ens-path-plots.js')
            else:
                hh.paragraph('Data for widget not available.  Generate with: make ... path-movie-N')
            # path movie
            hh.header('Path Movie of This Tour')
            if 'movie' in info:
                hh.paragraph('Path movie ' +
                    hh.link(WWW_DOC/'path-movie.html', 'format description', inner=True),
                    br=True)
                hh.video(info['movie'])
            else:
                hh.paragraph('No path movie available.  Generate with: make ... path-movie-N')
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

    def put_file_counts(self):
        # make an object
        info = dict()
        info['index_path_count'] = sum([len(g) for _,g in self.graphics.items()])
        info['index_path_gfx'] = len(self.path_graphics)
        info['index_gfx_count'] = sum([len(g) for _,g in self.path_graphics.items()])
        #info['experiment'] = "s_c1_0.01_c3_0.99"
        info['readme_info'] = self.readme_info
        info['readme_file'] = self.readme_file
        # dump the object
        fn = 'index-files.json'
        with open(fn, 'w') as fp:
            if False:
                print('Dumping to: {fn}')
            json.dump(info, fp, indent=2)


############################################################
#
# One-line Sim Summaries
#
############################################################

def exp_summary(d):
    r'''Summarize one Experiment/Family directory as an OrderedDict of strings.

    If d is None: instead, return the header corresponding to the summary.

    '''
    
    # TODONT: It is *not worthwhile* to attempt to roll up graphics file counts
    # from ensembles in to the tables covering Experiments/Families. Note, 
    # these counts cannot be known at the time of reduction, only after all 
    # graphics have been added, or at HTML generation time. 
    # Because HTML generation only *sometimes* recurses into subdirectories, 
    # the rolled-up count is in general not known and kept up-to-date. So 
    # the rolled-up count can only sometimes be known at HTML-generation time. 
    # Further, it mostly provides value at the ensemble-summary level, and this 
    # is where the count is tractable to glob for. Let. It. Go.
    # TODO: (2025) The tabulator.js version has obsoleted this warning.
    # We now generate summarizer files (index-files.json, index-files-byname.json).
    # The dataflow needs to be documented

    # Load the reduction summary file, reduce-info.csv
    rv = reduce_info_summary(d)
    return rv

def sim_summary(d):
    r'''Summarize one simulation directory as an OrderedDict of strings.

    If d is None: instead, return the header corresponding to the summary.'''
    # Load the reduction summary file, reduce-info.csv
    rv = reduce_info_summary(d)
    # the quick-out for header info
    if d is None:
        return rv
    # Performance note: the globs below run for each non-Family/Experiment
    # subdir below sims/ (i.e., if d itself has DRM's) when indexing the 
    # Sandbox. (2023/10: About 2/3 of the runtime is running these globs.)
    # Rather than being clever, the best solution is to put the old Ensembles
    # into a Family at the top-level, cutting runtime of many indexing
    # components including this one.
    # 
    # number of unique paths that were indexed in any way (fn = SEED-foo-bar.png)
    path_count = len(set(fn.split('-')[0] for fn in glob.glob(os.path.join(d, 'path/[0-9]*-*.*'))))
    # number of path graphics
    path_gfx = (len(glob.glob(os.path.join(d, 'path/[0-9]*.png'))) +
                len(glob.glob(os.path.join(d, 'path/[0-9]*.mp4'))))
    # total graphic file counts (non-path)
    gfx_count  = (len(glob.glob(os.path.join(d, 'gfx/*.png'))) +
                  len(glob.glob(os.path.join(d, 'path-ens/*.png'))))
    # insert the extra data, just retrieved above
    rv['path_count'] = path_count
    rv['path_gfx']   = path_gfx
    rv['gfx_count']  = gfx_count
    return rv

def sim_filecount_glob(d, prefix=''):
    r'''Count files, mostly graphics, in one simulation directory 

    If d is None: instead, return the header corresponding to the summary.'''
    # return value template
    rv = {
        f'{prefix}path_count':0,
        f'{prefix}path_gfx':0,
        f'{prefix}gfx_count':0}
    # the quick-out for header info
    if d is None:
        return rv
    # quick-out: glob will not work on nested structures - return 0
    if d.endswith(('.exp', '.fam')):
        return rv
    # Performance note: the globs below run for each non-Family/Experiment
    # subdir below sims/ (i.e., if d itself has DRM's) when indexing the 
    # Sandbox. (2023/10: About 2/3 of the runtime is running these globs.)
    # Rather than being clever, the best solution is to put the old Ensembles
    # into a Family at the top-level, cutting runtime of many indexing
    # components including this one.
    # 
    # number of unique paths that were indexed in any way
    # typical name: SEED-X-Y.(png|mp4), e.g.:
    #   SEED-obs-keepout-char.png
    #   SEED-obs-timeline.png
    #   SEED-final.png
    #   SEED.mp4
    #path_count = len(set(fn.split('-')[0] for fn in glob.glob(os.path.join(d, 'path/[0-9]*-*.*'))))
    # number of path graphics
    path_gfx = 0
    path_set = set()
    path_gfx_re = re.compile(r'(\d+).*\.(png|mp4)$')
    p = os.path.join(d, 'path')
    if os.path.isdir(p):
        for entry in os.scandir(p):
            if not entry.is_file(): continue
            m = path_gfx_re.match(entry.name)
            if not m: continue
            path_gfx += 1
            path_set.add(m.group(1)) # capture the seed part
    #path_gfx = (len(glob.glob(os.path.join(d, 'path/[0-9]*.png'))) +
    #            len(glob.glob(os.path.join(d, 'path/[0-9]*.mp4'))))
    path_count = len(path_set)

    # total graphic file counts (non-path, and path-ens is non-path)
    gfx_count = 0
    # count {d}/gfx/*.png
    p = os.path.join(d, 'gfx')
    if os.path.isdir(p):
        for entry in os.scandir(p):
            if entry.is_file() and entry.name.endswith('.png'):
                gfx_count += 1
    # add in {d}/path-ens/*.png
    p = os.path.join(d, 'path-ens')
    if os.path.isdir(p):
        for entry in os.scandir(p):
            if entry.is_file() and entry.name.endswith('.png'):
                gfx_count += 1

    #gfx_count  = (len(glob.glob(os.path.join(d, 'gfx/*.png'))) +
    #              len(glob.glob(os.path.join(d, 'path-ens/*.png'))))

    # insert the extra data, just retrieved above
    rv[f'{prefix}path_count'] = path_count
    rv[f'{prefix}path_gfx']   = path_gfx
    rv[f'{prefix}gfx_count']  = gfx_count
    return rv

def sim_filecount(d, prefix=''):
    r'''Obtain file-count, mostly graphics, in one simulation directory 

    If d is None: instead, return the header corresponding to the summary.'''
    # return value template
    rv = {
        f'{prefix}path_count':0,
        f'{prefix}path_gfx':0,
        f'{prefix}gfx_count':0}
    # the quick-out for header info
    if d is None:
        return rv

    fn = os.path.join(d, 'index-files.json')
    if not os.path.isfile(fn):
        # print('* Sad path')
        # needed when indexing hasn't been re-done by this code
        rv = sim_filecount_glob(d, prefix=prefix)
    else:
        #print('* Happy path')
        with open(os.path.join(d, 'index-files.json'), 'r') as fp:
            rv = json.load(fp)

    # insert the extra data, just retrieved above
    #rv[f'{prefix}path_count'] = path_count
    #rv[f'{prefix}path_gfx']   = path_gfx
    #rv[f'{prefix}gfx_count']  = gfx_count

    return rv


def reduce_info_summary(d):
    r'''Summarize one simulation directory as an OrderedDict of strings.

    If the argument is None, instead, return the header corresponding to the summary.'''
    # do not display floats to 12 digits of precision
    def fmt_float(x):
        return ('%.3f' % float(x))
    def fmt_int(x):
        return ('%d' % int(x))
    def fmt_str(x):
        return x
    # properties we want to maintain: key, label, formatter
    props = [
        ('ensemble_size',          'Ens. Size',       fmt_str),
        ('simtime',                'Last Sim. Date',  fmt_str),
        ('runtime',                'Reduction Date',  fmt_str),
        ('user',                   'User',            fmt_str),
        ('detections_earth_unique','Earths (Det.)',   fmt_float),
        ('chars_earth_unique',     'Earths (Char.)',  fmt_float),
        ('chars_earth_strict',     'Earths (Strict)', fmt_float),
        ('gfx_count',              'Ens. Graphs',     fmt_int),
        ('path_count',             "Path Summ's",     fmt_int),
        ('path_gfx',               'Path Graphs',     fmt_int),
        ]
    # Return just header columns if None was input
    if d is None:
        rv = OrderedDict()
        for dtag, dname, _ in props: rv[dtag] = dname
        return rv
    # grab the metadata for d from the CSV written by reduce_drms.py
    info_fn = os.path.join(d, 'reduce-info.csv')
    unknown_value = '(n/a)' # ()'s help sort order
    try:
        with open(info_fn) as f:
            info_items = csv.DictReader(f);
            info = next(info_items) # it is a 1-line csv
        # make the return value
        rv = OrderedDict()
        for dtag, dname, dfmt in props:
            if dtag in info:
                rv[dtag] = dfmt(info[dtag])
            else:
                rv[dtag] = unknown_value
    except IOError:
        # order does not matter
        rv = {dtag: unknown_value for dtag, _, _ in props}
    return rv


############################################################
#
# Handlers for Indexing Sims
#  (ensembles/experiments/families)
#
############################################################

def index_ensemble(args, path_sim, uplink):
    r'''Make an index.html for one ensemble.'''
    if not os.path.isdir(os.path.join(path_sim, 'drm')):
        return # not an ensemble -- e.g., a recursive call into an arb. dir
    outdir = 'html'
    with ChangeDir(path_sim):
        # ensure the output directory
        if not os.path.isdir(outdir):
            os.makedirs(outdir, 0o775)
        # gather information from all files present
        sim_info = SimSummary('.', name=path_sim)
        for root, dirs, files in os.walk('.'):
            level = root.replace(path_sim, '').count(os.sep)
            # delete some "dirs" to avoid walking them
            dirs[:] = [d for d in dirs if d not in sim_info.sim_dir_no_index]
            # examine all files
            for fn in files:
                fn_full = os.path.join(root, fn)
                if fn.startswith('.'):
                    continue # metadata files often of form ._FOO
                # print('Processing: ', fn_full)
                sim_info.add(fn_full)
        # render information
        fn = os.path.join(outdir, 'index.html')
        # index.html
        sim_info.render(fn, uplink)
        # render_paths handles all the seeds
        sim_info.render_paths(os.path.join(outdir, 'path-%s.html'))
        # drop a summary of file counts
        sim_info.put_file_counts()

def obtain_metadata_info(root, d):
    '''Find a one-line information block for a group, if it exists'''
    def meta_from_text(fp):
        '''Grab first nontrivial line as the "title"'''
        for l in fp.readlines():
            if l.strip():
                # final strip() removes whitespace after #
                meta = l.strip().replace('#', '').strip()
                return meta
        return ''
    def meta_from_html(fp):
        '''Extract the <title> element'''
        # (we don't want dependency on a parser)
        meta = ''
        for l in fp.readlines():
            if '<title>' in l:
                # assumes <title> is a one-liner
                lplus = ' '.join(l.split())
                # get what is between <title> and the next <
                meta = lplus.split('<title>')[1].split('<')[0]
                if not meta:
                    meta = 'HTML Title not understood'
            if '<body>' in l or '</head>' in l:
                # give up
                meta = 'HTML Title not found'
            if meta:
                return meta
        return 'HTML Title not found'
    # map from FN -> function, where
    #   FN = allowed filename (basename part)
    #   function = function that returns one-liner from file object
    meta_map = {
        'README.md': meta_from_text,
        'README.html': meta_from_html,
        }
    meta, meta_fn = '', ''
    # try allowed filenames
    for fn, mapper in meta_map.items():
        fn_full = os.path.join(root, d, fn)
        if os.path.isfile(fn_full):
            with open(fn_full) as fp:
                meta = mapper(fp) or f'Title not found in {fn_full}'
                meta_fn = os.path.join(d, fn)
            break
    return meta, meta_fn

def index_group(args, startpath, title, uplink):
    r'''Make an index.html with a table summarizing all ensembles below startpath.

    This function corresponds to Experiments and Families: groups of sims.
    If args.recurse, descend recursively and summarize all child directories
    as well. This may imply a recursive call to this routine.'''

    # optionally descend into startpath/*/
    # (-r => args.recurse, not usually done, make html-all is only current target)
    # TODO: perhaps put at top of function, so that the parents can use
    #   indexing results from children
    # TODO: lightly tested
    # TODO: will call index_one_sim on logdirs, etc.
    if args.recurse:
        if not os.path.isdir(startpath):
            print(f'{args.progname}: Given path {startpath} not present. Skipping.')
            return
        for entry in os.scandir(startpath):
            if entry.name.startswith('.'): continue
            if entry.is_dir():
                d = os.path.join(startpath, entry.name)
                # print(f'I_G:   recursive call into {d = }')
                index_one_sim(args, d)

    # output HTML file
    filename = os.path.join(startpath, 'index.html')
    with HTML_helper(filename, title) as hh:
        # title as H1
        hh.header(title, level=1)
        # navigation link
        if uplink:
            hh.link('../', 'Up to %s' % uplink)
        hh.paragraph('Sandbox ' + hh.link(WWW_RES/'doc_sandbox/', 'documentation', inner=True))
        # overall section
        hh.header('Ensembles')
        # summary
        hh.paragraph(f'Ensemble set: {startpath}')
        # link to the README.md/README.html, if present
        # (typically placed by reduce_drm_sets)
        for fn_stem in ('README.md', 'README.html'):
            fn = os.path.join(startpath, fn_stem)
            if os.path.isfile(fn):
                hh.paragraph('Experiment descriptive ' + hh.link(fn_stem, 'README', inner=True))
                break
        # table of individual sims
        vanilla_tables = False
        if vanilla_tables:
            # make the table be sortable so that the JS sorter knows about it
            hh.table_top(['Name'] + list(sim_summary(None).values()), elem_class='sortable')
            item_num = 0
            # FIXME: os.walk misleads (we only go down one level) - use os.listdir()
            for root, dirs, files in os.walk(startpath):
                # loop over dirs (i.e., sims enclosed by startpath)
                for d in sorted(dirs):
                    #print('I_G:   looping on {}'.format(d))
                    # 1: make table row with: link to the sim html + sim summary
                    #     no recursive descent in this block
                    if d.endswith('.exp') or d.endswith('.fam'):
                        # alink is an HTML fragment as a string
                        alink = hh.link(f'{d}/index.html', d, inner=True)
                        # augment alink (the row "name" field) with README.txt info ("extra")
                        extra, extra_fn = obtain_metadata_info(root, d)
                        if extra:
                            # tooltip HTML element classes are picked up by CSS styles (no JS needed)
                            dingbat = '&#9432;' # currently a circled "i"
                            extra_span = hh.span(html.escape(extra), **{'class': 'tooltiptext'})
                            alink += hh.link(extra_fn, f'&nbsp;{dingbat}{extra_span}',
                                             inner=True, **{'class': 'tooltip'})
                        properties = exp_summary(os.path.join(root, d))
                        hh.table_row([alink] + list(properties.values()))
                        item_num += 1
                    elif os.path.isdir(os.path.join(root, d, 'drm')):
                        alink = hh.link(f'{d}/html/index.html', d, inner=True)
                        properties = sim_summary(os.path.join(root, d))
                        # format the properties as a row
                        hh.table_row([alink] + list(properties.values()))
                        item_num += 1
                    else:
                        pass # dir name not recognized, skip it
                # finished all dirs below startpath -- do not descend further with os.walk()
                break
            # summary over the whole set of ensembles in the table (root/reduce-info.csv)
            properties = exp_summary(startpath)
            #properties = sim_summary(None)
            hh.table_foot()
            hh.table_row(['<b>SUMMARY</b> (%d items)' % item_num] + list(properties.values()))
            hh.table_end()
        else:
            hh.paragraph('Tabulated summary is below.', br=True)
            hh.div('<!-- tabulation goes here -->', id="scenario-table")
            # (presently, filenames are hard-coded, which is OK)
            ## # set up path to the data CSVs, and label mode, for the JS viewer
            ## hh.script('var ens_path_root = "../path-ens";', literal=True)
            ## hh.script([
            ##    'var ens_tab_root = "../path-ens";',
            ##    'var ens_tab_mode = "ensemble";'],
            ##          literal=True)
            hh.script(WWW_RES/'ensemble-tabulator.js')

        hh.paragraph('In the summary, ensemble size is cumulative.', br=True)
        hh.text('Simulation date reflects the most recent simulation run below this level.', br=True)
        hh.text('Reduction date and user reflect the most recent reduction below this level.', br=True)
        hh.text('Yields reflect the maximum over all ensembles below this level.', br=True)
        hh.text('Yield definitions:', br=True)
        hh.text('&nbsp;Earths (Det.) = Number of successful Earth detections, repeat visits not counted.', br=True)
        hh.text('&nbsp;Earths (Char.) = Number of successful Earth characterizations (any spectral band, status = &plusmn;1), repeat visits not counted.', br=True)
        hh.text('&nbsp;Earths (Strict) = Number of successful Earth characterizations (all spectral bands have status = +1), repeat visits not counted.', br=True)
        
    # make a summary of graphical files
    if not vanilla_tables:
        prop_all = []
        # TODO: sub os.scandir() for efficiency
        for entry in os.listdir(startpath):
            d = os.path.join(startpath, entry)
            if not os.path.isdir(d): continue
            if d.endswith(('.exp', '.fam')) or os.path.isdir(os.path.join(d, 'drm')):
                properties = sim_filecount(d, prefix='index_')
                properties['experiment'] = entry
                prop_all.append(properties)
        out_fn = os.path.join(startpath, 'index-files-byname.json')
        with open(out_fn, 'w') as fp:
            if args.DIAGNOSE:
                print('Dumping index to: %s' % out_fn)
            json.dump(prop_all, fp, indent=2)
        out_fn = os.path.join(startpath, 'index-files.json')
        with open(out_fn, 'w') as fp:
            if args.DIAGNOSE:
                print('Dumping index to: %s' % out_fn)
            rm_info, rm_file = obtain_metadata_info(startpath, '')
            prop_total = property_total(prop_all)
            prop_total['readme_info'] = rm_info
            prop_total['readme_file'] = rm_file
            json.dump(prop_total, fp, indent=2)
    return

def property_total(properties):
    skipped_keys = set(['experiment', 'readme_info', 'readme_file'])
    if not properties:
        return {}
    p_tot = defaultdict(int)
    for p1 in properties:
        for key in p1:
            if key in skipped_keys: continue
            p_tot[key] += p1[key]
    return dict(**p_tot)

############################################################
#
# Delegation of sims to handlers
#
############################################################

def index_one_sim(args, sim):
    r'''Index one simulation, with the specific action depending on its type.

    Note: the input argument is the simulation *directory* (with sims/ prefix).'''
    # remove trailing / if any (currently only happens with sim='sims/')
    sim = sim[:-1] if sim.endswith('/') else sim
    print('Indexing:', sim)
    # establish parent sim
    try:
        sim_up = sim.split(os.sep)[-2]
        sim_up_path = os.sep.join(sim.split(os.sep)[:-1])
    except IndexError:
        # sim = 'sims' -- these will be unused
        sim_up = ''
        sim_up_path = ''
    # text for HTML link to parent sim
    if sim_up.endswith('.exp'):
        uplink = 'Experiment'
    elif sim_up.endswith('.fam'):
        uplink = 'enclosing Family'
    elif sim_up == 'sims':
        uplink = 'Sandbox Root'
    else:
        uplink = 'Unknown'
    uplink_label = '{} ({})'.format(uplink, sim_up_path)
    # switch the index type depending on sim name
    if sim == 'sims':
        index_group(args, sim, 'Simulation Ensemble Root', '')
    elif sim.endswith('.exp'):
        index_group(args, sim, 'Experiment Root: ' + sim.split('/')[-1], uplink_label)
    elif sim.endswith('.fam'):
        index_group(args, sim, 'Family: '          + sim.split('/')[-1], uplink_label)
    else:
        index_ensemble(args, sim, uplink_label)


############################################################
#
# Main routine
#
############################################################

def sims_and_parents(sim_orig):
    '''Return a list of all parent paths of a list of input paths sim_orig.

    When indexing is turned on (-i option), then this code should index the given
    argument sims (sim_orig) as well as all parent sim directories. 
    That is, if sim_orig = ['test.fam/HabEx'], then test.fam/HabEx should be indexed,
    and then test.fam, and finally sims. 
    This code returns a list of all such parent sims, for a list of sims.
    The ordering in the returned list is by depth, so that the indexing 
    will proceed from leaves to root.'''

    def parent_sims(s):
        '''Return a list of (depth, path) for all parent paths of a single input path s.'''
        # normpath turns '' into '.', undesirable here
        subdirs = os.path.normpath(s).split(os.sep) if s else ''
        return [(depth+1, os.path.join(*subdirs[:depth+1])) for depth in range(len(subdirs))]

    # get a dict mapping sim_directories -> depth
    #   if sim_orig  = ['t.fam/f1.fam/f2.fam/HabEx']
    #   then sim_index = {
    #    't.fam': 1,
    #    't.fam/f1.fam': 2,
    #    't.fam/f1.fam/f2.fam': 3,
    #    't.fam/f1.fam/f2.fam/HabEx': 4}
    sim_index = dict()
    for s in sim_orig:
        for depth, s in parent_sims(s):
            sim_index[s] = depth # collisions are of same depth
    depth_max = max([0] + list(sim_index.values()))
    # compile a leaves-to-root list of sims
    sim_new = []
    # below loop runs depth_max...1, but not 0
    for d in range(depth_max, 0, -1):
        # extract and append all sims at depth == d
        sim_adds = [k for k in sim_index if sim_index[k] == d]
        sim_new.extend(sim_adds)
    # add the root
    sim_new.append('')
    return sim_new

def main(args):
    r'''Main routine: Summarize files below a given location to HTML.

    In our usage, a "sim" is the name of a simulation scenario, whether
    an ensemble, an Experiment, or a Family.
    Algorithm: If -i was supplied, add the chain of parents of each given SIM
    to the list of sims_to_visit. The ordering is such that branches/leaves 
    are visited first, and then the path from the branch/leaf up to the root.
    Conceptually -i augments the path upward from SIM to root.
    If -r was supplied, we also recurse downward from each supplied SIM to
    all child SIM's. Conceptually, -r augments the path downward from SIM to
    branches and leaves below SIM.
    Note that the SIM's added by giving the [-i] option are NOT recursed into
    due to the [-r] option. Only the original sims supplied on the command line
    are recursed into due to [-r].'''

    # a global for debugging issues deep within this code
    global DIAGNOSE
    DIAGNOSE = args.DIAGNOSE
    # sims_to_visit = args.sims + parents if indexing is on
    if args.index:
        sims_to_visit = sims_and_parents(args.sim)
    else:
        sims_to_visit = args.sim
    #print('Visiting: ' + ','.join(sims_to_visit))
    # copy the as-supplied argument
    can_recurse = args.recurse
    # visit all the desired sims, leaves-to-root
    for sim in sims_to_visit:
        # don't recurse beneath parent dirs that we're just indexing
        # (the below line re-uses the convenient god object, args.recurse)
        args.recurse = can_recurse and (sim in args.sim)
        # generate the desired index
        index_one_sim(args, os.path.join('sims', sim))

    
# not all options are currently used
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Summarize EXOSIMS sims to HTML.",
                                     epilog='')
    # default argument is the empty sim, a/k/a the root
    parser.add_argument('sim', metavar='SIM', nargs='*', default=[''],
                            help='''sim ensemble directory, no sims/ prefix (default = '' corresponds to sandbox root)''')
    parser.add_argument('-i', '--index', help='also generate index(es) within all sims enclosing each SIM',
                      dest='index', action='store_true', default=False)
    parser.add_argument('-r', '--recurse', help='recursive over all sims below each named SIM',
                      dest='recurse', action='store_true', default=False)
    parser.add_argument('-D', '--diagnose', help='diagnostic (debugging) output',
                      dest='DIAGNOSE', action='count', default=0)
    args = parser.parse_args()
    
    # measure elapsed time
    time_start = time.time()

    # set umask in hopes that files/dirs will be group-writable
    os.umask(0o002)

    # program name, for convenience
    args.progname = os.path.basename(sys.argv[0])
    
    # Announce updated version
    print(f'{args.progname}: tabulator.js summarizer version.')

    main(args)

    # report elapsed time
    time_delta = time.time() - time_start
    print(f'{args.progname}: Done. Elapsed time: {time_delta:.3f}s.')

    sys.exit(0)



