#
# Makefile for documentation for plots
#
# Note: requires "uv" for the python markdown package.
# 
# [1] To update documentation, edit the one corresponding *.md
# file here, and "make install".
#
# [2] To add a new documentation element, create a markdown file
# "plot.md" in this directory, add it to SOURCES below,
# and "make install"
# 
# If your documentation needs supplemental images, place
# them in "Media" below.
#
# The install can also be run from the main Makefile in the
# Sandbox root.
#
# turmon 11/2023

# accept bash-only to get proper failure reporting 
# for the python | sed line below
SHELL=/bin/bash -o pipefail

INDEX_HEAD=index-plot-descriptions-head.md
INDEX=index-plot-descriptions.md

# source files - markdown
SOURCES=obs-timeline.md star-obs-trace.md path-movie.md promotion-tabulation.md detection-tabulation.md per-star-metrics.md $(INDEX)

# grep '^Title: ' *.md | sed 's/Title: //' | awk -F : '{print "[" $2 "]" $1}'

# destination directory for above pages
DEST_DIR=../www-resources/doc

# final-destination webpages
PAGES=$(SOURCES:.md=.html)

####

# intermediate files
MARKED_DOWN=$(SOURCES:.md=.md_html)
MARKED_META=$(SOURCES:.md=.md_html.meta)

all: $(PAGES)

index: $(INDEX)

$(INDEX): $(INDEX_HEAD)
	@ echo "Remaking index file: $@"
	cp $(INDEX_HEAD) $(INDEX)
	util/compile-plot-index.sh >> $(INDEX)

# httpd has FollowSymlinks off, so we can't symlink the CSS
# simply rsync's the HTML and Media/* files to $DEST_DIR
install: all
	@ echo "Installing into $(DEST_DIR)"
	@ mkdir -p $(DEST_DIR)/Media
	cd $(DEST_DIR) && cp ../ensemble-reports.css ./
	rsync -q -av $(PAGES) $(DEST_DIR)/
	rsync -q -av Media $(DEST_DIR)

# render plot.md to plot.md_html
# the resulting file does not have an HTML header or footer
# the meta-extract.sh script simply extracts the Title:
# (and any other) field from the header of the plot.md file,
# for later use in generating the HTML title
# local UV cachedir because ./ is much more efficient than ~/
%.md_html: %.md
	UV_CACHE_DIR=./.uv_cache uv run --with markdown python -m markdown -x extra -x meta $< | sed 's/^/    /' > $@
	util/meta-extract.sh < $< > $@.meta

# For m4:
#   HTML_CONTENT is the plot.md_html filename from python markdown  
#                (it is just #included)
#   HTML-macros.m4 is some basic macros for document elements
#   $<.meta is an m4-readable file produced by meta-extract.sh
#   HTML-page-template.m4 is the header, content, footer template
#               that is macro-expanded
%.html: %.md_html
	m4 -D HTML_CONTENT=$< HTML-macros.m4 $<.meta HTML-page-template.m4 > $@

clean:
	rm -f $(MARKED_DOWN)
	rm -f $(MARKED_META)
	rm -f $(INDEX)

.PHONY: all install clean index

