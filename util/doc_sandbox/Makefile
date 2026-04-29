# 
# Makefile to preserve a few usage lines
#
# Other lines of interest:
#   uv run mkdocs build -- force a re-build
#   uv run mkdocs serve -- monitors this dir for changes, rebuilds, 
#                   if files change, and chatters to stdout/stderr
# 

# discard builtin rules
MAKEFLAGS += -r

PREFLIGHT_SENTINEL=docs/implementation/drm_tabulate.py

.PHONY: doc preflight default rebuild export

default: doc

preflight: $(PREFLIGHT_SENTINEL)

# this preflight routine makes the sentinel
# as well as some other files
$(PREFLIGHT_SENTINEL): ../drm_tabulate.py
	pre_mkdocs.py

# typical entry point
doc: preflight
	uv run mkdocs build

# one way to tell "mkdocs serve" to rebuild
rebuild:
	touch mkdocs.yml

export:
	rsync -av --delete --exclude=README ./site/ ../../Local/www-resources/doc_sandbox/


