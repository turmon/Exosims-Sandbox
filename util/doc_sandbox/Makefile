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

# remote holding the github pages branch; the JPL mirror is a separate remote,
# e.g. "make publish GH_REMOTE=jpl"
GH_REMOTE ?= origin

.PHONY: doc default rebuild export publish

default: doc

# typical entry point
doc:
	uv run mkdocs build

# one way to tell "mkdocs serve" to rebuild
rebuild:
	touch mkdocs.yml

export:
	rsync -q -av --delete --exclude=README.txt ./site/ ../../Local/www-resources/doc_sandbox/

# push the built site to the gh-pages branch of $(GH_REMOTE).
# NOTE: this pushes. Run it deliberately.
publish:
	uv run mkdocs gh-deploy --remote-name $(GH_REMOTE)
