# 
# Makefile to preserve a few usage lines
#
# Other lines of interest:
#   mkdocs build -- force a re-build
#   mkdocs serve -- monitors this dir for changes, rebuilds, 
#                   if files change, and chatters to stdout/stderr
# 

# venv containing some small Python dependencies
venv:
	@printf "Run this: \n  source ../../../Python-venvs/sandbox-docs/bin/activate\n"

preflight:
	pre_mkdocs.py

# one way to tell "mkdocs serve" to rebuild
rebuild:
	touch mkdocs.yml

export:
	rsync -av --delete --exclude=README ./site/ ../../Local/www-resources/doc_sandbox/


