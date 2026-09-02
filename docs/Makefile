# Minimal makefile for Sphinx documentation.
#
# `make html` builds into _build/html. The flags mirror what Read the Docs runs
# -- see .readthedocs.yaml, whose `fail_on_warning: true` becomes `-W
# --keep-going` there -- so a warning that builds here builds there too.

SPHINXOPTS    ?= -W --keep-going -j auto
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = .
BUILDDIR      = _build

.PHONY: help html linkcheck clean Makefile

help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

clean:
	rm -rf "$(BUILDDIR)"

%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
