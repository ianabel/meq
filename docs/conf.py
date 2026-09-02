# Configuration file for the Sphinx documentation builder.
#
# Build locally with:
#     python3 -m venv /tmp/docsvenv
#     /tmp/docsvenv/bin/pip install -r docs/requirements.txt
#     /tmp/docsvenv/bin/sphinx-build -W -j auto -b html docs docs/_build/html
#
# -W and -j auto both mirror what Read the Docs runs (it adds -W --keep-going
# for .readthedocs.yaml's fail_on_warning), so a warning here is a failed RTD
# build. docs/Makefile passes the same flags.
#
# There is deliberately no autodoc here, and no Doxygen/Breathe. meq's headers
# carry ordinary prose comments rather than Doxygen markup -- often long ones
# explaining a sign convention or an MFEM contract -- so Breathe would emit bare
# signatures with the reasoning stripped out. The C++ reference under api/ is
# therefore written by hand against the headers, using Sphinx's own C++ domain,
# which also lets the prose cross-reference individual methods.

project = "meq"
copyright = "2019-2026 Ian G. Abel; see the COPYRIGHT file"
author = "Ian G. Abel"

# meq carries no version number anywhere in the tree -- there is no VERSION
# file and no version macro -- so there is nothing to read one from. Leave it
# unset rather than invent a number that would immediately be wrong.
release = ""

extensions = [
    "sphinx.ext.mathjax",
    "sphinx.ext.intersphinx",
    "sphinxcontrib.bibtex",
]

# sphinx_material.get_html_context() returns a dict containing functions, which
# Sphinx cannot pickle into its environment cache, so it warns once per build:
#   cannot cache unpickable configuration value: 'html_context'
# It is harmless -- the context is rebuilt every run anyway -- but it is fatal
# under -W, and these docs are built with -W so that a broken cross-reference
# fails rather than scrolls past. Suppress just that one category.
suppress_warnings = ["config.cache"]

templates_path = ["_templates"]
exclude_patterns = ["_build", "manual", "requirements.txt", "README.md"]

# The default role, so `like this` renders as literal rather than as a broken
# cross-reference. Most inline markup in these pages is a configuration key, a
# file name or a C++ identifier.
default_role = "literal"

# The primary domain is C++: meq is a C++ library, so a bare `.. function::` or
# a :any: lookup should mean the C++ one.
primary_domain = "cpp"
highlight_language = "cpp"

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable", None),
}

# -- Bibliography ------------------------------------------------------------
#
# One .bib for the whole tree, rendered on the bibliography page. The papers
# are the reason meq is written the way it is, so citations appear throughout
# the mathematical chapters rather than being collected in a reading list.
#
# `refs/Refs.md` in the source tree is the same list with the DOIs to fetch each
# PDF by; docs/references.bib is generated from it by hand and the two are meant
# to agree.
bibtex_bibfiles = ["references.bib"]
bibtex_default_style = "unsrt"
bibtex_reference_style = "author_year"

# -- HTML output -------------------------------------------------------------

import sphinx_material  # noqa: E402

html_theme = "sphinx_material"
html_theme_path = sphinx_material.html_theme_path()
html_context = sphinx_material.get_html_context()

html_theme_options = {
    "nav_title": "meq",
    "color_primary": "deep-purple",
    "color_accent": "amber",
    "repo_url": "https://github.com/ianabel/meq",
    "repo_name": "meq",
    "repo_type": "github",
    "globaltoc_depth": 2,
    "globaltoc_collapse": True,
    "globaltoc_includehidden": True,
    "master_doc": False,
}

# sphinx-material needs these four sidebars; it ships the templates and renders
# an empty navigation column without them.
html_sidebars = {
    "**": ["logo-text.html", "globaltoc.html", "localtoc.html", "searchbox.html"]
}

html_static_path = ["_static"]
html_css_files = ["meq.css"]
html_title = "meq"
html_short_title = "meq"

# Where this is published, so generated pages carry a canonical link rather
# than letting the per-version URLs compete with each other in search results.
html_baseurl = "https://meq.readthedocs.io/en/latest/"

# -- MathJax -----------------------------------------------------------------
#
# meq's equations are full of the same few symbols: the poloidal flux, the
# barred gradient that is NOT the cylindrical gradient (see the discretisation
# chapter -- the distinction is the whole content of the 1/r and r weights),
# and the Grad-Shafranov operator. Define them once.
mathjax3_config = {
    "tex": {
        "macros": {
            "gradbar": r"\bar{\nabla}",
            "dstar": r"\Delta^{\!*}",
            "psiax": r"\psi_{\mathrm{ax}}",
            "psibnd": r"\psi_{\mathrm{bnd}}",
            "Th": r"\mathcal{T}_h",
        }
    }
}
