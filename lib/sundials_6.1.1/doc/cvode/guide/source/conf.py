# ------------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2022, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ------------------------------------------------------------------------------

import sys, os
sys.path.append(os.path.dirname(os.path.abspath('../../../shared/versions.py')))
from versions import *

# -- Create new object types --------------------------------------------------

from sphinx.application import Sphinx

def setup(app: Sphinx):
    app.add_object_type('cmakeoption', 'cmakeop', 'single: %s (CMake option)')
    app.add_config_value('package_name', '', 'env', types=[str])

# -- General configuration ----------------------------------------------------

# Set variable used to determine which package documentation this is
# Can be one of 'arkode', 'cvode', 'cvodes', 'ida', 'idas', 'kinsol' or 'super'
package_name = 'cvode'

# If your documentation needs a minimal Sphinx version, state it here.
needs_sphinx = '4.0'

# Add any Sphinx extension module names here, as strings. They can be extensions
# coming with Sphinx (named 'sphinx.ext.*') or your custom ones.
extensions = ['sphinx_rtd_theme', 'sphinx.ext.ifconfig', 'sphinx.ext.mathjax',
              'sphinxfortran.fortran_domain', 'sphinxcontrib.bibtex', 'sphinx_copybutton']

# References
bibtex_bibfiles = ['../../../shared/sundials.bib']

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix of source filenames.
source_suffix = '.rst'

# The encoding of source files.
#source_encoding = 'utf-8-sig'

# The master toctree document.
master_doc = 'index'

# General information about the project.
project = u'User Documentation for CVODE'
copyright = u'2012-2021, Alan C. Hindmarsh, Radu Serban, Cody J. Balos, David J. Gardner, Daniel R. Reynolds, and Carol S. Woodward. Release number UCRL-SM-208108'

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
#
# The short X.Y version.
version = '{cvode_version}'.format(cvode_version=cvode_version)
sun_version = '{sundials_version}'.format(sundials_version=sundials_version)

# Set the date format (full-month-name day, full-year)
today_fmt = '%B %d, %Y'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns = []

# The reST default role (used for this markup: `text`) to use for all documents.
#default_role = None

# If true, '()' will be appended to :func: etc. cross-reference text.
#add_function_parentheses = True

# If true, the current module name will be prepended to all description
# unit titles (such as .. function::).
#add_module_names = True

# If true, sectionauthor and moduleauthor directives will be shown in the
# output. They are ignored by default.
#show_authors = False

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'
highlight_language = "c"

# A list of ignored prefixes for module index sorting.
#modindex_common_prefix = []

# Number figures, tables, and code blocks (can reference by number with numref)
numfig = True

# Override format strings that numref/numfig uses
numfig_format = {
  'section': '§%s'
}

rst_epilog = """
.. |YEAR| replace:: 2021
.. |CVODE_VERSION| replace:: {cvode_version}
.. |CVODES_VERSION| replace:: {cvodes_version}
.. |ARKODE_VERSION| replace:: {arkode_version}
.. |IDA_VERSION| replace:: {ida_version}
.. |IDAS_VERSION| replace:: {idas_version}
.. |KINSOL_VERSION| replace:: {kinsol_version}
""".format(
cvode_version = cvode_version,
cvodes_version = cvodes_version,
arkode_version = arkode_version,
ida_version = ida_version,
idas_version = idas_version,
kinsol_version = kinsol_version
)

# -- Options for HTML output ---------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
html_theme = 'sphinx_rtd_theme'

# The name for this set of Sphinx documents.  If None, it defaults to
# "<project> v<release> documentation".
#html_title = None

# A shorter title for the navigation bar.  Default is the same as html_title.
#html_short_title = None

# The name of an image file (relative to this directory) to place at the top
# of the sidebar.
html_logo = 'figs/sundials_logo_blue.png'

# The name of an image file (within the static path) to use as favicon of the
# docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
#html_favicon = None

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['../../../shared/_static']

# These paths are either relative to html_static_path
# or fully qualified paths (eg. https://...)
html_css_files = [
  'css/custom.css'
]

# If not '', a 'Last updated on:' timestamp is inserted at every page bottom,
# using the given strftime format.
html_last_updated_fmt = '%b %d, %Y'

# If true, SmartyPants will be used to convert quotes and dashes to
# typographically correct entities.
#html_use_smartypants = True

# Custom sidebar templates, maps document names to template names.
#html_sidebars = {}

# Additional templates that should be rendered to pages, maps page names to
# template names.
#html_additional_pages = {}

# If false, no module index is generated.
#html_domain_indices = True

# If false, no index is generated.
#html_use_index = True

# If true, the index is split into individual pages for each letter.
#html_split_index = False

# If true, links to the reST sources are added to the pages.
html_show_sourcelink = False

# If true, "Created using Sphinx" is shown in the HTML footer. Default is True.
#html_show_sphinx = True

# If true, "(C) Copyright ..." is shown in the HTML footer. Default is True.
#html_show_copyright = True

# If true, an OpenSearch description file will be output, and all pages will
# contain a <link> tag referring to it.  The value of this option must be the
# base URL from which the finished HTML is served.
#html_use_opensearch = ''

# This is the file name suffix for HTML files (e.g. ".xhtml").
#html_file_suffix = None

# Output file base name for HTML help builder.
htmlhelp_basename = 'CVODEdoc'


# -- Options for LaTeX output --------------------------------------------------

# 1. the rst file name used to generate the LaTeX file
# 2. the name of the LaTeX file to generate (and the resulting PDF file name)
# 3. the document title
# 4. text for \author
# 5. the LaTeX theme
# 6. include the file from 1. in the output
tex_author = r'''
Alan C. Hindmarsh$^1$,
Radu Serban$^1$,
Cody J. Balos$^1$,
\\
David J. Gardner$^1$,
Daniel R. Reynolds$^2$, and
Carol S. Woodward$^1$
\\\\
{\em $^1$Center for Applied Scientific Computing, Lawrence Livermore National Laboratory}\\
{\em $^2$Department of Mathematics, Southern Methodist University}
'''

latex_documents = [('index', 'cv_guide.tex', project,
                    tex_author, 'manual', False)]

# make sure the doc logo gets copied to the build directory
latex_logo = 'figs/doc_logo_blue.pdf'

# LaTeX customizations
latex_elements = {
# paper size option of the document class
'papersize': 'letterpaper',
# font size option of the document class
'pointsize': '10pt',
# set the version number/release name
'releasename' : version,
# arguments to the sphinxsetup macro
'sphinxsetup':
    # the color for titles
    'TitleColor={RGB}{0,0,0},'+
    # disable frames around code-blocks
    'verbatimwithframe=false,'+
    # do not wrap long lines in code-blocks
    'verbatimwrapslines=false,'+
    # background color for code-blocks
    'VerbatimColor={RGB}{240.0,240.0,240.0},'+
    # font used by heading
    'HeaderFamily=\\rmfamily\\bfseries,' +
    # line breaks are allowed inside inline literals
    'inlineliteralwraps=true',
# disable the fncychap package
'fncychap':'',
# figure alignment options
'figure_align': 'htbp',
# additional preamble content
'preamble': r'''
% =====================================================
% Start custom preamble (see latex_elements in conf.py)
% =====================================================

% Use ragged-right for the whole document
% Commented out after adding hyphenat below to get better spacing
% \usepackage[document]{ragged2e}

% Hyphenate verbatim text e.g., text in double backticks, :c:func:, etc.
\usepackage[htt]{hyphenat}

% Specify depths for section numbering and table of contents
\setcounter{tocdepth}{1}
\setcounter{secnumdepth}{3}

% Link a footnote to its location in the text
\usepackage{footnotebackref}

% Replace macro to use same header font and color for tables
\makeatletter
\protected\def\sphinxstyletheadfamily{\py@HeaderFamily \py@TitleColor}
\makeatother

% Replace macro to use same font in continued tables
\makeatletter
\protected\def\sphinxtablecontinued#1{\textrm{#1}}
\makeatother

% Add new command for SUNDIALS version
'''
+
r'\newcommand{\sunreleasename}{' + sun_version + r'}'
+
r'''

% ===================================================
% End custom preamble (see latex_elements in conf.py)
% ===================================================
''',
# extra class options
'extraclassoptions': 'twoside,openright',
# custom maketitle
'maketitle': r'''
% =======================================================
% Start custom cover page (see latex_elements in conf.py)
% =======================================================

\makeatletter

% Start roman numbering
\pagenumbering{Roman}

% Title page
\begin{titlepage}
  \newpage
  \null
  \vskip 2em%
  \begin{center}%
    \let \footnote \thanks
    {\huge \rmfamily \@title \space \releasename \par}
    {\Large \rmfamily SUNDIALS \space \sunreleasename \par}
    \vskip 3.0em%
    {\large \lineskip .5em%
     \begin{tabular}[t]{c}%
       \@author
     \end{tabular}\par}%
    \vskip 1em%
    {\large \@date \par}%
    \vfill
    {\includegraphics[width=0.5\textwidth]{doc_logo_blue}}
    \vfill
    {\large \rmfamily UCRL-SM-208108}
    \vfill
  \end{center}
  \par
  \vskip 1.5em
\end{titlepage}

\makeatother

\clearpage

% Disclaimer
\thispagestyle{empty}% no number of this page
\vglue5\baselineskip
\begin{center}
  {\bf DISCLAIMER}
\end{center}
\noindent
This document was prepared as an account of work sponsored by an agency of
the United States government. Neither the United States government nor
Lawrence Livermore National Security, LLC, nor any of their employees makes
any warranty, expressed or implied, or assumes any legal liability or responsibility
for the accuracy, completeness, or usefulness of any information, apparatus, product,
or process disclosed, or represents that its use would not infringe privately owned rights.
Reference herein to any specific commercial product, process, or service by trade name,
trademark, manufacturer, or otherwise does not necessarily constitute or imply its endorsement,
recommendation, or favoring by the United States government or Lawrence Livermore National
Security, LLC. The views and opinions of authors expressed herein do not necessarily state
or reflect those of the United States government or Lawrence Livermore National Security, LLC,
and shall not be used for advertising or product endorsement purposes.

\vskip2\baselineskip
\noindent
This work was performed under the auspices of the U.S. Department of Energy by Lawrence
Livermore National Laboratory under Contract DE-AC52-07NA27344.
\vfill
\begin{center}
  Approved for public release; further dissemination unlimited
\end{center}

\clearpage

% Contributors
\thispagestyle{empty}% no number of this page
\vglue5\baselineskip
\begin{center}
{\bf CONTRIBUTORS}
\end{center}
\noindent
The SUNDIALS library has been developed over many years by a number of
contributors. The current SUNDIALS team consists of Cody J. Balos,
David J. Gardner, Alan C. Hindmarsh, Daniel R. Reynolds, and
Carol S. Woodward. We thank Radu Serban for significant and critical past
contributions.\\
\vskip 2em%
\noindent
Other contributors to SUNDIALS include: James Almgren-Bell, Lawrence E. Banks,
Peter N. Brown, George Byrne, Rujeko Chinomona, Scott D. Cohen, Aaron Collier,
Keith E. Grant, Steven L. Lee, Shelby L. Lockhart, John Loffeld, Daniel McGreer,
Slaven Peles, Cosmin Petra, H. Hunter Schwartz, Jean M. Sexton,
Dan Shumaker, Steve G. Smith, Allan G. Taylor, Hilari C. Tiedeman, Chris White,
Ting Yan, and Ulrike M. Yang.
\clearpage

% clear empty double page
\newpage{\pagestyle{empty}\cleardoublepage}

% Start arabic numbering
\pagenumbering{arabic}

% =====================================================
% End custom cover page (see latex_elements in conf.py)
% =====================================================
'''
}


# -- Options for manual page output --------------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    ('index', 'cvode', u'CVODE Documentation',
     [u'Alan C. Hindmarsh, Radu Serban, Cody J. Balos, David J. Gardner, Daniel R. Reynolds, and Carol S. Woodward'], 1)
]

# If true, show URL addresses after external links.
#man_show_urls = False


# -- Options for Texinfo output ------------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
  ('index', 'CVODE', u'CVODE Documentation',
   u'Alan C. Hindmarsh, Radu Serban, Cody J. Balos, David J. Gardner, Daniel R. Reynolds, and Carol S. Woodward', 'CVODE',
   'Time integration package for multi-rate systems of ordinary differntial equations.',
   'Miscellaneous'),
]

# Documents to append as an appendix to all manuals.
#texinfo_appendices = []

# If false, no module index is generated.
#texinfo_domain_indices = True

# How to display URL addresses: 'footnote', 'no', or 'inline'.
#texinfo_show_urls = 'footnote'
