# ===========================================================================
#  Copyright (c) 2024-2024 Barend Gehrels, Amsterdam, the Netherlands.  #
#  Use, modification and distribution is subject to the Boost Software License,
#  Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
#  http://www.boost.org/LICENSE_1_0.txt)
# ============================================================================

# Makes the documentation without calling b2, to be able to do it in steps
# and debug the individual steps.
#
# The following software must be in your path:
# - doxygen
# - xsltproc
# - python3
# - doxygen_xml2qbk (provided by Boost.Geometry)
# - quickbook (should be made, it is referred to in a relative path)
# 
# Execute the script from .../boost/libs/geometry/doc

echo "=== Generating documentation for Boost.Geometry"

base=`pwd`
boost_root=`realpath ${base}/../../..`

# Set the DocBook (and dtd and boost dtd, less relevant) versions. 

# 1.75.2 is the one setup within Boost for ages.
# It reports as debug message: "Computing chunks..."
# But it currently gives the error message 
#   "The called template 'id.attribute' was not found."
# For all libraries (not only Boost.Geometry)

# 1.79.1 works and gives correct documentation
# It does not output any extra debug message

# 1.79.2 (with namespaces, on mac) gives empty content for sub-toc pages for Boost.Geometry
# It reports as debug message: "Note: namesp. add : added namespace before processing   Geometry"
# 1.79.2 ("nons", no namespace) works and gives correct documentation
# It can be downloaded from https://www.linuxfromscratch.org/blfs/view/stable/pst/docbook-xsl.html

# Set the values with 4.2 / 1.75.2 - used by setup_boostbook.sh
docbook_dtd_version=4.2
docbook_xsl_version=1.75.2
docbook_dtd=${boost_root}/tools/boostbook/docbook-dtd-${docbook_dtd_version}
docbook_xsl=${boost_root}/tools/boostbook/docbook-xsl-${docbook_xsl_version}

# OVERRIDE the values with 4.5 / 1.79.1 - used by boostorg/release-tools/blob/develop/build_docs
docbook_dtd_version=4.5
docbook_xsl_version=1.79.1
docbook_dtd=~/data/docbook/docbook-dtd-${docbook_dtd_version}
docbook_xsl=~/data/docbook/docbook-xsl-${docbook_xsl_version}

# UNCOMMENT the values with 4.5 / 1.79.2 - for future usage, note the "nons", which is essential
# docbook_xsl_version=1.79.2
# docbook_xsl=~/data/docbook/docbook-xsl-nons-${docbook_xsl_version}

# Set other values of dtd, css, style sheets
boostbook_dtd=${boost_root}/tools/boostbook/dtd

file_css=${boost_root}/doc/src/boostbook.css

xsl_docbook=${boost_root}/tools/boostbook/xsl/docbook.xsl
xsl_html=${boost_root}/tools/boostbook/xsl/html.xsl

# Set the folder where to generate boostbook/docbook and some filenames
folder_generated=generated
file_boostbook=${folder_generated}/geometry.boostbook.xml
file_docbook=${folder_generated}/geometry.docbook.xml

# === step 1 ===

echo "=== Running Doxygen"
cd doxy
doxygen Doxyfile
cd ..

# === step 2 ===
# Note that doxygen_xml2qbk should be in your path!

echo "=== From Doxygen XML to QuickBook"
python3 make_qbk.py --skip_doxygen

# === step 3 ===

echo "=== From QuickBook to BoostBook"
${boost_root}/dist/bin/quickbook -I"${boost_root}" -I"." -D"enable_index" \
    --output-file=${file_boostbook} geometry.qbk

# === step 4 ===

# To avoid using b2 and user-config.jam, the catalog has to be written locally
filename_catalog=${folder_generated}/catalog.xml
echo "=== Write ${filename_catalog} for xstlproc with ${docbook_xsl}"

cat <<EOF >${filename_catalog}
<?xml version="1.0"?>
<!DOCTYPE catalog 
  PUBLIC "-//OASIS/DTD Entity Resolution XML Catalog V1.0//EN"
  "http://www.oasis-open.org/committees/entity/release/1.0/catalog.dtd">
<catalog xmlns="urn:oasis:names:tc:entity:xmlns:xml:catalog">
  <rewriteURI uriStartString="http://www.boost.org/tools/boostbook/dtd/" rewritePrefix="file:///${boostbook_dtd}/"/>
  <rewriteURI uriStartString="http://www.oasis-open.org/docbook/xml/${docbook_dtd_version}/" rewritePrefix="file:///${docbook_dtd}/"/>
  <rewriteURI uriStartString="http://docbook.sourceforge.net/release/xsl/current/" rewritePrefix="file:///${docbook_xsl}/"/>
</catalog>
EOF

# It is necessary to publish its location in an environment variable
export XML_CATALOG_FILES=${filename_catalog}

# Set parameters shared for Boostbook -> Docbook -> HTML
params="--xinclude"
params="${params} --stringparam boost.defaults Boost"

# Set boost.root (here, it should be relative, and go one above boost root itself)
params="${params} --stringparam boost.root ../../../.."

echo "=== From BoostBook to DocBook"
xsltproc ${params} --path ${folder_generated} -o ${file_docbook} ${xsl_docbook} ${file_boostbook}

# === step 5 ===

# These parameters define behaviour for DocBook.
# See also https://www.sagehill.net/docbookxsl/TOCcontrol.html
# The parameters values should correspond with the values in Jamfile.

# Their values are currently valid for docbook-xsl-1.75.2 (the version included in Boost).
# However, the sections.xsl file currently included boost do not work for 1.75.2
# If this is encountered, remove the part "<xsl:call-template name="id.attribute">"

# Their values are currently also valid for version 1.79.1
# But for version 1.79.2 is used, the chunking is broken.

chunk_section_depth=4
chunk_first_sections=1
toc_section_depth=3
toc_max_depth=2
generate_section_toc_level=4

params="${params} --stringparam chunk.section.depth ${chunk_section_depth}"
params="${params} --stringparam chunk.first.sections ${chunk_first_sections}"
params="${params} --stringparam toc.max.depth ${toc_max_depth}"
params="${params} --stringparam toc.section.depth ${toc_section_depth}"
params="${params} --stringparam generate.section.toc.level ${generate_section_toc_level}"

echo "=== From DocBook to html"
xsltproc ${params} --path ${folder_generated} -o "html/" ${xsl_html} ${file_docbook}

echo "=== Resulting html, generated with:"
cat html/index.html| grep -w -e generator

echo "=== Resulting contents of generated 'adapted.html', there should be sections"
cat html/geometry/reference/adapted.html | grep section

