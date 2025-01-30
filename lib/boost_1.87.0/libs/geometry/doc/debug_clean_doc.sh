# ===========================================================================
#  Copyright (c) 2024-2024 Barend Gehrels, Amsterdam, the Netherlands.
#
#  Use, modification and distribution is subject to the Boost Software License,
#  Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
#  http://www.boost.org/LICENSE_1_0.txt)
# ============================================================================

# Removes all generated output

b2 clean

rm -f index/generated/*.qbk
rm -f generated/*.qbk
rm -f generated/*.xml
rm -Rf html/geometry
rm -Rf doxy/doxygen_output/xml
rm -Rf doxy/doxygen_output/html_by_doxygen
rm -Rf index/xml
rm -Rf index/html_by_doxygen

git ls-files --others
