#
# Copyright (c) 2019 Mateusz Loskot <mateusz at loskot dot net>
#
# Distributed under the Boost Software License, Version 1.0.
# (See accompanying file LICENSE_1_0.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt)
#
add_subdirectory(core)

add_subdirectory(legacy)

add_subdirectory(extension)

# TODO: Split headers tests into core and extensions, see Jamfile-s
# if(GIL_BUILD_HEADERS_TESTS)
#   add_subdirectory(headers)
# endif()
