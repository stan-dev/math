#
# Copyright (c) 2017 Mateusz Loskot <mateusz at loskot dot net>
# All rights reserved.
#
# Distributed under the Boost Software License, Version 1.0.
# (See accompanying file LICENSE_1_0.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt)
#
message(STATUS "Boost.GIL: Configuring tests in test/extension/toolbox")

foreach(_name
  all)
  set(_test t_ext_toolbox_${_name})
  set(_target test_ext_toolbox_${_name})

  add_executable(${_target} "")
  target_sources(${_target}
    PRIVATE
      test.cpp
      channel_type.cpp
      channel_view.cpp
      cmyka.cpp
      get_num_bits.cpp
      get_pixel_type.cpp
      gray_alpha.cpp
      gray_to_rgba.cpp
      hsl_hsv_test.cpp
      indexed_image_test.cpp
      is_bit_aligned.cpp
      is_homogeneous.cpp
      is_similar.cpp
      lab_test.cpp
      pixel_bit_size.cpp
      rgb_to_luminance.cpp
      # TODO: Add subchroma_image.cpp after fixing run-time failure,
      #       for details see https://github.com/boostorg/gil/pull/164
      #subchroma_image.cpp
      xyz_test.cpp)
  target_link_libraries(${_target}
    PRIVATE
      gil_compile_options
      gil_include_directories
      gil_dependencies)
  add_test(NAME ${_test} COMMAND ${_target})

  unset(_name)
  unset(_target)
  unset(_test)
endforeach()
