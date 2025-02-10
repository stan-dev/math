# Boost.Geometry

# Copyright (c) 2024 Barend Gehrels, Amsterdam, the Netherlands.

# Use, modification and distribution is subject to the Boost Software License,
# Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt)

function(boost_geometry_add_example prefix item)
  set(example_name "boost_geometry_example_${prefix}_${item}")
  add_executable(${example_name} ${item}.cpp)

  # Add a dependency to Boost.Geometry
  target_link_libraries(${example_name}
    PRIVATE
      Boost::geometry
      )

  # Include the main Geometry test folder and the current folder
  target_include_directories(${example_name}
    PRIVATE
      .)

  # To compile with C++14
  target_compile_features(${example_name} PRIVATE cxx_std_14)

endfunction()

foreach(item IN ITEMS
    quick_start
)
  boost_geometry_add_example("example" ${item})
endforeach()

add_subdirectory(algorithms)
add_subdirectory(arithmetic)
add_subdirectory(core)
add_subdirectory(geometries)
add_subdirectory(io)
add_subdirectory(strategies)
add_subdirectory(views)
