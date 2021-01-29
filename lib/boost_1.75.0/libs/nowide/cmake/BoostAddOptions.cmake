# Add common options for the current library project
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR OR BOOST_SUPERPROJECT_SOURCE_DIR)
  set(def_INSTALL ON)
  set(def_WERROR ON)
else()
  set(def_INSTALL OFF)
  set(def_WERROR OFF)
endif()

string(TOUPPER ${PROJECT_NAME} name)
string(REPLACE BOOST_ Boost_ name ${name})

option(${name}_INSTALL "Install library" "${def_INSTALL}")
option(${name}_WERROR "Treat warnings as errors" "${def_WERROR}")
