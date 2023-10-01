# Copyright 2022 Alexander Grund
# Distributed under the Boost Software License, Version 1.0.
# https://www.boost.org/LICENSE_1_0.txt

include(BoostTestJamfile)

add_library(boost_locale_test INTERFACE)
# Add test folder to include directories, especially for systems
# where the current folder is not automatically added to the search path
target_include_directories(boost_locale_test INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}")

set(BOOST_TEST_COMPILE_DEFINITIONS "")
set(BOOST_TEST_LINK_LIBRARIES Boost::locale boost_locale_test)
set(BOOST_TEST_COMPILE_OPTIONS ${BOOST_LOCALE_WARNING_OPTIONS})
if(MSVC OR
  (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 11)
)
  list(APPEND BOOST_TEST_COMPILE_OPTIONS -Wall -Werror)
endif()

# Library config as determined by configure tests
if(BOOST_LOCALE_ENABLE_ICONV)
  list(APPEND BOOST_TEST_COMPILE_DEFINITIONS BOOST_LOCALE_WITH_ICONV=1)
endif()
if(BOOST_LOCALE_ENABLE_ICU)
  list(APPEND BOOST_TEST_COMPILE_DEFINITIONS BOOST_LOCALE_WITH_ICU=1)
  list(APPEND BOOST_TEST_LINK_LIBRARIES ICU::i18n ICU::uc)
endif()
if(NOT BOOST_LOCALE_ENABLE_STD)
  list(APPEND BOOST_TEST_COMPILE_DEFINITIONS BOOST_LOCALE_NO_STD_BACKEND=1)
endif()
if(NOT BOOST_LOCALE_ENABLE_WINAPI)
  list(APPEND BOOST_TEST_COMPILE_DEFINITIONS BOOST_LOCALE_NO_WINAPI_BACKEND=1)
endif()
if(NOT BOOST_LOCALE_ENABLE_POSIX)
  list(APPEND BOOST_TEST_COMPILE_DEFINITIONS BOOST_LOCALE_NO_POSIX_BACKEND=1)
endif()

boost_test_jamfile(FILE Jamfile.v2)

# Those require to be run in the test directory
foreach(name test_formatting test_message)
  boost_test(SOURCES ${name}.cpp WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
endforeach()
