#
# Copyright (c) 2019 Vinnie Falco (vinnie.falco@gmail.com)
# Copyright (c) 2021 Dmitry Arkhipov (grisumbras@gmail.com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#
# Official repository: https://github.com/boostorg/json
#

if(BOOST_JSON_IS_ROOT)
    set(EXCLUDE_TESTS_FROM_ALL)
else()
    set(EXCLUDE_TESTS_FROM_ALL EXCLUDE_FROM_ALL)
endif()

if(NOT TARGET tests)
    add_custom_target(tests ${EXCLUDE_TESTS_FROM_ALL})
    set_property(TARGET tests PROPERTY FOLDER _deps)
endif()


file(GLOB_RECURSE BOOST_JSON_TESTS_FILES CONFIGURE_DEPENDS Jamfile *.cpp *.hpp)
list(FILTER BOOST_JSON_TESTS_FILES EXCLUDE REGEX cmake_install_test/.*$)
list(FILTER BOOST_JSON_TESTS_FILES EXCLUDE REGEX cmake-subdir/.*$)
list(FILTER BOOST_JSON_TESTS_FILES EXCLUDE REGEX intrusive_macros\.cpp$)
list(FILTER BOOST_JSON_TESTS_FILES EXCLUDE REGEX limits\.cpp$)
list(FILTER BOOST_JSON_TESTS_FILES EXCLUDE REGEX no_exceptions\.cpp$)
list(FILTER BOOST_JSON_TESTS_FILES EXCLUDE REGEX printers\.cpp$)

source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} PREFIX "" FILES ${BOOST_JSON_TESTS_FILES})
add_executable(boost_json-tests ${EXCLUDE_TESTS_FROM_ALL} ${BOOST_JSON_TESTS_FILES})
target_include_directories(boost_json-tests PRIVATE .)
target_link_libraries(boost_json-tests PRIVATE Boost::json)
target_compile_options(boost_json-tests PRIVATE
    $<$<CXX_COMPILER_ID:MSVC>:/bigobj>
)
add_test(NAME boost_json-tests COMMAND boost_json-tests)
add_dependencies(tests boost_json-tests)


source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} PREFIX "" FILES limits.cpp main.cpp intrusive_macros.cpp)
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR}/../src PREFIX "" FILES ../src/src.cpp)
add_executable(boost_json-limits ${EXCLUDE_TESTS_FROM_ALL} limits.cpp main.cpp ../src/src.cpp Jamfile)
boost_json_setup_properties(boost_json-limits)

target_compile_definitions(boost_json-limits PRIVATE
    BOOST_JSON_MAX_STRING_SIZE=1000
    BOOST_JSON_MAX_STRUCTURED_SIZE=20
    BOOST_JSON_STACK_BUFFER_SIZE=256
    BOOST_JSON_NO_LIB=1
)

add_test(NAME boost_json-limits COMMAND boost_json-limits)
add_dependencies(tests boost_json-limits)

add_executable(boost_json-no_exceptions ${EXCLUDE_TESTS_FROM_ALL} no_exceptions.cpp main.cpp ../src/src.cpp Jamfile)
boost_json_setup_properties(boost_json-no_exceptions)

target_compile_definitions(boost_json-no_exceptions PRIVATE
    BOOST_JSON_NO_LIB=1
)
target_compile_options(boost_json-no_exceptions PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-fno-exceptions>
    $<$<CXX_COMPILER_ID:MSVC>:/wd4530>
    $<$<CXX_COMPILER_ID:MSVC>:/wd4577>
)

add_test(NAME boost_json-no_exceptions COMMAND boost_json-no_exceptions)
add_dependencies(tests boost_json-no_exceptions)

add_executable(boost_json-intrusive-macro-tests ${EXCLUDE_TESTS_FROM_ALL} intrusive_macros.cpp main.cpp Jamfile)
boost_json_setup_properties(boost_json-intrusive-macro-tests)
target_compile_definitions(boost_json-intrusive-macro-tests PRIVATE
    BOOST_JSON_NO_LIB=1
)
add_test(NAME boost_json-intrusive-macro-tests COMMAND boost_json-intrusive-macro-tests)
add_dependencies(tests boost_json-intrusive-macro-tests)

if(BoostPrettyPrinters_HAS_GDB)
    boost_pretty_printers_test_gdb_printers(
        TEST boost_json_test_gdb_printers
        SOURCES printers.cpp
        ${EXCLUDE_TESTS_FROM_ALL})
    boost_json_setup_properties(boost_json_test_gdb_printers)
    target_link_libraries(boost_json_test_gdb_printers PRIVATE Boost::json)
    add_dependencies(tests boost_json_test_gdb_printers)
endif()
