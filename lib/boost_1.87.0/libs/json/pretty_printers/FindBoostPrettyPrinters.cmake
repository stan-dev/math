#
# Copyright (c) 2024 Dmitry Arkhipov (grisumbras@yandex.ru)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#
# Official repository: https://github.com/boostorg/json
#

find_package(Python3 QUIET COMPONENTS Interpreter)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(BoostPrettyPrinters
    REQUIRED_VARS Python3_Interpreter_FOUND)

find_program(BoostPrettyPrinters_GDB gdb DOC "GDB executable tos use")
set(BoostPrettyPrinters_HAS_GDB "${BoostPrettyPrinters_GDB}")

set(BoostPrettyPrinters_GDB_HEADER_SCRIPT
    "${CMAKE_CURRENT_LIST_DIR}/generate-gdb-header.py")
set(BoostPrettyPrinters_GDB_TEST_SCRIPT
    "${CMAKE_CURRENT_LIST_DIR}/generate-gdb-test-runner.py")
set(BoostPrettyPrinters_INCLUDES "${CMAKE_CURRENT_LIST_DIR}/include")

function(boost_pretty_printers_gdb_python_header)
    set(options NO_DISABLE_MACRO EXCLUDE_FROM_ALL)
    set(oneValueArgs TARGET INPUT OUTPUT HEADER_GUARD DISABLE_MACRO)
    set(multiValueArgs)
    cmake_parse_arguments(BOOST_PPRINT_GDB_GEN
        "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
    foreach(kw TARGET INPUT OUTPUT)
        if(NOT DEFINED "BOOST_PPRINT_GDB_GEN_${kw}")
            message(FATAL_ERROR "Argument ${kw} is required for function \
                boost_pretty_printers_gdb_python_header.")
        endif()
    endforeach()

    if(DEFINED BOOST_PPRINT_GDB_GEN_HEADER_GUARD)
        set(BOOST_PPRINT_GDB_GEN_HEADER_GUARD
            "--header-guard=${BOOST_PPRINT_GDB_GEN_HEADER_GUARD}")
    endif()
    if(DEFINED BOOST_PPRINT_GDB_GEN_DISABLE_MACRO)
        set(BOOST_PPRINT_GDB_GEN_DISABLE_MACRO
            "--disable-macro=${BOOST_PPRINT_GDB_GEN_DISABLE_MACRO}")
    elseif(BOOST_PPRINT_GDB_GEN_NO_DISABLE_MACRO)
        set(BOOST_PPRINT_GDB_GEN_DISABLE_MACRO "--disable-macro=")
    endif()
    add_custom_command(
        OUTPUT "${BOOST_PPRINT_GDB_GEN_OUTPUT}"
        MAIN_DEPENDENCY "${BOOST_PPRINT_GDB_GEN_INPUT}"
        DEPENDS "${BoostPrettyPrinters_GDB_HEADER_SCRIPT}"
        COMMAND
            "${Python3_EXECUTABLE}"
            "${BoostPrettyPrinters_GDB_HEADER_SCRIPT}"
            "${CMAKE_CURRENT_SOURCE_DIR}/${BOOST_PPRINT_GDB_GEN_INPUT}"
            "${CMAKE_CURRENT_SOURCE_DIR}/${BOOST_PPRINT_GDB_GEN_OUTPUT}"
            ${BOOST_PPRINT_GDB_GEN_HEADER_GUARD}
            ${BOOST_PPRINT_GDB_GEN_DISABLE_MACRO}
        COMMENT "Regenerating ${BOOST_PPRINT_GDB_GEN_OUTPUT}")

    if(NOT BOOST_PPRINT_GDB_GEN_EXCLUDE_FROM_ALL)
        set(isInAll ALL)
    endif()
    add_custom_target(${BOOST_PPRINT_GDB_GEN_TARGET}
        ${isInAll}
        DEPENDS "${BOOST_PPRINT_GDB_GEN_OUTPUT}")
endfunction()


function(boost_pretty_printers_test_gdb_printers)
    set(options EXCLUDE_FROM_ALL)
    set(oneValueArgs TEST PROGRAM)
    set(multiValueArgs SOURCES)
    cmake_parse_arguments(BOOST_PPRINT_TEST_GDB
        "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
    foreach(kw TEST SOURCES)
        if(NOT DEFINED "BOOST_PPRINT_TEST_GDB_${kw}")
            message(FATAL_ERROR "Argument ${kw} is required for function \
                boost_pretty_printers_test_gdb_printers.")
        endif()
    endforeach()

    if(NOT DEFINED BOOST_PPRINT_TEST_GDB_PROGRAM)
        set(BOOST_PPRINT_TEST_GDB_PROGRAM ${BOOST_PPRINT_TEST_GDB_TEST})
    endif()
    if(BOOST_PPRINT_TEST_GDB_EXCLUDE_FROM_ALL)
        set(excludeFromAll EXCLUDE_FROM_ALL)
    else()
        set(includeInAll ALL)
    endif()

    LIST(GET BOOST_PPRINT_TEST_GDB_SOURCES 0 source0)
    add_custom_command(
        OUTPUT "${BOOST_PPRINT_TEST_GDB_TEST}.py"
        DEPENDS "${source0}"
        COMMAND
            "${Python3_EXECUTABLE}"
            "${BoostPrettyPrinters_GDB_TEST_SCRIPT}"
            "${CMAKE_CURRENT_SOURCE_DIR}/${source0}"
            "${BOOST_PPRINT_TEST_GDB_TEST}.py"
        COMMENT "Generating ${source0}")

    add_custom_target(${BOOST_PPRINT_TEST_GDB_TEST}_runner
        ${includeInAll}
        DEPENDS "${BOOST_PPRINT_TEST_GDB_TEST}.py")

    add_executable(${BOOST_PPRINT_TEST_GDB_PROGRAM}
        ${excludeFromAll}
        ${BOOST_PPRINT_TEST_GDB_SOURCES})
    add_dependencies(
        ${BOOST_PPRINT_TEST_GDB_PROGRAM}
        ${BOOST_PPRINT_TEST_GDB_TEST}_runner)

    add_test(
        NAME ${BOOST_PPRINT_TEST_GDB_TEST}
        COMMAND "${BoostPrettyPrinters_GDB}"
            --batch-silent
            -x "${BOOST_PPRINT_TEST_GDB_TEST}.py"
            $<TARGET_FILE:${BOOST_PPRINT_TEST_GDB_PROGRAM}>)
endfunction()
