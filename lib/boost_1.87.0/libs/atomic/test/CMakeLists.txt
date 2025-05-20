# Copyright 2024 Andrey Semashev
#
# Distributed under the Boost Software License, Version 1.0.
# See accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt

include(BoostTest OPTIONAL RESULT_VARIABLE HAVE_BOOST_TEST)

if (NOT HAVE_BOOST_TEST)
    return()
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

set(BOOST_TEST_COMPILE_FEATURES
    cxx_lambdas
    cxx_static_assert
    cxx_deleted_functions
    cxx_defaulted_functions
)

set(BOOST_TEST_LINK_LIBRARIES
    Boost::atomic
    Boost::config
    Boost::core
    Boost::type_traits
)

if (WIN32)
    set(BOOST_TEST_COMPILE_DEFINITIONS
        "_CRT_SECURE_NO_WARNINGS"
        "_CRT_SECURE_NO_DEPRECATE"
        "BOOST_USE_WINDOWS_H"
    )
    list(APPEND BOOST_TEST_LINK_LIBRARIES
        Boost::winapi
    )
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        list(APPEND BOOST_TEST_LINK_LIBRARIES
            kernel32
        )
    endif()
endif()

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(BOOST_TEST_COMPILE_OPTIONS "-Wall" "-Wextra" "-Werror")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    set(BOOST_TEST_COMPILE_OPTIONS "/W4" "/WX")
endif()

boost_test(TYPE run SOURCES atomic_api.cpp)
boost_test(TYPE run SOURCES atomic_ref_api.cpp)
boost_test(TYPE run SOURCES atomic_api.cpp COMPILE_DEFINITIONS BOOST_ATOMIC_FORCE_FALLBACK NAME fallback_atomic_api)
boost_test(TYPE run SOURCES atomic_ref_api.cpp COMPILE_DEFINITIONS BOOST_ATOMIC_FORCE_FALLBACK NAME fallback_atomic_ref_api)
boost_test(TYPE run SOURCES wait_api.cpp)
boost_test(TYPE run SOURCES wait_ref_api.cpp)
boost_test(TYPE run SOURCES wait_api.cpp COMPILE_DEFINITIONS BOOST_ATOMIC_FORCE_FALLBACK NAME fallback_wait_api)
boost_test(TYPE run SOURCES wait_ref_api.cpp COMPILE_DEFINITIONS BOOST_ATOMIC_FORCE_FALLBACK NAME fallback_wait_ref_api)
boost_test(TYPE run SOURCES wait_fuzz.cpp)
boost_test(TYPE run SOURCES wait_fuzz.cpp COMPILE_DEFINITIONS BOOST_ATOMIC_FORCE_FALLBACK NAME fallback_wait_fuzz)
boost_test(TYPE run SOURCES ipc_atomic_api.cpp)
boost_test(TYPE run SOURCES ipc_atomic_ref_api.cpp)
boost_test(TYPE run SOURCES ipc_wait_api.cpp)
boost_test(TYPE run SOURCES ipc_wait_ref_api.cpp)
boost_test(TYPE run SOURCES atomicity.cpp)
boost_test(TYPE run SOURCES atomicity_ref.cpp)
boost_test(TYPE run SOURCES ordering.cpp)
boost_test(TYPE run SOURCES ordering_ref.cpp)
boost_test(TYPE run SOURCES lockfree.cpp)

unset(BOOST_TEST_COMPILE_OPTIONS)

boost_test(TYPE compile-fail SOURCES cf_arith_void_ptr.cpp)
boost_test(TYPE compile-fail SOURCES cf_arith_func_ptr.cpp)
boost_test(TYPE compile-fail SOURCES cf_arith_mem_ptr.cpp)

boost_test(TYPE compile SOURCES c_implicit_ctor.cpp)
