# Test flags and precompiled headers shared by the test suite and the
# compile-only checks. A single set of flags and dependencies keeps PCH
# producers and consumers compatible. Do not enable IPO or add
# target-specific optimization flags.
set(STAN_EFFECTIVE_TEST_PCH ${STAN_TEST_PCH})
if(STAN_TEST_PCH AND CMAKE_CXX_COMPILER_ID STREQUAL "Clang"
    AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL "9"
    AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "11")
  # Clang 9/10 cannot deserialize Eigen's __builtin_bit_cast expressions.
  # Fixed in LLVM 95d7ccb70b9cbd53f1f137c0b2411852c42c122b (Clang 11).
  set(STAN_EFFECTIVE_TEST_PCH OFF)
  message(STATUS "Test PCH disabled: Clang 9/10 PCH bit_cast defect (fixed in Clang 11)")
endif()

add_library(stan_test_flags INTERFACE)
target_link_libraries(stan_test_flags INTERFACE stanmath)
if(MSVC)
  target_compile_options(stan_test_flags INTERFACE ${STAN_TEST_OPTIMIZATION_FLAGS} /W3)
else()
  target_compile_options(stan_test_flags INTERFACE
    ${STAN_TEST_OPTIMIZATION_FLAGS} -Wall -Wno-deprecated-declarations -Wno-sign-compare)
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(stan_test_flags INTERFACE
      -Werror=invalid-pch -Wno-misleading-indentation)
  endif()
endif()

# Producers and consumers of a PCH must link the same options target.
function(stan_add_pch target options)
  if(STAN_EFFECTIVE_TEST_PCH)
    add_library(${target} OBJECT EXCLUDE_FROM_ALL "${PROJECT_SOURCE_DIR}/cmake/empty.cpp")
    set_target_properties(${target} PROPERTIES POSITION_INDEPENDENT_CODE FALSE)
    target_link_libraries(${target} PRIVATE ${options})
    target_precompile_headers(${target} PRIVATE ${ARGN})
  endif()
endfunction()

# A PCH target exists only when precompiled headers are in effect.
function(stan_configure_test target options pch)
  set_target_properties(${target} PROPERTIES POSITION_INDEPENDENT_CODE FALSE)
  target_link_libraries(${target} PRIVATE ${options})
  if(TARGET "${pch}")
    target_precompile_headers(${target} REUSE_FROM ${pch})
  endif()
endfunction()
