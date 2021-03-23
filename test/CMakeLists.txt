add_library(stanmath_gtest_main STATIC ${googletest_SOURCE_DIR}/googletest/src/gtest-all.cc ${googletest_SOURCE_DIR}/googletest/src/gtest_main.cc)
target_compile_features(stanmath_gtest_main PRIVATE cxx_std_11)
target_include_directories(stanmath_gtest_main PRIVATE ${googletest_SOURCE_DIR}/googletest/include $${googletest_SOURCE_DIR}/googletest)

add_executable(test_ad
    ${CMAKE_CURRENT_SOURCE_DIR}/unit/math/test_ad_test.cpp
    )

if (NOT CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    target_compile_options(test_ad PRIVATE -Werror -Wextra)
endif()
target_compile_options(test_ad PRIVATE -g -Wall)
target_include_directories(test_ad PRIVATE
    ${googletest_SOURCE_DIR}/googletest/include
  ${CMAKE_CURRENT_SOURCE_DIR}/stan
${eigen_SOURCE_DIR}/Eigen)

target_link_libraries(test_ad stanmath_gtest_main)
if (NOT CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    target_link_libraries(test_ad pthread)
endif()
add_test(test_ad test_ad)
