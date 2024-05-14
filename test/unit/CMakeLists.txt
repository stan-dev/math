
add_library(test_ad_pch INTERFACE)
target_precompile_headers(test_ad_pch INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/unit/math/test_ad.hpp)
target_include_directories(test_ad_pch INTERFACE ${CMAKE_SOURCE_DIR} ${boost_SOURCE_DIR})
target_link_libraries(test_ad_pch INTERFACE benchmark::benchmark TBB::tbb Eigen3::Eigen sundials_kinsol sundials_cvodes)

add_gtest_grouped_test(${CMAKE_CURRENT_SOURCE_DIR} test_ad_pch)
