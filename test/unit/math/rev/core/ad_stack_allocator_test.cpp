#include <stan/math/rev/core.hpp>
#include <gtest/gtest.h>
#include <vector>

void ad_stack_allocator_test() {
  std::vector<int, stan::math::AD_stack_allocator<int>> v;
  v.reserve(3);
  for (int i = 0; i < 5; i++) {
    v.push_back(1);
  }

  std::vector<int, stan::math::AD_stack_allocator<int>> v2(2, 3);
  for (int i = 0; i < 5; i++) {
    v2.push_back(2);
  }

  std::vector<int, stan::math::AD_stack_allocator<int>> v3(v2);
  for (int i = 0; i < 15; i++) {
    v3.push_back(3);
  }
}

TEST(AgradRev, AD_stack_allocator_test) { EXPECT_NO_THROW(ad_stack_allocator_test()); }
