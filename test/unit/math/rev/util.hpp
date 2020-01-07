#ifndef TEST_UNIT_MATH_REV_UTIL_HPP
#define TEST_UNIT_MATH_REV_UTIL_HPP

#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <vector>

namespace test {

void check_varis_on_stack(const stan::math::var& x) {
  EXPECT_TRUE(stan::math::ChainableStack::instance_->memalloc_.in_stack(x.vi_))
      << "not on the stack";
}

void check_varis_on_stack(const std::vector<stan::math::var>& x) {
  for (size_t n = 0; n < x.size(); ++n)
    EXPECT_TRUE(
        stan::math::ChainableStack::instance_->memalloc_.in_stack(x[n].vi_))
        << n << " is not on the stack";
}

template <int R, int C>
void check_varis_on_stack(const Eigen::Matrix<stan::math::var, R, C>& x) {
  for (int j = 0; j < x.cols(); ++j)
    for (int i = 0; i < x.rows(); ++i)
      EXPECT_TRUE(stan::math::ChainableStack::instance_->memalloc_.in_stack(
          x(i, j).vi_))
          << i << ", " << j << " is not on the stack";
}

}  // namespace test
#endif
