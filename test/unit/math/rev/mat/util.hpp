#ifndef TEST_UNIT_MATH_REV_MAT_UTIL_HPP
#define TEST_UNIT_MATH_REV_MAT_UTIL_HPP

#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/arr/util.hpp>

namespace test {

template <typename T, stan::enable_if_eigen<T>* = nullptr,
          stan::enable_if_contains_var<T>* = nullptr>
void check_varis_on_stack(const T& x) {
  for (int j = 0; j < x.cols(); ++j)
    for (int i = 0; i < x.rows(); ++i)
      EXPECT_TRUE(stan::math::ChainableStack::instance_->memalloc_.in_stack(
          x(i, j).vi_))
          << i << ", " << j << " is not on the stack";
}

}  // namespace test
#endif
