#ifndef TEST_UNIT_MATH_REV_MAT_UTIL_HPP
#define TEST_UNIT_MATH_REV_MAT_UTIL_HPP

#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/arr/util.hpp>

namespace test {

  template <int R, int C>
  void check_varis_on_stack(const Eigen::Matrix<stan::math::var, R, C>& x) {
    for (int j = 0; j < x.cols(); ++j)
      for (int i = 0; i < x.rows(); ++i) 
        EXPECT_TRUE(stan::math::ChainableStack::memalloc_.in_stack(x(i, j).vi_))
          << i << ", " << j << " is not on the stack";
  }

}
#endif
