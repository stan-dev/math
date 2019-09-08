#ifndef TEST_UNIT_MATH_REV_MAT_UTIL_HPP
#define TEST_UNIT_MATH_REV_MAT_UTIL_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/rev/mat.hpp>
#include <stan/math/prim/meta.hpp>
#include <test/unit/math/rev/arr/util.hpp>
#include <gtest/gtest.h>

namespace test {

template <typename Mat, stan::require_eigen_t<stan::is_var, Mat>...>
void check_varis_on_stack(Mat&& x) {
  x.eval();
  for (int j = 0; j < x.cols(); ++j)
    for (int i = 0; i < x.rows(); ++i)
      EXPECT_TRUE(stan::math::ChainableStack::instance_->memalloc_.in_stack(
          x(i, j).vi_))
          << i << ", " << j << " is not on the stack";
}

}  // namespace test
#endif
