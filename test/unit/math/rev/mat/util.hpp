#ifndef TEST_UNIT_MATH_REV_MAT_UTIL_HPP
#define TEST_UNIT_MATH_REV_MAT_UTIL_HPP

#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/arr/util.hpp>

namespace test {

  void check_varis_on_stack(const stan::math::matrix_v& x) {
    for (int j = 0; j < x.cols(); ++j)
      for (int i = 0; i < x.rows(); ++i) 
        EXPECT_TRUE(stan::math::ChainableStack::memalloc_.in_stack(x(i, j).vi_))
          << i << ", " << j << " is not on the stack";
  }

  void check_varis_on_stack(const stan::math::vector_v& x) {
    for (int i = 0; i < x.rows(); ++i)
      EXPECT_TRUE(stan::math::ChainableStack::memalloc_.in_stack(x(i).vi_))
        << i << " is not on the stack";
  }

  void check_varis_on_stack(const stan::math::row_vector_v& x) {
    for (int j = 0; j < x.cols(); ++j)
      EXPECT_TRUE(stan::math::ChainableStack::memalloc_.in_stack(x(j).vi_))
        << j << ", " << j << " is not on the stack";
  }
  
}
#endif
