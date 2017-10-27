#ifndef TEST_UNIT_MATH_REV_SCAL_UTIL_HPP
#define TEST_UNIT_MATH_REV_SCAL_UTIL_HPP

#include <gtest/gtest.h>
#include <stan/math/rev/scal.hpp>

namespace test {

  void check_varis_on_stack(const stan::math::var& x) {
    EXPECT_TRUE(stan::math::ChainableStack::memalloc_.in_stack(x.vi_))
        << "not on the stack";
  }

}  // namespace test
#endif
