#ifndef TEST_UNIT_MATH_REV_ARR_UTIL_HPP
#define TEST_UNIT_MATH_REV_ARR_UTIL_HPP

#include <stan/math/rev/arr.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/scal/util.hpp>
#include <vector>

namespace test {

void check_varis_on_stack(const std::vector<stan::math::var>& x) {
  for (size_t n = 0; n < x.size(); ++n)
    EXPECT_TRUE(
        stan::math::ChainableStack::instance().memalloc_.in_stack(x[n].vi_))
        << n << " is not on the stack";
}

}  // namespace test
#endif
