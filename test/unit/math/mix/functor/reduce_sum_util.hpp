#ifndef TEST_UNIT_MATH_MIX_FUNCTOR_REDUCE_SUM_UTIL
#define TEST_UNIT_MATH_MIX_FUNCTOR_REDUCE_SUM_UTIL

#include <stan/math.hpp>
#include <test/unit/math/prim/functor/reduce_sum_util.hpp>
#include <test/unit/math/test_ad.hpp>
#include <gtest/gtest.h>
#include <algorithm>
#include <sstream>
#include <tuple>
#include <vector>

namespace stan {
namespace math {
namespace test {

template <typename T1, typename T2>
void expect_ad_reduce_sum_lpdf(T1&& data, T2&& arg) {
  using stan::math::test::reduce_sum_int_sum_lpdf;
  using stan::math::test::reduce_sum_static_int_sum_lpdf;
  using stan::math::test::reduce_sum_static_sum_lpdf;
  using stan::math::test::reduce_sum_sum_lpdf;
  stan::test::expect_ad(reduce_sum_static_int_sum_lpdf, arg);
  stan::test::expect_ad(reduce_sum_static_sum_lpdf, data, arg);
  stan::test::expect_ad(reduce_sum_int_sum_lpdf, arg);
  stan::test::expect_ad(reduce_sum_sum_lpdf, data, arg);
}

}  // namespace test
}  // namespace math
}  // namespace stan
#endif
