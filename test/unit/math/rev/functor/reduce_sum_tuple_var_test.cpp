#include <stan/math.hpp>
#include <test/unit/math/rev/util.hpp>
#include <gtest/gtest.h>
#include <tuple>
#include <vector>

namespace {

struct sum_tuple_var {
  template <typename Tuple>
  auto operator()(const std::vector<int>& slice, std::size_t start,
                  std::size_t end, std::ostream* msgs,
                  const Tuple& shared) const {
    return slice.size()
           * (stan::math::sum(std::get<0>(shared))
              + stan::math::sum(std::get<1>(shared)));
  }
};

// Regression test for https://github.com/stan-dev/math/issues/3041.
TEST_F(AgradRev, reduce_sum_accepts_tuple_containing_vars) {
  using stan::math::var;

  Eigen::Matrix<var, Eigen::Dynamic, 1> first(2);
  first << 1.0, 2.0;
  Eigen::Matrix<var, Eigen::Dynamic, 1> second(1);
  second << 3.0;
  auto shared = std::make_tuple(first, second);

  var result = stan::math::reduce_sum<sum_tuple_var>(std::vector<int>{0, 1}, 1,
                                                     nullptr, shared);

  EXPECT_FLOAT_EQ(12.0, result.val());
  result.grad();
  EXPECT_FLOAT_EQ(2.0, first(0).adj());
  EXPECT_FLOAT_EQ(2.0, first(1).adj());
  EXPECT_FLOAT_EQ(2.0, second(0).adj());
}

}  // namespace
