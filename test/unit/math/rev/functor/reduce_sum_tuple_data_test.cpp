#include <stan/math.hpp>
#include <test/unit/math/rev/util.hpp>
#include <gtest/gtest.h>
#include <tuple>
#include <vector>

namespace {

struct sum_with_laplace_options {
  template <typename Options>
  auto operator()(const std::vector<int>& slice, std::size_t start,
                  std::size_t end, std::ostream* msgs,
                  const stan::math::var& shared,
                  const Options& laplace_options) const {
    return slice.size()
           * (shared + static_cast<double>(std::get<2>(laplace_options)));
  }
};

// Regression test for https://github.com/stan-dev/math/issues/3359.
TEST_F(AgradRev, reduce_sum_accepts_const_laplace_options_tuple) {
  stan::math::var shared = 1.0;
  const auto laplace_options = stan::math::generate_laplace_options(1);

  stan::math::var result = stan::math::reduce_sum<sum_with_laplace_options>(
      std::vector<int>{0, 1}, 1, nullptr, shared, laplace_options);

  const double expected = 2.0 * (shared.val() + std::get<2>(laplace_options));
  EXPECT_FLOAT_EQ(expected, result.val());
  result.grad();
  EXPECT_FLOAT_EQ(2.0, shared.adj());
}

}  // namespace
