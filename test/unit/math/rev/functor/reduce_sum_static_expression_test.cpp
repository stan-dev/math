#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <gtest/gtest.h>

#include <vector>

namespace stan {
namespace math {
namespace test {

struct expression_sum_lpdf {
  template <typename VecT, typename EigExpr>
  inline auto operator()(const VecT& sub_slice, std::size_t start,
                         std::size_t end, std::ostream* msgs,
                         const EigExpr& shared) const {
    stan::math::var sum = 0;
    for (std::size_t i = 0; i < sub_slice.size(); ++i) {
      sum += sub_slice[i] * stan::math::sum(shared);
    }
    return sum;
  }
};

}  // namespace test
}  // namespace math
}  // namespace stan

// https://github.com/stan-dev/math/issues/3304
TEST_F(AgradRev, StanMathRev_reduce_sum_static_eigen_expression_arg_gradient) {
  using stan::math::var;
  using stan::math::test::expression_sum_lpdf;

  Eigen::Matrix<var, -1, -1> m = Eigen::MatrixXd::Ones(2, 2);
  std::vector<var> slice{1.0, 2.0, 3.0, 4.0};

  var out = stan::math::reduce_sum_static<expression_sum_lpdf>(
      slice, 1, nullptr,
      m.block(0, 0, 1, 1).unaryExpr([](const var& a) { return a; }));

  EXPECT_DOUBLE_EQ(10.0, out.val());

  out.grad();
  EXPECT_DOUBLE_EQ(10.0, m(0, 0).adj());
  EXPECT_DOUBLE_EQ(1.0, slice[0].adj());
}
