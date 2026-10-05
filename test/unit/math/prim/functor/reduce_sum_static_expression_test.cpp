#include <stan/math.hpp>
#include <gtest/gtest.h>

#include <atomic>
#include <vector>

namespace stan {
namespace math {
namespace test {

// Counts how many times the expression's coefficients are read.
template <typename T>
struct counting_op {
  std::atomic<int>* count_;

  explicit counting_op(std::atomic<int>* count) : count_(count) {}

  inline const T& operator()(const T& a) const {
    ++(*count_);
    return a;
  }
};

struct expression_sum_lpdf {
  template <typename VecT, typename EigExpr>
  inline auto operator()(const VecT& sub_slice, std::size_t start,
                         std::size_t end, std::ostream* msgs,
                         const EigExpr& shared) const {
    // Read the expression once per call so the counter tracks calls, not terms.
    stan::return_type_t<VecT, EigExpr> shared_sum = stan::math::sum(shared);
    stan::return_type_t<VecT, EigExpr> sum = 0;
    for (std::size_t i = 0; i < sub_slice.size(); ++i) {
      sum += sub_slice[i] * shared_sum;
    }
    return sum;
  }
};

}  // namespace test
}  // namespace math
}  // namespace stan

// https://github.com/stan-dev/math/issues/3304
TEST(StanMathPrim_reduce_sum_static, eigen_expression_arg_evaluated_once) {
  using stan::math::test::counting_op;
  using stan::math::test::expression_sum_lpdf;

  Eigen::MatrixXd m = Eigen::MatrixXd::Ones(2, 2);
  std::vector<double> slice{1.0, 2.0, 3.0, 4.0};

  std::atomic<int> reduce_sum_count{0};
  counting_op<double> op_dynamic(&reduce_sum_count);
  double from_reduce_sum = stan::math::reduce_sum<expression_sum_lpdf>(
      slice, 1, nullptr, m.block(0, 0, 1, 1).unaryExpr(op_dynamic));

  std::atomic<int> reduce_sum_static_count{0};
  counting_op<double> op_static(&reduce_sum_static_count);
  double from_reduce_sum_static
      = stan::math::reduce_sum_static<expression_sum_lpdf>(
          slice, 1, nullptr, m.block(0, 0, 1, 1).unaryExpr(op_static));

  EXPECT_DOUBLE_EQ(from_reduce_sum, 10.0);
  EXPECT_DOUBLE_EQ(from_reduce_sum_static, 10.0);

  EXPECT_EQ(1, reduce_sum_count.load());
  EXPECT_EQ(1, reduce_sum_static_count.load());
}
