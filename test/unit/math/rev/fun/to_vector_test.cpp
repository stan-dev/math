#include <stan/math/rev.hpp>
#include <test/unit/util.hpp>

#include <gtest/gtest.h>

TEST(MathFunRev, to_vector_var_value) {
  constexpr Eigen::Index n = 100;
  Eigen::MatrixXd a_val = Eigen::MatrixXd::Random(n, n);
  stan::math::var_value<Eigen::MatrixXd> a(a_val);
  const auto& tmp = stan::math::to_vector(a);
  EXPECT_TRUE((stan::is_col_vector<decltype(tmp.val())>::value));
  EXPECT_MATRIX_EQ(tmp.val(), Eigen::Map<Eigen::VectorXd>(a_val.data(), a.rows() * a.cols()));
}
