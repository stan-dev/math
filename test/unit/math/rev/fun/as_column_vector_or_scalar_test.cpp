#include <stan/math/rev.hpp>
#include <test/unit/util.hpp>

#include <gtest/gtest.h>

TEST(MathFunRev, as_column_vector_or_scalar_var_value_vector) {
  int n = 100;
  const Eigen::VectorXd a_val = Eigen::VectorXd::Random(n);
  stan::math::var_value<Eigen::VectorXd> a(a_val);
  const auto& tmp = stan::math::as_column_vector_or_scalar(a);
  EXPECT_TRUE((stan::is_col_vector<decltype(tmp.val())>::value));
  EXPECT_MATRIX_EQ(tmp.val(), a_val);
}

TEST(MathFunRev, as_column_vector_or_scalar_var_value_rowvector) {
  int n = 100;
  const Eigen::RowVectorXd a_val = Eigen::RowVectorXd::Random(n);
  stan::math::var_value<Eigen::RowVectorXd> a(a_val);
  const auto& tmp = stan::math::as_column_vector_or_scalar(a);
  EXPECT_TRUE((stan::is_col_vector<decltype(tmp.val())>::value));
  EXPECT_MATRIX_EQ(tmp.val(), a_val.transpose());
}
