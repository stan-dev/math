#include <stan/math/rev.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

// there is no agrad-defined version of promote_scalar, so this is
// just testing that it works with non-inter-convertible types (double
// can be assigned to var, but not vice-versa)

TEST(AgradRevFunctionsPromoteScalar, Mismatch) {
  using stan::math::promote_scalar;
  using stan::math::var;
  EXPECT_FLOAT_EQ(2.3, promote_scalar<var>(2.3).val());
}

TEST(AgradRevFunctionsPromoteScalar, var_matrix_promotions) {
  using stan::math::promote_scalar;
  using stan::math::promote_scalar_t;
  using stan::math::var_value;
  Eigen::MatrixXd tester = Eigen::MatrixXd::Random(2, 2);
  using var_mat = var_value<Eigen::MatrixXd>;
  using double_mat = promote_scalar_t<double, var_mat>;
  EXPECT_MATRIX_EQ(tester, (double_mat(promote_scalar<var_mat>(tester).val())));
}
