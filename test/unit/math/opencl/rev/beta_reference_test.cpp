#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>
#include <test/unit/math/prim/fun/beta_reference.hpp>
#include <cmath>
#include <limits>

namespace beta_reference_test {

Eigen::VectorXd column(int k) {
  const auto& rows = stan::test::beta_large_reference;
  Eigen::VectorXd v(rows.size());
  for (size_t i = 0; i < rows.size(); ++i) {
    v[i] = rows[i][k];
  }
  return v;
}

}  // namespace beta_reference_test

TEST(OpenCL_beta, prim_matches_reference) {
  using stan::math::matrix_cl;
  const auto& rows = stan::test::beta_large_reference;
  matrix_cl<double> a_cl(beta_reference_test::column(0));
  matrix_cl<double> b_cl(beta_reference_test::column(1));
  matrix_cl<double> res_cl = stan::math::beta(a_cl, b_cl);
  Eigen::VectorXd res = stan::math::from_matrix_cl<Eigen::VectorXd>(res_cl);
  for (size_t i = 0; i < rows.size(); ++i) {
    SCOPED_TRACE(::testing::Message()
                 << "a = " << rows[i][0] << ", b = " << rows[i][1]);
    EXPECT_NEAR(rows[i][2], res[i], 1e-13 * rows[i][2]);
  }
}

TEST(OpenCL_beta, rev_value_matches_reference) {
  using stan::math::matrix_cl;
  using stan::math::var_value;
  const auto& rows = stan::test::beta_large_reference;
  var_value<matrix_cl<double>> a{
      matrix_cl<double>(beta_reference_test::column(0))};
  var_value<matrix_cl<double>> b{
      matrix_cl<double>(beta_reference_test::column(1))};
  var_value<matrix_cl<double>> res = stan::math::beta(a, b);
  Eigen::VectorXd val = stan::math::from_matrix_cl<Eigen::VectorXd>(res.val());
  for (size_t i = 0; i < rows.size(); ++i) {
    SCOPED_TRACE(::testing::Message()
                 << "a = " << rows[i][0] << ", b = " << rows[i][1]);
    EXPECT_NEAR(rows[i][2], val[i], 1e-13 * rows[i][2]);
  }
  stan::math::recover_memory();
}

TEST(OpenCL_beta, limits) {
  using stan::math::matrix_cl;
  const double inf = std::numeric_limits<double>::infinity();
  const double nan = std::numeric_limits<double>::quiet_NaN();
  Eigen::VectorXd a(9);
  Eigen::VectorXd b(9);
  a << inf, 1.65, inf, 1e-300, 0.0, 0.0, 0.0, inf, nan;
  b << 2.16, inf, inf, inf, 2.16, 0.0, inf, 0.0, inf;
  matrix_cl<double> a_cl(a);
  matrix_cl<double> b_cl(b);
  matrix_cl<double> res_cl = stan::math::beta(a_cl, b_cl);
  Eigen::VectorXd res = stan::math::from_matrix_cl<Eigen::VectorXd>(res_cl);
  for (int i = 0; i < 4; ++i) {
    EXPECT_EQ(0.0, res[i]) << "i = " << i;
  }
  for (int i = 4; i < 8; ++i) {
    EXPECT_EQ(inf, res[i]) << "i = " << i;
  }
  EXPECT_TRUE(std::isnan(res[8]));
}
#endif
