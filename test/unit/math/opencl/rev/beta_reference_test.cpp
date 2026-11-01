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

TEST(OpenCL_beta, edge_values_match_cpu) {
  using stan::math::matrix_cl;
  const double inf = std::numeric_limits<double>::infinity();
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double edges[] = {nan, 0.0, -0.5, -2.5, inf, 1e-300, 0.5, 2.0};
  const int m = sizeof(edges) / sizeof(edges[0]);
  Eigen::VectorXd a(m * m);
  Eigen::VectorXd b(m * m);
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < m; ++j) {
      a[i * m + j] = edges[i];
      b[i * m + j] = edges[j];
    }
  }
  matrix_cl<double> a_cl(a);
  matrix_cl<double> b_cl(b);
  matrix_cl<double> res_cl = stan::math::beta(a_cl, b_cl);
  Eigen::VectorXd res = stan::math::from_matrix_cl<Eigen::VectorXd>(res_cl);
  for (int k = 0; k < m * m; ++k) {
    SCOPED_TRACE(::testing::Message() << "a = " << a[k] << ", b = " << b[k]);
    const double cpu = stan::math::beta(a[k], b[k]);
    if (std::isnan(cpu)) {
      EXPECT_TRUE(std::isnan(res[k]));
    } else if (std::isinf(cpu) || cpu == 0) {
      EXPECT_EQ(cpu, res[k]);
    } else {
      EXPECT_NEAR(cpu, res[k], 1e-11 * std::fabs(cpu));
    }
  }
}
#endif
