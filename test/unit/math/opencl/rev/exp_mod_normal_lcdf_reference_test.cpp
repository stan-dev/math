#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>
#include <test/unit/math/rev/prob/exp_mod_normal_lcdf_reference.hpp>
#include <cmath>
#include <vector>

namespace exp_mod_normal_lcdf_reference_test {

void expect_rel(double expected, double actual) {
  EXPECT_NEAR(expected, actual, 1e-12 * std::fabs(expected));
}

// All reference rows as one call on the device, every argument a vector, so
// each element of each adjoint is one partial of one row. Upper selects the
// log ccdf.
template <bool Upper>
void expect_device_matches_reference() {
  using stan::math::matrix_cl;
  using stan::math::var;
  using stan::math::var_value;
  const auto& rows = stan::test::exp_mod_normal_lcdf_reference;
  const int n = rows.size();
  Eigen::VectorXd args[4];
  for (auto& a : args) {
    a.resize(n);
  }
  double value = 0;
  double value_scale = 0;
  for (int i = 0; i < n; ++i) {
    for (int k = 0; k < 4; ++k) {
      args[k][i] = rows[i][k];
    }
    value += rows[i][Upper ? 9 : 4];
    value_scale += std::fabs(rows[i][Upper ? 9 : 4]);
  }
  var_value<matrix_cl<double>> y{matrix_cl<double>(args[0])};
  var_value<matrix_cl<double>> mu{matrix_cl<double>(args[1])};
  var_value<matrix_cl<double>> sigma{matrix_cl<double>(args[2])};
  var_value<matrix_cl<double>> lambda{matrix_cl<double>(args[3])};
  var lp = Upper ? stan::math::exp_mod_normal_lccdf(y, mu, sigma, lambda)
                 : stan::math::exp_mod_normal_lcdf(y, mu, sigma, lambda);
  lp.grad();
  EXPECT_NEAR(value, lp.val(), 1e-12 * value_scale);
  const Eigen::VectorXd adj[4] = {stan::math::from_matrix_cl(y.adj()),
                                  stan::math::from_matrix_cl(mu.adj()),
                                  stan::math::from_matrix_cl(sigma.adj()),
                                  stan::math::from_matrix_cl(lambda.adj())};
  for (int i = 0; i < n; ++i) {
    SCOPED_TRACE(i);
    for (int k = 0; k < 4; ++k) {
      expect_rel(rows[i][(Upper ? 10 : 5) + k], adj[k][i]);
    }
  }
  stan::math::recover_memory();
}

// Each row alone, as double arguments on the device: the value only.
template <bool Upper>
void expect_device_values_match_reference() {
  using stan::math::matrix_cl;
  for (const auto& r : stan::test::exp_mod_normal_lcdf_reference) {
    Eigen::VectorXd y(1);
    y << r[0];
    matrix_cl<double> y_cl(y);
    const double lp
        = Upper ? stan::math::exp_mod_normal_lccdf(y_cl, r[1], r[2], r[3])
                : stan::math::exp_mod_normal_lcdf(y_cl, r[1], r[2], r[3]);
    expect_rel(r[Upper ? 9 : 4], lp);
  }
}

TEST(ProbDistributionsExpModNormalLcdf, opencl_matches_reference) {
  expect_device_matches_reference<false>();
  expect_device_values_match_reference<false>();
}

TEST(ProbDistributionsExpModNormalLccdf, opencl_matches_reference) {
  expect_device_matches_reference<true>();
  expect_device_values_match_reference<true>();
}

}  // namespace exp_mod_normal_lcdf_reference_test
#endif
