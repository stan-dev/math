#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <test/unit/math/rev/util.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <type_traits>

using Eigen::VectorXd;
using stan::math::matrix_cl;
using stan::math::var;
using stan::math::var_value;
using stan::math::vector_v;
using stan::math::opencl::ScalarCl;
using stan::math::opencl::to_host;

TEST(ScalarClDistributions, normal_lpdf_return_types) {
  matrix_cl<double> y_cl(VectorXd::Random(5));
  var_value<matrix_cl<double>> y_v{matrix_cl<double>(VectorXd::Random(5))};
  EXPECT_TRUE((std::is_same<decltype(stan::math::normal_lpdf(y_cl, 0.5, 2.0)),
                            ScalarCl<double>>::value));
  EXPECT_TRUE((std::is_same<decltype(stan::math::normal_lpdf(
                                y_cl, ScalarCl<double>(0.5), 2.0)),
                            ScalarCl<double>>::value));
  EXPECT_TRUE((std::is_same<decltype(stan::math::normal_lpdf(y_v, 0.5, 2.0)),
                            ScalarCl<var>>::value));
  EXPECT_TRUE((std::is_same<decltype(stan::math::normal_lpdf(
                                y_cl, ScalarCl<var>(0.5), 2.0)),
                            ScalarCl<var>>::value));
  EXPECT_TRUE(
      (std::is_same<decltype(stan::math::normal_lpdf(y_cl, var(0.5), 2.0)),
                    ScalarCl<var>>::value));
  stan::math::recover_memory();
}

TEST_F(AgradRev, ScalarClDistributions_normal_lpdf_device_scalar_params) {
  VectorXd y_val = VectorXd::Random(1000);
  const double mu_val = 0.3;
  const double sigma_val = 1.7;

  // CPU reference
  vector_v y_cpu = y_val;
  var mu_cpu = mu_val;
  var sigma_cpu = sigma_val;
  var lp_cpu = stan::math::normal_lpdf(y_cpu, mu_cpu, sigma_cpu);
  lp_cpu.grad();
  VectorXd y_grad = y_cpu.adj();
  double mu_grad = mu_cpu.adj();
  double sigma_grad = sigma_cpu.adj();
  double lp_val = lp_cpu.val();
  stan::math::recover_memory();

  // everything on the device, scalars as device vars
  {
    var_value<matrix_cl<double>> y{matrix_cl<double>(y_val)};
    ScalarCl<var> mu(mu_val);
    ScalarCl<var> sigma(sigma_val);
    ScalarCl<var> lp = stan::math::normal_lpdf(y, mu, sigma);
    var lp_host = to_host(lp);
    EXPECT_NEAR(lp_host.val(), lp_val, 1e-9);
    lp_host.grad();
    EXPECT_MATRIX_NEAR(stan::math::from_matrix_cl(y.adj()), y_grad, 1e-9);
    EXPECT_NEAR(to_host(mu.adj()), mu_grad, 1e-9);
    EXPECT_NEAR(to_host(sigma.adj()), sigma_grad, 1e-9);
    stan::math::recover_memory();
  }
  // CPU var and device double scalars
  {
    matrix_cl<double> y(y_val);
    var mu = mu_val;
    ScalarCl<double> sigma(sigma_val);
    var lp_host = to_host(stan::math::normal_lpdf(y, mu, sigma));
    EXPECT_NEAR(lp_host.val(), lp_val, 1e-9);
    lp_host.grad();
    EXPECT_NEAR(mu.adj(), mu_grad, 1e-9);
    stan::math::recover_memory();
  }
  // device double location, device var scale
  {
    matrix_cl<double> y(y_val);
    ScalarCl<double> mu(mu_val);
    ScalarCl<var> sigma(sigma_val);
    var lp_host = to_host(stan::math::normal_lpdf(y, mu, sigma));
    EXPECT_NEAR(lp_host.val(), lp_val, 1e-9);
    lp_host.grad();
    EXPECT_NEAR(to_host(sigma.adj()), sigma_grad, 1e-9);
    stan::math::recover_memory();
  }
  // primitive device result
  {
    matrix_cl<double> y(y_val);
    ScalarCl<double> lp = stan::math::normal_lpdf(y, ScalarCl<double>(mu_val),
                                                  ScalarCl<double>(sigma_val));
    EXPECT_NEAR(to_host(lp), lp_val, 1e-9);
  }
}

TEST_F(AgradRev, ScalarClDistributions_normal_lpdf_checks_device_scalars) {
  matrix_cl<double> y(VectorXd::Random(10));
  EXPECT_THROW(stan::math::normal_lpdf(y, 0.0, ScalarCl<double>(-1.0)),
               std::domain_error);
  EXPECT_THROW(
      stan::math::normal_lpdf(y, ScalarCl<var>(stan::math::INFTY), 1.0),
      std::domain_error);
}

TEST_F(AgradRev, ScalarClDistributions_normal_lpdf_propto_and_empty) {
  matrix_cl<double> y(VectorXd::Random(10));
  matrix_cl<double> empty(VectorXd(0));
  EXPECT_EQ(to_host(stan::math::normal_lpdf<true>(y, 0.0, 1.0)), 0.0);
  EXPECT_EQ(to_host(stan::math::normal_lpdf(empty, 0.0, 1.0)), 0.0);
  ScalarCl<var> mu(0.1);
  var lp = to_host(stan::math::normal_lpdf(empty, mu, 1.0));
  EXPECT_EQ(lp.val(), 0.0);
  lp.grad();
  EXPECT_EQ(to_host(mu.adj()), 0.0);
}
#endif
