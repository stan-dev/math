#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>
#include <algorithm>
#include <cmath>
#include <vector>

TEST(ProbDistributionsStudentT, error_checking) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.3, -0.8, 1.0;
  Eigen::VectorXd y_size(N - 1);
  y_size << 0.3, 0.8;
  Eigen::VectorXd y_value(N);
  y_value << 0.3, NAN, 0.5;

  Eigen::VectorXd nu(N);
  nu << 0.3, 0.8, 1.5;
  Eigen::VectorXd nu_size(N - 1);
  nu_size << 0.3, 0.8;
  Eigen::VectorXd nu_value1(N);
  nu_value1 << 0.3, INFINITY, 0.5;
  Eigen::VectorXd nu_value2(N);
  nu_value2 << 0.3, -0.6, 0.5;

  Eigen::VectorXd mu(N);
  mu << 0.3, 0.8, 1.0;
  Eigen::VectorXd mu_size(N - 1);
  mu_size << 0.3, 0.8;
  Eigen::VectorXd mu_value(N);
  mu_value << 0.3, INFINITY, 0.5;

  Eigen::VectorXd sigma(N);
  sigma << 0.3, 0.8, 1.0;
  Eigen::VectorXd sigma_size(N - 1);
  sigma_size << 0.3, 0.8;
  Eigen::VectorXd sigma_value1(N);
  sigma_value1 << 0.3, -0.4, 0.5;
  Eigen::VectorXd sigma_value2(N);
  sigma_value2 << 0.3, INFINITY, 0.5;

  stan::math::matrix_cl<double> y_cl(y);
  stan::math::matrix_cl<double> y_size_cl(y_size);
  stan::math::matrix_cl<double> y_value_cl(y_value);
  stan::math::matrix_cl<double> nu_cl(nu);
  stan::math::matrix_cl<double> nu_size_cl(nu_size);
  stan::math::matrix_cl<double> nu_value1_cl(nu_value1);
  stan::math::matrix_cl<double> nu_value2_cl(nu_value2);
  stan::math::matrix_cl<double> mu_cl(mu);
  stan::math::matrix_cl<double> mu_size_cl(mu_size);
  stan::math::matrix_cl<double> mu_value_cl(mu_value);
  stan::math::matrix_cl<double> sigma_cl(sigma);
  stan::math::matrix_cl<double> sigma_size_cl(sigma_size);
  stan::math::matrix_cl<double> sigma_value1_cl(sigma_value1);
  stan::math::matrix_cl<double> sigma_value2_cl(sigma_value2);

  EXPECT_NO_THROW(stan::math::student_t_lpdf(y_cl, nu_cl, mu_cl, sigma_cl));

  EXPECT_THROW(stan::math::student_t_lpdf(y_size_cl, nu_cl, mu_cl, sigma_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::student_t_lpdf(y_cl, nu_size_cl, mu_cl, sigma_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::student_t_lpdf(y_cl, nu_cl, mu_size_cl, sigma_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::student_t_lpdf(y_cl, nu_cl, mu_cl, sigma_size_cl),
               std::invalid_argument);

  EXPECT_THROW(stan::math::student_t_lpdf(y_value_cl, nu_cl, mu_cl, sigma_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::student_t_lpdf(y_cl, nu_value1_cl, mu_cl, sigma_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::student_t_lpdf(y_cl, nu_value2_cl, mu_cl, sigma_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::student_t_lpdf(y_cl, nu_cl, mu_value_cl, sigma_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::student_t_lpdf(y_cl, nu_cl, mu_cl, sigma_value1_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::student_t_lpdf(y_cl, nu_cl, mu_cl, sigma_value2_cl),
               std::domain_error);
}

auto student_t_lpdf_functor
    = [](const auto& y, const auto& nu, const auto& mu, const auto& sigma) {
        return stan::math::student_t_lpdf(y, nu, mu, sigma);
      };
auto student_t_lpdf_functor_propto
    = [](const auto& y, const auto& nu, const auto& mu, const auto& sigma) {
        return stan::math::student_t_lpdf<true>(y, nu, mu, sigma);
      };

TEST(ProbDistributionsStudentT, opencl_matches_cpu_small) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.3, -0.8, 1.0;
  Eigen::VectorXd nu(N);
  nu << 0.3, 0.3, 1.5;
  Eigen::VectorXd mu(N);
  mu << 0.3, 0.8, -1.0;
  Eigen::VectorXd sigma(N);
  sigma << 0.3, 0.8, 1.0;

  stan::math::test::compare_cpu_opencl_prim_rev(student_t_lpdf_functor, y, nu,
                                                mu, sigma);
  stan::math::test::compare_cpu_opencl_prim_rev(student_t_lpdf_functor_propto,
                                                y, nu, mu, sigma);
  stan::math::test::compare_cpu_opencl_prim_rev(
      student_t_lpdf_functor, y.transpose().eval(), nu.transpose().eval(),
      mu.transpose().eval(), sigma.transpose().eval());
  stan::math::test::compare_cpu_opencl_prim_rev(
      student_t_lpdf_functor_propto, y.transpose().eval(),
      nu.transpose().eval(), mu.transpose().eval(), sigma.transpose().eval());
}

TEST(ProbDistributionsStudentT, opencl_broadcast_y) {
  int N = 3;

  double y_scal = 12.3;
  Eigen::VectorXd nu(N);
  nu << 0.5, 1.2, 1.0;
  Eigen::VectorXd mu(N);
  mu << 0.3, 0.8, -1.0;
  Eigen::VectorXd sigma(N);
  sigma << 0.3, 0.8, 1.0;

  stan::math::test::test_opencl_broadcasting_prim_rev<0>(student_t_lpdf_functor,
                                                         y_scal, nu, mu, sigma);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      student_t_lpdf_functor_propto, y_scal, nu, mu, sigma);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      student_t_lpdf_functor, y_scal, nu.transpose().eval(),
      mu.transpose().eval(), sigma);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      student_t_lpdf_functor_propto, y_scal, nu, mu.transpose().eval(),
      sigma.transpose().eval());
}

TEST(ProbDistributionsStudentT, opencl_matches_cpu_big) {
  int N = 153;

  Eigen::Matrix<double, Eigen::Dynamic, 1> y
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1);
  Eigen::Matrix<double, Eigen::Dynamic, 1> nu
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();
  Eigen::Matrix<double, Eigen::Dynamic, 1> mu
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1);
  Eigen::Matrix<double, Eigen::Dynamic, 1> sigma
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();

  stan::math::test::compare_cpu_opencl_prim_rev(student_t_lpdf_functor, y, nu,
                                                mu, sigma);
  stan::math::test::compare_cpu_opencl_prim_rev(student_t_lpdf_functor_propto,
                                                y, nu, mu, sigma);
  stan::math::test::compare_cpu_opencl_prim_rev(
      student_t_lpdf_functor, y.transpose().eval(), nu.transpose().eval(),
      mu.transpose().eval(), sigma);
  stan::math::test::compare_cpu_opencl_prim_rev(
      student_t_lpdf_functor_propto, y.transpose().eval(),
      nu.transpose().eval(), mu.transpose().eval(), sigma.transpose().eval());
}

TEST(ProbDistributionsStudentT, opencl_large_shapes_reference) {
  // The tests above compare OpenCL against the CPU and cannot see an error
  // that both share. These are absolute references from mpmath, as in
  // test/unit/math/rev/prob/large_shapes_test.cpp. The differences
  // lgamma(nu/2 + 1/2) - lgamma(nu/2) and digamma(nu/2 + 1/2) -
  // digamma(nu/2) keep no correct digits for large nu. The nu and sigma
  // partials are compared on the log scale. The last row is a small nu
  // where the plain differences are correct.
  using stan::math::var;
  struct TestValue {
    double y;
    double nu;
    double mu;
    double sigma;
    double value;
    double d_y;
    double d_nu;
    double d_mu;
    double d_sigma;
  };
  const std::vector<TestValue> test_values = {
      {0x1.0000000000000p+0, 0x1.6bcc41e900000p+46, 0x0.0p+0,
       0x1.0000000000000p+0, -1.4189385332046777, -1.0000000000000000,
       4.9999999999999833e-29, 1.0000000000000000, 2.6727647100921956e-51},
      {-0x1.cd1e504efb30cp+1, 0x1.7cc4b890abebfp+45, 0x1.64703afcf3380p-4,
       0x1.c43477d4ae376p+3, -3.6014210338578538, 1.8475570774651872e-2,
       1.0330563009181428e-28, -1.8475570774651872e-2, -6.5940664388226970e-2},
      {0x1.8000000000000p+1, 0x1.6bcc41e900000p+46, 0x0.0p+0,
       0x1.0000000000000p+1, -2.7370857137646191, -7.4999999999999062e-1,
       1.0937500000001266e-29, 7.4999999999999062e-1, 6.2499999999998594e-1},
      {-0x1.c8dd60359f470p-1, 0x1.da533b967d6d3p-4, -0x1.572a29ad6231cp-1,
       0x1.9d9944276e11dp-3, -1.6062891011707614, 4.5853946438927804,
       6.8870497923078714, -4.5853946438927804, 9.0518700772802596e-2},
  };
  auto near = [](double got, double expected, double tol) {
    return std::fabs(got - expected)
           <= tol * std::max(1.0, std::fabs(expected));
  };
  for (const auto& t : test_values) {
    Eigen::Matrix<var, Eigen::Dynamic, 1> y(1), nu(1), mu(1), sigma(1);
    y << t.y;
    nu << t.nu;
    mu << t.mu;
    sigma << t.sigma;
    auto y_cl = stan::math::to_matrix_cl(y);
    auto nu_cl = stan::math::to_matrix_cl(nu);
    auto mu_cl = stan::math::to_matrix_cl(mu);
    auto sigma_cl = stan::math::to_matrix_cl(sigma);
    var lp = stan::math::student_t_lpdf(y_cl, nu_cl, mu_cl, sigma_cl);
    lp.grad();
    EXPECT_TRUE(near(lp.val(), t.value, 1e-12)) << "value, nu = " << t.nu;
    EXPECT_TRUE(near(y(0).adj(), t.d_y, 1e-11)) << "d_y, nu = " << t.nu;
    EXPECT_TRUE(near(t.nu * nu(0).adj(), t.nu * t.d_nu, 1e-11))
        << "nu * d_nu, nu = " << t.nu;
    EXPECT_TRUE(near(mu(0).adj(), t.d_mu, 1e-11)) << "d_mu, nu = " << t.nu;
    EXPECT_TRUE(near(t.sigma * sigma(0).adj(), t.sigma * t.d_sigma, 1e-11))
        << "sigma * d_sigma, nu = " << t.nu;
    stan::math::recover_memory();
  }
}

#endif
