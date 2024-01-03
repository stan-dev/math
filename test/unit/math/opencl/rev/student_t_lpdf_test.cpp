#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>
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

#endif
