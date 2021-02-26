#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>
#include <vector>

auto student_t_lpdf2_functor
    = [](const auto& y, const auto& nu, const auto& mu, const auto& sigma) {
        return stan::math::student_t_lpdf(y, nu, mu, sigma);
      };
auto student_t_lpdf2_functor_propto
    = [](const auto& y, const auto& nu, const auto& mu, const auto& sigma) {
        return stan::math::student_t_lpdf<true>(y, nu, mu, sigma);
      };

TEST(ProbDistributionsStudentT, opencl_broadcast_nu) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.3, -0.8, 1.0;
  double nu_scal = 12.3;
  Eigen::VectorXd mu(N);
  mu << 0.3, 0.8, -1.0;
  Eigen::VectorXd sigma(N);
  sigma << 0.3, 0.8, 1.0;

  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      student_t_lpdf2_functor, y, nu_scal, mu, sigma);
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      student_t_lpdf2_functor_propto, y, nu_scal, mu, sigma);
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      student_t_lpdf2_functor, y, nu_scal, mu.transpose().eval(),
      sigma.transpose().eval());
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      student_t_lpdf2_functor_propto, y.transpose().eval(), nu_scal, mu,
      sigma.transpose().eval());
}

TEST(ProbDistributionsStudentT, opencl_broadcast_mu) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.3, -0.8, 1.0;
  Eigen::VectorXd nu(N);
  nu << 0.3, 0.3, 1.5;
  double mu_scal = 12.3;
  Eigen::VectorXd sigma(N);
  sigma << 0.3, 0.8, 1.0;

  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      student_t_lpdf2_functor, y, nu, mu_scal, sigma);
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      student_t_lpdf2_functor_propto, y, nu, mu_scal, sigma);
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      student_t_lpdf2_functor, y.transpose().eval(), nu, mu_scal,
      sigma.transpose().eval());
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      student_t_lpdf2_functor_propto, y.transpose().eval(),
      nu.transpose().eval(), mu_scal, sigma);
}

TEST(ProbDistributionsStudentT, opencl_broadcast_sigma) {
  int N = 3;

  Eigen::VectorXd y(N);
  y << 0.3, -0.8, 1.0;
  Eigen::VectorXd nu(N);
  nu << 0.3, 0.3, 1.5;
  Eigen::VectorXd mu(N);
  mu << 0.3, 0.8, -1.0;
  double sigma_scal = 12.3;

  stan::math::test::test_opencl_broadcasting_prim_rev<3>(
      student_t_lpdf2_functor, y, nu, mu, sigma_scal);
  stan::math::test::test_opencl_broadcasting_prim_rev<3>(
      student_t_lpdf2_functor_propto, y, nu, mu, sigma_scal);
  stan::math::test::test_opencl_broadcasting_prim_rev<3>(
      student_t_lpdf2_functor, y.transpose().eval(), nu.transpose().eval(), mu,
      sigma_scal);
  stan::math::test::test_opencl_broadcasting_prim_rev<3>(
      student_t_lpdf2_functor_propto, y, nu.transpose().eval(),
      mu.transpose().eval(), sigma_scal);
}

#endif
