#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>
#include <vector>

TEST(muProbDistributionsNegBinomial2, error_checking) {
  int N = 3;

  std::vector<int> n{1, 0, 12};
  std::vector<int> n_size{1, 0, 1, 0};
  std::vector<int> n_value{0, 1, -3};

  Eigen::VectorXd mu(N);
  mu << 0.3, 0.8, 1.3;
  Eigen::VectorXd mu_size(N - 1);
  mu_size << 0.3, 0.8;
  Eigen::VectorXd mu_value1(N);
  mu_value1 << 0.3, -0.3, 0.5;
  Eigen::VectorXd mu_value2(N);
  mu_value2 << 0.3, INFINITY, 0.5;

  Eigen::VectorXd phi(N);
  phi << 0.3, 0.8, 1.3;
  Eigen::VectorXd phi_size(N - 1);
  phi_size << 0.3, 0.8;
  Eigen::VectorXd phi_value1(N);
  phi_value1 << 0.3, -0.8, 0.5;
  Eigen::VectorXd phi_value2(N);
  phi_value2 << 0.3, INFINITY, 0.5;

  stan::math::matrix_cl<int> n_cl(n);
  stan::math::matrix_cl<int> n_size_cl(n_size);
  stan::math::matrix_cl<int> n_value_cl(n_value);
  stan::math::matrix_cl<double> mu_cl(mu);
  stan::math::matrix_cl<double> mu_size_cl(mu_size);
  stan::math::matrix_cl<double> mu_value1_cl(mu_value1);
  stan::math::matrix_cl<double> mu_value2_cl(mu_value2);
  stan::math::matrix_cl<double> phi_cl(phi);
  stan::math::matrix_cl<double> phi_size_cl(phi_size);
  stan::math::matrix_cl<double> phi_value1_cl(phi_value1);
  stan::math::matrix_cl<double> phi_value2_cl(phi_value2);

  EXPECT_NO_THROW(stan::math::neg_binomial_2_lpmf(n_cl, mu_cl, phi_cl));

  EXPECT_THROW(stan::math::neg_binomial_2_lpmf(n_size_cl, mu_cl, phi_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::neg_binomial_2_lpmf(n_cl, mu_size_cl, phi_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::neg_binomial_2_lpmf(n_cl, mu_cl, phi_size_cl),
               std::invalid_argument);

  EXPECT_THROW(stan::math::neg_binomial_2_lpmf(n_value_cl, mu_cl, phi_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::neg_binomial_2_lpmf(n_cl, mu_value1_cl, phi_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::neg_binomial_2_lpmf(n_cl, mu_value2_cl, phi_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::neg_binomial_2_lpmf(n_cl, mu_cl, phi_value1_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::neg_binomial_2_lpmf(n_cl, mu_cl, phi_value2_cl),
               std::domain_error);
}

auto neg_binomial_2_lpmf_functor
    = [](const auto& n, const auto& mu, const auto& phi) {
        return stan::math::neg_binomial_2_lpmf(n, mu, phi);
      };
auto neg_binomial_2_lpmf_functor_propto
    = [](const auto& n, const auto& mu, const auto& phi) {
        return stan::math::neg_binomial_2_lpmf<true>(n, mu, phi);
      };

TEST(muProbDistributionsNegBinomial2, opencl_matches_cpu_small) {
  int N = 3;
  int M = 2;

  std::vector<int> n{1, 0, 12};
  Eigen::VectorXd mu(N);
  mu << 0.3, 0.5, 1.8;
  Eigen::VectorXd phi(N);
  phi << 0.3, 0.8, 1.3;

  stan::math::test::compare_cpu_opencl_prim_rev(neg_binomial_2_lpmf_functor, n,
                                                mu, phi);
  stan::math::test::compare_cpu_opencl_prim_rev(
      neg_binomial_2_lpmf_functor_propto, n, mu, phi);
  stan::math::test::compare_cpu_opencl_prim_rev(neg_binomial_2_lpmf_functor, n,
                                                mu.transpose().eval(),
                                                phi.transpose().eval());
  stan::math::test::compare_cpu_opencl_prim_rev(
      neg_binomial_2_lpmf_functor_propto, n, mu.transpose().eval(),
      phi.transpose().eval());
}

TEST(muProbDistributionsNegBinomial2, opencl_broadcast_n) {
  int N = 3;

  int n = 2;
  Eigen::VectorXd mu(N);
  mu << 0.3, 0.5, 1.8;
  Eigen::VectorXd phi(N);
  phi << 0.3, 0.8, 1.3;

  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      neg_binomial_2_lpmf_functor, n, mu, phi);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      neg_binomial_2_lpmf_functor_propto, n, mu, phi);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      neg_binomial_2_lpmf_functor, n, mu.transpose().eval(), phi);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      neg_binomial_2_lpmf_functor_propto, n, mu, phi.transpose().eval());
}

TEST(muProbDistributionsNegBinomial2, opencl_broadcast_mu) {
  int N = 3;

  std::vector<int> n{1, 0, 12};
  double mu = 0.4;
  Eigen::VectorXd phi(N);
  phi << 0.3, 0.8, 1.3;

  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      neg_binomial_2_lpmf_functor, n, mu, phi);
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      neg_binomial_2_lpmf_functor_propto, n, mu, phi);
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      neg_binomial_2_lpmf_functor, n, mu, phi.transpose().eval());
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      neg_binomial_2_lpmf_functor_propto, n, mu, phi.transpose().eval());
}

TEST(muProbDistributionsNegBinomial2, opencl_broadcast_phi) {
  int N = 3;

  std::vector<int> n{1, 0, 12};
  Eigen::VectorXd mu(N);
  mu << 0.3, 0.5, 1.8;
  double phi = 0.4;

  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      neg_binomial_2_lpmf_functor, n, mu, phi);
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      neg_binomial_2_lpmf_functor_propto, n, mu, phi);
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      neg_binomial_2_lpmf_functor, n, mu.transpose().eval(), phi);
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      neg_binomial_2_lpmf_functor_propto, n, mu.transpose().eval(), phi);
}

TEST(muProbDistributionsNegBinomial2, opencl_matches_cpu_big) {
  int N = 153;

  std::vector<int> n(N);
  for (int i = 0; i < N; i++) {
    n[i] = Eigen::Array<int, Eigen::Dynamic, 1>::Random(1, 1).abs()(0) % 1000;
  }
  Eigen::Matrix<double, Eigen::Dynamic, 1> mu
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();
  Eigen::Matrix<double, Eigen::Dynamic, 1> phi
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();

  stan::math::test::compare_cpu_opencl_prim_rev(neg_binomial_2_lpmf_functor, n,
                                                mu, phi);
  stan::math::test::compare_cpu_opencl_prim_rev(
      neg_binomial_2_lpmf_functor_propto, n, mu, phi);
  stan::math::test::compare_cpu_opencl_prim_rev(neg_binomial_2_lpmf_functor, n,
                                                mu.transpose().eval(),
                                                phi.transpose().eval());
  stan::math::test::compare_cpu_opencl_prim_rev(
      neg_binomial_2_lpmf_functor_propto, n, mu.transpose().eval(),
      phi.transpose().eval());
}

TEST(ProbDistributionsNegBinomial2, opencl_scalar_n_mu) {
  int N = 3;
  int M = 2;

  int n = 1;
  double mu = 0.3;
  Eigen::VectorXd phi(N);
  phi << 0.3, 0.8, 1.3;

  stan::math::test::compare_cpu_opencl_prim_rev(neg_binomial_2_lpmf_functor, n,
                                                mu, phi);
  stan::math::test::compare_cpu_opencl_prim_rev(
      neg_binomial_2_lpmf_functor_propto, n, mu, phi);
}

#endif
