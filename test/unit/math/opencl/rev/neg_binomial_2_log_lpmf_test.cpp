#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>
#include <vector>

TEST(ProbDistributionsNegBinomial2Log, error_checking) {
  int N = 3;

  std::vector<int> n{1, 0, 12};
  std::vector<int> n_size{1, 0, 1, 0};
  std::vector<int> n_value{0, 1, -3};

  Eigen::VectorXd eta(N);
  eta << 0.3, 0.8, -1.3;
  Eigen::VectorXd eta_size(N - 1);
  eta_size << 0.3, 0.8;
  Eigen::VectorXd eta_value(N);
  eta_value << 0.3, INFINITY, 0.5;

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
  stan::math::matrix_cl<double> eta_cl(eta);
  stan::math::matrix_cl<double> eta_size_cl(eta_size);
  stan::math::matrix_cl<double> eta_value_cl(eta_value);
  stan::math::matrix_cl<double> phi_cl(phi);
  stan::math::matrix_cl<double> phi_size_cl(phi_size);
  stan::math::matrix_cl<double> phi_value1_cl(phi_value1);
  stan::math::matrix_cl<double> phi_value2_cl(phi_value2);

  EXPECT_NO_THROW(stan::math::neg_binomial_2_log_lpmf(n_cl, eta_cl, phi_cl));

  EXPECT_THROW(stan::math::neg_binomial_2_log_lpmf(n_size_cl, eta_cl, phi_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::neg_binomial_2_log_lpmf(n_cl, eta_size_cl, phi_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::neg_binomial_2_log_lpmf(n_cl, eta_cl, phi_size_cl),
               std::invalid_argument);

  EXPECT_THROW(stan::math::neg_binomial_2_log_lpmf(n_value_cl, eta_cl, phi_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::neg_binomial_2_log_lpmf(n_cl, eta_value_cl, phi_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::neg_binomial_2_log_lpmf(n_cl, eta_cl, phi_value1_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::neg_binomial_2_log_lpmf(n_cl, eta_cl, phi_value2_cl),
               std::domain_error);
}

auto neg_binomial_2_log_lpmf_functor
    = [](const auto& n, const auto& eta, const auto& phi) {
        return stan::math::neg_binomial_2_log_lpmf(n, eta, phi);
      };
auto neg_binomial_2_log_lpmf_functor_propto
    = [](const auto& n, const auto& eta, const auto& phi) {
        return stan::math::neg_binomial_2_log_lpmf<true>(n, eta, phi);
      };

TEST(ProbDistributionsNegBinomial2Log, opencl_matches_cpu_small) {
  int N = 3;
  int M = 2;

  std::vector<int> n{1, 0, 12};
  Eigen::VectorXd eta(N);
  eta << 0.3, 0.8, -1.3;
  Eigen::VectorXd phi(N);
  phi << 0.3, 0.8, 1.3;

  stan::math::test::compare_cpu_opencl_prim_rev(neg_binomial_2_log_lpmf_functor,
                                                n, eta, phi);
  stan::math::test::compare_cpu_opencl_prim_rev(
      neg_binomial_2_log_lpmf_functor_propto, n, eta, phi);
  stan::math::test::compare_cpu_opencl_prim_rev(neg_binomial_2_log_lpmf_functor,
                                                n, eta.transpose().eval(),
                                                phi.transpose().eval());
  stan::math::test::compare_cpu_opencl_prim_rev(
      neg_binomial_2_log_lpmf_functor_propto, n, eta.transpose().eval(),
      phi.transpose().eval());
}

TEST(ProbDistributionsNegBinomial2Log, opencl_broadcast_n) {
  int N = 3;

  int n = 2;
  Eigen::VectorXd eta(N);
  eta << 0.3, 0.8, 1.0;
  Eigen::VectorXd phi(N);
  phi << 0.3, 0.8, 1.3;

  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      neg_binomial_2_log_lpmf_functor, n, eta, phi);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      neg_binomial_2_log_lpmf_functor_propto, n, eta, phi);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      neg_binomial_2_log_lpmf_functor, n, eta.transpose().eval(), phi);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      neg_binomial_2_log_lpmf_functor_propto, n, eta, phi.transpose().eval());
}

TEST(ProbDistributionsNegBinomial2Log, opencl_broadcast_eta) {
  int N = 3;

  std::vector<int> n{1, 0, 12};
  double eta = 0.4;
  Eigen::VectorXd phi(N);
  phi << 0.3, 0.8, 1.3;

  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      neg_binomial_2_log_lpmf_functor, n, eta, phi);
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      neg_binomial_2_log_lpmf_functor_propto, n, eta, phi);
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      neg_binomial_2_log_lpmf_functor, n, eta, phi.transpose().eval());
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      neg_binomial_2_log_lpmf_functor_propto, n, eta, phi.transpose().eval());
}

TEST(ProbDistributionsNegBinomial2Log, opencl_broadcast_phi) {
  int N = 3;

  std::vector<int> n{1, 0, 12};
  Eigen::VectorXd eta(N);
  eta << 0.3, 0.8, 1.0;
  double phi = 0.4;

  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      neg_binomial_2_log_lpmf_functor, n, eta, phi);
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      neg_binomial_2_log_lpmf_functor_propto, n, eta, phi);
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      neg_binomial_2_log_lpmf_functor, n, eta.transpose().eval(), phi);
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      neg_binomial_2_log_lpmf_functor_propto, n, eta.transpose().eval(), phi);
}

TEST(ProbDistributionsNegBinomial2Log, opencl_matches_cpu_big) {
  int N = 153;

  std::vector<int> n(N);
  for (int i = 0; i < N; i++) {
    n[i] = Eigen::Array<int, Eigen::Dynamic, 1>::Random(1, 1).abs()(0) % 1000;
  }
  Eigen::Matrix<double, Eigen::Dynamic, 1> eta
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1);
  Eigen::Matrix<double, Eigen::Dynamic, 1> phi
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();

  stan::math::test::compare_cpu_opencl_prim_rev(neg_binomial_2_log_lpmf_functor,
                                                n, eta, phi);
  stan::math::test::compare_cpu_opencl_prim_rev(
      neg_binomial_2_log_lpmf_functor_propto, n, eta, phi);
  stan::math::test::compare_cpu_opencl_prim_rev(neg_binomial_2_log_lpmf_functor,
                                                n, eta.transpose().eval(),
                                                phi.transpose().eval());
  stan::math::test::compare_cpu_opencl_prim_rev(
      neg_binomial_2_log_lpmf_functor_propto, n, eta.transpose().eval(),
      phi.transpose().eval());
}

TEST(ProbDistributionsNegBinomial2Log, opencl_matches_cpu_eta_phi_scalar) {
  int N = 3;
  int M = 2;

  std::vector<int> n{1, 0, 12};
  double eta = 0.3;
  double phi = 0.8;

  stan::math::test::compare_cpu_opencl_prim_rev(neg_binomial_2_log_lpmf_functor,
                                                n, eta, phi);
  stan::math::test::compare_cpu_opencl_prim_rev(
      neg_binomial_2_log_lpmf_functor_propto, n, eta, phi);
}

#endif
