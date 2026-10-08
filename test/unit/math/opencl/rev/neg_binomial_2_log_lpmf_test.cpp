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

TEST(ProbDistributionsNegBinomial2Log, opencl_large_shapes_reference) {
  // The tests above compare OpenCL against the CPU and cannot see an error
  // that both share. These are absolute references from mpmath at 90 digits,
  // checked against 140 digits, as in
  // test/unit/math/rev/prob/large_shapes_test.cpp. The difference
  // digamma(n + phi) - digamma(phi) keeps no correct digits for large phi.
  // The phi partial is compared as phi * d/dphi. At phi = 1e15 the rounding
  // errors of the develop OpenCL formula cancel to 0 and it passes by
  // chance; the rows at phi = 1e10 and 1e12 show its error (3.6e-5 and
  // 1.4e-9 in phi * d/dphi). The last row is a small phi where the plain
  // difference is correct.
  using stan::math::var;
  struct TestValue {
    int n;
    double eta;
    double phi;
    double value;
    double d_eta;
    double d_phi;
  };
  const std::vector<TestValue> test_values = {
      {3, 0x1.193ea7aad030bp+0, 0x1.c6bf526340000p+49, -1.4959226032237274,
       -2.7213891705004510e-16, 1.4999999999999960e-30},
      {5, 0x1.f4bd2b7ac1bafp+1, 0x1.c6bf526340000p+49, -3.5227376715640301e+1,
       -4.4999999999997745e+1, -1.0099999999999289e-27},
      {2, -0x1.6d3c324e13f50p-2, 0x1.37807ed5e8000p+50, -2.1064970684374103,
       1.2999999999999994, 8.2582982577654707e-32},
      {3, 0x1.193ea7aad030bp+0, 0x1.2a05f20000000p+33, -1.4959226033737259,
       -2.7213891696840424e-16, 1.4999999996e-20},
      {5, 0x1.f4bd2b7ac1bafp+1, 0x1.2a05f20000000p+33, -3.5227376614641312e+1,
       -4.4999999774999996e+1, -1.0099999929136665e-17},
      {57, 0x1.193ea7aad030bp+0, 0x1.d1a94a2000000p+39, -1.1677494795148559e+2,
       5.3999999999838e+1, -1.429499999940379e-21},
      {0, 0x1.193ea7aad030bp+0, 0x1.9000000000000p+6, -2.9558802241544405,
       -2.9126213592233012, -4.3258864931139310e-4},
  };
  auto expect_reference = [](const TestValue& t, const var& lp, double eta_adj,
                             double phi_adj, const char* signature) {
    EXPECT_NEAR(lp.val(), t.value, 1e-12 * std::max(1.0, std::fabs(t.value)))
        << signature << ": n = " << t.n << ", eta = " << t.eta
        << ", phi = " << t.phi;
    EXPECT_NEAR(eta_adj, t.d_eta, 1e-11 * std::max(1.0, std::fabs(t.d_eta)))
        << signature << ": n = " << t.n << ", eta = " << t.eta
        << ", phi = " << t.phi;
    const double gp = t.phi * t.d_phi;
    EXPECT_NEAR(t.phi * phi_adj, gp, 1e-11 * std::max(1.0, std::fabs(gp)))
        << signature << ": n = " << t.n << ", eta = " << t.eta
        << ", phi = " << t.phi;
  };
  for (const auto& t : test_values) {
    const std::vector<int> n{t.n};
    stan::math::matrix_cl<int> n_cl(n);

    Eigen::Matrix<var, Eigen::Dynamic, 1> eta(1);
    eta << t.eta;
    Eigen::Matrix<var, Eigen::Dynamic, 1> phi(1);
    phi << t.phi;
    auto eta_cl = stan::math::to_matrix_cl(eta);
    auto phi_cl = stan::math::to_matrix_cl(phi);
    var lp = stan::math::neg_binomial_2_log_lpmf(n_cl, eta_cl, phi_cl);
    lp.grad();
    expect_reference(t, lp, eta(0).adj(), phi(0).adj(), "vector");
    stan::math::recover_memory();

    var eta_s = t.eta;
    var phi_s = t.phi;
    var lp_s = stan::math::neg_binomial_2_log_lpmf(n_cl, eta_s, phi_s);
    lp_s.grad();
    expect_reference(t, lp_s, eta_s.adj(), phi_s.adj(), "scalar");
    stan::math::recover_memory();
  }
}

#endif
