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

TEST(muProbDistributionsNegBinomial2, opencl_large_shapes_reference) {
  // The tests above compare OpenCL against the CPU and cannot see an error
  // that both share. These are absolute references from mpmath at 90 digits,
  // checked against 140 digits, as in
  // test/unit/math/rev/prob/large_shapes_test.cpp. The difference
  // digamma(n + phi) - digamma(phi) keeps no correct digits for large phi.
  // The phi partial is compared as phi * d/dphi. The last row is a small
  // phi where the plain difference is correct.
  using stan::math::var;
  struct TestValue {
    int n;
    double mu;
    double phi;
    double value;
    double d_mu;
    double d_phi;
  };
  const std::vector<TestValue> test_values = {
      {57, 0x1.8000000000000p+1, 0x1.2a05f20000000p+33, -1.1677494780996510e+2,
       1.7999999994600000e+1, -1.4294999940379000e-17},
      {5, 0x1.9000000000000p+5, 0x1.2a05f20000000p+33, -3.5227376614641316e+1,
       -8.9999999550000002e-1, -1.0099999929136667e-17},
      {1, 0x1.8000000000000p+1, 0x1.2a05f20000000p+33, -1.9013877111818903,
       -6.6666666646666667e-1, -1.4999999991000000e-20},
      {0, 0x1.8000000000000p+1, 0x1.9000000000000p+6, -2.9558802241544403,
       -9.7087378640776699e-1, -4.3258864931139302e-4},
  };
  auto expect_reference = [](const TestValue& t, const var& lp, double mu_adj,
                             double phi_adj, const char* signature) {
    EXPECT_NEAR(lp.val(), t.value, 1e-12 * std::max(1.0, std::fabs(t.value)))
        << signature << ": n = " << t.n << ", mu = " << t.mu
        << ", phi = " << t.phi;
    EXPECT_NEAR(mu_adj, t.d_mu, 1e-11 * std::max(1.0, std::fabs(t.d_mu)))
        << signature << ": n = " << t.n << ", mu = " << t.mu
        << ", phi = " << t.phi;
    const double gp = t.phi * t.d_phi;
    EXPECT_NEAR(t.phi * phi_adj, gp, 1e-11 * std::max(1.0, std::fabs(gp)))
        << signature << ": n = " << t.n << ", mu = " << t.mu
        << ", phi = " << t.phi;
  };
  for (const auto& t : test_values) {
    const std::vector<int> n{t.n};
    stan::math::matrix_cl<int> n_cl(n);

    Eigen::Matrix<var, Eigen::Dynamic, 1> mu(1);
    mu << t.mu;
    Eigen::Matrix<var, Eigen::Dynamic, 1> phi(1);
    phi << t.phi;
    auto mu_cl = stan::math::to_matrix_cl(mu);
    auto phi_cl = stan::math::to_matrix_cl(phi);
    var lp = stan::math::neg_binomial_2_lpmf(n_cl, mu_cl, phi_cl);
    lp.grad();
    expect_reference(t, lp, mu(0).adj(), phi(0).adj(), "vector");
    stan::math::recover_memory();

    var mu_s = t.mu;
    var phi_s = t.phi;
    var lp_s = stan::math::neg_binomial_2_lpmf(n_cl, mu_s, phi_s);
    lp_s.grad();
    expect_reference(t, lp_s, mu_s.adj(), phi_s.adj(), "scalar");
    stan::math::recover_memory();
  }
}

#endif
