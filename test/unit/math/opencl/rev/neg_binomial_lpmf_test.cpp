#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>
#include <vector>

TEST(muProbDistributionsNegBinomial, error_checking) {
  int N = 3;

  std::vector<int> n{1, 0, 12};
  std::vector<int> n_size{1, 0, 1, 0};
  std::vector<int> n_value{0, 1, -3};

  Eigen::VectorXd alpha(N);
  alpha << 0.3, 0.8, 1.3;
  Eigen::VectorXd alpha_size(N - 1);
  alpha_size << 0.3, 0.8;
  Eigen::VectorXd alpha_value1(N);
  alpha_value1 << 0.3, -0.3, 0.5;
  Eigen::VectorXd alpha_value2(N);
  alpha_value2 << 0.3, INFINITY, 0.5;

  Eigen::VectorXd beta(N);
  beta << 0.3, 0.8, 1.3;
  Eigen::VectorXd beta_size(N - 1);
  beta_size << 0.3, 0.8;
  Eigen::VectorXd beta_value1(N);
  beta_value1 << 0.3, -0.8, 0.5;
  Eigen::VectorXd beta_value2(N);
  beta_value2 << 0.3, INFINITY, 0.5;

  stan::math::matrix_cl<int> n_cl(n);
  stan::math::matrix_cl<int> n_size_cl(n_size);
  stan::math::matrix_cl<int> n_value_cl(n_value);
  stan::math::matrix_cl<double> alpha_cl(alpha);
  stan::math::matrix_cl<double> alpha_size_cl(alpha_size);
  stan::math::matrix_cl<double> alpha_value1_cl(alpha_value1);
  stan::math::matrix_cl<double> alpha_value2_cl(alpha_value2);
  stan::math::matrix_cl<double> beta_cl(beta);
  stan::math::matrix_cl<double> beta_size_cl(beta_size);
  stan::math::matrix_cl<double> beta_value1_cl(beta_value1);
  stan::math::matrix_cl<double> beta_value2_cl(beta_value2);

  EXPECT_NO_THROW(stan::math::neg_binomial_lpmf(n_cl, alpha_cl, beta_cl));

  EXPECT_THROW(stan::math::neg_binomial_lpmf(n_size_cl, alpha_cl, beta_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::neg_binomial_lpmf(n_cl, alpha_size_cl, beta_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::neg_binomial_lpmf(n_cl, alpha_cl, beta_size_cl),
               std::invalid_argument);

  EXPECT_THROW(stan::math::neg_binomial_lpmf(n_value_cl, alpha_cl, beta_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::neg_binomial_lpmf(n_cl, alpha_value1_cl, beta_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::neg_binomial_lpmf(n_cl, alpha_value2_cl, beta_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::neg_binomial_lpmf(n_cl, alpha_cl, beta_value1_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::neg_binomial_lpmf(n_cl, alpha_cl, beta_value2_cl),
               std::domain_error);
}

auto neg_binomial_lpmf_functor
    = [](const auto& n, const auto& alpha, const auto& beta) {
        return stan::math::neg_binomial_lpmf(n, alpha, beta);
      };
auto neg_binomial_lpmf_functor_propto
    = [](const auto& n, const auto& alpha, const auto& beta) {
        return stan::math::neg_binomial_lpmf<true>(n, alpha, beta);
      };

TEST(muProbDistributionsNegBinomial, opencl_matches_cpu_small) {
  int N = 3;
  int M = 2;

  std::vector<int> n{1, 0, 12};
  Eigen::VectorXd alpha(N);
  alpha << 0.3, 0.5, 1.8;
  Eigen::VectorXd beta(N);
  beta << 0.3, 0.8, 1.3;

  stan::math::test::compare_cpu_opencl_prim_rev(neg_binomial_lpmf_functor, n,
                                                alpha, beta);
  stan::math::test::compare_cpu_opencl_prim_rev(
      neg_binomial_lpmf_functor_propto, n, alpha, beta);
  stan::math::test::compare_cpu_opencl_prim_rev(neg_binomial_lpmf_functor, n,
                                                alpha.transpose().eval(),
                                                beta.transpose().eval());
  stan::math::test::compare_cpu_opencl_prim_rev(
      neg_binomial_lpmf_functor_propto, n, alpha.transpose().eval(),
      beta.transpose().eval());
}

TEST(muProbDistributionsNegBinomial, opencl_broadcast_n) {
  int N = 3;

  int n = 2;
  Eigen::VectorXd alpha(N);
  alpha << 0.3, 0.5, 1.8;
  Eigen::VectorXd beta(N);
  beta << 0.3, 0.8, 1.3;

  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      neg_binomial_lpmf_functor, n, alpha, beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      neg_binomial_lpmf_functor_propto, n, alpha, beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      neg_binomial_lpmf_functor, n, alpha.transpose().eval(), beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      neg_binomial_lpmf_functor_propto, n, alpha, beta.transpose().eval());
}

TEST(muProbDistributionsNegBinomial, opencl_broadcast_alpha) {
  int N = 3;

  std::vector<int> n{1, 0, 12};
  double alpha = 0.4;
  Eigen::VectorXd beta(N);
  beta << 0.3, 0.8, 1.3;

  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      neg_binomial_lpmf_functor, n, alpha, beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      neg_binomial_lpmf_functor_propto, n, alpha, beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      neg_binomial_lpmf_functor, n, alpha, beta.transpose().eval());
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      neg_binomial_lpmf_functor_propto, n, alpha, beta.transpose().eval());
}

TEST(muProbDistributionsNegBinomial, opencl_broadcast_beta) {
  int N = 3;

  std::vector<int> n{1, 0, 12};
  Eigen::VectorXd alpha(N);
  alpha << 0.3, 0.5, 1.8;
  double beta = 0.4;

  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      neg_binomial_lpmf_functor, n, alpha, beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      neg_binomial_lpmf_functor_propto, n, alpha, beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      neg_binomial_lpmf_functor, n, alpha.transpose().eval(), beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      neg_binomial_lpmf_functor_propto, n, alpha.transpose().eval(), beta);
}

TEST(muProbDistributionsNegBinomial, opencl_matches_cpu_big) {
  int N = 153;

  std::vector<int> n(N);
  for (int i = 0; i < N; i++) {
    n[i] = Eigen::Array<int, Eigen::Dynamic, 1>::Random(1, 1).abs()(0) % 1000;
  }
  Eigen::Matrix<double, Eigen::Dynamic, 1> alpha
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();
  Eigen::Matrix<double, Eigen::Dynamic, 1> beta
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();

  stan::math::test::compare_cpu_opencl_prim_rev(neg_binomial_lpmf_functor, n,
                                                alpha, beta);
  stan::math::test::compare_cpu_opencl_prim_rev(
      neg_binomial_lpmf_functor_propto, n, alpha, beta);
  stan::math::test::compare_cpu_opencl_prim_rev(neg_binomial_lpmf_functor, n,
                                                alpha.transpose().eval(),
                                                beta.transpose().eval());
  stan::math::test::compare_cpu_opencl_prim_rev(
      neg_binomial_lpmf_functor_propto, n, alpha.transpose().eval(),
      beta.transpose().eval());
}

TEST(muProbDistributionsNegBinomial, opencl_large_shapes_reference) {
  // The tests above compare OpenCL against the CPU and cannot see an error
  // that both share. These are absolute references from mpmath at 90 digits,
  // checked against 140 digits, as in
  // test/unit/math/rev/prob/large_shapes_test.cpp. The difference
  // digamma(alpha + n) - digamma(alpha) keeps no correct digits for alpha
  // near 1e15. The partials are compared as x * d/dx. The last row is a
  // small shape where the plain difference is correct.
  using stan::math::var;
  struct TestValue {
    int n;
    double alpha;
    double beta;
    double value;
    double d_alpha;
    double d_beta;
  };
  const std::vector<TestValue> test_values = {
      {3, 0x1.c6bf526340000p+49, 0x1.2f2a36ecd5555p+48, -1.4959226032237274,
       1.3124999999999966e-30, 5.6249999999999838e-31},
      {1, 0x1.c6bf526340000p+49, 0x1.2f2a36ecd5555p+48, -1.9013877113318889,
       -1.9999999999999957e-15, 5.9999999999999829e-15},
      {50, 0x1.c6bf526340000p+49, 0x1.2309ce5400000p+44, -2.8766166803657541,
       2.4999999999998758e-29, 0.0},
      {0, 0x1.9000000000000p+6, 0x1.0aaaaaaaaaaabp+5, -2.9558802241544401,
       -2.9558802241544401e-2, 8.7378640776699017e-2},
  };
  auto expect_reference = [](const TestValue& t, const var& lp,
                             double alpha_adj, double beta_adj,
                             const char* signature) {
    EXPECT_NEAR(lp.val(), t.value, 1e-12 * std::max(1.0, std::fabs(t.value)))
        << signature << ": n = " << t.n << ", alpha = " << t.alpha
        << ", beta = " << t.beta;
    const double ga = t.alpha * t.d_alpha;
    EXPECT_NEAR(t.alpha * alpha_adj, ga, 1e-11 * std::max(1.0, std::fabs(ga)))
        << signature << ": n = " << t.n << ", alpha = " << t.alpha
        << ", beta = " << t.beta;
    const double gb = t.beta * t.d_beta;
    EXPECT_NEAR(t.beta * beta_adj, gb, 1e-11 * std::max(1.0, std::fabs(gb)))
        << signature << ": n = " << t.n << ", alpha = " << t.alpha
        << ", beta = " << t.beta;
  };
  for (const auto& t : test_values) {
    const std::vector<int> n{t.n};
    stan::math::matrix_cl<int> n_cl(n);

    Eigen::Matrix<var, Eigen::Dynamic, 1> alpha(1);
    alpha << t.alpha;
    Eigen::Matrix<var, Eigen::Dynamic, 1> beta(1);
    beta << t.beta;
    auto alpha_cl = stan::math::to_matrix_cl(alpha);
    auto beta_cl = stan::math::to_matrix_cl(beta);
    var lp = stan::math::neg_binomial_lpmf(n_cl, alpha_cl, beta_cl);
    lp.grad();
    expect_reference(t, lp, alpha(0).adj(), beta(0).adj(), "vector");
    stan::math::recover_memory();

    var alpha_s = t.alpha;
    var beta_s = t.beta;
    var lp_s = stan::math::neg_binomial_lpmf(n_cl, alpha_s, beta_s);
    lp_s.grad();
    expect_reference(t, lp_s, alpha_s.adj(), beta_s.adj(), "scalar");
    stan::math::recover_memory();
  }
}

#endif
