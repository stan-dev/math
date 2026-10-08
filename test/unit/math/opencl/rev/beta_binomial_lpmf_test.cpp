#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

TEST(ProbDistributionsBetaBinomial, error_checking) {
  int N = 3;

  std::vector<int> n{2, 0, 12};
  std::vector<int> n_size{1, 0, 1, 0};
  std::vector<int> N_{2, 0, 123};
  std::vector<int> N_size{100, 100, 100, 100};
  std::vector<int> N_value{2, -1, 23};
  Eigen::VectorXd alpha(N);
  alpha << 0.3, 1.8, 1.3;
  Eigen::VectorXd alpha_size(N - 1);
  alpha_size << 0.3, 0.8;
  Eigen::VectorXd alpha_value1(N);
  alpha_value1 << 0, 0.4, 0.5;
  Eigen::VectorXd alpha_value2(N);
  alpha_value2 << 0.3, INFINITY, 0.5;
  Eigen::VectorXd beta(N);
  beta << 0.3, 1.8, 1.2;
  Eigen::VectorXd beta_size(N - 1);
  beta_size << 0.3, 0.8;
  Eigen::VectorXd beta_value1(N);
  beta_value1 << 0, 0.4, 0.5;
  Eigen::VectorXd beta_value2(N);
  beta_value2 << 0.3, INFINITY, 0.5;

  stan::math::matrix_cl<int> n_cl(n);
  stan::math::matrix_cl<int> n_size_cl(n_size);
  stan::math::matrix_cl<int> N_cl(N_);
  stan::math::matrix_cl<int> N_size_cl(N_size);
  stan::math::matrix_cl<int> N_value_cl(N_value);
  stan::math::matrix_cl<double> alpha_cl(alpha);
  stan::math::matrix_cl<double> alpha_size_cl(alpha_size);
  stan::math::matrix_cl<double> alpha_value1_cl(alpha_value1);
  stan::math::matrix_cl<double> alpha_value2_cl(alpha_value2);
  stan::math::matrix_cl<double> beta_cl(beta);
  stan::math::matrix_cl<double> beta_size_cl(beta_size);
  stan::math::matrix_cl<double> beta_value1_cl(beta_value1);
  stan::math::matrix_cl<double> beta_value2_cl(beta_value2);

  EXPECT_NO_THROW(
      stan::math::beta_binomial_lpmf(n_cl, N_cl, alpha_cl, beta_cl));

  EXPECT_THROW(
      stan::math::beta_binomial_lpmf(n_size_cl, N_cl, alpha_cl, beta_cl),
      std::invalid_argument);
  EXPECT_THROW(
      stan::math::beta_binomial_lpmf(n_cl, N_size_cl, alpha_cl, beta_cl),
      std::invalid_argument);
  EXPECT_THROW(
      stan::math::beta_binomial_lpmf(n_cl, N_cl, alpha_size_cl, beta_cl),
      std::invalid_argument);
  EXPECT_THROW(
      stan::math::beta_binomial_lpmf(n_cl, N_cl, alpha_cl, beta_size_cl),
      std::invalid_argument);

  EXPECT_THROW(
      stan::math::beta_binomial_lpmf(n_cl, N_value_cl, alpha_cl, beta_cl),
      std::domain_error);
  EXPECT_THROW(
      stan::math::beta_binomial_lpmf(n_cl, N_cl, alpha_value1_cl, beta_cl),
      std::domain_error);
  EXPECT_THROW(
      stan::math::beta_binomial_lpmf(n_cl, N_cl, alpha_value2_cl, beta_cl),
      std::domain_error);
  EXPECT_THROW(
      stan::math::beta_binomial_lpmf(n_cl, N_cl, alpha_cl, beta_value1_cl),
      std::domain_error);
  EXPECT_THROW(
      stan::math::beta_binomial_lpmf(n_cl, N_cl, alpha_cl, beta_value2_cl),
      std::domain_error);
}

auto beta_binomial_lpmf_functor
    = [](const auto& n, const auto& N, const auto& alpha, const auto& beta) {
        return stan::math::beta_binomial_lpmf(n, N, alpha, beta);
      };
auto beta_binomial_lpmf_functor_propto
    = [](const auto& n, const auto& N, const auto& alpha, const auto& beta) {
        return stan::math::beta_binomial_lpmf<true>(n, N, alpha, beta);
      };

TEST(ProbDistributionsBetaBinomial, opencl_matches_cpu_small) {
  int N_ = 3;

  std::vector<int> n{2, 0, 12};
  std::vector<int> N{2, 0, 123};
  Eigen::VectorXd alpha(N_);
  alpha << 0.3, 1.8, 1.3;
  Eigen::VectorXd beta(N_);
  beta << 0.3, 1.8, 1.2;

  stan::math::test::compare_cpu_opencl_prim_rev(beta_binomial_lpmf_functor, n,
                                                N, alpha, beta);
  stan::math::test::compare_cpu_opencl_prim_rev(
      beta_binomial_lpmf_functor_propto, n, N, alpha, beta);
  stan::math::test::compare_cpu_opencl_prim_rev(beta_binomial_lpmf_functor, n,
                                                N, alpha.transpose().eval(),
                                                beta.transpose().eval());
  stan::math::test::compare_cpu_opencl_prim_rev(
      beta_binomial_lpmf_functor_propto, n, N, alpha.transpose().eval(),
      beta.transpose().eval());
}

TEST(ProbDistributionsBetaBinomial, opencl_broadcast_n) {
  int N_ = 3;

  int n = 1;
  std::vector<int> N{2, 0, 123};
  Eigen::VectorXd alpha(N_);
  alpha << 0.3, 1.8, 1.3;
  Eigen::VectorXd beta(N_);
  beta << 0.3, 1.8, 1.2;

  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      beta_binomial_lpmf_functor, n, N, alpha, beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      beta_binomial_lpmf_functor_propto, n, N, alpha, beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      beta_binomial_lpmf_functor, n, N, alpha.transpose().eval(),
      beta.transpose().eval());
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      beta_binomial_lpmf_functor_propto, n, N, alpha.transpose().eval(),
      beta.transpose().eval());
}

TEST(ProbDistributionsBetaBinomial, opencl_broadcast_N) {
  int N_ = 3;

  std::vector<int> n{2, 0, 12};
  int N = 15;
  Eigen::VectorXd alpha(N_);
  alpha << 0.3, 1.8, 1.3;
  Eigen::VectorXd beta(N_);
  beta << 0.3, 1.8, 1.2;

  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      beta_binomial_lpmf_functor, n, N, alpha, beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      beta_binomial_lpmf_functor_propto, n, N, alpha, beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      beta_binomial_lpmf_functor, n, N, alpha.transpose().eval(),
      beta.transpose().eval());
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      beta_binomial_lpmf_functor_propto, n, N, alpha.transpose().eval(),
      beta.transpose().eval());
}

TEST(ProbDistributionsBetaBinomial, opencl_broadcast_alpha) {
  int N_ = 3;

  std::vector<int> n{2, 0, 12};
  std::vector<int> N{2, 0, 123};
  double alpha = 1.1;
  Eigen::VectorXd beta(N_);
  beta << 0.3, 1.8, 1.2;

  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      beta_binomial_lpmf_functor, n, N, alpha, beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      beta_binomial_lpmf_functor_propto, n, N, alpha, beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      beta_binomial_lpmf_functor, n, N, alpha, beta.transpose().eval());
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      beta_binomial_lpmf_functor_propto, n, N, alpha, beta.transpose().eval());
}

TEST(ProbDistributionsBetaBinomial, opencl_broadcast_beta) {
  int N_ = 3;

  std::vector<int> n{2, 0, 12};
  std::vector<int> N{2, 0, 123};
  Eigen::VectorXd alpha(N_);
  alpha << 0.3, 1.8, 1.3;
  double beta = 1.2;
  stan::math::test::test_opencl_broadcasting_prim_rev<3>(
      beta_binomial_lpmf_functor, n, N, alpha, beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<3>(
      beta_binomial_lpmf_functor_propto, n, N, alpha, beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<3>(
      beta_binomial_lpmf_functor, n, N, alpha.transpose().eval(), beta);
  stan::math::test::test_opencl_broadcasting_prim_rev<3>(
      beta_binomial_lpmf_functor_propto, n, N, alpha.transpose().eval(), beta);
}

TEST(ProbDistributionsBetaBinomial, opencl_matches_cpu_big) {
  int N_ = 153;

  std::vector<int> n(N_);
  std::vector<int> N(N_);
  for (int i = 0; i < N_; i++) {
    n[i] = Eigen::Array<int, Eigen::Dynamic, 1>::Random(1, 1).abs()(0) % 123;
    N[i] = Eigen::Array<int, Eigen::Dynamic, 1>::Random(1, 1).abs()(0) % 123
           + 123;
  }
  Eigen::Matrix<double, Eigen::Dynamic, 1> alpha
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N_, 1).abs();
  Eigen::Matrix<double, Eigen::Dynamic, 1> beta
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N_, 1).abs();

  stan::math::test::compare_cpu_opencl_prim_rev(beta_binomial_lpmf_functor, n,
                                                N, alpha, beta);
  stan::math::test::compare_cpu_opencl_prim_rev(
      beta_binomial_lpmf_functor_propto, n, N, alpha, beta);
  stan::math::test::compare_cpu_opencl_prim_rev(beta_binomial_lpmf_functor, n,
                                                N, alpha.transpose().eval(),
                                                beta.transpose().eval());
  stan::math::test::compare_cpu_opencl_prim_rev(
      beta_binomial_lpmf_functor_propto, n, N, alpha.transpose().eval(),
      beta.transpose().eval());
}

namespace beta_binomial_opencl_test_internal {
struct TestValue {
  int n;
  int N;
  double alpha;
  double beta;
  double value;
  double grad_log_alpha;  // alpha * d/dalpha
  double grad_log_beta;   // beta * d/dbeta
};

// The tests above compare OpenCL against the CPU and cannot see an error
// that both share. These are absolute references from mpmath at 80 digits:
// the value from mp.loggamma, the gradients from mp.digamma differences,
// both checked at 130 digits (the gradients with mp.diff). The first four
// rows are as in
// test/unit/math/rev/prob/beta_binomial_lpmf_test.cpp. The shapes are in
// hex so that they are exact. The plain differences of lbeta and digamma
// values keep no correct digits for shapes near 1e15.
// The last row is a small shape where the plain differences are correct.
const std::vector<TestValue> test_values = {
    {57, 117, 0x1.1f43fcc4b662cp+45, 0x1.1f43fcc4b662cp+45, -2.6471538352642870,
     -1.4999999999974545, 1.4999999999981384},
    {400, 1000, 0x1.b48eb57e00000p+44, 0x1.977420dc00000p+42,
     -4.1354072540579752e+2, -4.1081081080252486e+2, 4.1081081078769344e+2},
    {0, 117, 0x1.6345785d8a000p+56, 0x1.0a741a4627800p+58,
     -3.3658802476858363e+1, -2.9249999999999996e+1, 2.9249999999999990e+1},
    {117, 117, 0x1.6345785d8a000p+56, 0x1.0a741a4627800p+58,
     -1.6219644025102715e+2, 8.7749999999999936e+1, -8.7749999999999987e+1},
    // one shape below 10
    {5, 117, 0x1.c6bf526340000p+49, 0x1.35c28f5c28f5cp+2,
     -3.4143943462813216e+3, -1.1199999999999266e+2, 1.5906435196801570e+1},
    {30, 117, 0x1.4000000000000p+1, 0x1.6bcc41e900000p+46,
     -8.2341924150797106e+2, 6.9065498632679532, -2.9999999999966625e+1},
    {5, 20, 0x1.4000000000000p+3, 0x1.9000000000000p+4, -1.8540068216786033,
     -3.4626330704764243e-1, 5.0911144710080701e-1},
};

template <typename T_alpha, typename T_beta>
void expect_reference(const TestValue& t, const stan::math::var& lp,
                      const T_alpha& alpha_adj, const T_beta& beta_adj,
                      const std::string& signature) {
  EXPECT_NEAR(lp.val(), t.value, 1e-12 * std::max(1.0, std::fabs(t.value)))
      << signature << ": n = " << t.n << ", N = " << t.N
      << ", alpha = " << t.alpha << ", beta = " << t.beta;
  EXPECT_NEAR(t.alpha * alpha_adj, t.grad_log_alpha,
              1e-11 * std::max(1.0, std::fabs(t.grad_log_alpha)))
      << signature << ": n = " << t.n << ", N = " << t.N
      << ", alpha = " << t.alpha << ", beta = " << t.beta;
  EXPECT_NEAR(t.beta * beta_adj, t.grad_log_beta,
              1e-11 * std::max(1.0, std::fabs(t.grad_log_beta)))
      << signature << ": n = " << t.n << ", N = " << t.N
      << ", alpha = " << t.alpha << ", beta = " << t.beta;
}
}  // namespace beta_binomial_opencl_test_internal

TEST(ProbDistributionsBetaBinomial, opencl_large_shapes_reference) {
  using beta_binomial_opencl_test_internal::expect_reference;
  using beta_binomial_opencl_test_internal::test_values;
  using stan::math::var;
  for (const auto& t : test_values) {
    const std::vector<int> n{t.n};
    const std::vector<int> N{t.N};
    stan::math::matrix_cl<int> n_cl(n);
    stan::math::matrix_cl<int> N_cl(N);

    // shapes as vectors: the kernel generator functions run on the device
    Eigen::Matrix<var, Eigen::Dynamic, 1> alpha(1);
    alpha << t.alpha;
    Eigen::Matrix<var, Eigen::Dynamic, 1> beta(1);
    beta << t.beta;
    auto alpha_cl = stan::math::to_matrix_cl(alpha);
    auto beta_cl = stan::math::to_matrix_cl(beta);
    var lp = stan::math::beta_binomial_lpmf(n_cl, N_cl, alpha_cl, beta_cl);
    lp.grad();
    expect_reference(t, lp, alpha(0).adj(), beta(0).adj(), "vector shapes");
    stan::math::recover_memory();

    // shapes as scalars: alpha + beta is formed on the host
    var alpha_s = t.alpha;
    var beta_s = t.beta;
    var lp_s = stan::math::beta_binomial_lpmf(n_cl, N_cl, alpha_s, beta_s);
    lp_s.grad();
    expect_reference(t, lp_s, alpha_s.adj(), beta_s.adj(), "scalar shapes");
    stan::math::recover_memory();
  }
}

#endif
