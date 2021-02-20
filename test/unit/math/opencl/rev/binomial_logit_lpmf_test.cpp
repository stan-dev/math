#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>
#include <vector>

TEST(ProbDistributionsBinomialLogit, error_checking) {
  int N = 3;

  std::vector<int> n{1, 0, 4};
  std::vector<int> n_size{1, 0, 1, 0};
  std::vector<int> n_value{0, 1, 23};

  std::vector<int> m{1, 0, 34};
  std::vector<int> m_size{1, 0, 1, 0};
  std::vector<int> m_value{0, 1, 12};

  Eigen::VectorXd alpha(N);
  alpha << 0.0, -0.8, 1.0;
  Eigen::VectorXd alpha_size(N - 1);
  alpha_size << 0.3, 0.8;
  Eigen::VectorXd alpha_value(N);
  alpha_value << 0.3, 1.4, INFINITY;

  stan::math::matrix_cl<int> n_cl(n);
  stan::math::matrix_cl<int> n_size_cl(n_size);
  stan::math::matrix_cl<int> n_value_cl(n_value);
  stan::math::matrix_cl<int> m_cl(n);
  stan::math::matrix_cl<int> m_size_cl(m_size);
  stan::math::matrix_cl<int> m_value_cl(m_value);
  stan::math::matrix_cl<double> alpha_cl(alpha);
  stan::math::matrix_cl<double> alpha_size_cl(alpha_size);
  stan::math::matrix_cl<double> alpha_value_cl(alpha_value);

  EXPECT_NO_THROW(stan::math::binomial_logit_lpmf(n_cl, m_cl, alpha_cl));

  EXPECT_THROW(stan::math::binomial_logit_lpmf(n_size_cl, m_cl, alpha_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::binomial_logit_lpmf(n_cl, m_size_cl, alpha_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::binomial_logit_lpmf(n_cl, m_cl, alpha_size_cl),
               std::invalid_argument);

  EXPECT_THROW(stan::math::binomial_logit_lpmf(n_value_cl, m_cl, alpha_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::binomial_logit_lpmf(n_cl, m_value_cl, alpha_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::binomial_logit_lpmf(n_cl, m_cl, alpha_value_cl),
               std::domain_error);
}

auto binomial_logit_lpmf_functor
    = [](const auto& n, const auto N, const auto& alpha) {
        return stan::math::binomial_logit_lpmf(n, N, alpha);
      };
auto binomial_logit_lpmf_functor_propto
    = [](const auto& n, const auto N, const auto& alpha) {
        return stan::math::binomial_logit_lpmf<true>(n, N, alpha);
      };

TEST(ProbDistributionsBinomialLogit, opencl_matches_cpu_small) {
  int N = 3;

  std::vector<int> n{0, 1, 3};
  std::vector<int> m{0, 1, 5};
  Eigen::VectorXd alpha(N);
  alpha << 0.0, -0.8, 1.0;

  stan::math::test::compare_cpu_opencl_prim_rev(binomial_logit_lpmf_functor, n,
                                                m, alpha);
  stan::math::test::compare_cpu_opencl_prim_rev(
      binomial_logit_lpmf_functor_propto, n, m, alpha);
  stan::math::test::compare_cpu_opencl_prim_rev(binomial_logit_lpmf_functor, n,
                                                m, alpha.transpose().eval());
  stan::math::test::compare_cpu_opencl_prim_rev(
      binomial_logit_lpmf_functor_propto, n, m, alpha.transpose().eval());
}

TEST(ProbDistributionsBinomialLogit, opencl_broadcast_n) {
  int N = 3;

  int n = 1;
  std::vector<int> m{2, 1, 123};
  Eigen::VectorXd alpha(N);
  alpha << 0.0, -0.8, 1.0;

  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      binomial_logit_lpmf_functor, n, m, alpha);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      binomial_logit_lpmf_functor_propto, n, m, alpha);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      binomial_logit_lpmf_functor, n, m, alpha.transpose().eval());
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      binomial_logit_lpmf_functor_propto, n, m, alpha.transpose().eval());
}

TEST(ProbDistributionsBinomialLogit, opencl_broadcast_N) {
  int N = 3;

  std::vector<int> n{0, 1, 12};
  int m = 14;
  Eigen::VectorXd alpha(N);
  alpha << 0.0, -0.8, 1.0;

  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      binomial_logit_lpmf_functor, n, m, alpha);
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      binomial_logit_lpmf_functor_propto, n, m, alpha);
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      binomial_logit_lpmf_functor, n, m, alpha.transpose().eval());
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      binomial_logit_lpmf_functor_propto, n, m, alpha.transpose().eval());
}

TEST(ProbDistributionsBinomialLogit, opencl_broadcast_alpha) {
  int N = 3;

  std::vector<int> n{0, 1, 12};
  std::vector<int> m{0, 1, 123};
  double alpha_scal = 0.4;

  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      binomial_logit_lpmf_functor, n, m, alpha_scal);
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      binomial_logit_lpmf_functor_propto, n, m, alpha_scal);
}

TEST(ProbDistributionsBinomialLogit, opencl_matches_cpu_big) {
  int N = 153;

  std::vector<int> n(N);
  std::vector<int> m(N);
  for (int i = 0; i < N; i++) {
    n[i] = Eigen::Array<int, Eigen::Dynamic, 1>::Random(1, 1).abs()(0) % 1000;
    m[i] = Eigen::Array<int, Eigen::Dynamic, 1>::Random(1, 1).abs()(0) % 1000
           + n[i];
  }
  Eigen::Matrix<double, Eigen::Dynamic, 1> alpha
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1);

  stan::math::test::compare_cpu_opencl_prim_rev(binomial_logit_lpmf_functor, n,
                                                m, alpha);
  stan::math::test::compare_cpu_opencl_prim_rev(
      binomial_logit_lpmf_functor_propto, n, m, alpha);
  stan::math::test::compare_cpu_opencl_prim_rev(binomial_logit_lpmf_functor, n,
                                                m, alpha.transpose().eval());
  stan::math::test::compare_cpu_opencl_prim_rev(
      binomial_logit_lpmf_functor_propto, n, m, alpha.transpose().eval());
}

TEST(ProbDistributionsBinomialLogit, opencl_n_N_scalar) {
  int N = 3;

  int n = 1;
  int m = 5;
  Eigen::VectorXd alpha(N);
  alpha << 0.3, 0.8, 0.9;

  stan::math::test::compare_cpu_opencl_prim_rev(binomial_logit_lpmf_functor, n,
                                                m, alpha);
  stan::math::test::compare_cpu_opencl_prim_rev(
      binomial_logit_lpmf_functor_propto, n, m, alpha);
}

#endif
