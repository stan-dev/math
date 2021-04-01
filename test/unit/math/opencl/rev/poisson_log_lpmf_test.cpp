#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>
#include <vector>

TEST(ProbDistributionsPoissonLog, error_checking) {
  int N = 3;

  std::vector<int> n{1, 0, 5};
  std::vector<int> n_size{1, 6, 1, 0};
  std::vector<int> n_value{0, -1, 23};
  Eigen::VectorXd alpha(N);
  alpha << 0.3, 0.8, 1.5;
  Eigen::VectorXd alpha_size(N - 1);
  alpha_size << 0.3, 0.8;
  Eigen::VectorXd alpha_value(N);
  alpha_value << 0.3, NAN, 0.5;

  stan::math::matrix_cl<int> n_cl(n);
  stan::math::matrix_cl<int> n_size_cl(n_size);
  stan::math::matrix_cl<int> n_value_cl(n_value);
  stan::math::matrix_cl<double> alpha_cl(alpha);
  stan::math::matrix_cl<double> alpha_size_cl(alpha_size);
  stan::math::matrix_cl<double> alpha_value_cl(alpha_value);

  EXPECT_NO_THROW(stan::math::poisson_log_lpmf(n_cl, alpha_cl));

  EXPECT_THROW(stan::math::poisson_log_lpmf(n_size_cl, alpha_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::poisson_log_lpmf(n_cl, alpha_size_cl),
               std::invalid_argument);

  EXPECT_THROW(stan::math::poisson_log_lpmf(n_value_cl, alpha_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::poisson_log_lpmf(n_cl, alpha_value_cl),
               std::domain_error);
}

auto poisson_log_lpmf_functor = [](const auto& n, const auto& alpha) {
  return stan::math::poisson_log_lpmf(n, alpha);
};
auto poisson_log_lpmf_functor_propto = [](const auto& n, const auto& alpha) {
  return stan::math::poisson_log_lpmf<true>(n, alpha);
};

TEST(ProbDistributionsPoissonLog, opencl_matches_cpu_small) {
  int N = 3;
  int M = 2;

  std::vector<int> n{0, 1, 6};
  Eigen::VectorXd alpha(N);
  alpha << 0.3, 0.8, 2.0;

  stan::math::test::compare_cpu_opencl_prim_rev(poisson_log_lpmf_functor, n,
                                                alpha);
  stan::math::test::compare_cpu_opencl_prim_rev(poisson_log_lpmf_functor_propto,
                                                n, alpha);
  stan::math::test::compare_cpu_opencl_prim_rev(poisson_log_lpmf_functor, n,
                                                alpha.transpose().eval());
  stan::math::test::compare_cpu_opencl_prim_rev(poisson_log_lpmf_functor_propto,
                                                n, alpha.transpose().eval());
}

TEST(ProbDistributionsPoissonLog, opencl_broadcast_n) {
  int N = 3;

  int n_scal = 1;
  Eigen::VectorXd alpha(N);
  alpha << 0.3, 0.8, 1.0;

  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      poisson_log_lpmf_functor, n_scal, alpha);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      poisson_log_lpmf_functor_propto, n_scal, alpha);
}

TEST(ProbDistributionsPoissonLog, opencl_broadcast_alpha) {
  int N = 3;

  std::vector<int> n{0, 1, 5};
  double alpha_scal = 0.4;

  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      poisson_log_lpmf_functor, n, alpha_scal);
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      poisson_log_lpmf_functor_propto, n, alpha_scal);
}

TEST(ProbDistributionsPoissonLog, opencl_matches_cpu_big) {
  int N = 153;

  std::vector<int> n(N);
  for (int i = 0; i < N; i++) {
    n[i] = Eigen::Array<int, Eigen::Dynamic, 1>::Random(1, 1).abs()(0);
  }
  Eigen::Matrix<double, Eigen::Dynamic, 1> alpha
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1);

  stan::math::test::compare_cpu_opencl_prim_rev(poisson_log_lpmf_functor, n,
                                                alpha);
  stan::math::test::compare_cpu_opencl_prim_rev(poisson_log_lpmf_functor_propto,
                                                n, alpha);
  stan::math::test::compare_cpu_opencl_prim_rev(poisson_log_lpmf_functor, n,
                                                alpha.transpose().eval());
  stan::math::test::compare_cpu_opencl_prim_rev(poisson_log_lpmf_functor_propto,
                                                n, alpha.transpose().eval());
}

#endif
