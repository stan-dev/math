#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>
#include <vector>

TEST(ProbDistributionsBernoulli, error_checking) {
  int N = 3;

  std::vector<int> n{1, 0, 1};
  std::vector<int> n_size{1, 0, 1, 0};
  std::vector<int> n_value{0, 1, 23};
  Eigen::VectorXd theta(N);
  theta << 0.3, 0.8, 1.0;
  Eigen::VectorXd theta_size(N - 1);
  theta_size << 0.3, 0.8;
  Eigen::VectorXd theta_value(N);
  theta_value << 0.3, -0.8, 0.5;

  stan::math::matrix_cl<int> n_cl(n);
  stan::math::matrix_cl<int> n_size_cl(n_size);
  stan::math::matrix_cl<int> n_value_cl(n_value);
  stan::math::matrix_cl<double> theta_cl(theta);
  stan::math::matrix_cl<double> theta_size_cl(theta_size);
  stan::math::matrix_cl<double> theta_value_cl(theta_value);

  EXPECT_NO_THROW(stan::math::bernoulli_lpmf(n_cl, theta_cl));

  EXPECT_THROW(stan::math::bernoulli_lpmf(n_size_cl, theta_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::bernoulli_lpmf(n_cl, theta_size_cl),
               std::invalid_argument);

  EXPECT_THROW(stan::math::bernoulli_lpmf(n_value_cl, theta_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::bernoulli_lpmf(n_cl, theta_value_cl),
               std::domain_error);
}

auto bernoulli_lpmf_functor = [](const auto& n, const auto& theta) {
  return stan::math::bernoulli_lpmf(n, theta);
};
auto bernoulli_lpmf_functor_propto = [](const auto& n, const auto& theta) {
  return stan::math::bernoulli_lpmf<true>(n, theta);
};

TEST(ProbDistributionsBernoulli, opencl_matches_cpu_small) {
  int N = 3;
  int M = 2;

  std::vector<int> n{0, 1, 0};
  Eigen::VectorXd theta(N);
  theta << 0.3, 0.8, 1.0;

  stan::math::test::compare_cpu_opencl_prim_rev(bernoulli_lpmf_functor, n,
                                                theta);
  stan::math::test::compare_cpu_opencl_prim_rev(bernoulli_lpmf_functor_propto,
                                                n, theta);
  stan::math::test::compare_cpu_opencl_prim_rev(bernoulli_lpmf_functor, n,
                                                theta.transpose().eval());
  stan::math::test::compare_cpu_opencl_prim_rev(bernoulli_lpmf_functor_propto,
                                                n, theta.transpose().eval());
}

TEST(ProbDistributionsBernoulli, opencl_broadcast_n) {
  int N = 3;

  int n_scal = 1;
  Eigen::VectorXd theta(N);
  theta << 0.3, 0.8, 1.0;

  stan::math::test::test_opencl_broadcasting_prim_rev<0>(bernoulli_lpmf_functor,
                                                         n_scal, theta);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      bernoulli_lpmf_functor_propto, n_scal, theta);
}

TEST(ProbDistributionsBernoulli, opencl_broadcast_theta) {
  int N = 3;

  std::vector<int> n{0, 1, 0};
  double theta_scal = 0.4;

  stan::math::test::test_opencl_broadcasting_prim_rev<1>(bernoulli_lpmf_functor,
                                                         n, theta_scal);
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      bernoulli_lpmf_functor_propto, n, theta_scal);
}

TEST(ProbDistributionsBernoulli, opencl_matches_cpu_big) {
  int N = 153;

  std::vector<int> n(N);
  for (int i = 0; i < N; i++) {
    n[i] = Eigen::Array<int, Eigen::Dynamic, 1>::Random(1, 1).abs()(0) % 2;
  }
  Eigen::Matrix<double, Eigen::Dynamic, 1> theta
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();

  stan::math::test::compare_cpu_opencl_prim_rev(bernoulli_lpmf_functor, n,
                                                theta);
  stan::math::test::compare_cpu_opencl_prim_rev(bernoulli_lpmf_functor_propto,
                                                n, theta);
  stan::math::test::compare_cpu_opencl_prim_rev(bernoulli_lpmf_functor, n,
                                                theta.transpose().eval());
  stan::math::test::compare_cpu_opencl_prim_rev(bernoulli_lpmf_functor_propto,
                                                n, theta.transpose().eval());
}

#endif
