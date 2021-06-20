#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/opencl/util.hpp>
#include <vector>

TEST(ProbDistributionsBinomial, error_checking) {
  int N = 3;

  std::vector<int> n{1, 0, 4};
  std::vector<int> n_size{1, 0, 1, 0};
  std::vector<int> n_value{0, 1, 23};

  std::vector<int> m{1, 0, 34};
  std::vector<int> m_size{1, 0, 1, 0};
  std::vector<int> m_value{0, 1, 12};

  Eigen::VectorXd theta(N);
  theta << 0.0, 0.8, 1.0;
  Eigen::VectorXd theta_size(N - 1);
  theta_size << 0.3, 0.8;
  Eigen::VectorXd theta_value1(N);
  theta_value1 << 0.3, -0.8, 0.5;
  Eigen::VectorXd theta_value2(N);
  theta_value2 << 0.3, 1.4, 0.5;

  stan::math::matrix_cl<int> n_cl(n);
  stan::math::matrix_cl<int> n_size_cl(n_size);
  stan::math::matrix_cl<int> n_value_cl(n_value);
  stan::math::matrix_cl<int> m_cl(n);
  stan::math::matrix_cl<int> m_size_cl(m_size);
  stan::math::matrix_cl<int> m_value_cl(m_value);
  stan::math::matrix_cl<double> theta_cl(theta);
  stan::math::matrix_cl<double> theta_size_cl(theta_size);
  stan::math::matrix_cl<double> theta_value1_cl(theta_value1);
  stan::math::matrix_cl<double> theta_value2_cl(theta_value2);

  EXPECT_NO_THROW(stan::math::binomial_lpmf(n_cl, m_cl, theta_cl));

  EXPECT_THROW(stan::math::binomial_lpmf(n_size_cl, m_cl, theta_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::binomial_lpmf(n_cl, m_size_cl, theta_cl),
               std::invalid_argument);
  EXPECT_THROW(stan::math::binomial_lpmf(n_cl, m_cl, theta_size_cl),
               std::invalid_argument);

  EXPECT_THROW(stan::math::binomial_lpmf(n_value_cl, m_cl, theta_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::binomial_lpmf(n_cl, m_value_cl, theta_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::binomial_lpmf(n_cl, m_cl, theta_value1_cl),
               std::domain_error);
  EXPECT_THROW(stan::math::binomial_lpmf(n_cl, m_cl, theta_value2_cl),
               std::domain_error);
}

auto binomial_lpmf_functor
    = [](const auto& n, const auto N, const auto& theta) {
        return stan::math::binomial_lpmf(n, N, theta);
      };
auto binomial_lpmf_functor_propto
    = [](const auto& n, const auto N, const auto& theta) {
        return stan::math::binomial_lpmf<true>(n, N, theta);
      };

TEST(ProbDistributionsBinomial, opencl_matches_cpu_small) {
  int N = 3;

  std::vector<int> n{0, 1, 3};
  std::vector<int> m{0, 1, 5};
  Eigen::VectorXd theta(N);
  theta << 0.3, 0.8, 0.9;

  stan::math::test::compare_cpu_opencl_prim_rev(binomial_lpmf_functor, n, m,
                                                theta);
  stan::math::test::compare_cpu_opencl_prim_rev(binomial_lpmf_functor_propto, n,
                                                m, theta);
  stan::math::test::compare_cpu_opencl_prim_rev(binomial_lpmf_functor, n, m,
                                                theta.transpose().eval());
  stan::math::test::compare_cpu_opencl_prim_rev(binomial_lpmf_functor_propto, n,
                                                m, theta.transpose().eval());
}

TEST(ProbDistributionsBinomial, opencl_broadcast_n) {
  int N = 3;

  int n = 1;
  std::vector<int> m{2, 1, 123};
  Eigen::VectorXd theta(N);
  theta << 0.3, 0.8, 0.9;

  stan::math::test::test_opencl_broadcasting_prim_rev<0>(binomial_lpmf_functor,
                                                         n, m, theta);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      binomial_lpmf_functor_propto, n, m, theta);
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      binomial_lpmf_functor, n, m, theta.transpose().eval());
  stan::math::test::test_opencl_broadcasting_prim_rev<0>(
      binomial_lpmf_functor_propto, n, m, theta.transpose().eval());
}

TEST(ProbDistributionsBinomial, opencl_broadcast_N) {
  int N = 3;

  std::vector<int> n{0, 1, 12};
  int m = 14;
  Eigen::VectorXd theta(N);
  theta << 0.3, 0.8, 0.9;

  stan::math::test::test_opencl_broadcasting_prim_rev<1>(binomial_lpmf_functor,
                                                         n, m, theta);
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      binomial_lpmf_functor_propto, n, m, theta);
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      binomial_lpmf_functor, n, m, theta.transpose().eval());
  stan::math::test::test_opencl_broadcasting_prim_rev<1>(
      binomial_lpmf_functor_propto, n, m, theta.transpose().eval());
}

TEST(ProbDistributionsBinomial, opencl_broadcast_theta) {
  int N = 3;

  std::vector<int> n{0, 1, 12};
  std::vector<int> m{0, 1, 123};
  double theta_scal = 0.4;

  stan::math::test::test_opencl_broadcasting_prim_rev<2>(binomial_lpmf_functor,
                                                         n, m, theta_scal);
  stan::math::test::test_opencl_broadcasting_prim_rev<2>(
      binomial_lpmf_functor_propto, n, m, theta_scal);
}

TEST(ProbDistributionsBinomial, opencl_matches_cpu_big) {
  int N = 153;

  std::vector<int> n(N);
  std::vector<int> m(N);
  for (int i = 0; i < N; i++) {
    n[i] = Eigen::Array<int, Eigen::Dynamic, 1>::Random(1, 1).abs()(0) % 1000;
    m[i] = Eigen::Array<int, Eigen::Dynamic, 1>::Random(1, 1).abs()(0) % 1000
           + n[i];
  }
  Eigen::Matrix<double, Eigen::Dynamic, 1> theta
      = Eigen::Array<double, Eigen::Dynamic, 1>::Random(N, 1).abs();

  stan::math::test::compare_cpu_opencl_prim_rev(binomial_lpmf_functor, n, m,
                                                theta);
  stan::math::test::compare_cpu_opencl_prim_rev(binomial_lpmf_functor_propto, n,
                                                m, theta);
  stan::math::test::compare_cpu_opencl_prim_rev(binomial_lpmf_functor, n, m,
                                                theta.transpose().eval());
  stan::math::test::compare_cpu_opencl_prim_rev(binomial_lpmf_functor_propto, n,
                                                m, theta.transpose().eval());
}

TEST(ProbDistributionsBinomial, opencl_n_N_scalar) {
  int N = 3;

  int n = 1;
  int m = 5;
  Eigen::VectorXd theta(N);
  theta << 0.3, 0.8, 0.9;

  stan::math::test::compare_cpu_opencl_prim_rev(binomial_lpmf_functor, n, m,
                                                theta);
  stan::math::test::compare_cpu_opencl_prim_rev(binomial_lpmf_functor_propto, n,
                                                m, theta);
}

#endif
