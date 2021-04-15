#ifdef STAN_OPENCL
#include <stan/math.hpp>
#include <test/unit/math/opencl/util.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <algorithm>

TEST(OpenCLPrimGpExponentialCov, exceptions) {
  Eigen::VectorXd a(3);
  a << 1, 2, 3;
  Eigen::VectorXd b(2);
  b << -3, 4;
  std::vector<Eigen::VectorXd> x1{a, a};
  std::vector<Eigen::VectorXd> x2{b, b};
  stan::math::matrix_cl<double> x1_cl(x1);
  stan::math::matrix_cl<double> x2_cl(x2);

  double sigma = 1.3;
  double l1 = 1.4;
  std::vector<double> l2 = {1.2, 0.7, 2.3};
  stan::math::matrix_cl<double> l2_cl(l2);

  EXPECT_THROW(stan::math::gp_exponential_cov(x1_cl, x2_cl, sigma, l1),
               std::invalid_argument);
  EXPECT_THROW(stan::math::gp_exponential_cov(x1_cl, x2_cl, sigma, l2_cl),
               std::invalid_argument);
}

auto gp_exponential_cov1 = [](const auto x, const auto sigma, const auto l) {
  return stan::math::gp_exponential_cov(x, sigma, l);
};
auto gp_exponential_cov2
    = [](const auto x1, const auto x2, const auto sigma, const auto l) {
        return stan::math::gp_exponential_cov(x1, x2, sigma, l);
      };

TEST(OpenCLPrimGpExponentialCov, small) {
  Eigen::VectorXd a(3);
  a << 1, 2, 3;
  Eigen::VectorXd b(3);
  b << -3, 4, -1;
  Eigen::VectorXd c(3);
  c << 4, -5, 3;
  Eigen::VectorXd d(3);
  d << -4, 5, 5;
  std::vector<Eigen::VectorXd> x1{a, b, c};
  std::vector<Eigen::VectorXd> x2{c, d, d, d};

  double sigma = 1.3;
  double l1 = 1.4;
  std::vector<double> l2 = {1.2, 0.7, 2.3};

  stan::math::test::compare_cpu_opencl_prim(gp_exponential_cov1, x1, sigma, l1);
  stan::math::test::compare_cpu_opencl_prim(gp_exponential_cov2, x1, x2, sigma,
                                            l1);
  stan::math::test::compare_cpu_opencl_prim(gp_exponential_cov1, x1, sigma, l2);
  stan::math::test::compare_cpu_opencl_prim(gp_exponential_cov2, x1, x2, sigma,
                                            l2);
}

TEST(OpenCLPrimGpExponentialCov, large) {
  int N1 = 67;
  int N2 = 73;
  std::vector<Eigen::VectorXd> x1;
  std::vector<double> l2;
  for (int i = 0; i < N1; i++) {
    x1.push_back(Eigen::VectorXd::Random(N1));
    l2.push_back(abs(Eigen::VectorXd::Random(1)[0]));
  }
  std::vector<Eigen::VectorXd> x2;
  for (int i = 0; i < N2; i++) {
    x2.push_back(Eigen::VectorXd::Random(N1));
  }

  double sigma = 1.3;
  double l1 = 1.4;

  stan::math::test::compare_cpu_opencl_prim(gp_exponential_cov1, x1, sigma, l1);
  stan::math::test::compare_cpu_opencl_prim(gp_exponential_cov2, x1, x2, sigma,
                                            l1);
  stan::math::test::compare_cpu_opencl_prim(gp_exponential_cov1, x1, sigma, l2);
  stan::math::test::compare_cpu_opencl_prim(gp_exponential_cov2, x1, x2, sigma,
                                            l2);
}

#endif
