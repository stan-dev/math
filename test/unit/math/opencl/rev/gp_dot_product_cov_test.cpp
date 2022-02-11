#ifdef STAN_OPENCL
#include <stan/math.hpp>
#include <test/unit/math/opencl/util.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <algorithm>

TEST(OpenCLRevGpDotProdCov, exceptions) {
  Eigen::VectorXd a(3);
  a << 1, 2, 3;
  Eigen::VectorXd b(3);
  b << -3, 4, 3;
  Eigen::VectorXd c(3);
  c << -3, 4, NAN;
  std::vector<Eigen::VectorXd> x{a, a, b};
  std::vector<Eigen::VectorXd> x_val{a, c, b};
  stan::math::matrix_cl<double> x_cl(x);
  stan::math::matrix_cl<double> x_val_cl(x_val);

  stan::math::var sigma = 1.3;
  stan::math::var sigma_val = -1.3;

  EXPECT_NO_THROW(stan::math::gp_dot_prod_cov(x_cl, sigma));
  EXPECT_THROW(stan::math::gp_dot_prod_cov(x_val_cl, sigma), std::domain_error);
  EXPECT_THROW(stan::math::gp_dot_prod_cov(x_cl, sigma_val), std::domain_error);
}

auto gp_dot_prod_cov_functor = [](const auto& x, const auto sigma) {
  return stan::math::gp_dot_prod_cov(x, sigma);
};
auto gp_dot_prod_cov_functor2
    = [](const auto& x, const auto& y, const auto sigma) {
        return stan::math::gp_dot_prod_cov(x, y, sigma);
      };

TEST(OpenCLRevGpDotProdCov, small) {
  Eigen::VectorXd a(3);
  a << 1, 2, 3;
  Eigen::VectorXd b(3);
  b << -3, 4, -1;
  Eigen::VectorXd c(3);
  c << 4, -5, 3;
  Eigen::VectorXd d(3);
  d << -4, 5, 5;
  std::vector<Eigen::VectorXd> x{a, b};
  std::vector<Eigen::VectorXd> y{b, c, d, d};

  double sigma = 1.3;

  stan::math::test::compare_cpu_opencl_prim_rev(gp_dot_prod_cov_functor, x,
                                                sigma);
  stan::math::test::compare_cpu_opencl_prim_rev(gp_dot_prod_cov_functor2, x, y,
                                                sigma);
}

TEST(OpenCLRevGpDotProdCov, large) {
  int N1 = 67;
  int N2 = 73;
  std::vector<Eigen::VectorXd> x;
  std::vector<Eigen::VectorXd> y;
  for (int i = 0; i < N1; i++) {
    x.push_back(Eigen::VectorXd::Random(N1));
  }
  for (int i = 0; i < N2; i++) {
    y.push_back(Eigen::VectorXd::Random(N1));
  }

  double sigma = 1.3;

  stan::math::test::compare_cpu_opencl_prim_rev(gp_dot_prod_cov_functor, x,
                                                sigma);
  stan::math::test::compare_cpu_opencl_prim_rev(gp_dot_prod_cov_functor2, x, y,
                                                sigma);
}

#endif
