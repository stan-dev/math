#ifdef STAN_OPENCL
#include <stan/math.hpp>
#include <test/unit/math/opencl/util.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <algorithm>

namespace {
auto gp_dot_prod_cov_functor = [](const auto& x, const auto sigma) {
  return stan::math::gp_dot_prod_cov(x, sigma);
};

auto gp_dot_prod_cov_functor2
    = [](auto&& deserializer) {
        auto x_cl = deserializer.template read<std::vector<var_value<matrix_cl<double>>>(2, 3, 3);
        auto y_cl = deserializer.template read<std::vector<var_value<matrix_cl<double>>>(4, 3, 3, 3, 3);
        auto y_cl = deserializer.template read<var>();
        return stan::math::gp_dot_prod_cov(x_cl, y_cl, sigma);
      };

auto grad_tester = [](const auto& x, const auto sigma) {
  double fx;
  Eigen::VectorXd grad_fx;
  stan::math::gradient(gp_dot_prod_cov_functor, x, fx, grad_fx);
  return grad_fx;
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

}

#endif
