#ifdef STAN_OPENCL
#include <stan/math/opencl/prim.hpp>
#include <stan/math.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <cmath>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using stan::math::from_matrix_cl;
using stan::math::matrix_cl;
using stan::math::opencl::ScalarCl;

namespace {
template <typename F>
void expect_gp_scalar_cl(const F& f) {
  MatrixXd x = MatrixXd::Random(3, 20);
  MatrixXd y = MatrixXd::Random(3, 7);
  matrix_cl<double> x_cl(x);
  matrix_cl<double> y_cl(y);
  const double sigma = 1.3;
  const double l = 0.7;
  MatrixXd expected = from_matrix_cl(f(x_cl, sigma, l));
  EXPECT_MATRIX_NEAR(
      from_matrix_cl(f(x_cl, ScalarCl<double>(sigma), ScalarCl<double>(l))),
      expected, 1e-12);
  EXPECT_MATRIX_NEAR(from_matrix_cl(f(x_cl, ScalarCl<double>(sigma), l)),
                     expected, 1e-12);
  EXPECT_MATRIX_NEAR(from_matrix_cl(f(x_cl, sigma, ScalarCl<double>(l))),
                     expected, 1e-12);
}
template <typename F>
void expect_gp_cross_scalar_cl(const F& f) {
  MatrixXd x = MatrixXd::Random(3, 20);
  MatrixXd y = MatrixXd::Random(3, 7);
  matrix_cl<double> x_cl(x);
  matrix_cl<double> y_cl(y);
  const double sigma = 1.3;
  const double l = 0.7;
  MatrixXd expected = from_matrix_cl(f(x_cl, y_cl, sigma, l));
  EXPECT_MATRIX_NEAR(from_matrix_cl(f(x_cl, y_cl, ScalarCl<double>(sigma),
                                      ScalarCl<double>(l))),
                     expected, 1e-12);
}
}  // namespace

TEST(ScalarClGp, scalar_hyperparameters) {
  expect_gp_scalar_cl([](const auto& x, const auto& s, const auto& l) {
    return stan::math::gp_exp_quad_cov(x, s, l);
  });
  expect_gp_scalar_cl([](const auto& x, const auto& s, const auto& l) {
    return stan::math::gp_exponential_cov(x, s, l);
  });
  expect_gp_scalar_cl([](const auto& x, const auto& s, const auto& l) {
    return stan::math::gp_matern32_cov(x, s, l);
  });
  expect_gp_scalar_cl([](const auto& x, const auto& s, const auto& l) {
    return stan::math::gp_matern52_cov(x, s, l);
  });
}

TEST(ScalarClGp, cross_scalar_hyperparameters) {
  expect_gp_cross_scalar_cl(
      [](const auto& x, const auto& y, const auto& s, const auto& l) {
        return stan::math::gp_exp_quad_cov(x, y, s, l);
      });
  expect_gp_cross_scalar_cl(
      [](const auto& x, const auto& y, const auto& s, const auto& l) {
        return stan::math::gp_exponential_cov(x, y, s, l);
      });
  expect_gp_cross_scalar_cl(
      [](const auto& x, const auto& y, const auto& s, const auto& l) {
        return stan::math::gp_matern32_cov(x, y, s, l);
      });
  expect_gp_cross_scalar_cl(
      [](const auto& x, const auto& y, const auto& s, const auto& l) {
        return stan::math::gp_matern52_cov(x, y, s, l);
      });
}

TEST(ScalarClGp, vector_length_scale_scalar_sigma) {
  MatrixXd x = MatrixXd::Random(3, 20);
  VectorXd l(3);
  l << 0.5, 1.0, 2.0;
  matrix_cl<double> x_cl(x);
  matrix_cl<double> l_cl(l);
  EXPECT_MATRIX_NEAR(
      from_matrix_cl(
          stan::math::gp_matern32_cov(x_cl, ScalarCl<double>(1.3), l_cl)),
      from_matrix_cl(stan::math::gp_matern32_cov(x_cl, 1.3, l_cl)), 1e-12);
  EXPECT_MATRIX_NEAR(
      from_matrix_cl(
          stan::math::gp_exponential_cov(x_cl, ScalarCl<double>(1.3), l_cl)),
      from_matrix_cl(stan::math::gp_exponential_cov(x_cl, 1.3, l_cl)), 1e-12);
}
TEST(ScalarClGp, exp_quad_cov_matches_cpu) {
  MatrixXd x = MatrixXd::Random(3, 15);
  const double sigma = 1.3;
  const double l = 0.7;
  MatrixXd expected(15, 15);
  for (int i = 0; i < 15; ++i) {
    for (int j = 0; j < 15; ++j) {
      expected(i, j)
          = sigma * sigma
            * std::exp(-0.5 / (l * l) * (x.col(i) - x.col(j)).squaredNorm());
    }
  }
  matrix_cl<double> x_cl(x);
  EXPECT_MATRIX_NEAR(
      from_matrix_cl(stan::math::gp_exp_quad_cov(x_cl, sigma, l)), expected,
      1e-12);
  EXPECT_MATRIX_NEAR(from_matrix_cl(stan::math::gp_exp_quad_cov(
                         x_cl, ScalarCl<double>(sigma), ScalarCl<double>(l))),
                     expected, 1e-12);
}

#endif
