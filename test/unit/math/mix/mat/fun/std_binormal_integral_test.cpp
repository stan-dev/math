#include <stan/math/mix/scal.hpp>
#include <stan/math/mix/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/scal/fun/util.hpp>
#include <test/unit/math/mix/scal/fun/nan_util.hpp>
#include <cmath>
#include <typeinfo>

struct binorm_functor {
  template <typename T>
  inline T operator()(
      const Eigen::Matrix<T, Eigen::Dynamic, 1>& inp_vec) const {
    return stan::math::std_binormal_integral(inp_vec(0), inp_vec(1),
                                             inp_vec(2));
  }
};

TEST(MathFunctions, binormal_integral_using) {
  using stan::math::std_binormal_integral;
}

TEST(MathFunctions, binormal_integral_throw_RV_1_nan_fv) {
  using stan::math::fvar;
  using stan::math::var;
  fvar<var> nan(std::numeric_limits<var>::quiet_NaN(), 1.0);
  fvar<var> rho(0.3, 1.0);
  fvar<var> a = nan;
  fvar<var> b(2, 1.0);
  EXPECT_THROW(stan::math::std_binormal_integral(a, b, rho), std::domain_error);
}
TEST(MathFunctions, binormal_integral_throw_RV_2_nan_fv) {
  using stan::math::fvar;
  using stan::math::var;
  fvar<var> nan(std::numeric_limits<var>::quiet_NaN(), 1.0);
  fvar<var> rho(0.3, 1.0);
  fvar<var> a(2.0, 1.0);
  fvar<var> b = nan;
  EXPECT_THROW(stan::math::std_binormal_integral(a, b, rho), std::domain_error);
}
TEST(MathFunctions, binormal_integral_throw_rho_nan_fv) {
  using stan::math::fvar;
  using stan::math::var;
  fvar<var> nan(std::numeric_limits<var>::quiet_NaN(), 1.0);
  fvar<var> rho = nan;
  fvar<var> a(2.0, 1.0);
  fvar<var> b(2.0, 1.0);
  EXPECT_THROW(stan::math::std_binormal_integral(a, b, rho), std::domain_error);
}
TEST(MathFunctions, binormal_integral_throw_corr_coef_neg_fv) {
  using stan::math::fvar;
  using stan::math::var;
  fvar<var> rho(-1.3, 1.0);
  fvar<var> a(2.0, 1.0);
  fvar<var> b(2.0, 1.0);
  EXPECT_THROW(stan::math::std_binormal_integral(a, b, rho), std::domain_error);
}
TEST(MathFunctions, binormal_integral_throw_corr_coef_gt_one_fv) {
  using stan::math::fvar;
  using stan::math::var;
  fvar<var> rho(1.3, 1.0);
  fvar<var> a(2.0, 1.0);
  fvar<var> b(2.0, 1.0);
  EXPECT_THROW(stan::math::std_binormal_integral(a, b, rho), std::domain_error);
}
TEST(MathFunctions, binormal_integral_no_throw_fd) {
  using stan::math::fvar;
  using stan::math::var;
  fvar<var> rho(0.3, 1.0);
  fvar<var> a(2.0, 1.0);
  fvar<var> b(2.0, 1.0);
  EXPECT_NO_THROW(stan::math::std_binormal_integral(a, b, rho));
}
// ffv
TEST(MathFunctions, binormal_integral_throw_RV_1_nan_ffv) {
  using stan::math::fvar;
  using stan::math::var;
  fvar<fvar<var>> nan = std::numeric_limits<double>::quiet_NaN();
  fvar<fvar<var>> rho = 0.3;
  fvar<fvar<var>> a = nan;
  fvar<fvar<var>> b = 2;
  EXPECT_THROW(stan::math::std_binormal_integral(a, b, rho), std::domain_error);
}
TEST(MathFunctions, binormal_integral_throw_RV_2_nan_ffd) {
  using stan::math::fvar;
  fvar<fvar<double>> nan = std::numeric_limits<double>::quiet_NaN();
  fvar<fvar<double>> rho = 0.3;
  fvar<fvar<double>> a = 2;
  fvar<fvar<double>> b = nan;
  EXPECT_THROW(stan::math::std_binormal_integral(a, b, rho), std::domain_error);
}
TEST(MathFunctions, binormal_integral_throw_rho_nan_ffd) {
  using stan::math::fvar;
  fvar<fvar<double>> nan = std::numeric_limits<double>::quiet_NaN();
  fvar<fvar<double>> rho = nan;
  fvar<fvar<double>> a = 2;
  fvar<fvar<double>> b = 2;
  EXPECT_THROW(stan::math::std_binormal_integral(a, b, rho), std::domain_error);
}
TEST(MathFunctions, binormal_integral_throw_corr_coef_neg_ffd) {
  using stan::math::fvar;
  fvar<fvar<double>> rho = -1.3;
  fvar<fvar<double>> a = 2;
  fvar<fvar<double>> b = 1;
  EXPECT_THROW(stan::math::std_binormal_integral(a, b, rho), std::domain_error);
}
TEST(MathFunctions, binormal_integral_throw_corr_coef_gt_one_ffd) {
  using stan::math::fvar;
  fvar<fvar<double>> rho = 1.3;
  fvar<fvar<double>> a = 2;
  fvar<fvar<double>> b = 1;
  EXPECT_THROW(stan::math::std_binormal_integral(a, b, rho), std::domain_error);
}
TEST(MathFunctions, binormal_integral_no_throw_ffd) {
  using stan::math::fvar;
  fvar<fvar<double>> rho = 0.3;
  fvar<fvar<double>> a = 2;
  fvar<fvar<double>> b = 1;
  EXPECT_NO_THROW(stan::math::std_binormal_integral(a, b, rho));
}
TEST(MathFunctions, binormal_hessian_v_finite_diff) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  binorm_functor norm;
  Matrix<double, Dynamic, 1> norm_vec(3);

  norm_vec << 0.5, 0.3, 0.7;

  double norm_eval(0);
  double norm_fin_diff_eval(0);

  Matrix<double, Dynamic, 1> grad_norm;
  Matrix<double, Dynamic, Dynamic> H_norm;
  Matrix<double, Dynamic, Dynamic> fin_diff_H_norm;
  Matrix<double, Dynamic, 1> fin_diff_grad_norm;
  Matrix<double, Dynamic, Dynamic> fin_diff_auto_H_norm;

  stan::math::hessian(norm, norm_vec, norm_eval, grad_norm, H_norm);
  stan::math::finite_diff_hessian(norm, norm_vec, norm_fin_diff_eval,
                                  fin_diff_grad_norm, fin_diff_H_norm);

  EXPECT_FLOAT_EQ(norm_eval, norm_fin_diff_eval);

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      EXPECT_NEAR(H_norm(i, j), fin_diff_H_norm(i, j), 1e-09)
          << "i: " << i << " j: " << j;
    }
    EXPECT_NEAR(grad_norm(i), fin_diff_grad_norm(i), 1e-10);
  }
}
TEST(MathFunctions, binormal_grad_hessian_finite_diff) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  binorm_functor norm;
  Matrix<double, Dynamic, 1> norm_vec(3);
  norm_vec << 0.5, 0.3, 0.7;

  double norm_eval;
  Matrix<double, Dynamic, Dynamic> H_norm;
  std::vector<Matrix<double, Dynamic, Dynamic>> grad_H_norm;
  stan::math::grad_hessian(norm, norm_vec, norm_eval, H_norm, grad_H_norm);

  double norm_fin_diff_eval;
  Matrix<double, Dynamic, Dynamic> fin_diff_H_norm;
  std::vector<Matrix<double, Dynamic, Dynamic>> fin_diff_grad_H_norm;
  stan::math::finite_diff_grad_hessian(norm, norm_vec, norm_fin_diff_eval,
                                       fin_diff_H_norm, fin_diff_grad_H_norm);

  EXPECT_FLOAT_EQ(norm_eval, norm_fin_diff_eval);

  for (size_t i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < 3; ++k) {
        EXPECT_NEAR(grad_H_norm[i](j, k), fin_diff_grad_H_norm[i](j, k), 1e-10)
            << " i: " << i << " j: " << j << " k: " << k;
        EXPECT_FLOAT_EQ(H_norm(j, k), fin_diff_H_norm(j, k));
      }
    }
  }
}
