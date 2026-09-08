#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <tuple>

TEST(MathPrimFun, interp1CubicHermiteSetup) {
  Eigen::VectorXd xk(4);
  Eigen::VectorXd yk(4);
  Eigen::VectorXd dydxk(4);
  xk << 0.0, 1.0, 2.5, 4.0;
  yk << 10.0, 12.0, 9.0, 11.0;
  dydxk << 1.0, 0.5, -0.5, 1.5;

  auto W = stan::math::interp1_cubic_hermite_setup(xk, yk, dydxk);
  const auto& coef = std::get<0>(W);
  const auto& knots = std::get<1>(W);

  EXPECT_EQ(3, coef.rows());
  EXPECT_EQ(4, coef.cols());
  EXPECT_EQ(xk.size(), knots.size());

  for (Eigen::Index i = 0; i < coef.size(); ++i) {
    EXPECT_FLOAT_EQ(0.0, coef(i));
  }
  for (Eigen::Index i = 0; i < xk.size(); ++i) {
    EXPECT_FLOAT_EQ(xk(i), knots(i));
  }
}

TEST(MathPrimFun, interp1CubicHermiteEvalIdentityPlaceholder) {
  Eigen::VectorXd xk(3);
  Eigen::VectorXd yk(3);
  Eigen::VectorXd dydxk(3);
  xk << 0.0, 1.0, 3.0;
  yk << 2.0, 4.0, 8.0;
  dydxk << 1.0, 2.0, 1.5;
  auto W = stan::math::interp1_cubic_hermite_setup(xk, yk, dydxk);

  // TODO: Once cubic Hermite evaluation is implemented, fill in expected values
  // using the cubic Hermite spline formula from boost:
  // y(t) = (1-t)^2 * (y0*(1+2*t) + s0*(x-x0)) + t^2 * (y1*(3-2*t) + dx*s1*(t-1))
  // where t = (x - x0) / (x1 - x0), and s_i = dydx_i
  EXPECT_FLOAT_EQ(-1.5, stan::math::interp1_cubic_hermite_eval(-1.5, W));
  EXPECT_FLOAT_EQ(0.0, stan::math::interp1_cubic_hermite_eval(0.0, W));
  EXPECT_FLOAT_EQ(1.25, stan::math::interp1_cubic_hermite_eval(1.25, W));
  EXPECT_FLOAT_EQ(7.0, stan::math::interp1_cubic_hermite_eval(7.0, W));
}

TEST(MathPrimFun, interp1CubicHermiteSetupThrows) {
  Eigen::VectorXd xk2(2);
  Eigen::VectorXd yk2(2);
  Eigen::VectorXd dydxk2(2);
  xk2 << 0.0, 1.0;
  yk2 << 2.0, 3.0;
  dydxk2 << 1.0, 0.5;
  EXPECT_NO_THROW(stan::math::interp1_cubic_hermite_setup(xk2, yk2, dydxk2));

  Eigen::VectorXd xk1(1);
  Eigen::VectorXd yk1(1);
  Eigen::VectorXd dydxk1(1);
  xk1 << 0.0;
  yk1 << 1.0;
  dydxk1 << 0.5;
  EXPECT_THROW(stan::math::interp1_cubic_hermite_setup(xk1, yk1, dydxk1),
               std::domain_error);

  Eigen::VectorXd yk3(3);
  yk3 << 2.0, 3.0, 4.0;
  EXPECT_THROW(stan::math::interp1_cubic_hermite_setup(xk2, yk3, dydxk2),
               std::invalid_argument);

  Eigen::VectorXd dydxk3(3);
  dydxk3 << 1.0, 2.0, 3.0;
  EXPECT_THROW(stan::math::interp1_cubic_hermite_setup(xk2, yk2, dydxk3),
               std::invalid_argument);

  Eigen::VectorXd xk_nan = xk2;
  xk_nan(1) = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(stan::math::interp1_cubic_hermite_setup(xk_nan, yk2, dydxk2),
               std::domain_error);

  Eigen::VectorXd yk_inf = yk2;
  yk_inf(1) = std::numeric_limits<double>::infinity();
  EXPECT_THROW(stan::math::interp1_cubic_hermite_setup(xk2, yk_inf, dydxk2),
               std::domain_error);

  Eigen::VectorXd dydxk_nan = dydxk2;
  dydxk_nan(0) = std::numeric_limits<double>::quiet_NaN();
  EXPECT_THROW(stan::math::interp1_cubic_hermite_setup(xk2, yk2, dydxk_nan),
               std::domain_error);

  Eigen::VectorXd xk_duplicate(3);
  Eigen::VectorXd yk_duplicate(3);
  Eigen::VectorXd dydxk_duplicate(3);
  xk_duplicate << 0.0, 1.0, 1.0;
  yk_duplicate << 2.0, 3.0, 4.0;
  dydxk_duplicate << 1.0, 2.0, 3.0;
  EXPECT_THROW(stan::math::interp1_cubic_hermite_setup(xk_duplicate, yk_duplicate,
                                                       dydxk_duplicate),
               std::domain_error);

  Eigen::VectorXd xk_decreasing(3);
  xk_decreasing << 0.0, 2.0, 1.0;
  EXPECT_THROW(stan::math::interp1_cubic_hermite_setup(xk_decreasing, yk_duplicate,
                                                       dydxk_duplicate),
               std::domain_error);
}

TEST(MathPrimFun, interp1CubicHermiteEvalThrowsForInvalidSetupShape) {
  stan::math::matrix_d bad_coef_cols(2, 3);
  stan::math::vector_d knots3(3);
  bad_coef_cols.setZero();
  knots3 << 0.0, 1.0, 2.0;
  auto bad_cols = std::make_tuple(bad_coef_cols, knots3);
  EXPECT_THROW(stan::math::interp1_cubic_hermite_eval(1.0, bad_cols),
               std::invalid_argument);

  stan::math::matrix_d coef2(2, 4);
  stan::math::vector_d bad_knots(2);
  coef2.setZero();
  bad_knots << 0.0, 1.0;
  auto bad_rows = std::make_tuple(coef2, bad_knots);
  EXPECT_THROW(stan::math::interp1_cubic_hermite_eval(1.0, bad_rows),
               std::invalid_argument);
}
