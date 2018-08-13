#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/fun/jacobian.hpp>
#include <test/unit/math/rev/mat/util.hpp>
#include <vector>
#include <random>

TEST(prob_transform, ordered_jacobian_ad) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::determinant;
  using stan::math::ordered_constrain;
  using stan::math::var;

  Matrix<double, Dynamic, 1> x(3);
  x << -12.0, 3.0, -1.9;
  double lp = 0.0;
  Matrix<double, Dynamic, 1> y = ordered_constrain(x, lp);

  Matrix<var, Dynamic, 1> xv(3);
  xv << -12.0, 3.0, -1.9;

  std::vector<var> xvec(3);
  for (int i = 0; i < 3; ++i)
    xvec[i] = xv[i];

  Matrix<var, Dynamic, 1> yv = ordered_constrain(xv);

  EXPECT_EQ(y.size(), yv.size());
  for (int i = 0; i < y.size(); ++i)
    EXPECT_FLOAT_EQ(y(i), yv(i).val());

  std::vector<var> yvec(3);
  for (unsigned int i = 0; i < 3; ++i)
    yvec[i] = yv[i];

  std::vector<std::vector<double> > j;
  stan::math::jacobian(yvec, xvec, j);

  Matrix<double, Dynamic, Dynamic> J(3, 3);
  for (int m = 0; m < 3; ++m)
    for (int n = 0; n < 3; ++n)
      J(m, n) = j[m][n];

  double log_abs_jacobian_det = log(fabs(determinant(J)));
  EXPECT_FLOAT_EQ(log_abs_jacobian_det, lp);
}

TEST(prob_transform, ordered_constrain_length_zero_no_segfault) {
  using stan::math::var;
  Eigen::Matrix<var, Eigen::Dynamic, 1> xv(0);
  var out = sum(ordered_constrain(xv));

  out.grad();
}

TEST(prob_transform, ordered_constrain_length_one_no_segfault) {
  using stan::math::var;
  Eigen::Matrix<var, Eigen::Dynamic, 1> xv(1);
  xv << 1.0;
  var out = sum(ordered_constrain(xv));

  out.grad();

  EXPECT_FLOAT_EQ(xv(0).adj(), 1.0);
}

TEST(prob_transform, ordered_constrain_analytical_gradients) {
  using stan::math::var;
  Eigen::Matrix<var, Eigen::Dynamic, 1> xv(4);
  xv << -1.0, 0.0, -1.1, 0.5;
  var out = sum(ordered_constrain(xv));

  out.grad();

  EXPECT_FLOAT_EQ(xv(0).adj(), 4.0);
  EXPECT_FLOAT_EQ(xv(1).adj(), 3.0 * exp(xv(1).val()));
  EXPECT_FLOAT_EQ(xv(2).adj(), 2.0 * exp(xv(2).val()));
  EXPECT_FLOAT_EQ(xv(3).adj(), exp(xv(3).val()));
}

TEST(prob_transform, ordered_constrain_analytical_grads_rng) {
  using stan::math::var;

  std::mt19937 rng;
  std::uniform_real_distribution<> uniform(-5.0, 5.0);

  for (int i = 0; i < 5; ++i) {
    Eigen::Matrix<var, Eigen::Dynamic, 1> xv(4);
    for (int j = 0; j < xv.size(); j++) {
      xv(j) = uniform(rng);
    }
    var out = sum(ordered_constrain(xv));

    out.grad();

    EXPECT_FLOAT_EQ(xv(0).adj(), xv.size());
    for (int j = 1; j < xv.size(); ++j) {
      EXPECT_FLOAT_EQ(xv(j).adj(), (xv.size() - j) * exp(xv(j).val()));
    }
  }
}

TEST(AgradRevMatrix, check_varis_on_stack) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> x(3);

  x << -12.0, 3.0, -1.9;
  stan::math::var lp = 0.0;

  test::check_varis_on_stack(stan::math::ordered_constrain(x, lp));
  test::check_varis_on_stack(stan::math::ordered_constrain(x));
}
