#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/rev/mat/util.hpp>

std::vector<double>
unit_vector_grad(Eigen::Matrix<double,Eigen::Dynamic,1>& y_dbl,
                 int k) {
  using Eigen::Matrix;
  using Eigen::Dynamic;
  using stan::math::var;
  Matrix<var,Dynamic,1> y(y_dbl.size());
  for (int i = 0; i < y.size(); ++i)
    y(i) = y_dbl(i);

  std::vector<var> x(y.size());
  for (size_t i = 0; i < x.size(); ++i)
    x[i] = y(i);

  var fx_k = stan::math::unit_vector_constrain(y)[k];
  std::vector<double> grad(y.size());
  fx_k.grad(x,grad);
  return grad;
}
TEST(AgradRevUnitVectorConstrain, Grad) {
  using stan::math::unit_vector_constrain;
  using stan::math::var;
  using Eigen::Matrix;
  using Eigen::Dynamic;
  for (int k = 0; k < 3; ++k) {
    Matrix<AVAR,Dynamic,1> y(3);
    y << 0.0, 3.0, -1.0;
    Matrix<double,Dynamic,1> y_dbl(3);
    y_dbl << 0.0, 3.0, -1.0;

    AVEC x(3);
    for (int i = 0; i < 3; ++i)
      x[i] = y(i);
    Matrix<AVAR,Dynamic,1> theta = unit_vector_constrain(y);
    AVAR fx_k = theta(k);
    std::vector<double> grad;
    fx_k.grad(x,grad);

    std::vector<double> grad_expected = unit_vector_grad(y_dbl,k);
    EXPECT_EQ(grad_expected.size(), grad.size());
    for (size_t i = 0; i < grad_expected.size(); ++i)
      EXPECT_FLOAT_EQ(grad_expected[i], grad[i]);
  }
}
TEST(AgradRevUnitVectorConstrain, exceptions) {
  using stan::math::unit_vector_constrain;
  using stan::math::var;
  Eigen::Matrix<var, Eigen::Dynamic, 1> x(3);
  x.setZero();
  EXPECT_THROW(unit_vector_constrain(x),std::domain_error);
  x.setOnes();
  x(0) = std::numeric_limits<var>::quiet_NaN();
  EXPECT_THROW(unit_vector_constrain(x),std::domain_error);
  x(0) = std::numeric_limits<var>::infinity();
  EXPECT_THROW(unit_vector_constrain(x),std::domain_error);
}

TEST(AgradRevMatrix, check_varis_on_stack) {
  using stan::math::var;
  using stan::math::to_var;
  Eigen::Matrix<var, Eigen::Dynamic, 1> y(3);
  y << 0.0, 3.0, -1.0;
  var lp(0);
  
  test::check_varis_on_stack(stan::math::unit_vector_constrain(y, lp));
  test::check_varis_on_stack(stan::math::unit_vector_constrain(y));
}
