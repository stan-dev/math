#include <stan/math/fwd/mat.hpp>
#include <gtest/gtest.h>

TEST(AgradFwdMatrixUnitVectorConstrain, fd) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;
  using stan::math::unit_vector_constrain;
  using stan::math::vector_fd;
  using std::sqrt;

  EXPECT_THROW(unit_vector_constrain(vector_fd()), std::invalid_argument);

  Matrix<fvar<double>, Dynamic, 1> x(1);
  x << 0.7;
  x(0).d_ = 1.0;

  Matrix<fvar<double>, Dynamic, 1> theta = unit_vector_constrain(x);
  EXPECT_EQ(1, theta.size());
  EXPECT_FLOAT_EQ(1.0, theta[0].val_);
  EXPECT_NEAR(0.0, theta[0].d_, 3e-16);

  Matrix<fvar<double>, Dynamic, 1> x3(3);
  x3.setRandom();
  for (int i = 0; i < x3.size(); ++i)
    x3(i).d_ = 1.0;

  Matrix<fvar<double>, Dynamic, 1> theta3 = unit_vector_constrain(x3);
  EXPECT_EQ(3, theta3.size());
  const double eps = 10e-20;
  for (int i = 0; i < x3.size(); ++i) {
    Eigen::VectorXcd cx3(x3.size());
    for (int j = 0; j < x3.size(); ++j)
      cx3.real().coeffRef(j) = x3(j).val_;
    cx3.imag().setZero();
    cx3.imag().coeffRef(2) = eps;
    // should be cx3.squaredNorm() but Eigen has a bug?
    std::complex<double> SN(0.0);
    for (int j = 0; j < x3.size(); ++j)
      SN += cx3(j) * cx3(j);
    Matrix<double, Dynamic, 1> d = ((cx3 / sqrt(SN)) / eps).imag();
    EXPECT_FLOAT_EQ(d.coeff(i), theta3[i].d_);
  }
}

TEST(AgradFwdMatrixSoftmax, ffd) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::fvar;
  using stan::math::unit_vector_constrain;
  using stan::math::vector_ffd;

  EXPECT_THROW(unit_vector_constrain(vector_ffd()), std::invalid_argument);

  Matrix<fvar<fvar<double> >, Dynamic, 1> x(1);
  x << 0.7;
  x(0).d_ = 1.0;

  Matrix<fvar<fvar<double> >, Dynamic, 1> theta = unit_vector_constrain(x);
  EXPECT_EQ(1, theta.size());
  EXPECT_FLOAT_EQ(1.0, theta[0].val_.val());
  EXPECT_NEAR(0.0, theta[0].d_.val(), 3e-16);

  Matrix<fvar<fvar<double> >, Dynamic, 1> x3(3);
  x3.setRandom();
  for (int i = 0; i < x3.size(); ++i)
    x3(i).d_ = 1.0;

  Matrix<fvar<fvar<double> >, Dynamic, 1> theta3 = unit_vector_constrain(x3);
  EXPECT_EQ(3, theta3.size());
  const double eps = 10e-20;
  for (int i = 0; i < x3.size(); ++i) {
    Eigen::VectorXcd cx3(x3.size());
    for (int j = 0; j < x3.size(); ++j)
      cx3.real().coeffRef(j) = x3(j).val_.val();
    cx3.imag().setZero();
    cx3.imag().coeffRef(2) = eps;
    // should be cx3.squaredNorm() but Eigen has a bug?
    std::complex<double> SN(0.0);
    for (int j = 0; j < x3.size(); ++j)
      SN += cx3(j) * cx3(j);
    Matrix<double, Dynamic, 1> d = ((cx3 / sqrt(SN)) / eps).imag();
    EXPECT_FLOAT_EQ(d.coeff(i), theta3[i].d_.val());
  }
}
