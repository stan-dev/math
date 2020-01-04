#include <stan/math/fwd.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathMetaFwd, OperandsAndPartialsFvar) {
  using stan::math::fvar;
  using stan::math::operands_and_partials;

  fvar<double> x1 = 2.0;
  fvar<double> x2 = 3.0;
  fvar<double> x3 = 5.0;
  x1.d_ = 2.0;
  x2.d_ = -1.0;
  x3.d_ = 4.0;

  operands_and_partials<fvar<double>, fvar<double>, fvar<double>> o(x1, x2, x3);
  o.edge1_.partials_[0] += 17.0;
  o.edge2_.partials_[0] += 19.0;
  o.edge3_.partials_[0] += 23.0;

  fvar<double> y = o.build(-1.0);
  EXPECT_FLOAT_EQ(107, y.d_);
  EXPECT_FLOAT_EQ(-1, y.val_);
}

TEST(MathMetaFwd, OperandsAndPartialsFvarScal) {
  using stan::math::fvar;
  using stan::math::operands_and_partials;

  fvar<double> x3 = 5.0;
  x3.d_ = 4.0;

  Eigen::VectorXd dx1(2);
  dx1 << 17.0, 13.0;

  operands_and_partials<fvar<double>> o(x3);
  o.edge1_.partials_[0] += 23.0;
  o.edge1_.partials_[0] += 23.0;
  fvar<double> y = o.build(-1.0);

  EXPECT_FLOAT_EQ(2 * 4 * 23, y.d_);
  EXPECT_FLOAT_EQ(-1, y.val_);
}

TEST(MathMetaFwd, OperandsAndPartialsFvarVec) {
  using stan::math::fvar;
  using stan::math::operands_and_partials;

  std::vector<fvar<double>> x1;
  x1.push_back(fvar<double>(2.0, 2.0));
  x1.push_back(fvar<double>(1.0, 3.0));

  fvar<double> x2 = 3.0;
  fvar<double> x3 = 5.0;
  x2.d_ = -1.0;
  x3.d_ = 4.0;

  Eigen::VectorXd dx1(2);
  dx1 << 17.0, 13.0;

  operands_and_partials<std::vector<fvar<double>>, fvar<double>, fvar<double>>
      o(x1, x2, x3);
  o.edge1_.partials_vec_[0] += dx1;
  o.edge2_.partials_[0] += 19.0;
  o.edge2_.partials_[0] += 19.0;
  o.edge3_.partials_[0] += 23.0;
  o.edge3_.partials_[0] += 23.0;
  fvar<double> y = o.build(-1.0);

  EXPECT_FLOAT_EQ(2 * 17 + 3 * 13 - 2 * 19 + 2 * 4 * 23, y.d_);
  EXPECT_FLOAT_EQ(-1, y.val_);
}

TEST(MathMetaFwd, OperandsAndPartialsFvarMat) {
  using stan::math::fvar;
  using stan::math::operands_and_partials;

  Eigen::Matrix<fvar<double>, -1, -1> x1(2, 2);
  x1 << (fvar<double>(2.0, 2.0)), (fvar<double>(1.0, 3.0)),
      (fvar<double>(4.0, 5.0)), (fvar<double>(6.0, 7.0));

  Eigen::MatrixXd dx1(2, 2);
  dx1 << 17.0, 13.0, 23.0, 32.0;

  operands_and_partials<Eigen::Matrix<fvar<double>, -1, -1>> o(x1);
  o.edge1_.partials_vec_[0] += dx1;

  fvar<double> y = o.build(-1.0);

  EXPECT_FLOAT_EQ(2 * 17 + 3 * 13 + 5 * 23 + 7 * 32, y.d_);
  EXPECT_FLOAT_EQ(-1, y.val_);
}
