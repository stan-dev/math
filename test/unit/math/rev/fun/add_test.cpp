#include <stan/math/rev.hpp>
#include <gtest/gtest.h>

TEST(MathRevAdd, add_v_d) {
  using stan::math::add;
  using stan::math::matrix_d;
  using stan::math::matrix_v;
  using stan::math::var;
  matrix_d Ad = Eigen::MatrixXd::Random(4, 6);
  matrix_d Bd = Eigen::MatrixXd::Random(4, 6);
  var c = 3.2;
  var d = -2.73;
  matrix_v Av = c * Ad;
  matrix_v Bv = d * Bd;
  matrix_v Cvd = add(Av, Bd);
  matrix_v Cdv = add(Ad, Bv);
  matrix_v Cvv = add(Av, Bv);
  matrix_v Cdvd = add(Ad, add(Bv, Eigen::MatrixXd::Random(4, 6)));
  matrix_v Adcv = add(c, Ad);
  matrix_v Bvcv = add(c, Bv);
  matrix_v Bvscal = add(0.7, Bv);
  Cvd(0, 1).grad();
  EXPECT_FLOAT_EQ(c.adj(), Ad(0, 1));
  EXPECT_FLOAT_EQ(d.adj(), 0.0);

  stan::math::set_zero_all_adjoints();
  Cdv(2, 3).grad();
  EXPECT_FLOAT_EQ(d.adj(), Bd(2, 3));
  EXPECT_FLOAT_EQ(c.adj(), 0.0);

  stan::math::set_zero_all_adjoints();
  Cvv(3, 2).grad();
  EXPECT_FLOAT_EQ(d.adj(), Bd(3, 2));
  EXPECT_FLOAT_EQ(c.adj(), Ad(3, 2));

  stan::math::set_zero_all_adjoints();
  Cdvd(3, 4).grad();
  EXPECT_FLOAT_EQ(d.adj(), Bd(3, 4));
  EXPECT_FLOAT_EQ(c.adj(), 0.0);

  stan::math::set_zero_all_adjoints();
  Adcv(1, 1).grad();
  EXPECT_FLOAT_EQ(c.adj(), 1.0);
  EXPECT_FLOAT_EQ(d.adj(), 0.0);

  stan::math::set_zero_all_adjoints();
  Bvcv(2, 1).grad();
  EXPECT_FLOAT_EQ(c.adj(), 1.0);
  EXPECT_FLOAT_EQ(d.adj(), Bd(2, 1));

  stan::math::set_zero_all_adjoints();
  Bvscal(2, 0).grad();
  EXPECT_FLOAT_EQ(d.adj(), Bd(2, 0));
  EXPECT_FLOAT_EQ(c.adj(), 0.0);
}

TEST(MathRevAdd, add_check_dim) {
  using stan::math::add;
  using stan::math::matrix_d;
  using stan::math::matrix_v;
  matrix_v A, B;

  A.resize(3, 2);
  B.resize(3, 2);
  EXPECT_NO_THROW(add(A, B));

  A.resize(1, 0);
  B.resize(1, 0);
  EXPECT_NO_THROW(add(A, B));

  A.resize(0, 0);
  B.resize(0, 0);
  EXPECT_NO_THROW(add(A, B));

  A.resize(2, 3);
  B.resize(3, 2);
  EXPECT_THROW(add(A, B), std::invalid_argument);
};
