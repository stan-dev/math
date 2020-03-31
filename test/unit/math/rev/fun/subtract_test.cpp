#include <stan/math/rev.hpp>
#include <gtest/gtest.h>

TEST(MathRevSubtract, subtract_v_d) {
  using stan::math::matrix_d;
  using stan::math::matrix_v;
  using stan::math::subtract;
  using stan::math::var;
  matrix_d Ad = Eigen::MatrixXd::Random(4, 6);
  matrix_d Bd = Eigen::MatrixXd::Random(4, 6);
  var c(3.2);
  var d(-2.73);
  matrix_v Av = c * Ad;
  matrix_v Bv = d * Bd;
  matrix_v Cvd = subtract(Av, Bd);
  matrix_v Cdv = subtract(Ad, Bv);
  matrix_v Cvv = subtract(Av, Bv);
  matrix_v Cdvd = subtract(Ad, subtract(Bv, Eigen::MatrixXd::Random(4, 6)));
  matrix_v Cddv = subtract(Ad, subtract(Eigen::MatrixXd::Random(4, 6), Bv));
  matrix_v Adcv = subtract(Ad, c);
  matrix_v cvAd = subtract(c, Ad);
  matrix_v Bvcv = subtract(Bv, c);
  matrix_v cvBv = subtract(c, Bv);
  matrix_v Bvscal = subtract(Bv, 0.7);
  matrix_v scalBv = subtract(0.7, Bv);

  Cvd(0, 1).grad();
  EXPECT_FLOAT_EQ(c.adj(), Ad(0, 1));

  stan::math::set_zero_all_adjoints();
  Cdv(2, 3).grad();
  EXPECT_FLOAT_EQ(d.adj(), -Bd(2, 3));

  stan::math::set_zero_all_adjoints();
  Cvv(3, 2).grad();
  EXPECT_FLOAT_EQ(d.adj(), -Bd(3, 2));
  EXPECT_FLOAT_EQ(c.adj(), Ad(3, 2));

  stan::math::set_zero_all_adjoints();
  Cdvd(3, 4).grad();
  EXPECT_FLOAT_EQ(d.adj(), -Bd(3, 4));

  stan::math::set_zero_all_adjoints();
  Cddv(2, 5).grad();
  EXPECT_FLOAT_EQ(d.adj(), Bd(2, 5));

  stan::math::set_zero_all_adjoints();
  Adcv(1, 1).grad();
  EXPECT_FLOAT_EQ(c.adj(), -1.0);
  EXPECT_FLOAT_EQ(d.adj(), 0.0);

  stan::math::set_zero_all_adjoints();
  cvAd(1, 1).grad();
  EXPECT_FLOAT_EQ(c.adj(), 1.0);
  EXPECT_FLOAT_EQ(d.adj(), 0.0);

  stan::math::set_zero_all_adjoints();
  Bvcv(2, 1).grad();
  EXPECT_FLOAT_EQ(c.adj(), -1.0);
  EXPECT_FLOAT_EQ(d.adj(), Bd(2, 1));

  stan::math::set_zero_all_adjoints();
  cvBv(2, 1).grad();
  EXPECT_FLOAT_EQ(c.adj(), 1.0);
  EXPECT_FLOAT_EQ(d.adj(), -Bd(2, 1));

  stan::math::set_zero_all_adjoints();
  Bvscal(2, 0).grad();
  EXPECT_FLOAT_EQ(d.adj(), Bd(2, 0));
  EXPECT_FLOAT_EQ(c.adj(), 0.0);

  stan::math::set_zero_all_adjoints();
  scalBv(2, 0).grad();
  EXPECT_FLOAT_EQ(d.adj(), -Bd(2, 0));
  EXPECT_FLOAT_EQ(c.adj(), 0.0);
}

TEST(MathRevSubtract, subtract_check_dim) {
  using stan::math::matrix_d;
  using stan::math::matrix_v;
  using stan::math::subtract;
  matrix_v A, B;

  A.resize(3, 2);
  B.resize(3, 2);
  EXPECT_NO_THROW(subtract(A, B));

  A.resize(1, 0);
  B.resize(1, 0);
  EXPECT_NO_THROW(subtract(A, B));

  A.resize(0, 0);
  B.resize(0, 0);
  EXPECT_NO_THROW(subtract(A, B));

  A.resize(2, 3);
  B.resize(3, 2);
  EXPECT_THROW(subtract(A, B), std::invalid_argument);
};
