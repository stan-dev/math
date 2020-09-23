#include <stan/math/prim.hpp>
#include <gtest/gtest.h>

TEST(MathMatrixPrim, quad_form_sym_mat) {
  using stan::math::matrix_d;
  using stan::math::quad_form_sym;

  matrix_d resd;
  matrix_d m0;

  matrix_d m02(0, 2);
  matrix_d m22(2, 2);
  m22 << 1, 2, 2, 4;
  EXPECT_THROW(quad_form_sym(m02, m22), std::invalid_argument);
  EXPECT_THROW(quad_form_sym(m22, m02), std::invalid_argument);

  resd = quad_form_sym(m0, m0);
  EXPECT_EQ(0, resd.rows());
  EXPECT_EQ(0, resd.cols());

  resd = quad_form_sym(m0, m02);
  EXPECT_EQ(2, resd.rows());
  EXPECT_EQ(2, resd.cols());

  matrix_d m1(1, 1);
  m1 << 2;
  resd = quad_form_sym(m1, m1);
  EXPECT_FLOAT_EQ(8, resd(0));

  matrix_d m2(1, 2);
  m2 << 1, 2;
  resd = quad_form_sym(m1, m2);
  EXPECT_FLOAT_EQ(2, resd(0, 0));
  EXPECT_FLOAT_EQ(4, resd(0, 1));
  EXPECT_FLOAT_EQ(4, resd(1, 0));
  EXPECT_FLOAT_EQ(8, resd(1, 1));

  matrix_d ad(4, 4);
  matrix_d bd(4, 2);

  bd << 100, 10, 0, 1, -3, -3, 5, 2;
  ad << 2.0, 3.0, 4.0, 5.0, 3.0, 10.0, 2.0, 2.0, 4.0, 2.0, 7.0, 1.0, 5.0, 2.0,
      1.0, 112.0;

  resd = quad_form_sym(ad, bd);
  EXPECT_FLOAT_EQ(25433, resd(0, 0));
  EXPECT_FLOAT_EQ(3396, resd(0, 1));
  EXPECT_FLOAT_EQ(3396, resd(1, 0));
  EXPECT_FLOAT_EQ(725, resd(1, 1));

  bd.resize(4, 3);
  bd << 100, 10, 11, 0, 1, 12, -3, -3, 34, 5, 2, 44;
  resd = quad_form_sym(ad, bd);
  EXPECT_EQ(resd(1, 0), resd(0, 1));
  EXPECT_EQ(resd(2, 0), resd(0, 2));
  EXPECT_EQ(resd(2, 1), resd(1, 2));
}

TEST(MathMatrixPrim, quad_form_sym_vec) {
  using stan::math::matrix_d;
  using stan::math::quad_form_sym;
  using stan::math::vector_d;

  matrix_d m0;
  vector_d v0;
  EXPECT_EQ(0, quad_form_sym(m0, v0));

  matrix_d m1(1, 1);
  m1 << 2;
  vector_d v1(1);
  v1 << 2;
  EXPECT_FLOAT_EQ(8, quad_form_sym(m1, v1));

  matrix_d ad(4, 4);
  vector_d bd(4);
  double res;

  bd << 100, 0, -3, 5;
  ad << 2.0, 3.0, 4.0, 5.0, 3.0, 10.0, 2.0, 2.0, 4.0, 2.0, 7.0, 1.0, 5.0, 2.0,
      1.0, 112.0;

  res = quad_form_sym(ad, bd);
  EXPECT_FLOAT_EQ(25433, res);
}

TEST(MathMatrixPrim, quad_form_sym_asymmetric) {
  using stan::math::matrix_d;

  matrix_d ad(4, 4);
  matrix_d bd(4, 2);

  bd << 100, 10, 0, 1, -3, -3, 5, 2;
  ad << 2.0, 3.0, 4.0, 5.0, 6.0, 10.0, 2.0, 2.0, 7.0, 2.0, 7.0, 1.0, 8.0, 2.0,
      1.0, 112.0;

  EXPECT_THROW(stan::math::quad_form_sym(ad, bd), std::domain_error);
}

TEST(MathMatrixPrim, quad_form_sym_2095) {
  using stan::math::matrix_d;

  matrix_d ad(2, 2);
  matrix_d bd(2, 2);

  ad << 1.25882993696386514, -0.03325949401909023, -0.03325949605426004,
      1.85523447220884385;

  bd << 2.55474619740069508, 0.66362927717720988, -1.92917223922387349,
      2.01853255256721731;

  matrix_d cd = stan::math::quad_form_sym(ad, bd);
  EXPECT_FLOAT_EQ(0, cd(1, 0) - cd(0, 1));
}
