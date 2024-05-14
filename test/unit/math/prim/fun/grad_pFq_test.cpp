#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/fun/util.hpp>

TEST(PrimMath, grad_2F2) {
  using stan::math::grad_pFq;
  using stan::math::hypergeometric_pFq;
  using stan::math::vector_d;

  vector_d a_d(2);
  a_d << 4, 2;

  vector_d b_d(2);
  b_d << 6, 3;

  double z_d = 4;

  auto pfq_val = hypergeometric_pFq(a_d, b_d, z_d);
  auto grad_tuple = grad_pFq<true, true, true>(pfq_val, a_d, b_d, z_d);

  EXPECT_FLOAT_EQ(3.924636646666071, std::get<0>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(6.897245961898751, std::get<0>(grad_tuple)[1]);

  EXPECT_FLOAT_EQ(-2.775051002566842, std::get<1>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(-4.980095849781222, std::get<1>(grad_tuple)[1]);

  EXPECT_FLOAT_EQ(4.916522138006060, std::get<2>(grad_tuple));

  grad_tuple = grad_pFq<true, false, false>(pfq_val, a_d, b_d, z_d);

  EXPECT_FLOAT_EQ(3.924636646666071, std::get<0>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(6.897245961898751, std::get<0>(grad_tuple)[1]);

  grad_tuple = grad_pFq<false, true, false>(pfq_val, a_d, b_d, z_d);

  EXPECT_FLOAT_EQ(-2.775051002566842, std::get<1>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(-4.980095849781222, std::get<1>(grad_tuple)[1]);

  grad_tuple = grad_pFq<false, false, true>(pfq_val, a_d, b_d, z_d);

  EXPECT_FLOAT_EQ(4.916522138006060, std::get<2>(grad_tuple));

  grad_tuple = grad_pFq<true, true, false>(pfq_val, a_d, b_d, z_d);

  EXPECT_FLOAT_EQ(3.924636646666071, std::get<0>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(6.897245961898751, std::get<0>(grad_tuple)[1]);
  EXPECT_FLOAT_EQ(-2.775051002566842, std::get<1>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(-4.980095849781222, std::get<1>(grad_tuple)[1]);

  grad_tuple = grad_pFq<true, false, true>(pfq_val, a_d, b_d, z_d);

  EXPECT_FLOAT_EQ(3.924636646666071, std::get<0>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(6.897245961898751, std::get<0>(grad_tuple)[1]);
  EXPECT_FLOAT_EQ(4.916522138006060, std::get<2>(grad_tuple));

  grad_tuple = grad_pFq<false, true, true>(pfq_val, a_d, b_d, z_d);

  EXPECT_FLOAT_EQ(-2.775051002566842, std::get<1>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(-4.980095849781222, std::get<1>(grad_tuple)[1]);
  EXPECT_FLOAT_EQ(4.916522138006060, std::get<2>(grad_tuple));
}

TEST(PrimMath, grad_2F3) {
  using stan::math::grad_pFq;
  using stan::math::hypergeometric_pFq;
  using stan::math::vector_d;

  vector_d a_d(2);
  a_d << 2, 3;
  vector_d b_d(3);
  b_d << 2, 4, 5;
  double z_d = 1;

  auto pfq_val = hypergeometric_pFq(a_d, b_d, z_d);
  auto grad_tuple = grad_pFq(pfq_val, a_d, b_d, z_d);

  EXPECT_FLOAT_EQ(0.08377717301140296, std::get<0>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(0.05615450733193106, std::get<0>(grad_tuple)[1]);

  EXPECT_FLOAT_EQ(-0.08377717301140296, std::get<1>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(-0.04225296806591615, std::get<1>(grad_tuple)[1]);
  EXPECT_FLOAT_EQ(-0.03387575989873739, std::get<1>(grad_tuple)[2]);

  EXPECT_FLOAT_EQ(0.1712340452215524, std::get<2>(grad_tuple));
}

TEST(PrimMath, grad_4F3) {
  using stan::math::grad_pFq;
  using stan::math::hypergeometric_pFq;
  using stan::math::vector_d;

  vector_d a_d(4);
  a_d << 1, 2, 3, 4;
  vector_d b_d(3);
  b_d << 5, 6, 7;
  double z_d = 1;

  auto pfq_val = hypergeometric_pFq(a_d, b_d, z_d);
  auto grad_tuple = grad_pFq(pfq_val, a_d, b_d, z_d);

  EXPECT_FLOAT_EQ(0.1587098625610631, std::get<0>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(0.08249338029396848, std::get<0>(grad_tuple)[1]);
  EXPECT_FLOAT_EQ(0.05611368752226367, std::get<0>(grad_tuple)[2]);
  EXPECT_FLOAT_EQ(0.04261209968272329, std::get<0>(grad_tuple)[3]);

  EXPECT_FLOAT_EQ(-0.03438035893346993, std::get<1>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(-0.02882791253333995, std::get<1>(grad_tuple)[1]);
  EXPECT_FLOAT_EQ(-0.02482622713079612, std::get<1>(grad_tuple)[2]);

  EXPECT_FLOAT_EQ(0.1800529055890911, std::get<2>(grad_tuple));
}

TEST(PrimMath, grad_2F1_derivs_match) {
  using stan::math::grad_2F1;
  using stan::math::grad_pFq;
  using stan::math::hypergeometric_pFq;
  using stan::math::vector_d;

  vector_d a_d(2);
  a_d << 1, 1;
  vector_d b_d(1);
  b_d << 1;
  double z_d = 0.6;

  auto grad_2F1_tuple = grad_2F1<true>(a_d[0], a_d[1], b_d[0], z_d);
  auto pfq_val = hypergeometric_pFq(a_d, b_d, z_d);
  auto grad_tuple = grad_pFq(pfq_val, a_d, b_d, z_d);

  EXPECT_FLOAT_EQ(std::get<0>(grad_2F1_tuple), std::get<0>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(std::get<1>(grad_2F1_tuple), std::get<0>(grad_tuple)[1]);
  EXPECT_FLOAT_EQ(std::get<2>(grad_2F1_tuple), std::get<1>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(std::get<3>(grad_2F1_tuple), std::get<2>(grad_tuple));
}

TEST(PrimMath, grad2F1_2) {
  using stan::math::grad_2F1;
  using stan::math::grad_pFq;
  using stan::math::hypergeometric_pFq;
  using stan::math::vector_d;

  vector_d a_d(2);
  a_d << 1, 31;
  vector_d b_d(1);
  b_d << 41;
  double z_d = 1;

  auto grad_2F1_tuple = grad_2F1<true>(a_d[0], a_d[1], b_d[0], z_d);
  auto pfq_val = hypergeometric_pFq(a_d, b_d, z_d);
  auto grad_tuple = grad_pFq(pfq_val, a_d, b_d, z_d);

  EXPECT_FLOAT_EQ(std::get<0>(grad_2F1_tuple), std::get<0>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(std::get<1>(grad_2F1_tuple), std::get<0>(grad_tuple)[1]);
  EXPECT_FLOAT_EQ(std::get<2>(grad_2F1_tuple), std::get<1>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(std::get<3>(grad_2F1_tuple), std::get<2>(grad_tuple));
}

TEST(PrimMath, grad2F1_3) {
  using stan::math::grad_2F1;
  using stan::math::grad_pFq;
  using stan::math::hypergeometric_pFq;
  using stan::math::vector_d;

  vector_d a_d(2);
  a_d << 1.0, -2.1;
  vector_d b_d(1);
  b_d << 41.0;
  double z_d = 1.0;

  auto grad_2F1_tuple = grad_2F1<true>(a_d[0], a_d[1], b_d[0], z_d);
  auto pfq_val = hypergeometric_pFq(a_d, b_d, z_d);
  auto grad_tuple = grad_pFq(pfq_val, a_d, b_d, z_d);

  EXPECT_FLOAT_EQ(std::get<0>(grad_2F1_tuple), std::get<0>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(std::get<1>(grad_2F1_tuple), std::get<0>(grad_tuple)[1]);
  EXPECT_FLOAT_EQ(std::get<2>(grad_2F1_tuple), std::get<1>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(std::get<3>(grad_2F1_tuple), std::get<2>(grad_tuple));
}

TEST(PrimMath, grad2F1_6) {
  using stan::math::grad_2F1;
  using stan::math::grad_pFq;
  using stan::math::hypergeometric_pFq;
  using stan::math::vector_d;

  vector_d a_d(2);
  a_d << 1, -0.5;
  vector_d b_d(1);
  b_d << 10.6;
  double z_d = 0.3;

  auto grad_2F1_tuple = grad_2F1<true>(a_d[0], a_d[1], b_d[0], z_d);
  auto pfq_val = hypergeometric_pFq(a_d, b_d, z_d);
  auto grad_tuple = grad_pFq(pfq_val, a_d, b_d, z_d);

  EXPECT_FLOAT_EQ(std::get<0>(grad_2F1_tuple), std::get<0>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(std::get<1>(grad_2F1_tuple), std::get<0>(grad_tuple)[1]);
  EXPECT_FLOAT_EQ(std::get<2>(grad_2F1_tuple), std::get<1>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(std::get<3>(grad_2F1_tuple), std::get<2>(grad_tuple));
}

TEST(PrimMath, grad2F1_7) {
  using stan::math::grad_2F1;
  using stan::math::grad_pFq;
  using stan::math::hypergeometric_pFq;
  using stan::math::vector_d;

  vector_d a_d(2);
  a_d << 1, -0.5;
  vector_d b_d(1);
  b_d << 10;
  double z_d = 0.3;

  auto grad_2F1_tuple = grad_2F1<true>(a_d[0], a_d[1], b_d[0], z_d);
  auto pfq_val = hypergeometric_pFq(a_d, b_d, z_d);
  auto grad_tuple = grad_pFq(pfq_val, a_d, b_d, z_d);

  EXPECT_FLOAT_EQ(std::get<0>(grad_2F1_tuple), std::get<0>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(std::get<1>(grad_2F1_tuple), std::get<0>(grad_tuple)[1]);
  EXPECT_FLOAT_EQ(std::get<2>(grad_2F1_tuple), std::get<1>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(std::get<3>(grad_2F1_tuple), std::get<2>(grad_tuple));
}

TEST(PrimMath, grad2F1_8) {
  using stan::math::grad_2F1;
  using stan::math::grad_pFq;
  using stan::math::hypergeometric_pFq;
  using stan::math::vector_d;

  vector_d a_d(2);
  a_d << -0.5, -4.5;
  vector_d b_d(1);
  b_d << 11;
  double z_d = 0.3;

  auto grad_2F1_tuple = grad_2F1<true>(a_d[0], a_d[1], b_d[0], z_d);
  auto pfq_val = hypergeometric_pFq(a_d, b_d, z_d);
  auto grad_tuple = grad_pFq(pfq_val, a_d, b_d, z_d);

  EXPECT_FLOAT_EQ(std::get<0>(grad_2F1_tuple), std::get<0>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(std::get<1>(grad_2F1_tuple), std::get<0>(grad_tuple)[1]);
  EXPECT_FLOAT_EQ(std::get<2>(grad_2F1_tuple), std::get<1>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(std::get<3>(grad_2F1_tuple), std::get<2>(grad_tuple));
}

TEST(PrimMath, grad2F1_9) {
  using stan::math::grad_2F1;
  using stan::math::grad_pFq;
  using stan::math::hypergeometric_pFq;
  using stan::math::vector_d;

  vector_d a_d(2);
  a_d << -0.5, -4.5;
  vector_d b_d(1);
  b_d << -3.2;
  double z_d = 0.9;

  auto grad_2F1_tuple = grad_2F1<true>(a_d[0], a_d[1], b_d[0], z_d);
  auto pfq_val = hypergeometric_pFq(a_d, b_d, z_d);
  auto grad_tuple = grad_pFq(pfq_val, a_d, b_d, z_d);

  EXPECT_FLOAT_EQ(std::get<0>(grad_2F1_tuple), std::get<0>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(std::get<1>(grad_2F1_tuple), std::get<0>(grad_tuple)[1]);
  EXPECT_FLOAT_EQ(std::get<2>(grad_2F1_tuple), std::get<1>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(std::get<3>(grad_2F1_tuple), std::get<2>(grad_tuple));
}

TEST(PrimMath, grad2F1_10) {
  using stan::math::grad_2F1;
  using stan::math::grad_pFq;
  using stan::math::hypergeometric_pFq;
  using stan::math::vector_d;

  vector_d a_d(2);
  a_d << 2, 1;
  vector_d b_d(1);
  b_d << 2;
  double z_d = 0.4;

  auto grad_2F1_tuple = grad_2F1<true>(a_d[0], a_d[1], b_d[0], z_d);
  auto pfq_val = hypergeometric_pFq(a_d, b_d, z_d);
  auto grad_tuple = grad_pFq(pfq_val, a_d, b_d, z_d);

  EXPECT_FLOAT_EQ(std::get<0>(grad_2F1_tuple), std::get<0>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(std::get<1>(grad_2F1_tuple), std::get<0>(grad_tuple)[1]);
  EXPECT_FLOAT_EQ(std::get<2>(grad_2F1_tuple), std::get<1>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(std::get<3>(grad_2F1_tuple), std::get<2>(grad_tuple));
}

TEST(PrimMath, grad2F1_11) {
  using stan::math::grad_2F1;
  using stan::math::grad_pFq;
  using stan::math::hypergeometric_pFq;
  using stan::math::vector_d;

  vector_d a_d(2);
  a_d << 3.70975, 1;
  vector_d b_d(1);
  b_d << 2.70975;
  double z_d = 0.999696;

  auto grad_2F1_tuple = grad_2F1<true>(a_d[0], a_d[1], b_d[0], z_d);
  auto pfq_val = hypergeometric_pFq(a_d, b_d, z_d);
  auto grad_tuple = grad_pFq(pfq_val, a_d, b_d, z_d);

  EXPECT_FLOAT_EQ(std::get<0>(grad_2F1_tuple), std::get<0>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(std::get<1>(grad_2F1_tuple), std::get<0>(grad_tuple)[1]);
  EXPECT_FLOAT_EQ(std::get<2>(grad_2F1_tuple), std::get<1>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(std::get<3>(grad_2F1_tuple), std::get<2>(grad_tuple));
}

TEST(PrimMath, F32_converges_by_z) {
  using stan::math::grad_F32;
  using stan::math::grad_pFq;
  using stan::math::hypergeometric_pFq;
  using stan::math::vector_d;

  vector_d a_d(3);
  a_d << 1.0, 1.0, 1.0;
  vector_d b_d(2);
  b_d << 1.0, 1.0;
  double z_d = 0.6;

  double g_calc[6];

  grad_F32(g_calc, 1.0, 1.0, 1.0, 1.0, 1.0, 0.6, 1e-10);
  auto pfq_val = hypergeometric_pFq(a_d, b_d, z_d);
  auto grad_tuple = grad_pFq(pfq_val, a_d, b_d, z_d);

  EXPECT_FLOAT_EQ(g_calc[0], std::get<0>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(g_calc[1], std::get<0>(grad_tuple)[1]);
  EXPECT_FLOAT_EQ(g_calc[2], std::get<0>(grad_tuple)[2]);
  EXPECT_FLOAT_EQ(g_calc[3], std::get<1>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(g_calc[4], std::get<1>(grad_tuple)[1]);
  EXPECT_FLOAT_EQ(g_calc[5], std::get<2>(grad_tuple));
}

TEST(PrimMath, grad_F32_double_sign_flip_1) {
  using stan::math::grad_F32;
  using stan::math::grad_pFq;
  using stan::math::hypergeometric_pFq;
  using stan::math::vector_d;

  vector_d a_d(3);
  a_d << 1.0, -0.5, -2.5;
  vector_d b_d(2);
  b_d << 10.0, 1.0;
  double z_d = 0.3;

  double g_calc[6];

  grad_F32(g_calc, 1.0, -.5, -2.5, 10.0, 1.0, 0.3, 1e-10);
  auto pfq_val = hypergeometric_pFq(a_d, b_d, z_d);
  auto grad_tuple = grad_pFq(pfq_val, a_d, b_d, z_d);

  EXPECT_FLOAT_EQ(g_calc[0], std::get<0>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(g_calc[1], std::get<0>(grad_tuple)[1]);
  EXPECT_FLOAT_EQ(g_calc[2], std::get<0>(grad_tuple)[2]);
  EXPECT_FLOAT_EQ(g_calc[3], std::get<1>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(g_calc[4], std::get<1>(grad_tuple)[1]);
  EXPECT_FLOAT_EQ(g_calc[5], std::get<2>(grad_tuple));
}

TEST(PrimMath, grad_F32_double_sign_flip_2) {
  using stan::math::grad_F32;
  using stan::math::grad_pFq;
  using stan::math::hypergeometric_pFq;
  using stan::math::vector_d;

  vector_d a_d(3);
  a_d << 1.0, -0.5, -4.5;
  vector_d b_d(2);
  b_d << 10.0, 1.0;
  double z_d = 0.3;

  double g_calc[6];

  grad_F32(g_calc, 1.0, -.5, -4.5, 10.0, 1.0, 0.3, 1e-10);
  auto pfq_val = hypergeometric_pFq(a_d, b_d, z_d);
  auto grad_tuple = grad_pFq(pfq_val, a_d, b_d, z_d);

  EXPECT_FLOAT_EQ(g_calc[0], std::get<0>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(g_calc[1], std::get<0>(grad_tuple)[1]);
  EXPECT_FLOAT_EQ(g_calc[2], std::get<0>(grad_tuple)[2]);
  EXPECT_FLOAT_EQ(g_calc[3], std::get<1>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(g_calc[4], std::get<1>(grad_tuple)[1]);
  EXPECT_FLOAT_EQ(g_calc[5], std::get<2>(grad_tuple));
}

TEST(PrimMath, grad_2F1_negative_z) {
  using stan::math::grad_pFq;
  using stan::math::hypergeometric_pFq;
  using stan::math::vector_d;

  vector_d a_d(2);
  a_d << 3.70975, 1.0;
  vector_d b_d(1);
  b_d << 2.70975;
  double z_d = -0.2;

  auto pfq_val = hypergeometric_pFq(a_d, b_d, z_d);
  auto grad_tuple = grad_pFq(pfq_val, a_d, b_d, z_d);

  EXPECT_FLOAT_EQ(-0.0488658806159776, std::get<0>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(-0.193844936204681, std::get<0>(grad_tuple)[1]);
  EXPECT_FLOAT_EQ(0.0677809985598383, std::get<1>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(0.865295247272367, std::get<2>(grad_tuple));
}

TEST(PrimMath, grad_3F2_cross_zero) {
  using stan::math::grad_pFq;
  using stan::math::hypergeometric_pFq;
  using stan::math::vector_d;

  vector_d a_d(3);
  a_d << 1, 1, -1;

  vector_d b_d(2);
  b_d << 2, 2;

  double z_d = 0.292893;

  auto pfq_val = hypergeometric_pFq(a_d, b_d, z_d);
  auto grad_tuple = grad_pFq(pfq_val, a_d, b_d, z_d);

  // Values from Mathematica
  EXPECT_FLOAT_EQ(-0.0732232500000000, std::get<0>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(-0.0732232500000000, std::get<0>(grad_tuple)[1]);
  EXPECT_FLOAT_EQ(0.0681675749282880, std::get<0>(grad_tuple)[2]);
  EXPECT_FLOAT_EQ(0.0366116250000000, std::get<1>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(0.0366116250000000, std::get<1>(grad_tuple)[1]);
  EXPECT_FLOAT_EQ(-0.250000000000000, std::get<2>(grad_tuple));
}
