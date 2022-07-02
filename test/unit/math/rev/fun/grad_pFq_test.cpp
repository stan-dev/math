#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/fun/util.hpp>

TEST(RevMath, grad_2F2) {
  using stan::math::grad_pFq;
  using stan::math::var;
  using stan::math::vector_d;
  using stan::math::vector_v;

  vector_v a_v(2);
  vector_d a_d(2);
  a_v << 4, 2;
  a_d << 4, 2;

  vector_v b_v(2);
  vector_d b_d(2);
  b_v << 6, 3;
  b_d << 6, 3;

  var z_v = 4;
  double z_d = 4;

  auto grad_tuple = grad_pFq(a_v, b_v, z_v);

  EXPECT_FLOAT_EQ(3.924636646666071, std::get<0>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(6.897245961898751, std::get<0>(grad_tuple)[1]);

  EXPECT_FLOAT_EQ(-2.775051002566842, std::get<1>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(-4.980095849781222, std::get<1>(grad_tuple)[1]);

  EXPECT_FLOAT_EQ(4.916522138006060, std::get<2>(grad_tuple));

  grad_tuple = grad_pFq(a_v, b_d, z_d);

  EXPECT_FLOAT_EQ(3.924636646666071, std::get<0>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(6.897245961898751, std::get<0>(grad_tuple)[1]);

  grad_tuple = grad_pFq(a_d, b_v, z_d);

  EXPECT_FLOAT_EQ(-2.775051002566842, std::get<1>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(-4.980095849781222, std::get<1>(grad_tuple)[1]);

  grad_tuple = grad_pFq(a_d, b_d, z_v);

  EXPECT_FLOAT_EQ(4.916522138006060, std::get<2>(grad_tuple));

  grad_tuple = grad_pFq(a_v, b_v, z_d);

  EXPECT_FLOAT_EQ(3.924636646666071, std::get<0>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(6.897245961898751, std::get<0>(grad_tuple)[1]);
  EXPECT_FLOAT_EQ(-2.775051002566842, std::get<1>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(-4.980095849781222, std::get<1>(grad_tuple)[1]);

  grad_tuple = grad_pFq(a_v, b_d, z_v);

  EXPECT_FLOAT_EQ(3.924636646666071, std::get<0>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(6.897245961898751, std::get<0>(grad_tuple)[1]);
  EXPECT_FLOAT_EQ(4.916522138006060, std::get<2>(grad_tuple));

  grad_tuple = grad_pFq(a_d, b_v, z_v);

  EXPECT_FLOAT_EQ(-2.775051002566842, std::get<1>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(-4.980095849781222, std::get<1>(grad_tuple)[1]);
  EXPECT_FLOAT_EQ(4.916522138006060, std::get<2>(grad_tuple));
}
TEST(RevMath, grad_2F3) {
  using stan::math::grad_pFq;
  using stan::math::var;
  using stan::math::vector_v;

  vector_v a_v(2);
  a_v << 2, 3;
  vector_v b_v(3);
  b_v << 2, 4, 5;
  var z_v = 1;

  auto grad_tuple = grad_pFq(a_v, b_v, z_v);

  EXPECT_FLOAT_EQ(0.08377717301140296, std::get<0>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(0.05615450733193106, std::get<0>(grad_tuple)[1]);

  EXPECT_FLOAT_EQ(-0.08377717301140296, std::get<1>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(-0.04225296806591615, std::get<1>(grad_tuple)[1]);
  EXPECT_FLOAT_EQ(-0.03387575989873739, std::get<1>(grad_tuple)[2]);

  EXPECT_FLOAT_EQ(0.1712340452215524, std::get<2>(grad_tuple));
}

TEST(RevMath, grad_4F3) {
  using stan::math::grad_pFq;
  using stan::math::var;
  using stan::math::vector_v;

  vector_v a_v(4);
  a_v << 1, 2, 3, 4;
  vector_v b_v(3);
  b_v << 5, 6, 7;
  var z_v = 1;

  auto grad_tuple = grad_pFq(a_v, b_v, z_v);

  EXPECT_FLOAT_EQ(0.1587098625610631, std::get<0>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(0.08249338029396848, std::get<0>(grad_tuple)[1]);
  EXPECT_FLOAT_EQ(0.05611368752226367, std::get<0>(grad_tuple)[2]);
  EXPECT_FLOAT_EQ(0.04261209968272329, std::get<0>(grad_tuple)[3]);

  EXPECT_FLOAT_EQ(-0.03438035893346993, std::get<1>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(-0.02882791253333995, std::get<1>(grad_tuple)[1]);
  EXPECT_FLOAT_EQ(-0.02482622713079612, std::get<1>(grad_tuple)[2]);

  EXPECT_FLOAT_EQ(0.1800529055890911, std::get<2>(grad_tuple));
}

TEST(RevMath, grad_2F1_derivs_match) {
  using stan::math::grad_2F1;
  using stan::math::grad_pFq;
  using stan::math::var;
  using stan::math::vector_v;

  vector_v a_v(2);
  a_v << 1, 1;
  vector_v b_v(1);
  b_v << 1;
  var z_v = 0.6;

  double g_a1;
  double g_a2;
  double g_b1;

  grad_2F1(g_a1, g_a2, g_b1, a_v[0].val(), a_v[1].val(), b_v[0].val(),
           z_v.val());
  auto grad_tuple = grad_pFq(a_v, b_v, z_v);

  EXPECT_FLOAT_EQ(g_a1, std::get<0>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(g_a2, std::get<0>(grad_tuple)[1]);
  EXPECT_FLOAT_EQ(g_b1, std::get<1>(grad_tuple)[0]);
}

TEST(RevMath, grad2F1_2) {
  using stan::math::grad_2F1;
  using stan::math::grad_pFq;
  using stan::math::var;
  using stan::math::vector_v;

  vector_v a_v(2);
  a_v << 1, 31;
  vector_v b_v(1);
  b_v << 41;
  var z_v = 1;

  double g_a1;
  double g_a2;
  double g_b1;

  grad_2F1(g_a1, g_a2, g_b1, a_v[0].val(), a_v[1].val(), b_v[0].val(),
           z_v.val());
  auto grad_tuple = grad_pFq(a_v, b_v, z_v);

  EXPECT_FLOAT_EQ(g_a1, std::get<0>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(g_a2, std::get<0>(grad_tuple)[1]);
  EXPECT_FLOAT_EQ(g_b1, std::get<1>(grad_tuple)[0]);
}

TEST(RevMath, grad2F1_3) {
  using stan::math::grad_2F1;
  using stan::math::grad_pFq;
  using stan::math::var;
  using stan::math::vector_v;

  vector_v a_v(2);
  a_v << 1.0, -2.1;
  vector_v b_v(1);
  b_v << 41.0;
  var z_v = 1.0;

  double g_a1;
  double g_a2;
  double g_b1;

  grad_2F1(g_a1, g_a2, g_b1, a_v[0].val(), a_v[1].val(), b_v[0].val(),
           z_v.val());
  auto grad_tuple = grad_pFq(a_v, b_v, z_v);

  EXPECT_FLOAT_EQ(g_a1, std::get<0>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(g_a2, std::get<0>(grad_tuple)[1]);
  EXPECT_FLOAT_EQ(g_b1, std::get<1>(grad_tuple)[0]);
}

TEST(RevMath, grad2F1_6) {
  using stan::math::grad_2F1;
  using stan::math::grad_pFq;
  using stan::math::var;
  using stan::math::vector_v;

  vector_v a_v(2);
  a_v << 1, -0.5;
  vector_v b_v(1);
  b_v << 10.6;
  var z_v = 0.3;

  double g_a1;
  double g_a2;
  double g_b1;

  grad_2F1(g_a1, g_a2, g_b1, a_v[0].val(), a_v[1].val(), b_v[0].val(),
           z_v.val());
  auto grad_tuple = grad_pFq(a_v, b_v, z_v);

  EXPECT_FLOAT_EQ(g_a1, std::get<0>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(g_a2, std::get<0>(grad_tuple)[1]);
  EXPECT_FLOAT_EQ(g_b1, std::get<1>(grad_tuple)[0]);
}

TEST(RevMath, grad2F1_7) {
  using stan::math::grad_2F1;
  using stan::math::grad_pFq;
  using stan::math::var;
  using stan::math::vector_v;

  vector_v a_v(2);
  a_v << 1, -0.5;
  vector_v b_v(1);
  b_v << 10;
  var z_v = 0.3;

  double g_a1;
  double g_a2;
  double g_b1;

  grad_2F1(g_a1, g_a2, g_b1, a_v[0].val(), a_v[1].val(), b_v[0].val(),
           z_v.val());
  auto grad_tuple = grad_pFq(a_v, b_v, z_v);

  EXPECT_FLOAT_EQ(g_a1, std::get<0>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(g_a2, std::get<0>(grad_tuple)[1]);
  EXPECT_FLOAT_EQ(g_b1, std::get<1>(grad_tuple)[0]);
}

TEST(RevMath, grad2F1_8) {
  using stan::math::grad_2F1;
  using stan::math::grad_pFq;
  using stan::math::var;
  using stan::math::vector_v;

  vector_v a_v(2);
  a_v << -0.5, -4.5;
  vector_v b_v(1);
  b_v << 11;
  var z_v = 0.3;

  double g_a1;
  double g_a2;
  double g_b1;

  grad_2F1(g_a1, g_a2, g_b1, a_v[0].val(), a_v[1].val(), b_v[0].val(),
           z_v.val());
  auto grad_tuple = grad_pFq(a_v, b_v, z_v);

  EXPECT_FLOAT_EQ(g_a1, std::get<0>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(g_a2, std::get<0>(grad_tuple)[1]);
  EXPECT_FLOAT_EQ(g_b1, std::get<1>(grad_tuple)[0]);
}

TEST(RevMath, grad2F1_9) {
  using stan::math::grad_2F1;
  using stan::math::grad_pFq;
  using stan::math::var;
  using stan::math::vector_v;

  vector_v a_v(2);
  a_v << -0.5, -4.5;
  vector_v b_v(1);
  b_v << -3.2;
  var z_v = 0.9;

  double g_a1;
  double g_a2;
  double g_b1;

  grad_2F1(g_a1, g_a2, g_b1, a_v[0].val(), a_v[1].val(), b_v[0].val(),
           z_v.val());
  auto grad_tuple = grad_pFq(a_v, b_v, z_v);

  EXPECT_FLOAT_EQ(g_a1, std::get<0>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(g_a2, std::get<0>(grad_tuple)[1]);
  EXPECT_FLOAT_EQ(g_b1, std::get<1>(grad_tuple)[0]);
}

TEST(RevMath, grad2F1_10) {
  using stan::math::grad_2F1;
  using stan::math::grad_pFq;
  using stan::math::var;
  using stan::math::vector_v;

  vector_v a_v(2);
  a_v << 2, 1;
  vector_v b_v(1);
  b_v << 2;
  var z_v = 0.4;

  double g_a1;
  double g_a2;
  double g_b1;

  grad_2F1(g_a1, g_a2, g_b1, a_v[0].val(), a_v[1].val(), b_v[0].val(),
           z_v.val());
  auto grad_tuple = grad_pFq(a_v, b_v, z_v);

  EXPECT_FLOAT_EQ(g_a1, std::get<0>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(g_a2, std::get<0>(grad_tuple)[1]);
  EXPECT_FLOAT_EQ(g_b1, std::get<1>(grad_tuple)[0]);
}
/*
TEST(RevMath, grad2F1_11) {
  using stan::math::grad_2F1;
  using stan::math::grad_pFq;
  using stan::math::var;
  using stan::math::vector_v;

  vector_v a_v(2);
  a_v << 3.70975, 1;
  vector_v b_v(1);
  b_v << 2.70975;
  var z_v = 0.999696;

  double g_a1;
  double g_b1;

  grad_2F1(g_a1, g_b1, a_v[0].val(), a_v[1].val(), b_v[0].val(), z_v.val());
  auto grad_tuple = grad_pFq(a_v, b_v, z_v);

  EXPECT_FLOAT_EQ(g_a1, std::get<0>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(g_b1, std::get<1>(grad_tuple)[0]);
}*/

TEST(RevMath, F32_converges_by_z) {
  using stan::math::grad_F32;
  using stan::math::grad_pFq;
  using stan::math::var;
  using stan::math::vector_v;

  vector_v a_v(3);
  a_v << 1.0, 1.0, 1.0;
  vector_v b_v(2);
  b_v << 1.0, 1.0;
  var z_v = 0.6;

  double g_calc[6];

  grad_F32(g_calc, 1.0, 1.0, 1.0, 1.0, 1.0, 0.6, 1e-10);
  auto grad_tuple = grad_pFq(a_v, b_v, z_v);

  EXPECT_FLOAT_EQ(g_calc[0], std::get<0>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(g_calc[1], std::get<0>(grad_tuple)[1]);
  EXPECT_FLOAT_EQ(g_calc[2], std::get<0>(grad_tuple)[2]);
  EXPECT_FLOAT_EQ(g_calc[3], std::get<1>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(g_calc[4], std::get<1>(grad_tuple)[1]);
  EXPECT_FLOAT_EQ(g_calc[5], std::get<2>(grad_tuple));
}

TEST(RevMath, grad_F32_double_sign_flip_1) {
  using stan::math::grad_F32;
  using stan::math::grad_pFq;
  using stan::math::var;
  using stan::math::vector_v;

  vector_v a_v(3);
  a_v << 1.0, -0.5, -2.5;
  vector_v b_v(2);
  b_v << 10.0, 1.0;
  var z_v = 0.3;

  double g_calc[6];

  grad_F32(g_calc, 1.0, -.5, -2.5, 10.0, 1.0, 0.3, 1e-10);
  auto grad_tuple = grad_pFq(a_v, b_v, z_v);

  EXPECT_FLOAT_EQ(g_calc[0], std::get<0>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(g_calc[1], std::get<0>(grad_tuple)[1]);
  EXPECT_FLOAT_EQ(g_calc[2], std::get<0>(grad_tuple)[2]);
  EXPECT_FLOAT_EQ(g_calc[3], std::get<1>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(g_calc[4], std::get<1>(grad_tuple)[1]);
  EXPECT_FLOAT_EQ(g_calc[5], std::get<2>(grad_tuple));
}

TEST(RevMath, grad_F32_double_sign_flip_2) {
  using stan::math::grad_F32;
  using stan::math::grad_pFq;
  using stan::math::var;
  using stan::math::vector_v;

  vector_v a_v(3);
  a_v << 1.0, -0.5, -4.5;
  vector_v b_v(2);
  b_v << 10.0, 1.0;
  var z_v = 0.3;

  double g_calc[6];

  grad_F32(g_calc, 1.0, -.5, -4.5, 10.0, 1.0, 0.3, 1e-10);
  auto grad_tuple = grad_pFq(a_v, b_v, z_v);

  EXPECT_FLOAT_EQ(g_calc[0], std::get<0>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(g_calc[1], std::get<0>(grad_tuple)[1]);
  EXPECT_FLOAT_EQ(g_calc[2], std::get<0>(grad_tuple)[2]);
  EXPECT_FLOAT_EQ(g_calc[3], std::get<1>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(g_calc[4], std::get<1>(grad_tuple)[1]);
  EXPECT_FLOAT_EQ(g_calc[5], std::get<2>(grad_tuple));
}

TEST(RevMath, grad_2F1_negative_z) {
  using stan::math::grad_pFq;
  using stan::math::var;
  using stan::math::vector_v;

  vector_v a_v(2);
  a_v << 3.70975, 1.0;
  vector_v b_v(1);
  b_v << 2.70975;
  var z_v = -0.2;

  auto grad_tuple = grad_pFq(a_v, b_v, z_v);

  EXPECT_FLOAT_EQ(-0.0488658806159776, std::get<0>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(-0.193844936204681, std::get<0>(grad_tuple)[1]);
  EXPECT_FLOAT_EQ(0.0677809985598383, std::get<1>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(0.865295247272367, std::get<2>(grad_tuple));
}
