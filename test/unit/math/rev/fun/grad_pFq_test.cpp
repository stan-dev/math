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
  double g_b1;

  grad_2F1(g_a1, g_b1, a_v[0].val(), a_v[1].val(), b_v[0].val(), z_v.val());
  auto grad_tuple = grad_pFq(a_v, b_v, z_v);

  EXPECT_FLOAT_EQ(g_a1, std::get<0>(grad_tuple)[0]);
  EXPECT_FLOAT_EQ(g_b1, std::get<1>(grad_tuple)[0]);
}
