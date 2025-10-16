#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/fun/util.hpp>
#include <vector>

TEST(AgradRevMatrix, mv_squaredNorm) {
  using stan::math::matrix_v;

  matrix_v a(2, 2);
  a << -1.0, 2.0, 5.0, 10.0;

  std::vector<stan::math::var> x{a(0, 0), a(0, 1), a(1, 0), a(1, 1)};

  stan::math::var s = a.squaredNorm();
  EXPECT_FLOAT_EQ(130.0, s.val());

  std::vector<double> g;
  s.grad(x, g);
  EXPECT_FLOAT_EQ(-2.0, g[0]);
  EXPECT_FLOAT_EQ(4.0, g[1]);
  EXPECT_FLOAT_EQ(10.0, g[2]);
  EXPECT_FLOAT_EQ(20.0, g[3]);
}
TEST(AgradRevMatrix, mv_norm) {
  using stan::math::matrix_v;

  matrix_v a(2, 1);
  a << -3.0, 4.0;

  std::vector<stan::math::var> x{a(0, 0), a(1, 0)};

  stan::math::var s = a.norm();
  EXPECT_FLOAT_EQ(5.0, s.val());

  // (see hypot in special_functions_test)
  std::vector<double> g;
  s.grad(x, g);
  EXPECT_FLOAT_EQ(-3.0 / 5.0, g[0]);
  EXPECT_FLOAT_EQ(4.0 / 5.0, g[1]);
}
TEST(AgradRevMatrix, mv_lp_norm) {
  using stan::math::matrix_v;

  matrix_v a(2, 2);
  a << -1.0, 2.0, 5.0, 0.0;

  std::vector<stan::math::var> x{a(0, 0), a(0, 1), a(1, 0), a(1, 1)};

  stan::math::var s = a.lpNorm<1>();
  EXPECT_FLOAT_EQ(8.0, s.val());

  std::vector<double> g;
  s.grad(x, g);
  EXPECT_FLOAT_EQ(-1.0, g[0]);
  EXPECT_FLOAT_EQ(1.0, g[1]);
  EXPECT_FLOAT_EQ(1.0, g[2]);
  // ? depends on impl here, could be -1 or 1
  EXPECT_FLOAT_EQ(0.0, g[3]);
}
TEST(AgradRevMatrix, mv_lp_norm_inf) {
  using stan::math::matrix_v;

  matrix_v a(2, 2);
  a << -1.0, 2.0, -5.0, 0.0;

  std::vector<stan::math::var> x = {a(0, 0), a(0, 1), a(1, 0), a(1, 1)};

  stan::math::var s = a.lpNorm<Eigen::Infinity>();
  EXPECT_FLOAT_EQ(5.0, s.val());

  std::vector<double> g;
  s.grad(x, g);
  EXPECT_FLOAT_EQ(0.0, g[0]);
  EXPECT_FLOAT_EQ(0.0, g[1]);
  EXPECT_FLOAT_EQ(-1.0, g[2]);
  EXPECT_FLOAT_EQ(0.0, g[3]);
}

TEST(AgradRevMatrix, UserCase1) {
  using stan::math::assign;
  using stan::math::dot_product;
  using stan::math::get_base1;
  using stan::math::matrix_v;
  using stan::math::multiply;
  using stan::math::subtract;
  using stan::math::transpose;
  using stan::math::vector_d;
  using stan::math::vector_v;
  using std::vector;

  // also tried DpKm1 > H
  size_t H = 3;
  size_t DpKm1 = 3;

  vector_v vk(DpKm1);
  for (size_t k = 0; k < DpKm1; ++k)
    vk[k] = (k + 1) * (k + 2);

  matrix_v L_etaprec(DpKm1, DpKm1);
  for (size_t m = 0; m < DpKm1; ++m)
    for (size_t n = 0; n < DpKm1; ++n)
      L_etaprec(m, n) = (m + 1) * (n + 1);

  vector_d etamu(DpKm1);
  for (size_t k = 0; k < DpKm1; ++k)
    etamu[k] = 10 + (k * k);

  vector<vector_d> eta(H, vector_d(DpKm1));
  for (size_t h = 0; h < H; ++h)
    for (size_t k = 0; k < DpKm1; ++k)
      eta[h][k] = (h + 1) * (k + 10);

  stan::math::var lp__ = 0.0;

  for (size_t h = 1; h <= H; ++h) {
    assign(vk, multiply(transpose(L_etaprec),
                        subtract(get_base1(eta, h, "eta", 1), etamu)));
    assign(lp__, (lp__ - (0.5 * dot_product(vk, vk))));
  }

  EXPECT_NE(lp__.val(), 0);
}
