#include <stan/math/rev/mat.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/rev/mat/util.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(ProbDistributionsOrdLog, vv) {
  using stan::math::ordered_logistic_lpmf;
  using stan::math::var;
  using stan::math::vector_d;
  using stan::math::vector_v;

  int y = 1;

  AVAR lam_v = -1.32;

  vector_v c_v(3);
  c_v << -0.95, -0.10, 0.95;

  AVAR out_v = ordered_logistic_lpmf(y, lam_v, c_v);
  out_v.grad();

  EXPECT_FLOAT_EQ(out_v.val(), -0.52516294973063);
  EXPECT_FLOAT_EQ(lam_v.adj(), -0.40854102156722);
  EXPECT_FLOAT_EQ(c_v[0].adj(), -lam_v.adj());
  EXPECT_FLOAT_EQ(c_v[1].adj(), 0.0);
  EXPECT_FLOAT_EQ(c_v[2].adj(), 0.0);
}

TEST(ProbDistributionsOrdLog, vd) {
  using stan::math::ordered_logistic_lpmf;
  using stan::math::var;
  using stan::math::vector_d;
  using stan::math::vector_v;

  int y = 1;

  AVAR lam_v = -1.32;
  double lam_d = -1.32;

  vector_v c_v(3);
  c_v << -0.95, -0.10, 0.95;

  vector_d c_d(3);
  c_d << -0.95, -0.10, 0.95;

  AVAR out_v = ordered_logistic_lpmf(y, lam_v, c_d);
  out_v.grad();

  EXPECT_FLOAT_EQ(out_v.val(), -0.52516294973063);
  EXPECT_FLOAT_EQ(lam_v.adj(), -0.40854102156722);

  out_v = ordered_logistic_lpmf(y, lam_d, c_v);
  out_v.grad();

  EXPECT_FLOAT_EQ(out_v.val(), -0.52516294973063);
  EXPECT_FLOAT_EQ(c_v[0].adj(), 0.40854102156722);
  EXPECT_FLOAT_EQ(c_v[1].adj(), 0.0);
  EXPECT_FLOAT_EQ(c_v[2].adj(), 0.0);
}

TEST(ProbDistributionsOrdLog, vv_vec) {
  using stan::math::ordered_logistic_lpmf;
  using stan::math::var;
  using stan::math::vector_d;
  using stan::math::vector_v;

  std::vector<int> y{1, 2, 3, 4};

  vector_v lam_v(4);
  lam_v << -1.32, -0.05, 0.56, 1.55;

  vector_v c_v(3);
  c_v << -0.95, -0.10, 0.95;

  AVAR out_v = ordered_logistic_lpmf(y, lam_v, c_v);
  out_v.grad();

  EXPECT_FLOAT_EQ(out_v.val(), -3.9442226297351447171906);

  EXPECT_FLOAT_EQ(lam_v[0].adj(), -0.40854102156722);
  EXPECT_FLOAT_EQ(lam_v[1].adj(), -0.223446899109214);
  EXPECT_FLOAT_EQ(lam_v[2].adj(), -0.062977689154598);
  EXPECT_FLOAT_EQ(lam_v[3].adj(), 0.354343693774205);

  EXPECT_FLOAT_EQ(c_v[0].adj(), -0.62697485850346);
  EXPECT_FLOAT_EQ(c_v[1].adj(), 0.37990895444219);
  EXPECT_FLOAT_EQ(c_v[2].adj(), 0.587687820118095);
}

TEST(ProbDistributionsOrdLog, vd_vec) {
  using stan::math::ordered_logistic_lpmf;
  using stan::math::var;
  using stan::math::vector_d;
  using stan::math::vector_v;

  std::vector<int> y{1, 2, 3, 4};

  vector_v lam_v(4);
  lam_v << -1.32, -0.05, 0.56, 1.55;

  vector_d lam_d(4);
  lam_d << -1.32, -0.05, 0.56, 1.55;

  vector_v c_v(3);
  c_v << -0.95, -0.10, 0.95;

  vector_d c_d(3);
  c_d << -0.95, -0.10, 0.95;

  AVAR out_v = ordered_logistic_lpmf(y, lam_v, c_d);
  out_v.grad();

  EXPECT_FLOAT_EQ(out_v.val(), -3.9442226297351447171906);

  EXPECT_FLOAT_EQ(lam_v[0].adj(), -0.40854102156722);
  EXPECT_FLOAT_EQ(lam_v[1].adj(), -0.223446899109214);
  EXPECT_FLOAT_EQ(lam_v[2].adj(), -0.062977689154598);
  EXPECT_FLOAT_EQ(lam_v[3].adj(), 0.354343693774205);

  out_v = ordered_logistic_lpmf(y, lam_d, c_v);
  out_v.grad();

  EXPECT_FLOAT_EQ(out_v.val(), -3.9442226297351447171906);

  EXPECT_FLOAT_EQ(c_v[0].adj(), -0.62697485850346);
  EXPECT_FLOAT_EQ(c_v[1].adj(), 0.37990895444219);
  EXPECT_FLOAT_EQ(c_v[2].adj(), 0.587687820118095);
}

TEST(ProbDistributionsOrdLog, vv_stvec) {
  using stan::math::ordered_logistic_lpmf;
  using stan::math::var;
  using stan::math::vector_d;
  using stan::math::vector_v;

  std::vector<int> y{1, 2, 3, 4};

  vector_v lam_v(4);
  lam_v << -1.32, -0.05, 0.56, 1.55;

  vector_v c1_v(3);
  c1_v << -0.95, -0.10, 0.95;

  vector_v c2_v(3);
  c2_v << -0.95, -0.10, 0.95;

  vector_v c3_v(3);
  c3_v << -0.95, -0.10, 0.95;

  vector_v c4_v(3);
  c4_v << -0.95, -0.10, 0.95;

  std::vector<vector_v> std_c_v{c1_v, c2_v, c3_v, c4_v};

  AVAR out_v = ordered_logistic_lpmf(y, lam_v, std_c_v);
  out_v.grad();

  EXPECT_FLOAT_EQ(out_v.val(), -3.9442226297351447171906);

  EXPECT_FLOAT_EQ(lam_v[0].adj(), -0.40854102156722);
  EXPECT_FLOAT_EQ(lam_v[1].adj(), -0.223446899109214);
  EXPECT_FLOAT_EQ(lam_v[2].adj(), -0.062977689154598);
  EXPECT_FLOAT_EQ(lam_v[3].adj(), 0.354343693774205);

  EXPECT_FLOAT_EQ(std_c_v[0][0].adj(), 0.40854102156722);
  EXPECT_FLOAT_EQ(std_c_v[0][1].adj(), 0.0);
  EXPECT_FLOAT_EQ(std_c_v[0][2].adj(), 0.0);

  EXPECT_FLOAT_EQ(std_c_v[1][0].adj(), -1.03551588007068);
  EXPECT_FLOAT_EQ(std_c_v[1][1].adj(), 1.25896277917989);
  EXPECT_FLOAT_EQ(std_c_v[1][2].adj(), 0.0);

  EXPECT_FLOAT_EQ(std_c_v[2][0].adj(), 0.0);
  EXPECT_FLOAT_EQ(std_c_v[2][1].adj(), -0.879053824737702);
  EXPECT_FLOAT_EQ(std_c_v[2][2].adj(), 0.942031513892299);

  EXPECT_FLOAT_EQ(std_c_v[3][0].adj(), 0.0);
  EXPECT_FLOAT_EQ(std_c_v[3][1].adj(), 0.0);
  EXPECT_FLOAT_EQ(std_c_v[3][2].adj(), -0.354343693774205);
}

TEST(ProbDistributionsOrdLog, vd_stvec) {
  using stan::math::ordered_logistic_lpmf;
  using stan::math::var;
  using stan::math::vector_d;
  using stan::math::vector_v;

  std::vector<int> y{1, 2, 3, 4};

  vector_v lam_v(4);
  lam_v << -1.32, -0.05, 0.56, 1.55;

  vector_v c1_v(3);
  c1_v << -0.95, -0.10, 0.95;

  vector_v c2_v(3);
  c2_v << -0.95, -0.10, 0.95;

  vector_v c3_v(3);
  c3_v << -0.95, -0.10, 0.95;

  vector_v c4_v(3);
  c4_v << -0.95, -0.10, 0.95;

  vector_d lam_d(4);
  lam_d << -1.32, -0.05, 0.56, 1.55;

  vector_d c1_d(3);
  c1_d << -0.95, -0.10, 0.95;

  vector_d c2_d(3);
  c2_d << -0.95, -0.10, 0.95;

  vector_d c3_d(3);
  c3_d << -0.95, -0.10, 0.95;

  vector_d c4_d(3);
  c4_d << -0.95, -0.10, 0.95;

  std::vector<vector_v> std_c_v{c1_v, c2_v, c3_v, c4_v};
  std::vector<vector_d> std_c_d{c1_d, c2_d, c3_d, c4_d};

  AVAR out_v = ordered_logistic_lpmf(y, lam_v, std_c_d);
  out_v.grad();

  EXPECT_FLOAT_EQ(out_v.val(), -3.9442226297351447171906);

  EXPECT_FLOAT_EQ(lam_v[0].adj(), -0.40854102156722);
  EXPECT_FLOAT_EQ(lam_v[1].adj(), -0.223446899109214);
  EXPECT_FLOAT_EQ(lam_v[2].adj(), -0.062977689154598);
  EXPECT_FLOAT_EQ(lam_v[3].adj(), 0.354343693774205);

  out_v = ordered_logistic_lpmf(y, lam_d, std_c_v);
  out_v.grad();

  EXPECT_FLOAT_EQ(out_v.val(), -3.9442226297351447171906);

  EXPECT_FLOAT_EQ(std_c_v[0][0].adj(), 0.40854102156722);
  EXPECT_FLOAT_EQ(std_c_v[0][1].adj(), 0.0);
  EXPECT_FLOAT_EQ(std_c_v[0][2].adj(), 0.0);

  EXPECT_FLOAT_EQ(std_c_v[1][0].adj(), -1.03551588007068);
  EXPECT_FLOAT_EQ(std_c_v[1][1].adj(), 1.25896277917989);
  EXPECT_FLOAT_EQ(std_c_v[1][2].adj(), 0.0);

  EXPECT_FLOAT_EQ(std_c_v[2][0].adj(), 0.0);
  EXPECT_FLOAT_EQ(std_c_v[2][1].adj(), -0.879053824737702);
  EXPECT_FLOAT_EQ(std_c_v[2][2].adj(), 0.942031513892299);

  EXPECT_FLOAT_EQ(std_c_v[3][0].adj(), 0.0);
  EXPECT_FLOAT_EQ(std_c_v[3][1].adj(), 0.0);
  EXPECT_FLOAT_EQ(std_c_v[3][2].adj(), -0.354343693774205);
}
