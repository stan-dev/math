#include <stan/math/mix/mat.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(ProbDistributionsOrdLog, fv_fv) {
  using stan::math::ordered_logistic_lpmf;
  using stan::math::fvar;
  using stan::math::var;
  using stan::math::vector_d;
  using stan::math::vector_ffv;
  using stan::math::vector_fv;
  using stan::math::vector_d;

  int y = 1;

  fvar<var> lam_fv = -1.32;
  lam_fv.d_ = 1.0;

  vector_fv c_fv(3);
  c_fv << -0.95, -0.10, 0.95;
  for (int i = 0; i < 3; i++)
    c_fv[i].d_ = 1.0;

  fvar<fvar<var>> lam_ffv;
  lam_ffv.val_ = -1.32;
  lam_ffv.d_ = 1.0;
  lam_ffv.val_.d_ = 1.0;

  vector_ffv c_ffv(3);
  c_ffv << -0.95, -0.10, 0.95;
  for (int i = 0; i < 3; i++) {
    c_ffv[i].d_ = 1.0;
    c_ffv[i].val_.d_ = 1.0;
  }

  fvar<var> out_fv = ordered_logistic_lpmf(y, lam_fv, c_fv);
  out_fv.d_.grad();

  EXPECT_FLOAT_EQ(out_fv.val_.val(), -0.52516294973063);
  EXPECT_FLOAT_EQ(out_fv.d_.val() + 1, 0.0 + 1);
  EXPECT_FLOAT_EQ(lam_fv.d_.adj(), -0.40854102156722);
  EXPECT_FLOAT_EQ(c_fv[0].d_.adj(), -lam_fv.d_.adj());
  EXPECT_FLOAT_EQ(c_fv[1].d_.adj(), 0.0);
  EXPECT_FLOAT_EQ(c_fv[2].d_.adj(), 0.0);

  fvar<fvar<var>> out_ffv = ordered_logistic_lpmf(y, lam_ffv, c_ffv);
  out_ffv.d_.val_.grad();

  EXPECT_FLOAT_EQ(out_ffv.val_.val_.val(), -0.52516294973063);
  EXPECT_FLOAT_EQ(out_ffv.d_.val_.val() + 1, 0.0 + 1);
  EXPECT_FLOAT_EQ(lam_ffv.d_.val_.adj(), -0.40854102156722);
  EXPECT_FLOAT_EQ(c_ffv[0].d_.val_.adj(), -lam_ffv.d_.val_.adj());
  EXPECT_FLOAT_EQ(c_ffv[1].d_.val_.adj(), 0.0);
  EXPECT_FLOAT_EQ(c_ffv[2].d_.val_.adj(), 0.0);
}

TEST(ProbDistributionsOrdLog, fv_d) {
  using stan::math::ordered_logistic_lpmf;
  using stan::math::fvar;
  using stan::math::var;
  using stan::math::vector_d;
  using stan::math::vector_ffv;
  using stan::math::vector_fv;
  using stan::math::vector_d;

  int y = 1;

  fvar<var> lam_fv = -1.32;
  lam_fv.d_ = 1.0;

  double lam_d = -1.32;

  vector_fv c_fv(3);
  c_fv << -0.95, -0.10, 0.95;
  for (int i = 0; i < 3; i++)
    c_fv[i].d_ = 1.0;

  vector_d c_d(3);
  c_d << -0.95, -0.10, 0.95;

  fvar<fvar<var>> lam_ffv;
  lam_ffv.val_ = -1.32;
  lam_ffv.d_ = 1.0;
  lam_ffv.val_.d_ = 1.0;

  vector_ffv c_ffv(3);
  c_ffv << -0.95, -0.10, 0.95;
  for (int i = 0; i < 3; i++) {
    c_ffv[i].d_ = 1.0;
    c_ffv[i].val_.d_ = 1.0;
  }

  fvar<var> out = ordered_logistic_lpmf(y, lam_fv, c_d);
  out.d_.grad();

  EXPECT_FLOAT_EQ(out.val_.val(), -0.52516294973063);
  EXPECT_FLOAT_EQ(out.d_.val(), -0.40854102156722);
  EXPECT_FLOAT_EQ(lam_fv.d_.adj(), -0.40854102156722);

  out = ordered_logistic_lpmf(y, lam_d, c_fv);
  out.d_.grad();

  EXPECT_FLOAT_EQ(out.val_.val(), -0.52516294973063);
  EXPECT_FLOAT_EQ(out.d_.val(), 0.40854102156722);
  EXPECT_FLOAT_EQ(c_fv[0].d_.adj(), 0.40854102156722);
  EXPECT_FLOAT_EQ(c_fv[1].d_.adj(), 0.0);
  EXPECT_FLOAT_EQ(c_fv[2].d_.adj(), 0.0);

  fvar<fvar<var>> out_ffv = ordered_logistic_lpmf(y, lam_ffv, c_d);
  out_ffv.d_.val_.grad();

  EXPECT_FLOAT_EQ(out_ffv.val_.val_.val(), -0.52516294973063);
  EXPECT_FLOAT_EQ(out_ffv.d_.val_.val(), -0.40854102156722);
  EXPECT_FLOAT_EQ(lam_ffv.d_.val_.adj(), -0.40854102156722);

  out_ffv = ordered_logistic_lpmf(y, lam_d, c_ffv);
  out_ffv.d_.val_.grad();

  EXPECT_FLOAT_EQ(out_ffv.val_.val_.val(), -0.52516294973063);
  EXPECT_FLOAT_EQ(out_ffv.d_.val_.val(), 0.40854102156722);
  EXPECT_FLOAT_EQ(c_ffv[0].d_.val_.adj(), 0.40854102156722);
  EXPECT_FLOAT_EQ(c_ffv[1].d_.val_.adj(), 0.0);
  EXPECT_FLOAT_EQ(c_ffv[2].d_.val_.adj(), 0.0);
}

TEST(ProbDistributionsOrdLog, fv_fv_vec) {
  using stan::math::ordered_logistic_lpmf;
  using stan::math::fvar;
  using stan::math::var;
  using stan::math::vector_d;
  using stan::math::vector_ffv;
  using stan::math::vector_fv;
  using stan::math::vector_d;

  std::vector<int> y{1, 2, 3, 4};

  vector_fv lam_fv(4);
  lam_fv << -1.32, -0.05, 0.56, 1.55;
  for (int i = 0; i < 4; i++)
    lam_fv[i].d_ = 1.0;

  vector_fv c_fv(3);
  c_fv << -0.95, -0.10, 0.95;
  for (int i = 0; i < 3; i++)
    c_fv[i].d_ = 1.0;

  vector_ffv lam_ffv(4);
  lam_ffv << -1.32, -0.05, 0.56, 1.55;
  for (int i = 0; i < 4; i++) {
    lam_ffv[i].d_ = 1.0;
    lam_ffv[i].val_.d_ = 1.0;
  }

  vector_ffv c_ffv(3);
  c_ffv << -0.95, -0.10, 0.95;
  for (int i = 0; i < 3; i++) {
    c_ffv[i].d_ = 1.0;
    c_ffv[i].val_.d_ = 1.0;
  }

  fvar<var> out_fv = ordered_logistic_lpmf(y, lam_fv, c_fv);
  out_fv.d_.grad();

  EXPECT_FLOAT_EQ(out_fv.val_.val(), -3.9442226297351447171906);
  EXPECT_FLOAT_EQ(out_fv.d_.val() + 1, 0.0 + 1);

  EXPECT_FLOAT_EQ(lam_fv[0].d_.adj(), -0.40854102156722);
  EXPECT_FLOAT_EQ(lam_fv[1].d_.adj(), -0.223446899109214);
  EXPECT_FLOAT_EQ(lam_fv[2].d_.adj(), -0.062977689154598);
  EXPECT_FLOAT_EQ(lam_fv[3].d_.adj(), 0.354343693774205);

  EXPECT_FLOAT_EQ(c_fv[0].d_.adj(), -0.62697485850346);
  EXPECT_FLOAT_EQ(c_fv[1].d_.adj(), 0.37990895444219);
  EXPECT_FLOAT_EQ(c_fv[2].d_.adj(), 0.587687820118095);

  fvar<fvar<var>> out_ffv = ordered_logistic_lpmf(y, lam_ffv, c_ffv);
  out_ffv.d_.val_.grad();

  EXPECT_FLOAT_EQ(out_ffv.val_.val_.val(), -3.9442226297351447171906);
  EXPECT_FLOAT_EQ(out_ffv.d_.val_.val() + 1, 0.0 + 1);

  EXPECT_FLOAT_EQ(lam_ffv[0].d_.val_.adj(), -0.40854102156722);
  EXPECT_FLOAT_EQ(lam_ffv[1].d_.val_.adj(), -0.223446899109214);
  EXPECT_FLOAT_EQ(lam_ffv[2].d_.val_.adj(), -0.062977689154598);
  EXPECT_FLOAT_EQ(lam_ffv[3].d_.val_.adj(), 0.354343693774205);

  EXPECT_FLOAT_EQ(c_ffv[0].d_.val_.adj(), -0.62697485850346);
  EXPECT_FLOAT_EQ(c_ffv[1].d_.val_.adj(), 0.37990895444219);
  EXPECT_FLOAT_EQ(c_ffv[2].d_.val_.adj(), 0.587687820118095);
}

TEST(ProbDistributionsOrdLog, fv_d_vec) {
  using stan::math::ordered_logistic_lpmf;
  using stan::math::fvar;
  using stan::math::var;
  using stan::math::vector_d;
  using stan::math::vector_ffv;
  using stan::math::vector_fv;
  using stan::math::vector_d;

  std::vector<int> y{1, 2, 3, 4};

  vector_fv lam_fv(4);
  lam_fv << -1.32, -0.05, 0.56, 1.55;
  for (int i = 0; i < 4; i++)
    lam_fv[i].d_ = 1.0;

  vector_fv c_fv(3);
  c_fv << -0.95, -0.10, 0.95;
  for (int i = 0; i < 3; i++)
    c_fv[i].d_ = 1.0;

  vector_ffv lam_ffv(4);
  lam_ffv << -1.32, -0.05, 0.56, 1.55;
  for (int i = 0; i < 4; i++) {
    lam_ffv[i].d_ = 1.0;
    lam_ffv[i].val_.d_ = 1.0;
  }

  vector_ffv c_ffv(3);
  c_ffv << -0.95, -0.10, 0.95;
  for (int i = 0; i < 3; i++) {
    c_ffv[i].d_ = 1.0;
    c_ffv[i].val_.d_ = 1.0;
  }

  vector_d lam_d(4);
  lam_d << -1.32, -0.05, 0.56, 1.55;

  vector_d c_d(3);
  c_d << -0.95, -0.10, 0.95;

  fvar<var> out = ordered_logistic_lpmf(y, lam_fv, c_d);
  out.d_.grad();

  EXPECT_FLOAT_EQ(out.val_.val(), -3.9442226297351447171906);
  EXPECT_FLOAT_EQ(out.d_.val(), -0.340621916056826);

  EXPECT_FLOAT_EQ(lam_fv[0].d_.adj(), -0.40854102156722);
  EXPECT_FLOAT_EQ(lam_fv[1].d_.adj(), -0.223446899109214);
  EXPECT_FLOAT_EQ(lam_fv[2].d_.adj(), -0.062977689154598);
  EXPECT_FLOAT_EQ(lam_fv[3].d_.adj(), 0.354343693774205);

  out = ordered_logistic_lpmf(y, lam_d, c_fv);
  out.d_.grad();

  EXPECT_FLOAT_EQ(out.val_.val(), -3.9442226297351447171906);
  EXPECT_FLOAT_EQ(out.d_.val(), 0.340621916056825);

  EXPECT_FLOAT_EQ(c_fv[0].d_.adj(), -0.62697485850346);
  EXPECT_FLOAT_EQ(c_fv[1].d_.adj(), 0.37990895444219);
  EXPECT_FLOAT_EQ(c_fv[2].d_.adj(), 0.587687820118095);

  fvar<fvar<var>> out_ffv = ordered_logistic_lpmf(y, lam_ffv, c_d);
  out_ffv.d_.val_.grad();

  EXPECT_FLOAT_EQ(out_ffv.val_.val_.val(), -3.9442226297351447171906);
  EXPECT_FLOAT_EQ(out_ffv.d_.val_.val(), -0.340621916056826);

  EXPECT_FLOAT_EQ(lam_ffv[0].d_.val_.adj(), -0.40854102156722);
  EXPECT_FLOAT_EQ(lam_ffv[1].d_.val_.adj(), -0.223446899109214);
  EXPECT_FLOAT_EQ(lam_ffv[2].d_.val_.adj(), -0.062977689154598);
  EXPECT_FLOAT_EQ(lam_ffv[3].d_.val_.adj(), 0.354343693774205);

  out_ffv = ordered_logistic_lpmf(y, lam_d, c_ffv);
  out_ffv.d_.val_.grad();

  EXPECT_FLOAT_EQ(out_ffv.val_.val_.val(), -3.9442226297351447171906);
  EXPECT_FLOAT_EQ(out_ffv.d_.val_.val(), 0.340621916056826);

  EXPECT_FLOAT_EQ(c_ffv[0].d_.val_.adj(), -0.62697485850346);
  EXPECT_FLOAT_EQ(c_ffv[1].d_.val_.adj(), 0.37990895444219);
  EXPECT_FLOAT_EQ(c_ffv[2].d_.val_.adj(), 0.587687820118095);
}

TEST(ProbDistributionsOrdLog, fv_fv_stvec) {
  using stan::math::ordered_logistic_lpmf;
  using stan::math::fvar;
  using stan::math::var;
  using stan::math::vector_d;
  using stan::math::vector_ffv;
  using stan::math::vector_fv;
  using stan::math::vector_d;

  std::vector<int> y{1, 2, 3, 4};

  vector_fv lam_fv(4);
  lam_fv << -1.32, -0.05, 0.56, 1.55;
  for (int i = 0; i < 4; i++)
    lam_fv[i].d_ = 1.0;

  vector_fv c1_fv(3);
  c1_fv << -0.95, -0.10, 0.95;
  for (int i = 0; i < 3; i++)
    c1_fv[i].d_ = 1.0;

  vector_fv c2_fv(3);
  c2_fv << -0.95, -0.10, 0.95;
  for (int i = 0; i < 3; i++)
    c2_fv[i].d_ = 1.0;

  vector_fv c3_fv(3);
  c3_fv << -0.95, -0.10, 0.95;
  for (int i = 0; i < 3; i++)
    c3_fv[i].d_ = 1.0;

  vector_fv c4_fv(3);
  c4_fv << -0.95, -0.10, 0.95;
  for (int i = 0; i < 3; i++)
    c4_fv[i].d_ = 1.0;

  vector_ffv lam_ffv(4);
  lam_ffv << -1.32, -0.05, 0.56, 1.55;
  for (int i = 0; i < 4; i++) {
    lam_ffv[i].d_ = 1.0;
    lam_ffv[i].val_.d_ = 1.0;
  }

  vector_ffv c1_ffv(3);
  c1_ffv << -0.95, -0.10, 0.95;
  for (int i = 0; i < 3; i++) {
    c1_ffv[i].d_ = 1.0;
    c1_ffv[i].val_.d_ = 1.0;
  }

  vector_ffv c2_ffv(3);
  c2_ffv << -0.95, -0.10, 0.95;
  for (int i = 0; i < 3; i++) {
    c2_ffv[i].d_ = 1.0;
    c2_ffv[i].val_.d_ = 1.0;
  }

  vector_ffv c3_ffv(3);
  c3_ffv << -0.95, -0.10, 0.95;
  for (int i = 0; i < 3; i++) {
    c3_ffv[i].d_ = 1.0;
    c3_ffv[i].val_.d_ = 1.0;
  }

  vector_ffv c4_ffv(3);
  c4_ffv << -0.95, -0.10, 0.95;
  for (int i = 0; i < 3; i++) {
    c4_ffv[i].d_ = 1.0;
    c4_ffv[i].val_.d_ = 1.0;
  }

  std::vector<vector_fv> std_c_fv{c1_fv, c2_fv, c3_fv, c4_fv};
  std::vector<vector_ffv> std_c_ffv{c1_ffv, c2_ffv, c3_ffv, c4_ffv};

  fvar<var> out_fv = ordered_logistic_lpmf(y, lam_fv, std_c_fv);
  out_fv.d_.grad();

  EXPECT_FLOAT_EQ(out_fv.val_.val(), -3.9442226297351447171906);
  EXPECT_FLOAT_EQ(out_fv.d_.val() + 1, 0.0 + 1);

  EXPECT_FLOAT_EQ(lam_fv[0].d_.adj(), -0.40854102156722);
  EXPECT_FLOAT_EQ(lam_fv[1].d_.adj(), -0.223446899109214);
  EXPECT_FLOAT_EQ(lam_fv[2].d_.adj(), -0.062977689154598);
  EXPECT_FLOAT_EQ(lam_fv[3].d_.adj(), 0.354343693774205);

  EXPECT_FLOAT_EQ(std_c_fv[0][0].d_.adj(), 0.40854102156722);
  EXPECT_FLOAT_EQ(std_c_fv[0][1].d_.adj(), 0.0);
  EXPECT_FLOAT_EQ(std_c_fv[0][2].d_.adj(), 0.0);

  EXPECT_FLOAT_EQ(std_c_fv[1][0].d_.adj(), -1.03551588007068);
  EXPECT_FLOAT_EQ(std_c_fv[1][1].d_.adj(), 1.25896277917989);
  EXPECT_FLOAT_EQ(std_c_fv[1][2].d_.adj(), 0.0);

  EXPECT_FLOAT_EQ(std_c_fv[2][0].d_.adj(), 0.0);
  EXPECT_FLOAT_EQ(std_c_fv[2][1].d_.adj(), -0.879053824737702);
  EXPECT_FLOAT_EQ(std_c_fv[2][2].d_.adj(), 0.942031513892299);

  EXPECT_FLOAT_EQ(std_c_fv[3][0].d_.adj(), 0.0);
  EXPECT_FLOAT_EQ(std_c_fv[3][1].d_.adj(), 0.0);
  EXPECT_FLOAT_EQ(std_c_fv[3][2].d_.adj(), -0.354343693774205);

  fvar<fvar<var>> out_ffv = ordered_logistic_lpmf(y, lam_ffv, std_c_ffv);
  out_ffv.d_.val_.grad();

  EXPECT_FLOAT_EQ(out_ffv.val_.val_.val(), -3.9442226297351447171906);
  EXPECT_FLOAT_EQ(out_ffv.d_.val_.val() + 1, 0.0 + 1);

  EXPECT_FLOAT_EQ(lam_ffv[0].d_.val_.adj(), -0.40854102156722);
  EXPECT_FLOAT_EQ(lam_ffv[1].d_.val_.adj(), -0.223446899109214);
  EXPECT_FLOAT_EQ(lam_ffv[2].d_.val_.adj(), -0.062977689154598);
  EXPECT_FLOAT_EQ(lam_ffv[3].d_.val_.adj(), 0.354343693774205);

  EXPECT_FLOAT_EQ(std_c_ffv[0][0].d_.val_.adj(), 0.40854102156722);
  EXPECT_FLOAT_EQ(std_c_ffv[0][1].d_.val_.adj(), 0.0);
  EXPECT_FLOAT_EQ(std_c_ffv[0][2].d_.val_.adj(), 0.0);

  EXPECT_FLOAT_EQ(std_c_ffv[1][0].d_.val_.adj(), -1.03551588007068);
  EXPECT_FLOAT_EQ(std_c_ffv[1][1].d_.val_.adj(), 1.25896277917989);
  EXPECT_FLOAT_EQ(std_c_ffv[1][2].d_.val_.adj(), 0.0);

  EXPECT_FLOAT_EQ(std_c_ffv[2][0].d_.val_.adj(), 0.0);
  EXPECT_FLOAT_EQ(std_c_ffv[2][1].d_.val_.adj(), -0.879053824737702);
  EXPECT_FLOAT_EQ(std_c_ffv[2][2].d_.val_.adj(), 0.942031513892299);

  EXPECT_FLOAT_EQ(std_c_ffv[3][0].d_.val_.adj(), 0.0);
  EXPECT_FLOAT_EQ(std_c_ffv[3][1].d_.val_.adj(), 0.0);
  EXPECT_FLOAT_EQ(std_c_ffv[3][2].d_.val_.adj(), -0.354343693774205);
}

TEST(ProbDistributionsOrdLog, fv_d_stvec) {
  using stan::math::ordered_logistic_lpmf;
  using stan::math::fvar;
  using stan::math::var;
  using stan::math::vector_d;
  using stan::math::vector_ffv;
  using stan::math::vector_fv;
  using stan::math::vector_d;

  std::vector<int> y{1, 2, 3, 4};

  vector_fv lam_fv(4);
  lam_fv << -1.32, -0.05, 0.56, 1.55;
  for (int i = 0; i < 4; i++)
    lam_fv[i].d_ = 1.0;

  vector_fv c1_fv(3);
  c1_fv << -0.95, -0.10, 0.95;
  for (int i = 0; i < 3; i++)
    c1_fv[i].d_ = 1.0;

  vector_fv c2_fv(3);
  c2_fv << -0.95, -0.10, 0.95;
  for (int i = 0; i < 3; i++)
    c2_fv[i].d_ = 1.0;

  vector_fv c3_fv(3);
  c3_fv << -0.95, -0.10, 0.95;
  for (int i = 0; i < 3; i++)
    c3_fv[i].d_ = 1.0;

  vector_fv c4_fv(3);
  c4_fv << -0.95, -0.10, 0.95;
  for (int i = 0; i < 3; i++)
    c4_fv[i].d_ = 1.0;

  vector_ffv lam_ffv(4);
  lam_ffv << -1.32, -0.05, 0.56, 1.55;
  for (int i = 0; i < 4; i++) {
    lam_ffv[i].d_ = 1.0;
    lam_ffv[i].val_.d_ = 1.0;
  }

  vector_ffv c1_ffv(3);
  c1_ffv << -0.95, -0.10, 0.95;
  for (int i = 0; i < 3; i++) {
    c1_ffv[i].d_ = 1.0;
    c1_ffv[i].val_.d_ = 1.0;
  }

  vector_ffv c2_ffv(3);
  c2_ffv << -0.95, -0.10, 0.95;
  for (int i = 0; i < 3; i++) {
    c2_ffv[i].d_ = 1.0;
    c2_ffv[i].val_.d_ = 1.0;
  }

  vector_ffv c3_ffv(3);
  c3_ffv << -0.95, -0.10, 0.95;
  for (int i = 0; i < 3; i++) {
    c3_ffv[i].d_ = 1.0;
    c3_ffv[i].val_.d_ = 1.0;
  }

  vector_ffv c4_ffv(3);
  c4_ffv << -0.95, -0.10, 0.95;
  for (int i = 0; i < 3; i++) {
    c4_ffv[i].d_ = 1.0;
    c4_ffv[i].val_.d_ = 1.0;
  }

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

  std::vector<vector_fv> std_c_fv{c1_fv, c2_fv, c3_fv, c4_fv};
  std::vector<vector_ffv> std_c_ffv{c1_ffv, c2_ffv, c3_ffv, c4_ffv};
  std::vector<vector_d> std_c_d{c1_d, c2_d, c3_d, c4_d};

  fvar<var> out_fv = ordered_logistic_lpmf(y, lam_fv, std_c_d);
  out_fv.d_.grad();

  EXPECT_FLOAT_EQ(out_fv.val_.val(), -3.9442226297351447171906);
  EXPECT_FLOAT_EQ(out_fv.d_.val(), -0.340621916056822);

  EXPECT_FLOAT_EQ(lam_fv[0].d_.adj(), -0.40854102156722);
  EXPECT_FLOAT_EQ(lam_fv[1].d_.adj(), -0.223446899109214);
  EXPECT_FLOAT_EQ(lam_fv[2].d_.adj(), -0.062977689154598);
  EXPECT_FLOAT_EQ(lam_fv[3].d_.adj(), 0.354343693774205);

  out_fv = ordered_logistic_lpmf(y, lam_d, std_c_fv);
  out_fv.d_.grad();

  EXPECT_FLOAT_EQ(out_fv.val_.val(), -3.9442226297351447171906);
  EXPECT_FLOAT_EQ(out_fv.d_.val(), 0.340621916056822);

  EXPECT_FLOAT_EQ(std_c_fv[0][0].d_.adj(), 0.40854102156722);
  EXPECT_FLOAT_EQ(std_c_fv[0][1].d_.adj(), 0.0);
  EXPECT_FLOAT_EQ(std_c_fv[0][2].d_.adj(), 0.0);

  EXPECT_FLOAT_EQ(std_c_fv[1][0].d_.adj(), -1.03551588007068);
  EXPECT_FLOAT_EQ(std_c_fv[1][1].d_.adj(), 1.25896277917989);
  EXPECT_FLOAT_EQ(std_c_fv[1][2].d_.adj(), 0.0);

  EXPECT_FLOAT_EQ(std_c_fv[2][0].d_.adj(), 0.0);
  EXPECT_FLOAT_EQ(std_c_fv[2][1].d_.adj(), -0.879053824737702);
  EXPECT_FLOAT_EQ(std_c_fv[2][2].d_.adj(), 0.942031513892299);

  EXPECT_FLOAT_EQ(std_c_fv[3][0].d_.adj(), 0.0);
  EXPECT_FLOAT_EQ(std_c_fv[3][1].d_.adj(), 0.0);
  EXPECT_FLOAT_EQ(std_c_fv[3][2].d_.adj(), -0.354343693774205);

  fvar<fvar<var>> out_ffv = ordered_logistic_lpmf(y, lam_ffv, std_c_d);
  out_ffv.d_.val_.grad();

  EXPECT_FLOAT_EQ(out_ffv.val_.val_.val(), -3.9442226297351447171906);
  EXPECT_FLOAT_EQ(out_ffv.d_.val_.val(), -0.340621916056822);

  EXPECT_FLOAT_EQ(lam_ffv[0].d_.val_.adj(), -0.40854102156722);
  EXPECT_FLOAT_EQ(lam_ffv[1].d_.val_.adj(), -0.223446899109214);
  EXPECT_FLOAT_EQ(lam_ffv[2].d_.val_.adj(), -0.062977689154598);
  EXPECT_FLOAT_EQ(lam_ffv[3].d_.val_.adj(), 0.354343693774205);

  out_ffv = ordered_logistic_lpmf(y, lam_d, std_c_ffv);
  out_ffv.d_.val_.grad();

  EXPECT_FLOAT_EQ(out_ffv.val_.val_.val(), -3.9442226297351447171906);
  EXPECT_FLOAT_EQ(out_ffv.d_.val_.val(), 0.340621916056822);

  EXPECT_FLOAT_EQ(std_c_ffv[0][0].d_.val_.adj(), 0.40854102156722);
  EXPECT_FLOAT_EQ(std_c_ffv[0][1].d_.val_.adj(), 0.0);
  EXPECT_FLOAT_EQ(std_c_ffv[0][2].d_.val_.adj(), 0.0);

  EXPECT_FLOAT_EQ(std_c_ffv[1][0].d_.val_.adj(), -1.03551588007068);
  EXPECT_FLOAT_EQ(std_c_ffv[1][1].d_.val_.adj(), 1.25896277917989);
  EXPECT_FLOAT_EQ(std_c_ffv[1][2].d_.val_.adj(), 0.0);

  EXPECT_FLOAT_EQ(std_c_ffv[2][0].d_.val_.adj(), 0.0);
  EXPECT_FLOAT_EQ(std_c_ffv[2][1].d_.val_.adj(), -0.879053824737702);
  EXPECT_FLOAT_EQ(std_c_ffv[2][2].d_.val_.adj(), 0.942031513892299);

  EXPECT_FLOAT_EQ(std_c_ffv[3][0].d_.val_.adj(), 0.0);
  EXPECT_FLOAT_EQ(std_c_ffv[3][1].d_.val_.adj(), 0.0);
  EXPECT_FLOAT_EQ(std_c_ffv[3][2].d_.val_.adj(), -0.354343693774205);
}
