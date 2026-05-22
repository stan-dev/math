#include <stan/math/mix.hpp>
#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/rev/fun/util.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST_F(AgradRev, ProbDistributionsOrdLog_fv_fv) {
  using stan::math::fvar;
  using stan::math::ordered_logistic_lpmf;
  using stan::math::var;
  using stan::math::vector_d;
  using stan::math::vector_ffv;
  using stan::math::vector_fv;

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

TEST_F(AgradRev, ProbDistributionsOrdLog_fv_d) {
  using stan::math::fvar;
  using stan::math::ordered_logistic_lpmf;
  using stan::math::var;
  using stan::math::vector_d;
  using stan::math::vector_ffv;
  using stan::math::vector_fv;

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

TEST_F(AgradRev, ProbDistributionsOrdLog_fv_fv_vec) {
  using stan::math::fvar;
  using stan::math::ordered_logistic_lpmf;
  using stan::math::var;
  using stan::math::vector_d;
  using stan::math::vector_ffv;
  using stan::math::vector_fv;

  std::vector<int> y{1, 2, 3, 4};

  vector_fv lam_fv(4);
  lam_fv << 1.25, -0.33, 1.36, 2.11;
  for (int i = 0; i < 4; i++)
    lam_fv[i].d_ = 1.0;

  vector_fv c_fv(3);
  c_fv << -1.21, -1.01, 0.90;
  for (int i = 0; i < 3; i++)
    c_fv[i].d_ = 1.0;

  vector_ffv lam_ffv(4);
  lam_ffv << 1.25, -0.33, 1.36, 2.11;
  for (int i = 0; i < 4; i++) {
    lam_ffv[i].d_ = 1.0;
    lam_ffv[i].val_.d_ = 1.0;
  }

  vector_ffv c_ffv(3);
  c_ffv << -1.21, -1.01, 0.90;
  for (int i = 0; i < 3; i++) {
    c_ffv[i].d_ = 1.0;
    c_ffv[i].val_.d_ = 1.0;
  }

  fvar<var> out_fv = ordered_logistic_lpmf(y, lam_fv, c_fv);
  out_fv.d_.grad();

  EXPECT_FLOAT_EQ(out_fv.val_.val(), -7.14656827285528);
  EXPECT_FLOAT_EQ(out_fv.d_.val() + 1, 0.0 + 1);

  EXPECT_FLOAT_EQ(lam_fv[0].d_.adj(), -0.921289662829465);
  EXPECT_FLOAT_EQ(lam_fv[1].d_.adj(), -0.370560918497922);
  EXPECT_FLOAT_EQ(lam_fv[2].d_.adj(), -0.527525036704529);
  EXPECT_FLOAT_EQ(lam_fv[3].d_.adj(), 0.229701050953398);

  EXPECT_FLOAT_EQ(c_fv[0].d_.adj(), -3.88854368220396);
  EXPECT_FLOAT_EQ(c_fv[1].d_.adj(), 4.92108545347799);
  EXPECT_FLOAT_EQ(c_fv[2].d_.adj(), 0.557132795804491);

  fvar<fvar<var>> out_ffv = ordered_logistic_lpmf(y, lam_ffv, c_ffv);
  out_ffv.d_.val_.grad();

  EXPECT_FLOAT_EQ(out_ffv.val_.val_.val(), -7.14656827285528);
  EXPECT_FLOAT_EQ(out_ffv.d_.val_.val() + 1, 0.0 + 1);

  EXPECT_FLOAT_EQ(lam_ffv[0].d_.val_.adj(), -0.921289662829465);
  EXPECT_FLOAT_EQ(lam_ffv[1].d_.val_.adj(), -0.370560918497922);
  EXPECT_FLOAT_EQ(lam_ffv[2].d_.val_.adj(), -0.527525036704529);
  EXPECT_FLOAT_EQ(lam_ffv[3].d_.val_.adj(), 0.229701050953398);

  EXPECT_FLOAT_EQ(c_ffv[0].d_.val_.adj(), -3.88854368220396);
  EXPECT_FLOAT_EQ(c_ffv[1].d_.val_.adj(), 4.92108545347799);
  EXPECT_FLOAT_EQ(c_ffv[2].d_.val_.adj(), 0.557132795804491);
}

TEST_F(AgradRev, ProbDistributionsOrdLog_fv_d_vec) {
  using stan::math::fvar;
  using stan::math::ordered_logistic_lpmf;
  using stan::math::var;
  using stan::math::vector_d;
  using stan::math::vector_ffv;
  using stan::math::vector_fv;

  std::vector<int> y{1, 2, 3, 4};

  vector_fv lam_fv(4);
  lam_fv << 1.25, -0.33, 1.36, 2.11;
  for (int i = 0; i < 4; i++)
    lam_fv[i].d_ = 1.0;

  vector_fv c_fv(3);
  c_fv << -2.22, -1.55, -0.36;
  for (int i = 0; i < 3; i++)
    c_fv[i].d_ = 1.0;

  vector_ffv lam_ffv(4);
  lam_ffv << 1.25, -0.33, 1.36, 2.11;
  for (int i = 0; i < 4; i++) {
    lam_ffv[i].d_ = 1.0;
    lam_ffv[i].val_.d_ = 1.0;
  }

  vector_ffv c_ffv(3);
  c_ffv << -2.22, -1.55, -0.36;
  for (int i = 0; i < 3; i++) {
    c_ffv[i].d_ = 1.0;
    c_ffv[i].val_.d_ = 1.0;
  }

  vector_d lam_d(4);
  lam_d << 1.25, -0.33, 1.36, 2.11;

  vector_d c_d(3);
  c_d << -2.22, -1.55, -0.36;

  fvar<var> out = ordered_logistic_lpmf(y, lam_fv, c_d);
  out.d_.grad();

  EXPECT_FLOAT_EQ(out.val_.val(), -8.21855481819114);
  EXPECT_FLOAT_EQ(out.d_.val(), -2.32912026425027);

  EXPECT_FLOAT_EQ(lam_fv[0].d_.adj(), -0.969822018514124);
  EXPECT_FLOAT_EQ(lam_fv[1].d_.adj(), -0.640819079988261);
  EXPECT_FLOAT_EQ(lam_fv[2].d_.adj(), -0.796467400877256);
  EXPECT_FLOAT_EQ(lam_fv[3].d_.adj(), 0.077988235129366);

  out = ordered_logistic_lpmf(y, lam_d, c_fv);
  out.d_.grad();

  EXPECT_FLOAT_EQ(out.val_.val(), -8.21855481819114);
  EXPECT_FLOAT_EQ(out.d_.val(), 2.32912026425027);

  EXPECT_FLOAT_EQ(c_fv[0].d_.adj(), -0.209379786457621);
  EXPECT_FLOAT_EQ(c_fv[1].d_.adj(), 1.33112093047159);
  EXPECT_FLOAT_EQ(c_fv[2].d_.adj(), 1.20737912023631);

  fvar<fvar<var>> out_ffv = ordered_logistic_lpmf(y, lam_ffv, c_d);
  out_ffv.d_.val_.grad();

  EXPECT_FLOAT_EQ(out_ffv.val_.val_.val(), -8.21855481819114);
  EXPECT_FLOAT_EQ(out_ffv.d_.val_.val(), -2.32912026425027);

  EXPECT_FLOAT_EQ(lam_ffv[0].d_.val_.adj(), -0.969822018514124);
  EXPECT_FLOAT_EQ(lam_ffv[1].d_.val_.adj(), -0.640819079988261);
  EXPECT_FLOAT_EQ(lam_ffv[2].d_.val_.adj(), -0.796467400877256);
  EXPECT_FLOAT_EQ(lam_ffv[3].d_.val_.adj(), 0.077988235129366);

  out_ffv = ordered_logistic_lpmf(y, lam_d, c_ffv);
  out_ffv.d_.val_.grad();

  EXPECT_FLOAT_EQ(out_ffv.val_.val_.val(), -8.21855481819114);
  EXPECT_FLOAT_EQ(out_ffv.d_.val_.val(), 2.32912026425027);

  EXPECT_FLOAT_EQ(c_ffv[0].d_.val_.adj(), -0.209379786457621);
  EXPECT_FLOAT_EQ(c_ffv[1].d_.val_.adj(), 1.33112093047159);
  EXPECT_FLOAT_EQ(c_ffv[2].d_.val_.adj(), 1.20737912023631);
}

TEST_F(AgradRev, ProbDistributionsOrdLog_fv_fv_stvec) {
  using stan::math::fvar;
  using stan::math::ordered_logistic_lpmf;
  using stan::math::var;
  using stan::math::vector_d;
  using stan::math::vector_ffv;
  using stan::math::vector_fv;

  std::vector<int> y{1, 2, 3, 4};

  vector_fv lam_fv(4);
  lam_fv << 0.61, 2.63, -0.06, 1.04;
  for (int i = 0; i < 4; i++)
    lam_fv[i].d_ = 1.0;

  vector_fv c1_fv(3);
  c1_fv << -2.58, -1.66, -0.64;
  for (int i = 0; i < 3; i++)
    c1_fv[i].d_ = 1.0;

  vector_fv c2_fv(3);
  c2_fv << -1.20, 0.22, 1.34;
  for (int i = 0; i < 3; i++)
    c2_fv[i].d_ = 1.0;

  vector_fv c3_fv(3);
  c3_fv << -1.68, -0.28, 1.33;
  for (int i = 0; i < 3; i++)
    c3_fv[i].d_ = 1.0;

  vector_fv c4_fv(3);
  c4_fv << -2.51, -0.64, 1.03;
  for (int i = 0; i < 3; i++)
    c4_fv[i].d_ = 1.0;

  vector_ffv lam_ffv(4);
  lam_ffv << 0.61, 2.63, -0.06, 1.04;
  for (int i = 0; i < 4; i++) {
    lam_ffv[i].d_ = 1.0;
    lam_ffv[i].val_.d_ = 1.0;
  }

  vector_ffv c1_ffv(3);
  c1_ffv << -2.58, -1.66, -0.64;
  for (int i = 0; i < 3; i++) {
    c1_ffv[i].d_ = 1.0;
    c1_ffv[i].val_.d_ = 1.0;
  }

  vector_ffv c2_ffv(3);
  c2_ffv << -1.20, 0.22, 1.34;
  for (int i = 0; i < 3; i++) {
    c2_ffv[i].d_ = 1.0;
    c2_ffv[i].val_.d_ = 1.0;
  }

  vector_ffv c3_ffv(3);
  c3_ffv << -1.68, -0.28, 1.33;
  for (int i = 0; i < 3; i++) {
    c3_ffv[i].d_ = 1.0;
    c3_ffv[i].val_.d_ = 1.0;
  }

  vector_ffv c4_ffv(3);
  c4_ffv << -2.51, -0.64, 1.03;
  for (int i = 0; i < 3; i++) {
    c4_ffv[i].d_ = 1.0;
    c4_ffv[i].val_.d_ = 1.0;
  }

  std::vector<vector_fv> std_c_fv{c1_fv, c2_fv, c3_fv, c4_fv};
  std::vector<vector_ffv> std_c_ffv{c1_ffv, c2_ffv, c3_ffv, c4_ffv};

  fvar<var> out_fv = ordered_logistic_lpmf(y, lam_fv, std_c_fv);
  out_fv.d_.grad();

  EXPECT_FLOAT_EQ(out_fv.val_.val(), -7.74727840068321);
  EXPECT_FLOAT_EQ(out_fv.d_.val() + 1, 0.0 + 1);

  EXPECT_FLOAT_EQ(lam_fv[0].d_.adj(), -0.960456220449047);
  EXPECT_FLOAT_EQ(lam_fv[1].d_.adj(), -0.896338359161074);
  EXPECT_FLOAT_EQ(lam_fv[2].d_.adj(), 0.245813008044117);
  EXPECT_FLOAT_EQ(lam_fv[3].d_.adj(), 0.497500020833125);

  EXPECT_FLOAT_EQ(std_c_fv[0][0].d_.adj(), 0.960456220449047);
  EXPECT_FLOAT_EQ(std_c_fv[0][1].d_.adj(), 0.0);
  EXPECT_FLOAT_EQ(std_c_fv[0][2].d_.adj(), 0.0);

  EXPECT_FLOAT_EQ(std_c_fv[1][0].d_.adj(), -0.340011984816337);
  EXPECT_FLOAT_EQ(std_c_fv[1][1].d_.adj(), 1.23635034397741);
  EXPECT_FLOAT_EQ(std_c_fv[1][2].d_.adj(), 0.0);

  EXPECT_FLOAT_EQ(std_c_fv[2][0].d_.adj(), 0.0);
  EXPECT_FLOAT_EQ(std_c_fv[2][1].d_.adj(), -0.695045186550866);
  EXPECT_FLOAT_EQ(std_c_fv[2][2].d_.adj(), 0.44923217850675);

  EXPECT_FLOAT_EQ(std_c_fv[3][0].d_.adj(), 0.0);
  EXPECT_FLOAT_EQ(std_c_fv[3][1].d_.adj(), 0.0);
  EXPECT_FLOAT_EQ(std_c_fv[3][2].d_.adj(), -0.497500020833125);

  fvar<fvar<var>> out_ffv = ordered_logistic_lpmf(y, lam_ffv, std_c_ffv);
  out_ffv.d_.val_.grad();

  EXPECT_FLOAT_EQ(out_ffv.val_.val_.val(), -7.74727840068321);
  EXPECT_FLOAT_EQ(out_ffv.d_.val_.val() + 1, 0.0 + 1);

  EXPECT_FLOAT_EQ(lam_ffv[0].d_.val_.adj(), -0.960456220449047);
  EXPECT_FLOAT_EQ(lam_ffv[1].d_.val_.adj(), -0.896338359161074);
  EXPECT_FLOAT_EQ(lam_ffv[2].d_.val_.adj(), 0.245813008044117);
  EXPECT_FLOAT_EQ(lam_ffv[3].d_.val_.adj(), 0.497500020833125);

  EXPECT_FLOAT_EQ(std_c_ffv[0][0].d_.val_.adj(), 0.960456220449047);
  EXPECT_FLOAT_EQ(std_c_ffv[0][1].d_.val_.adj(), 0.0);
  EXPECT_FLOAT_EQ(std_c_ffv[0][2].d_.val_.adj(), 0.0);

  EXPECT_FLOAT_EQ(std_c_ffv[1][0].d_.val_.adj(), -0.340011984816337);
  EXPECT_FLOAT_EQ(std_c_ffv[1][1].d_.val_.adj(), 1.23635034397741);
  EXPECT_FLOAT_EQ(std_c_ffv[1][2].d_.val_.adj(), 0.0);

  EXPECT_FLOAT_EQ(std_c_ffv[2][0].d_.val_.adj(), 0.0);
  EXPECT_FLOAT_EQ(std_c_ffv[2][1].d_.val_.adj(), -0.695045186550866);
  EXPECT_FLOAT_EQ(std_c_ffv[2][2].d_.val_.adj(), 0.44923217850675);

  EXPECT_FLOAT_EQ(std_c_ffv[3][0].d_.val_.adj(), 0.0);
  EXPECT_FLOAT_EQ(std_c_ffv[3][1].d_.val_.adj(), 0.0);
  EXPECT_FLOAT_EQ(std_c_ffv[3][2].d_.val_.adj(), -0.497500020833125);
}

TEST_F(AgradRev, ProbDistributionsOrdLog_fv_d_stvec) {
  using stan::math::fvar;
  using stan::math::ordered_logistic_lpmf;
  using stan::math::var;
  using stan::math::vector_d;
  using stan::math::vector_ffv;
  using stan::math::vector_fv;

  std::vector<int> y{1, 2, 3, 4};

  vector_fv lam_fv(4);
  lam_fv << -3.18, -2.06, 0.52, 1.82;
  for (int i = 0; i < 4; i++)
    lam_fv[i].d_ = 1.0;

  vector_fv c1_fv(3);
  c1_fv << -1.02, -0.13, 0.86;
  for (int i = 0; i < 3; i++)
    c1_fv[i].d_ = 1.0;

  vector_fv c2_fv(3);
  c2_fv << -2.38, -1.80, -0.60;
  for (int i = 0; i < 3; i++)
    c2_fv[i].d_ = 1.0;

  vector_fv c3_fv(3);
  c3_fv << -0.61, 0.25, 1.36;
  for (int i = 0; i < 3; i++)
    c3_fv[i].d_ = 1.0;

  vector_fv c4_fv(3);
  c4_fv << -1.07, -0.37, 2.69;
  for (int i = 0; i < 3; i++)
    c4_fv[i].d_ = 1.0;

  vector_ffv lam_ffv(4);
  lam_ffv << -3.18, -2.06, 0.52, 1.82;
  for (int i = 0; i < 4; i++) {
    lam_ffv[i].d_ = 1.0;
    lam_ffv[i].val_.d_ = 1.0;
  }

  vector_ffv c1_ffv(3);
  c1_ffv << -1.02, -0.13, 0.86;
  for (int i = 0; i < 3; i++) {
    c1_ffv[i].d_ = 1.0;
    c1_ffv[i].val_.d_ = 1.0;
  }

  vector_ffv c2_ffv(3);
  c2_ffv << -2.38, -1.80, -0.60;
  for (int i = 0; i < 3; i++) {
    c2_ffv[i].d_ = 1.0;
    c2_ffv[i].val_.d_ = 1.0;
  }

  vector_ffv c3_ffv(3);
  c3_ffv << -0.61, 0.25, 1.36;
  for (int i = 0; i < 3; i++) {
    c3_ffv[i].d_ = 1.0;
    c3_ffv[i].val_.d_ = 1.0;
  }

  vector_ffv c4_ffv(3);
  c4_ffv << -1.07, -0.37, 2.69;
  for (int i = 0; i < 3; i++) {
    c4_ffv[i].d_ = 1.0;
    c4_ffv[i].val_.d_ = 1.0;
  }

  vector_d lam_d(4);
  lam_d << -3.18, -2.06, 0.52, 1.82;

  vector_d c1_d(3);
  c1_d << -1.02, -0.13, 0.86;

  vector_d c2_d(3);
  c2_d << -2.38, -1.80, -0.60;

  vector_d c3_d(3);
  c3_d << -0.61, 0.25, 1.36;

  vector_d c4_d(3);
  c4_d << -1.07, -0.37, 2.69;

  std::vector<vector_fv> std_c_fv{c1_fv, c2_fv, c3_fv, c4_fv};
  std::vector<vector_ffv> std_c_ffv{c1_ffv, c2_ffv, c3_ffv, c4_ffv};
  std::vector<vector_d> std_c_d{c1_d, c2_d, c3_d, c4_d};

  fvar<var> out_fv = ordered_logistic_lpmf(y, lam_fv, std_c_d);
  out_fv.d_.grad();

  EXPECT_FLOAT_EQ(out_fv.val_.val(), -4.59320177226145);
  EXPECT_FLOAT_EQ(out_fv.d_.val(), 0.718029597231206);

  EXPECT_FLOAT_EQ(lam_fv[0].d_.adj(), -0.10340045145825);
  EXPECT_FLOAT_EQ(lam_fv[1].d_.adj(), -0.01468796034572);
  EXPECT_FLOAT_EQ(lam_fv[2].d_.adj(), 0.131372311037084);
  EXPECT_FLOAT_EQ(lam_fv[3].d_.adj(), 0.704745697998091);

  out_fv = ordered_logistic_lpmf(y, lam_d, std_c_fv);
  out_fv.d_.grad();

  EXPECT_FLOAT_EQ(out_fv.val_.val(), -4.59320177226145);
  EXPECT_FLOAT_EQ(out_fv.d_.val(), -0.718029597231206);

  EXPECT_FLOAT_EQ(std_c_fv[0][0].d_.adj(), 0.10340045145825);
  EXPECT_FLOAT_EQ(std_c_fv[0][1].d_.adj(), 0.0);
  EXPECT_FLOAT_EQ(std_c_fv[0][2].d_.adj(), 0.0);

  EXPECT_FLOAT_EQ(std_c_fv[1][0].d_.adj(), -1.69287817572205);
  EXPECT_FLOAT_EQ(std_c_fv[1][1].d_.adj(), 1.70756613606778);
  EXPECT_FLOAT_EQ(std_c_fv[1][2].d_.adj(), 0.0);

  EXPECT_FLOAT_EQ(std_c_fv[2][0].d_.adj(), 0.0);
  EXPECT_FLOAT_EQ(std_c_fv[2][1].d_.adj(), -0.924462566644255);
  EXPECT_FLOAT_EQ(std_c_fv[2][2].d_.adj(), 0.793090255607171);

  EXPECT_FLOAT_EQ(std_c_fv[3][0].d_.adj(), 0.0);
  EXPECT_FLOAT_EQ(std_c_fv[3][1].d_.adj(), 0.0);
  EXPECT_FLOAT_EQ(std_c_fv[3][2].d_.adj(), -0.704745697998091);

  fvar<fvar<var>> out_ffv = ordered_logistic_lpmf(y, lam_ffv, std_c_d);
  out_ffv.d_.val_.grad();

  EXPECT_FLOAT_EQ(out_ffv.val_.val_.val(), -4.59320177226145);
  EXPECT_FLOAT_EQ(out_ffv.d_.val_.val(), 0.718029597231206);

  EXPECT_FLOAT_EQ(lam_ffv[0].d_.val_.adj(), -0.10340045145825);
  EXPECT_FLOAT_EQ(lam_ffv[1].d_.val_.adj(), -0.01468796034572);
  EXPECT_FLOAT_EQ(lam_ffv[2].d_.val_.adj(), 0.131372311037084);
  EXPECT_FLOAT_EQ(lam_ffv[3].d_.val_.adj(), 0.704745697998091);

  out_ffv = ordered_logistic_lpmf(y, lam_d, std_c_ffv);
  out_ffv.d_.val_.grad();

  EXPECT_FLOAT_EQ(out_ffv.val_.val_.val(), -4.59320177226145);
  EXPECT_FLOAT_EQ(out_ffv.d_.val_.val(), -0.718029597231206);

  EXPECT_FLOAT_EQ(std_c_ffv[0][0].d_.val_.adj(), 0.10340045145825);
  EXPECT_FLOAT_EQ(std_c_ffv[0][1].d_.val_.adj(), 0.0);
  EXPECT_FLOAT_EQ(std_c_ffv[0][2].d_.val_.adj(), 0.0);

  EXPECT_FLOAT_EQ(std_c_ffv[1][0].d_.val_.adj(), -1.69287817572205);
  EXPECT_FLOAT_EQ(std_c_ffv[1][1].d_.val_.adj(), 1.70756613606778);
  EXPECT_FLOAT_EQ(std_c_ffv[1][2].d_.val_.adj(), 0.0);

  EXPECT_FLOAT_EQ(std_c_ffv[2][0].d_.val_.adj(), 0.0);
  EXPECT_FLOAT_EQ(std_c_ffv[2][1].d_.val_.adj(), -0.924462566644255);
  EXPECT_FLOAT_EQ(std_c_ffv[2][2].d_.val_.adj(), 0.793090255607171);

  EXPECT_FLOAT_EQ(std_c_ffv[3][0].d_.val_.adj(), 0.0);
  EXPECT_FLOAT_EQ(std_c_ffv[3][1].d_.val_.adj(), 0.0);
  EXPECT_FLOAT_EQ(std_c_ffv[3][2].d_.val_.adj(), -0.704745697998091);
}

// Regression test for the higher-order-AD NaN bug reported 2026-05-21:
// `ordered_logistic_lpmf` returned NaN under `laplace_marginal_tol` when y
// contained the top category K. Root cause: the partials block materialized
// `exp(-cut1)` and `exp(cut1)` (in an unselected branch of Eigen .select),
// which for y == K (cut1 = -INF) and y == 1 (cut2 = +INF) pushed an orphan
// `vari` with val == +INFINITY onto the autodiff stack. `stan::math::grad()`
// then ran the standard `exp_vari::chain` body
//   a.adj() += vi.adj() * vi.val();
// with vi.adj() == 0 for the orphan (e.g. in laplace_marginal_density's
// compute_s2 where most tangent seeds are zero), producing 0 * INF = NaN on
// the val-side adjoint chain of the input. The pre-fix mix tests only
// asserted `.d_.adj()` (tangent-side) with all `.d_` seeds equal to 1, so
// they hit 1 * INF = +INF (finite-ish, escaping the EXPECT_FLOAT_EQ) and
// never inspected the val-side adjoint where the NaN actually lands.
//
// The tests below specifically exercise:
//   1. y contains K (top category): triggers cut1 = -INF.
//   2. y contains 1 (bottom category): triggers cut2 = +INF.
//   3. Some .d_ tangent seeds are zero: triggers 0 * INF = NaN in chain().
//   4. Val-side adjoints (`.val_.val_.adj()` for fvar<fvar<var>>) are
//      asserted finite \u2014 these are the adjoints that Laplace reads.
TEST_F(AgradRev, ProbDistributionsOrdLog_top_category_higher_order_ad) {
  using stan::math::fvar;
  using stan::math::ordered_logistic_lpmf;
  using stan::math::var;
  using stan::math::vector_ffv;

  // K = 3 (cutpoints has size 2). y contains both 1 (bottom) and K=3 (top).
  std::vector<int> y{1, 2, 3, 2, 1};

  vector_ffv lam_ffv(5);
  lam_ffv << -0.5, 0.2, 0.7, 0.1, -0.3;
  // Mixed tangent seeds: some zero, some non-zero. This matches Laplace's
  // compute_s2 pattern where v(j) = 0 for most indices when
  // hessian_block_size > 1. Zero seeds combined with an orphan +INF vari
  // produce 0 * INF = NaN in the buggy version.
  lam_ffv[0].d_ = 1.0;
  lam_ffv[0].val_.d_ = 1.0;
  lam_ffv[1].d_ = 0.0;  // zero outer tangent
  lam_ffv[1].val_.d_ = 1.0;
  lam_ffv[2].d_ = 1.0;
  lam_ffv[2].val_.d_ = 0.0;  // zero inner tangent
  lam_ffv[3].d_ = 0.0;
  lam_ffv[3].val_.d_ = 0.0;  // both zero
  lam_ffv[4].d_ = 1.0;
  lam_ffv[4].val_.d_ = 1.0;

  vector_ffv c_ffv(2);
  c_ffv << 0.0, 1.0;
  for (int i = 0; i < 2; i++) {
    c_ffv[i].d_ = 1.0;
    c_ffv[i].val_.d_ = 1.0;
  }

  fvar<fvar<var>> out_ffv = ordered_logistic_lpmf(y, lam_ffv, c_ffv);

  // Value must be finite.
  EXPECT_TRUE(std::isfinite(out_ffv.val_.val_.val()))
      << "Forward value is non-finite when y contains top category K";
  EXPECT_TRUE(std::isfinite(out_ffv.d_.val_.val()))
      << "First-order tangent value is non-finite";
  EXPECT_TRUE(std::isfinite(out_ffv.val_.d_.val()))
      << "Inner tangent value is non-finite";
  EXPECT_TRUE(std::isfinite(out_ffv.d_.d_.val()))
      << "Second-order tangent value is non-finite";

  // Hessian-style grad: differentiate the second-order tangent. This is the
  // exact code path used by laplace_likelihood::compute_s2 /
  // laplace_likelihood::third_diff.
  out_ffv.d_.d_.grad();

  for (int i = 0; i < 5; i++) {
    // VAL-SIDE adjoints (the ones the pre-fix tests never checked):
    EXPECT_TRUE(std::isfinite(lam_ffv[i].val_.val_.adj()))
        << "lam_ffv[" << i << "].val_.val_.adj() is NaN/INF -- "
        << "regression: orphan exp(+INF) vari injected NaN into val-side "
        << "chain. y[i]=" << y[i] << ", lam.d_=" << lam_ffv[i].d_.val_.val()
        << ", lam.val_.d_=" << lam_ffv[i].val_.d_.val();
    EXPECT_TRUE(std::isfinite(lam_ffv[i].d_.val_.adj()));
    EXPECT_TRUE(std::isfinite(lam_ffv[i].val_.d_.adj()));
    EXPECT_TRUE(std::isfinite(lam_ffv[i].d_.d_.adj()));
  }
  for (int i = 0; i < 2; i++) {
    EXPECT_TRUE(std::isfinite(c_ffv[i].val_.val_.adj()))
        << "c_ffv[" << i << "].val_.val_.adj() is NaN/INF";
    EXPECT_TRUE(std::isfinite(c_ffv[i].d_.val_.adj()));
    EXPECT_TRUE(std::isfinite(c_ffv[i].val_.d_.adj()));
    EXPECT_TRUE(std::isfinite(c_ffv[i].d_.d_.adj()));
  }
}

// Same regression test but for scalar y, scalar lambda \u2014 mirrors the
// production failure pattern in mixtureis/test_ordered_logistic_laplace_bug
// where ordered_logistic_lpmf is called once per observation with y == K.
TEST_F(AgradRev, ProbDistributionsOrdLog_scalar_top_category_higher_order_ad) {
  using stan::math::fvar;
  using stan::math::ordered_logistic_lpmf;
  using stan::math::var;
  using stan::math::vector_ffv;

  // K = 3, y = K (top category) \u2014 the case that produced NaN.
  int y_top = 3;
  int y_bot = 1;

  fvar<fvar<var>> lam_ffv;
  lam_ffv.val_ = 0.7;
  lam_ffv.d_ = 0.0;          // zero outer tangent (compute_s2 pattern)
  lam_ffv.val_.d_ = 1.0;

  vector_ffv c_ffv(2);
  c_ffv << 0.0, 1.0;
  for (int i = 0; i < 2; i++) {
    c_ffv[i].d_ = 0.0;
    c_ffv[i].val_.d_ = 1.0;
  }

  fvar<fvar<var>> out_top = ordered_logistic_lpmf(y_top, lam_ffv, c_ffv);
  EXPECT_TRUE(std::isfinite(out_top.val_.val_.val()));
  EXPECT_TRUE(std::isfinite(out_top.d_.val_.val()));
  EXPECT_TRUE(std::isfinite(out_top.val_.d_.val()));
  EXPECT_TRUE(std::isfinite(out_top.d_.d_.val()));

  out_top.d_.d_.grad();
  EXPECT_TRUE(std::isfinite(lam_ffv.val_.val_.adj()))
      << "scalar y == K: val-side adjoint of lambda is non-finite";
  EXPECT_TRUE(std::isfinite(c_ffv[0].val_.val_.adj()));
  EXPECT_TRUE(std::isfinite(c_ffv[1].val_.val_.adj()));

  // Same for y == 1 (bottom category), which triggers cut2 = +INF.
  stan::math::set_zero_all_adjoints();
  fvar<fvar<var>> out_bot = ordered_logistic_lpmf(y_bot, lam_ffv, c_ffv);
  EXPECT_TRUE(std::isfinite(out_bot.val_.val_.val()));
  out_bot.d_.d_.grad();
  EXPECT_TRUE(std::isfinite(lam_ffv.val_.val_.adj()))
      << "scalar y == 1: val-side adjoint of lambda is non-finite";
  EXPECT_TRUE(std::isfinite(c_ffv[0].val_.val_.adj()));
  EXPECT_TRUE(std::isfinite(c_ffv[1].val_.val_.adj()));
}

// Numerical-value regression: verify that the value AND first-order
// gradients are unchanged by the fix when y contains K. The reference
// values are computed from the algebraically-equivalent hand-rolled form
//   log p(y = K | lambda, c) = log_inv_logit(lambda - c[K-2])
// and its gradient
//   dlog p/dlambda = inv_logit(c[K-2] - lambda)
//   dlog p/dc[K-2] = -inv_logit(c[K-2] - lambda)
TEST_F(AgradRev, ProbDistributionsOrdLog_top_category_value_and_gradient) {
  using stan::math::ordered_logistic_lpmf;
  using stan::math::var;

  // K = 3, y = K. lambda = 0.7, c = (0.0, 1.0). c[K-2] = c[1] = 1.0.
  // sigmoid(c[K-2] - lambda) = sigmoid(0.3) = 0.574442516811659
  // log p = log_inv_logit(lambda - c[K-2]) = log_inv_logit(-0.3)
  //       = -log(1 + exp(0.3)) = -0.854355244468526
  int y = 3;
  var lam = 0.7;
  Eigen::Matrix<var, Eigen::Dynamic, 1> c(2);
  c << 0.0, 1.0;

  var out = ordered_logistic_lpmf(y, lam, c);
  EXPECT_FLOAT_EQ(out.val(), -0.854355244468526);

  out.grad();
  EXPECT_FLOAT_EQ(lam.adj(), 0.574442516811659);
  EXPECT_FLOAT_EQ(c[0].adj(), 0.0);
  EXPECT_FLOAT_EQ(c[1].adj(), -0.574442516811659);
}
