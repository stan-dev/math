#include <stan/math/mix.hpp>
#include <test/unit/math/rev/fun/util.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(ProbDistributionsOrdLog, fv_fv) {
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

TEST(ProbDistributionsOrdLog, fv_d) {
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

TEST(ProbDistributionsOrdLog, fv_fv_vec) {
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

TEST(ProbDistributionsOrdLog, fv_d_vec) {
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

TEST(ProbDistributionsOrdLog, fv_fv_stvec) {
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

TEST(ProbDistributionsOrdLog, fv_d_stvec) {
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
