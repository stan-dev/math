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
  lam_v << -2.95, -1.68, 0.96, 2.68;

  vector_v c_v(3);
  c_v << -2.68, -1.53, 0.46;

  AVAR out_v = ordered_logistic_lpmf(y, lam_v, c_v);
  out_v.grad();

  EXPECT_FLOAT_EQ(out_v.val(), -3.18600232948363);

  EXPECT_FLOAT_EQ(lam_v[0].adj(), -0.432907095034546);
  EXPECT_FLOAT_EQ(lam_v[1].adj(), -0.193628733286255);
  EXPECT_FLOAT_EQ(lam_v[2].adj(), -0.545897133850042);
  EXPECT_FLOAT_EQ(lam_v[3].adj(), 0.097968804297554);

  EXPECT_FLOAT_EQ(c_v[0].adj(), -0.299384935778303);
  EXPECT_FLOAT_EQ(c_v[1].adj(), 0.6910188226123);
  EXPECT_FLOAT_EQ(c_v[2].adj(), 0.682830271039293);
}

TEST(ProbDistributionsOrdLog, vd_vec) {
  using stan::math::ordered_logistic_lpmf;
  using stan::math::var;
  using stan::math::vector_d;
  using stan::math::vector_v;

  std::vector<int> y{1, 2, 3, 4};

  vector_v lam_v(4);
  lam_v << -2.95, -1.68, 0.96, 2.68;

  vector_d lam_d(4);
  lam_d << -2.95, -1.68, 0.96, 2.68;

  vector_v c_v(3);
  c_v << -3.62, -2.89, -1.57;

  vector_d c_d(3);
  c_d << -3.62, -2.89, -1.57;

  AVAR out_v = ordered_logistic_lpmf(y, lam_v, c_d);
  out_v.grad();

  EXPECT_FLOAT_EQ(out_v.val(), -6.29875293667739);

  EXPECT_FLOAT_EQ(lam_v[0].adj(), -0.661503159202952);
  EXPECT_FLOAT_EQ(lam_v[1].adj(), -0.644651092531256);
  EXPECT_FLOAT_EQ(lam_v[2].adj(), -0.905382008878854);
  EXPECT_FLOAT_EQ(lam_v[3].adj(), 0.014063627043246);

  out_v = ordered_logistic_lpmf(y, lam_d, c_v);
  out_v.grad();

  EXPECT_FLOAT_EQ(out_v.val(), -6.29875293667739);

  EXPECT_FLOAT_EQ(c_v[0].adj(), -0.394307508232632);
  EXPECT_FLOAT_EQ(c_v[1].adj(), 1.3151170652951);
  EXPECT_FLOAT_EQ(c_v[2].adj(), 1.27666307650735);
}

TEST(ProbDistributionsOrdLog, vv_stvec) {
  using stan::math::ordered_logistic_lpmf;
  using stan::math::var;
  using stan::math::vector_d;
  using stan::math::vector_v;

  std::vector<int> y{1, 2, 3, 4};

  vector_v lam_v(4);
  lam_v << -2.95, -1.68, 0.96, 2.68;

  vector_v c1_v(3);
  c1_v << -2.68, -1.53, 0.46;

  vector_v c2_v(3);
  c2_v << -3.62, -2.89, -1.57;

  vector_v c3_v(3);
  c3_v << -1.21, -0.83, 1.87;

  vector_v c4_v(3);
  c4_v << -0.14, 1.53, 3.87;

  std::vector<vector_v> std_c_v{c1_v, c2_v, c3_v, c4_v};

  AVAR out_v = ordered_logistic_lpmf(y, lam_v, std_c_v);
  out_v.grad();

  EXPECT_FLOAT_EQ(out_v.val(), -4.84793751666795);

  EXPECT_FLOAT_EQ(lam_v[0].adj(), -0.432907095034546);
  EXPECT_FLOAT_EQ(lam_v[1].adj(), -0.644651092531256);
  EXPECT_FLOAT_EQ(lam_v[2].adj(), -0.143927113767638);
  EXPECT_FLOAT_EQ(lam_v[3].adj(), 0.766741064228543);

  EXPECT_FLOAT_EQ(std_c_v[0][0].adj(), 0.432907095034546);
  EXPECT_FLOAT_EQ(std_c_v[0][1].adj(), 0.0);
  EXPECT_FLOAT_EQ(std_c_v[0][2].adj(), 0.0);

  EXPECT_FLOAT_EQ(std_c_v[1][0].adj(), -1.05581066743558);
  EXPECT_FLOAT_EQ(std_c_v[1][1].adj(), 1.70046175996684);
  EXPECT_FLOAT_EQ(std_c_v[1][2].adj(), 0.0);

  EXPECT_FLOAT_EQ(std_c_v[2][0].adj(), 0.0);
  EXPECT_FLOAT_EQ(std_c_v[2][1].adj(), -0.215120225537197);
  EXPECT_FLOAT_EQ(std_c_v[2][2].adj(), 0.359047339304835);

  EXPECT_FLOAT_EQ(std_c_v[3][0].adj(), 0.0);
  EXPECT_FLOAT_EQ(std_c_v[3][1].adj(), 0.0);
  EXPECT_FLOAT_EQ(std_c_v[3][2].adj(), -0.766741064228543);
}

TEST(ProbDistributionsOrdLog, vd_stvec) {
  using stan::math::ordered_logistic_lpmf;
  using stan::math::var;
  using stan::math::vector_d;
  using stan::math::vector_v;

  std::vector<int> y{1, 2, 3, 4};

  vector_v lam_v(4);
  lam_v << -1.52, -3.51, -0.56, 1.55;

  vector_v c1_v(3);
  c1_v << -3.18, -1.84, 0.58;

  vector_v c2_v(3);
  c2_v << -5.14, -2.81, -0.93;

  vector_v c3_v(3);
  c3_v << -3.24, -1.62, 2.61;

  vector_v c4_v(3);
  c4_v << -2.17, -0.24, 2.89;

  vector_d lam_d(4);
  lam_d << -1.52, -3.51, -0.56, 1.55;

  vector_d c1_d(3);
  c1_d << -3.18, -1.84, 0.58;

  vector_d c2_d(3);
  c2_d << -5.14, -2.81, -0.93;

  vector_d c3_d(3);
  c3_d << -3.24, -1.62, 2.61;

  vector_d c4_d(3);
  c4_d << -2.17, -0.24, 2.89;

  std::vector<vector_v> std_c_v{c1_v, c2_v, c3_v, c4_v};
  std::vector<vector_d> std_c_d{c1_d, c2_d, c3_d, c4_d};

  AVAR out_v = ordered_logistic_lpmf(y, lam_v, std_c_d);
  out_v.grad();

  EXPECT_FLOAT_EQ(out_v.val(), -4.44439619529986);

  EXPECT_FLOAT_EQ(lam_v[0].adj(), -0.840238003056331);
  EXPECT_FLOAT_EQ(lam_v[1].adj(), -0.167981866608475);
  EXPECT_FLOAT_EQ(lam_v[2].adj(), 0.216999039274605);
  EXPECT_FLOAT_EQ(lam_v[3].adj(), 0.792489941440364);

  out_v = ordered_logistic_lpmf(y, lam_d, std_c_v);
  out_v.grad();

  EXPECT_FLOAT_EQ(out_v.val(), -4.44439619529986);

  EXPECT_FLOAT_EQ(std_c_v[0][0].adj(), 0.840238003056331);
  EXPECT_FLOAT_EQ(std_c_v[0][1].adj(), 0.0);
  EXPECT_FLOAT_EQ(std_c_v[0][2].adj(), 0.0);

  EXPECT_FLOAT_EQ(std_c_v[1][0].adj(), -0.271612889970549);
  EXPECT_FLOAT_EQ(std_c_v[1][1].adj(), 0.439594756579024);
  EXPECT_FLOAT_EQ(std_c_v[1][2].adj(), 0.0);

  EXPECT_FLOAT_EQ(std_c_v[2][0].adj(), 0.0);
  EXPECT_FLOAT_EQ(std_c_v[2][1].adj(), -0.272076744615968);
  EXPECT_FLOAT_EQ(std_c_v[2][2].adj(), 0.055077705341364);

  EXPECT_FLOAT_EQ(std_c_v[3][0].adj(), 0.0);
  EXPECT_FLOAT_EQ(std_c_v[3][1].adj(), 0.0);
  EXPECT_FLOAT_EQ(std_c_v[3][2].adj(), -0.792489941440364);
}

TEST(ProbDistributionsOrdLog, intErrors) {
  using stan::math::ordered_logistic_lpmf;
  using stan::math::vector_v;

  std::vector<int> y{1, 2, 3, 4};

  std::vector<int> lam_i{-1, -3, -1, 2};

  vector_v c1_v(3);
  c1_v << -3.18, -1.84, 0.58;

  vector_v c2_v(3);
  c2_v << -5.14, -2.81, -0.93;

  vector_v c3_v(3);
  c3_v << -3.24, -1.62, 2.61;

  vector_v c4_v(3);
  c4_v << -2.17, -0.24, 2.89;

  std::vector<vector_v> std_c_v{c1_v, c2_v, c3_v, c4_v};

  AVAR out_v = ordered_logistic_lpmf(y, lam_i, std_c_v);
  out_v.grad();

  EXPECT_FLOAT_EQ(out_v.val(), -4.80919533214772);

  EXPECT_FLOAT_EQ(std_c_v[0][0].adj(), 0.898439072102363);
  EXPECT_FLOAT_EQ(std_c_v[0][1].adj(), 0.0);
  EXPECT_FLOAT_EQ(std_c_v[0][2].adj(), 0.0);

  EXPECT_FLOAT_EQ(std_c_v[1][0].adj(), -0.213051918269699);
  EXPECT_FLOAT_EQ(std_c_v[1][1].adj(), 0.5604249106041);
  EXPECT_FLOAT_EQ(std_c_v[1][2].adj(), 0.0);

  EXPECT_FLOAT_EQ(std_c_v[2][0].adj(), 0.0);
  EXPECT_FLOAT_EQ(std_c_v[2][1].adj(), -0.364548741344827);
  EXPECT_FLOAT_EQ(std_c_v[2][2].adj(), 0.041106609543937);

  EXPECT_FLOAT_EQ(std_c_v[3][0].adj(), 0.0);
  EXPECT_FLOAT_EQ(std_c_v[3][1].adj(), 0.0);
  EXPECT_FLOAT_EQ(std_c_v[3][2].adj(), -0.70889017256612);
}
