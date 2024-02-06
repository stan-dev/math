#include <stan/math/prim/functor.hpp>
#include <stan/math/prim/fun.hpp>
#include <stan/math/prim/prob.hpp>
#include <stan/math/rev.hpp>

#include <gtest/gtest.h>
#include <vector>

TEST(mathPrimScalProbWienerFullScal, valid) {
  using stan::math::INFTY;
  using stan::math::wiener_lpdf;
  double rt = 1, a = 1, v = -1, w = 0.5, t0 = 0.1, sv = 0.2, sw = 0.2,
         st0 = 0.1;
  EXPECT_NO_THROW(wiener_lpdf(rt, a, t0, w, v, sv, sw, st0));
  rt = 5;
  a = 1;
  v = 1;
  w = 0.5;
  t0 = 0.0;
  sv = 0.0;
  sw = 0.0;
  st0 = 0.0;
  EXPECT_NO_THROW(wiener_lpdf(rt, a, t0, w, v, sv, sw, st0));
}

// rt
TEST(mathPrimScalProbWienerFullScal, invalid_rt) {
  using stan::math::INFTY;
  using stan::math::wiener_lpdf;
  double a = 1, v = -1, w = 0.5, t0 = 0.1, sv = 0.2, sw = 0.2, st0 = 0.1;
  EXPECT_THROW(wiener_lpdf(0, a, t0, w, v, sv, sw, st0), std::domain_error);
  EXPECT_THROW(wiener_lpdf(-1, a, t0, w, v, sv, sw, st0), std::domain_error);
  EXPECT_THROW(wiener_lpdf(INFTY, a, t0, w, v, sv, sw, st0), std::domain_error);
  EXPECT_THROW(wiener_lpdf(-INFTY, a, t0, w, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lpdf(NAN, a, t0, w, v, sv, sw, st0), std::domain_error);
}

// a
TEST(mathPrimScalProbWienerFullScal, invalid_a) {
  using stan::math::INFTY;
  using stan::math::wiener_lpdf;
  double rt = 1, v = -1, w = 0.5, t0 = 0.1, sv = 0.2, sw = 0.2, st0 = 0.1;
  EXPECT_THROW(wiener_lpdf(rt, 0, t0, w, v, sv, sw, st0), std::domain_error);
  EXPECT_THROW(wiener_lpdf(rt, -1, t0, w, v, sv, sw, st0), std::domain_error);
  EXPECT_THROW(wiener_lpdf(rt, INFTY, t0, w, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lpdf(rt, -INFTY, t0, w, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lpdf(rt, NAN, t0, w, v, sv, sw, st0), std::domain_error);
}

// v
TEST(mathPrimScalProbWienerFullScal, invalid_v) {
  using stan::math::INFTY;
  using stan::math::wiener_lpdf;
  double rt = 1, a = 1, v = -1, w = 0.5, t0 = 0.1, sv = 0.2, sw = 0.2,
         st0 = 0.1;
  EXPECT_THROW(wiener_lpdf(rt, a, t0, w, INFTY, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lpdf(rt, a, t0, w, -INFTY, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lpdf(rt, a, t0, w, NAN, sv, sw, st0), std::domain_error);
}

// w
TEST(mathPrimScalProbWienerFullScal, invalid_w) {
  using stan::math::INFTY;
  using stan::math::wiener_lpdf;
  double rt = 1, a = 1, v = -1, t0 = 0.1, sv = 0.2, sw = 0.2, st0 = 0.1;
  EXPECT_THROW(wiener_lpdf(rt, a, t0, -0.1, v, sv, sw, st0), std::domain_error);
  EXPECT_THROW(wiener_lpdf(rt, a, t0, 0, v, sv, sw, st0), std::domain_error);
  EXPECT_THROW(wiener_lpdf(rt, a, t0, 1, v, sv, sw, st0), std::domain_error);
  EXPECT_THROW(wiener_lpdf(rt, a, t0, 1.1, v, sv, sw, st0), std::domain_error);
  EXPECT_THROW(wiener_lpdf(rt, a, t0, INFTY, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lpdf(rt, a, t0, -INFTY, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lpdf(rt, a, t0, NAN, v, sv, sw, st0), std::domain_error);
}

// t0
TEST(mathPrimScalProbWienerFullScal, invalid_t0) {
  using stan::math::INFTY;
  using stan::math::wiener_lpdf;
  double rt = 1, a = 1, v = -1, w = 0.5, sv = 0.2, sw = 0.2, st0 = 0.1;
  EXPECT_THROW(wiener_lpdf(rt, a, 2, w, v, sv, sw, st0),
               std::domain_error);  // rt must be greater than t0
  EXPECT_THROW(wiener_lpdf(rt, a, -1, w, v, sv, sw, st0), std::domain_error);
  EXPECT_THROW(wiener_lpdf(rt, a, INFTY, w, v, sv, sw, st0), std::domain_error);
  EXPECT_THROW(wiener_lpdf(rt, a, -INFTY, w, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lpdf(rt, a, NAN, w, v, sv, sw, st0), std::domain_error);
}

// sv
TEST(mathPrimScalProbWienerFullScal, invalid_sv) {
  using stan::math::INFTY;
  using stan::math::wiener_lpdf;
  double rt = 1, a = 1, v = -1, w = 0.5, t0 = 0.1, sw = 0.2, st0 = 0.1;
  EXPECT_THROW(wiener_lpdf(rt, a, t0, w, v, -1, sw, st0), std::domain_error);
  EXPECT_THROW(wiener_lpdf(rt, a, t0, w, v, INFTY, sw, st0), std::domain_error);
  EXPECT_THROW(wiener_lpdf(rt, a, t0, w, v, -INFTY, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lpdf(rt, a, t0, w, v, NAN, sw, st0), std::domain_error);
}

// sw
TEST(mathPrimScalProbWienerFullScal, invalid_sw) {
  using stan::math::INFTY;
  using stan::math::wiener_lpdf;
  double rt = 1, a = 1, v = -1, w = 0.5, t0 = 0.1, sv = 0.2, st0 = 0.1;
  EXPECT_THROW(wiener_lpdf(rt, a, t0, w, v, sv, -1, st0), std::domain_error);
  EXPECT_THROW(wiener_lpdf(rt, a, t0, 0.8, v, sv, 0.5, st0),
               std::domain_error);  // sw must be smaller than 2*(1-w)
  EXPECT_THROW(wiener_lpdf(rt, a, t0, 0.3, v, sv, 0.7, st0),
               std::domain_error);  // sw must be smaller than 2*w
  EXPECT_THROW(wiener_lpdf(rt, a, t0, w, v, sv, INFTY, st0), std::domain_error);
  EXPECT_THROW(wiener_lpdf(rt, a, t0, w, v, sv, -INFTY, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lpdf(rt, a, t0, w, v, sv, NAN, st0), std::domain_error);
}

// st0
TEST(mathPrimScalProbWienerFullScal, invalid_st0) {
  using stan::math::INFTY;
  using stan::math::wiener_lpdf;
  double rt = 1, a = 1, v = -1, w = 0.5, t0 = 0.1, sv = 0.2, sw = 0.2;
  EXPECT_THROW(wiener_lpdf(rt, a, t0, w, v, sv, sw, -1), std::domain_error);
  EXPECT_THROW(wiener_lpdf(rt, a, t0, w, v, sv, sw, INFTY), std::domain_error);
  EXPECT_THROW(wiener_lpdf(rt, a, t0, w, v, sv, sw, -INFTY), std::domain_error);
  EXPECT_THROW(wiener_lpdf(rt, a, t0, w, v, sv, sw, NAN), std::domain_error);
}

TEST(mathPrimScalProbWienerFullPrecScal, valid) {
  using stan::math::INFTY;
  using stan::math::wiener_lpdf;
  double rt = 1, a = 1, v = -1, w = 0.5, t0 = 0.1, sv = 0.2, sw = 0.2,
         st0 = 0.1;
  EXPECT_NO_THROW(wiener_lpdf(rt, a, t0, w, v, sv, sw, st0, 1e-4));
  rt = 5;
  a = 1;
  v = 1;
  w = 0.5;
  t0 = 0.0;
  sv = 0.0;
  sw = 0.0;
  st0 = 0.0;
  EXPECT_NO_THROW(wiener_lpdf(rt, a, t0, w, v, sv, sw, st0, 1e-4));
}

// rt
TEST(mathPrimScalProbWienerFullPrecScal, invalid_rt) {
  using stan::math::INFTY;
  using stan::math::wiener_lpdf;
  double a = 1, v = -1, w = 0.5, t0 = 0.1, sv = 0.2, sw = 0.2, st0 = 0.1;
  EXPECT_THROW(wiener_lpdf(0, a, t0, w, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lpdf(-1, a, t0, w, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lpdf(INFTY, a, t0, w, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lpdf(-INFTY, a, t0, w, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lpdf(NAN, a, t0, w, v, sv, sw, st0, 1e-4),
               std::domain_error);
}

// a
TEST(mathPrimScalProbWienerFullPrecScal, invalid_a) {
  using stan::math::INFTY;
  using stan::math::wiener_lpdf;
  double rt = 1, v = -1, w = 0.5, t0 = 0.1, sv = 0.2, sw = 0.2, st0 = 0.1;
  EXPECT_THROW(wiener_lpdf(rt, 0, t0, w, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lpdf(rt, -1, t0, w, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lpdf(rt, INFTY, t0, w, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lpdf(rt, -INFTY, t0, w, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lpdf(rt, NAN, t0, w, v, sv, sw, st0, 1e-4),
               std::domain_error);
}

// v
TEST(mathPrimScalProbWienerFullPrecScal, invalid_v) {
  using stan::math::INFTY;
  using stan::math::wiener_lpdf;
  double rt = 1, a = 1, v = -1, w = 0.5, t0 = 0.1, sv = 0.2, sw = 0.2,
         st0 = 0.1;
  EXPECT_THROW(wiener_lpdf(rt, a, t0, w, INFTY, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lpdf(rt, a, t0, w, -INFTY, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lpdf(rt, a, t0, w, NAN, sv, sw, st0, 1e-4),
               std::domain_error);
}

// w
TEST(mathPrimScalProbWienerFullPrecScal, invalid_w) {
  using stan::math::INFTY;
  using stan::math::wiener_lpdf;
  double rt = 1, a = 1, v = -1, t0 = 0.1, sv = 0.2, sw = 0.2, st0 = 0.1;
  EXPECT_THROW(wiener_lpdf(rt, a, t0, -0.1, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lpdf(rt, a, t0, 0, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lpdf(rt, a, t0, 1, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lpdf(rt, a, t0, 1.1, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lpdf(rt, a, t0, INFTY, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lpdf(rt, a, t0, -INFTY, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lpdf(rt, a, t0, NAN, v, sv, sw, st0, 1e-4),
               std::domain_error);
}

// t0
TEST(mathPrimScalProbWienerFullPrecScal, invalid_t0) {
  using stan::math::INFTY;
  using stan::math::wiener_lpdf;
  double rt = 1, a = 1, v = -1, w = 0.5, sv = 0.2, sw = 0.2, st0 = 0.1;
  EXPECT_THROW(wiener_lpdf(rt, a, 2, w, v, sv, sw, st0, 1e-4),
               std::domain_error);  // rt must be greater than t0
  EXPECT_THROW(wiener_lpdf(rt, a, -1, w, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lpdf(rt, a, INFTY, w, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lpdf(rt, a, -INFTY, w, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lpdf(rt, a, NAN, w, v, sv, sw, st0, 1e-4),
               std::domain_error);
}

// sv
TEST(mathPrimScalProbWienerFullPrecScal, invalid_sv) {
  using stan::math::INFTY;
  using stan::math::wiener_lpdf;
  double rt = 1, a = 1, v = -1, w = 0.5, t0 = 0.1, sw = 0.2, st0 = 0.1;
  EXPECT_THROW(wiener_lpdf(rt, a, t0, w, v, -1, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lpdf(rt, a, t0, w, v, INFTY, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lpdf(rt, a, t0, w, v, -INFTY, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lpdf(rt, a, t0, w, v, NAN, sw, st0, 1e-4),
               std::domain_error);
}

// sw
TEST(mathPrimScalProbWienerFullPrecScal, invalid_sw) {
  using stan::math::INFTY;
  using stan::math::wiener_lpdf;
  double rt = 1, a = 1, v = -1, w = 0.5, t0 = 0.1, sv = 0.2, st0 = 0.1;
  EXPECT_THROW(wiener_lpdf(rt, a, t0, w, v, sv, -1, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lpdf(rt, a, t0, 0.8, v, sv, 0.5, st0, 1e-4),
               std::domain_error);  // sw must be smaller than 2*(1-w)
  EXPECT_THROW(wiener_lpdf(rt, a, t0, 0.3, v, sv, 0.7, st0, 1e-4),
               std::domain_error);  // sw must be smaller than 2*w
  EXPECT_THROW(wiener_lpdf(rt, a, t0, w, v, sv, INFTY, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lpdf(rt, a, t0, w, v, sv, -INFTY, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lpdf(rt, a, t0, w, v, sv, NAN, st0, 1e-4),
               std::domain_error);
}

// st0
TEST(mathPrimScalProbWienerFullPrecScal, invalid_st0) {
  using stan::math::INFTY;
  using stan::math::wiener_lpdf;
  double rt = 1, a = 1, v = -1, w = 0.5, t0 = 0.1, sv = 0.2, sw = 0.2;
  EXPECT_THROW(wiener_lpdf(rt, a, t0, w, v, sv, sw, -1, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lpdf(rt, a, t0, w, v, sv, sw, INFTY, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lpdf(rt, a, t0, w, v, sv, sw, -INFTY, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lpdf(rt, a, t0, w, v, sv, sw, NAN, 1e-4),
               std::domain_error);
}

TEST(mathPrimCorrectValues, wiener_lpdf) {
  // Test concrete values. True values are computed in R using the R-package
  // WienR and the function WienerPDF.
  std::vector<double> y_vec = {2, 3, 4, 5, 6, 7, 8, 8.85, 8.9, 9, 1};
  std::vector<double> a_vec
      = {2.0, 2.0, 10.0, 4.0, 10.0, 1.0, 3.0, 1.7, 2.4, 11.0, 1.5};
  std::vector<double> v_vec
      = {2.0, 2.0, 4.0, 3.0, -3.0, 1.0, -1.0, -7.3, -4.9, 4.5, 3};
  std::vector<double> w_vec = {.1, 0.5, .8, 0.7, .1, .9, .7, .92, .9, .12, 0.5};
  std::vector<double> t0_vec
      = {1e-9, 0.01, .01, .01, .01, .01, .01, .01, .01, .01, 0.1};
  std::vector<double> sv_vec = {0, 0.2, 0, 0, .2, .2, 0, .7, 0, .7, 0.5};
  std::vector<double> sw_vec = {0, 0, .1, 0, .1, 0, .1, .01, 0, .1, 0.2};
  std::vector<double> st0_vec
      = {0, 0, 0, 0.007, 0, .007, .007, .009, .009, .009, 0};

  std::vector<double> true_dens
      = {-4.28564747866615, -7.52379235146909, -26.1551056209248,
         -22.1939134892089, -50.0587553794834, -37.2817263586318,
         -10.5428662079438, -61.5915905674246, -117.238967959795,
         -12.5788594249676, -3.1448097740735};
  std::vector<double> true_grad_y
      = {-3.22509339523307, -2.91155058614589, -8.21331631900955,
         -4.82948967379739, -1.50069056428102, -5.25831601347426,
         -1.04831896413742, -2.67457492096193, -12.8617364931501,
         -1.12047317491985, -5.68799957241344};
  std::vector<double> true_grad_a = {
      3.25018678924105,  3.59980430191399,  0.876602303160642, 1.2215517888504,
      -3.02928674030948, 67.0322498959921,  1.95334514374631,  16.4642201959135,
      5.02038145619773,  0.688439187670968, 2.63200041459657};
  std::vector<double> true_grad_t0
      = {3.22509339523307, 2.91155058614589, 8.21331631900955, 4.82948967379739,
         1.50069056428102, 5.25831601347426, 1.04831896413742, 2.67457492096193,
         12.8617364931501, 1.12047317491985, 5.68799957241344};
  std::vector<double> true_grad_w
      = {5.67120184517318,  -3.64396221090076, -38.7775057146792,
         -14.1837930137393, -34.5869239580708, -10.4535345681946,
         0.679597983582904, -9.93144540834201, 2.09117200953597,
         -6.0858540417876,  -3.74870310978083};
  std::vector<double> true_grad_v = {
      -2.199999998,     -4.44801714898178, -13.6940602985224, -13.7593709622169,
      21.5540563802381, -5.38233555673517, 8.88475440789056,  12.1280680728793,
      43.7785246930371, -5.68143495684294, -1.57639220567218};
  std::vector<double> true_grad_sv = {0,
                                      3.42285198319565,
                                      0,
                                      0,
                                      91.9551438876654,
                                      4.70180879974639,
                                      0,
                                      101.80250964211,
                                      0,
                                      21.4332628706595,
                                      0.877556017134384};
  std::vector<double> true_grad_sw = {0,
                                      0,
                                      10.1052188867058,
                                      0,
                                      8.72398,
                                      0,
                                      -0.122807217815892,
                                      -0.0506322723373748,
                                      0,
                                      -0.0704990526706635,
                                      0.0827817310725268};
  std::vector<double> true_grad_st0 = {0,
                                       0,
                                       0,
                                       2.42836139121338,
                                       0,
                                       2.64529825657625,
                                       0.524800556172613,
                                       1.34278261179603,
                                       6.55490874737353,
                                       0.561295838843035,
                                       0};

  using stan::math::var;
  double err_tol_dens = 1e-6;
  double err_tol = 1e-4;
  for (int i = 0; i < y_vec.size(); i++) {
    var y = y_vec[i];
    var a = a_vec[i];
    var t0 = t0_vec[i];
    var w = w_vec[i];
    var v = v_vec[i];
    var sv = sv_vec[i];
    var sw = sw_vec[i];
    var st0 = st0_vec[i];
    var dens = stan::math::wiener_lpdf(y, a, t0, w, v, sv, sw, st0);
    dens.grad();
    EXPECT_NEAR(dens.val(), true_dens[i], err_tol_dens);
    EXPECT_NEAR(y.adj(), true_grad_y[i], err_tol);
    EXPECT_NEAR(a.adj(), true_grad_a[i], err_tol);
    EXPECT_NEAR(t0.adj(), true_grad_t0[i], err_tol);
    EXPECT_NEAR(w.adj(), true_grad_w[i], err_tol);
    EXPECT_NEAR(v.adj(), true_grad_v[i], err_tol);
    EXPECT_NEAR(sv.adj(), true_grad_sv[i], err_tol);
    EXPECT_NEAR(sw.adj(), true_grad_sw[i], err_tol);
    EXPECT_NEAR(st0.adj(), true_grad_st0[i], err_tol);
  }
}

TEST(mathPrimCorrectValuesFourParamModel, wiener_lpdf) {
  // Test concrete values. True values are computed in R using the R-package
  // WienR and the function WienerPDF.
  std::vector<double> y_vec = {2.0};
  std::vector<double> a_vec = {2.0};
  std::vector<double> v_vec = {2.0};
  std::vector<double> w_vec = {.1};
  std::vector<double> t0_vec = {1e-9};

  std::vector<double> true_dens = {-4.28564747866615};
  std::vector<double> true_grad_y = {-3.22509339523307};
  std::vector<double> true_grad_a = {3.25018678924105};
  std::vector<double> true_grad_t0 = {3.22509339523307};
  std::vector<double> true_grad_w = {5.67120184517318};
  std::vector<double> true_grad_v = {-2.199999998};

  using stan::math::var;
  double err_tol_dens = 1e-6;
  double err_tol = 1e-4;
  for (int i = 0; i < y_vec.size(); i++) {
    var y = y_vec[i];
    var a = a_vec[i];
    var t0 = t0_vec[i];
    var w = w_vec[i];
    var v = v_vec[i];
    var dens = stan::math::wiener_lpdf(y, a, t0, w, v);
    dens.grad();
    EXPECT_NEAR(dens.val(), true_dens[i], err_tol_dens);
    EXPECT_NEAR(y.adj(), true_grad_y[i], err_tol);
    EXPECT_NEAR(a.adj(), true_grad_a[i], err_tol);
    EXPECT_NEAR(t0.adj(), true_grad_t0[i], err_tol);
    EXPECT_NEAR(w.adj(), true_grad_w[i], err_tol);
    EXPECT_NEAR(v.adj(), true_grad_v[i], err_tol);
  }
}

TEST(mathPrimCorrectValuesFiveParameterModel, wiener_lpdf) {
  // Test concrete values. True values are computed in R using the R-package
  // WienR and the function WienerPDF.
  std::vector<double> y_vec = {3.0};
  std::vector<double> a_vec = {2.0};
  std::vector<double> v_vec = {2.0};
  std::vector<double> w_vec = {0.5};
  std::vector<double> t0_vec = {0.01};
  std::vector<double> sv_vec = {0.2};

  std::vector<double> true_dens = {-7.52379235146909};
  std::vector<double> true_grad_y = {-2.91155058614589};
  std::vector<double> true_grad_a = {3.59980430191399};
  std::vector<double> true_grad_t0 = {2.91155058614589};
  std::vector<double> true_grad_w = {-3.64396221090076};
  std::vector<double> true_grad_v = {-4.44801714898178};
  std::vector<double> true_grad_sv = {3.42285198319565};

  using stan::math::var;
  double err_tol_dens = 1e-6;
  double err_tol = 1e-4;
  for (int i = 0; i < y_vec.size(); i++) {
    var y = y_vec[i];
    var a = a_vec[i];
    var t0 = t0_vec[i];
    var w = w_vec[i];
    var v = v_vec[i];
    var sv = sv_vec[i];
    var dens = stan::math::wiener_lpdf(y, a, t0, w, v, sv);
    dens.grad();
    EXPECT_NEAR(dens.val(), true_dens[i], err_tol_dens);
    EXPECT_NEAR(y.adj(), true_grad_y[i], err_tol);
    EXPECT_NEAR(a.adj(), true_grad_a[i], err_tol);
    EXPECT_NEAR(t0.adj(), true_grad_t0[i], err_tol);
    EXPECT_NEAR(w.adj(), true_grad_w[i], err_tol);
    EXPECT_NEAR(v.adj(), true_grad_v[i], err_tol);
    EXPECT_NEAR(sv.adj(), true_grad_sv[i], err_tol);
  }
}
