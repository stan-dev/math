#include <stan/math/prim/functor.hpp>
#include <stan/math/prim/fun.hpp>
#include <stan/math/prim/prob.hpp>
#include <stan/math/rev.hpp>

#include <gtest/gtest.h>
#include <vector>

TEST(mathPrimScalProbWienerFullLcdfScal, valid) {
  using stan::math::INFTY;
  using stan::math::wiener_lcdf_unnorm;
  double rt = 1, a = 1, v = -1, w = 0.5, t0 = 0.1, sv = 0.2, sw = 0.2,
         st0 = 0.1;
  EXPECT_NO_THROW(wiener_lcdf_unnorm(rt, a, t0, w, v, sv, sw, st0));
  rt = 5;
  a = 1;
  v = 1;
  w = 0.5;
  t0 = 0.0;
  sv = 0.0;
  sw = 0.0;
  st0 = 0.0;
  EXPECT_NO_THROW(wiener_lcdf_unnorm(rt, a, t0, w, v, sv, sw, st0));
}

// rt
TEST(mathPrimScalProbWienerFullLcdfScal, invalid_rt) {
  using stan::math::INFTY;
  using stan::math::wiener_lcdf_unnorm;
  double a = 1, v = -1, w = 0.5, t0 = 0.1, sv = 0.2, sw = 0.2, st0 = 0.1;
  EXPECT_THROW(wiener_lcdf_unnorm(0, a, t0, w, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lcdf_unnorm(-1, a, t0, w, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lcdf_unnorm(INFTY, a, t0, w, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lcdf_unnorm(-INFTY, a, t0, w, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lcdf_unnorm(NAN, a, t0, w, v, sv, sw, st0),
               std::domain_error);
}

// a
TEST(mathPrimScalProbWienerFullLcdfScal, invalid_a) {
  using stan::math::INFTY;
  using stan::math::wiener_lcdf_unnorm;
  double rt = 1, v = -1, w = 0.5, t0 = 0.1, sv = 0.2, sw = 0.2, st0 = 0.1;
  EXPECT_THROW(wiener_lcdf_unnorm(rt, 0, t0, w, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lcdf_unnorm(rt, -1, t0, w, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lcdf_unnorm(rt, INFTY, t0, w, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lcdf_unnorm(rt, -INFTY, t0, w, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lcdf_unnorm(rt, NAN, t0, w, v, sv, sw, st0),
               std::domain_error);
}

// v
TEST(mathPrimScalProbWienerFullLcdfScal, invalid_v) {
  using stan::math::INFTY;
  using stan::math::wiener_lcdf_unnorm;
  double rt = 1, a = 1, v = -1, w = 0.5, t0 = 0.1, sv = 0.2, sw = 0.2,
         st0 = 0.1;
  EXPECT_THROW(wiener_lcdf_unnorm(rt, a, t0, w, INFTY, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lcdf_unnorm(rt, a, t0, w, -INFTY, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lcdf_unnorm(rt, a, t0, w, NAN, sv, sw, st0),
               std::domain_error);
}

// w
TEST(mathPrimScalProbWienerFullLcdfScal, invalid_w) {
  using stan::math::INFTY;
  using stan::math::wiener_lcdf_unnorm;
  double rt = 1, a = 1, v = -1, t0 = 0.1, sv = 0.2, sw = 0.2, st0 = 0.1;
  EXPECT_THROW(wiener_lcdf_unnorm(rt, a, t0, -0.1, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lcdf_unnorm(rt, a, t0, 0, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lcdf_unnorm(rt, a, t0, 1, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lcdf_unnorm(rt, a, t0, 1.1, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lcdf_unnorm(rt, a, t0, INFTY, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lcdf_unnorm(rt, a, t0, -INFTY, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lcdf_unnorm(rt, a, t0, NAN, v, sv, sw, st0),
               std::domain_error);
}

// t0
TEST(mathPrimScalProbWienerFullLcdfScal, invalid_t0) {
  using stan::math::INFTY;
  using stan::math::wiener_lcdf_unnorm;
  double rt = 1, a = 1, v = -1, w = 0.5, sv = 0.2, sw = 0.2, st0 = 0.1;
  EXPECT_THROW(wiener_lcdf_unnorm(rt, a, 2, w, v, sv, sw, st0),
               std::domain_error);  // rt must be greater than t0
  EXPECT_THROW(wiener_lcdf_unnorm(rt, a, -1, w, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lcdf_unnorm(rt, a, INFTY, w, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lcdf_unnorm(rt, a, -INFTY, w, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lcdf_unnorm(rt, a, NAN, w, v, sv, sw, st0),
               std::domain_error);
}

// sv
TEST(mathPrimScalProbWienerFullLcdfScal, invalid_sv) {
  using stan::math::INFTY;
  using stan::math::wiener_lcdf_unnorm;
  double rt = 1, a = 1, v = -1, w = 0.5, t0 = 0.1, sw = 0.2, st0 = 0.1;
  EXPECT_THROW(wiener_lcdf_unnorm(rt, a, t0, w, v, -1, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lcdf_unnorm(rt, a, t0, w, v, INFTY, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lcdf_unnorm(rt, a, t0, w, v, -INFTY, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lcdf_unnorm(rt, a, t0, w, v, NAN, sw, st0),
               std::domain_error);
}

// sw
TEST(mathPrimScalProbWienerFullLcdfScal, invalid_sw) {
  using stan::math::INFTY;
  using stan::math::wiener_lcdf_unnorm;
  double rt = 1, a = 1, v = -1, w = 0.5, t0 = 0.1, sv = 0.2, st0 = 0.1;
  EXPECT_THROW(wiener_lcdf_unnorm(rt, a, t0, w, v, sv, -1, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lcdf_unnorm(rt, a, t0, 0.8, v, sv, 0.5, st0),
               std::domain_error);  // sw must be smaller than 2*(1-w)
  EXPECT_THROW(wiener_lcdf_unnorm(rt, a, t0, 0.3, v, sv, 0.7, st0),
               std::domain_error);  // sw must be smaller than 2*w
  EXPECT_THROW(wiener_lcdf_unnorm(rt, a, t0, w, v, sv, INFTY, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lcdf_unnorm(rt, a, t0, w, v, sv, -INFTY, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lcdf_unnorm(rt, a, t0, w, v, sv, NAN, st0),
               std::domain_error);
}

// st0
TEST(mathPrimScalProbWienerFullLcdfScal, invalid_st0) {
  using stan::math::INFTY;
  using stan::math::wiener_lcdf_unnorm;
  double rt = 1, a = 1, v = -1, w = 0.5, t0 = 0.1, sv = 0.2, sw = 0.2;
  EXPECT_THROW(wiener_lcdf_unnorm(rt, a, t0, w, v, sv, sw, -1),
               std::domain_error);
  EXPECT_THROW(wiener_lcdf_unnorm(rt, a, t0, w, v, sv, sw, INFTY),
               std::domain_error);
  EXPECT_THROW(wiener_lcdf_unnorm(rt, a, t0, w, v, sv, sw, -INFTY),
               std::domain_error);
  EXPECT_THROW(wiener_lcdf_unnorm(rt, a, t0, w, v, sv, sw, NAN),
               std::domain_error);
}

TEST(mathPrimScalProbWienerFullLcdfPrecScal, valid) {
  using stan::math::INFTY;
  using stan::math::wiener_lcdf_unnorm;
  double rt = 1, a = 1, v = -1, w = 0.5, t0 = 0.1, sv = 0.2, sw = 0.2,
         st0 = 0.1;
  EXPECT_NO_THROW(wiener_lcdf_unnorm(rt, a, t0, w, v, sv, sw, st0, 1e-4));
  rt = 5;
  a = 1;
  v = 1;
  w = 0.5;
  t0 = 0.0;
  sv = 0.0;
  sw = 0.0;
  st0 = 0.0;
  EXPECT_NO_THROW(wiener_lcdf_unnorm(rt, a, t0, w, v, sv, sw, st0, 1e-4));
}

// rt
TEST(mathPrimScalProbWienerFullLcdfPrecScal, invalid_rt) {
  using stan::math::INFTY;
  using stan::math::wiener_lcdf_unnorm;
  double a = 1, v = -1, w = 0.5, t0 = 0.1, sv = 0.2, sw = 0.2, st0 = 0.1;
  EXPECT_THROW(wiener_lcdf_unnorm(0, a, t0, w, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lcdf_unnorm(-1, a, t0, w, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lcdf_unnorm(INFTY, a, t0, w, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lcdf_unnorm(-INFTY, a, t0, w, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lcdf_unnorm(NAN, a, t0, w, v, sv, sw, st0, 1e-4),
               std::domain_error);
}

// a
TEST(mathPrimScalProbWienerFullLcdfPrecScal, invalid_a) {
  using stan::math::INFTY;
  using stan::math::wiener_lcdf_unnorm;
  double rt = 1, v = -1, w = 0.5, t0 = 0.1, sv = 0.2, sw = 0.2, st0 = 0.1;
  EXPECT_THROW(wiener_lcdf_unnorm(rt, 0, t0, w, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lcdf_unnorm(rt, -1, t0, w, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lcdf_unnorm(rt, INFTY, t0, w, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lcdf_unnorm(rt, -INFTY, t0, w, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lcdf_unnorm(rt, NAN, t0, w, v, sv, sw, st0, 1e-4),
               std::domain_error);
}

// v
TEST(mathPrimScalProbWienerFullLcdfPrecScal, invalid_v) {
  using stan::math::INFTY;
  using stan::math::wiener_lcdf_unnorm;
  double rt = 1, a = 1, v = -1, w = 0.5, t0 = 0.1, sv = 0.2, sw = 0.2,
         st0 = 0.1;
  EXPECT_THROW(wiener_lcdf_unnorm(rt, a, t0, w, INFTY, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lcdf_unnorm(rt, a, t0, w, -INFTY, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lcdf_unnorm(rt, a, t0, w, NAN, sv, sw, st0, 1e-4),
               std::domain_error);
}

// w
TEST(mathPrimScalProbWienerFullLcdfPrecScal, invalid_w) {
  using stan::math::INFTY;
  using stan::math::wiener_lcdf_unnorm;
  double rt = 1, a = 1, v = -1, t0 = 0.1, sv = 0.2, sw = 0.2, st0 = 0.1;
  EXPECT_THROW(wiener_lcdf_unnorm(rt, a, t0, -0.1, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lcdf_unnorm(rt, a, t0, 0, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lcdf_unnorm(rt, a, t0, 1, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lcdf_unnorm(rt, a, t0, 1.1, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lcdf_unnorm(rt, a, t0, INFTY, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lcdf_unnorm(rt, a, t0, -INFTY, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lcdf_unnorm(rt, a, t0, NAN, v, sv, sw, st0, 1e-4),
               std::domain_error);
}

// t0
TEST(mathPrimScalProbWienerFullLcdfPrecScal, invalid_t0) {
  using stan::math::INFTY;
  using stan::math::wiener_lcdf_unnorm;
  double rt = 1, a = 1, v = -1, w = 0.5, sv = 0.2, sw = 0.2, st0 = 0.1;
  EXPECT_THROW(wiener_lcdf_unnorm(rt, a, 2, w, v, sv, sw, st0, 1e-4),
               std::domain_error);  // rt must be greater than t0
  EXPECT_THROW(wiener_lcdf_unnorm(rt, a, -1, w, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lcdf_unnorm(rt, a, INFTY, w, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lcdf_unnorm(rt, a, -INFTY, w, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lcdf_unnorm(rt, a, NAN, w, v, sv, sw, st0, 1e-4),
               std::domain_error);
}

// sv
TEST(mathPrimScalProbWienerFullLcdfPrecScal, invalid_sv) {
  using stan::math::INFTY;
  using stan::math::wiener_lcdf_unnorm;
  double rt = 1, a = 1, v = -1, w = 0.5, t0 = 0.1, sw = 0.2, st0 = 0.1;
  EXPECT_THROW(wiener_lcdf_unnorm(rt, a, t0, w, v, -1, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lcdf_unnorm(rt, a, t0, w, v, INFTY, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lcdf_unnorm(rt, a, t0, w, v, -INFTY, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lcdf_unnorm(rt, a, t0, w, v, NAN, sw, st0, 1e-4),
               std::domain_error);
}

// sw
TEST(mathPrimScalProbWienerFullLcdfPrecScal, invalid_sw) {
  using stan::math::INFTY;
  using stan::math::wiener_lcdf_unnorm;
  double rt = 1, a = 1, v = -1, w = 0.5, t0 = 0.1, sv = 0.2, st0 = 0.1;
  EXPECT_THROW(wiener_lcdf_unnorm(rt, a, t0, w, v, sv, -1, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lcdf_unnorm(rt, a, t0, 0.8, v, sv, 0.5, st0, 1e-4),
               std::domain_error);  // sw must be smaller than 2*(1-w)
  EXPECT_THROW(wiener_lcdf_unnorm(rt, a, t0, 0.3, v, sv, 0.7, st0, 1e-4),
               std::domain_error);  // sw must be smaller than 2*w
  EXPECT_THROW(wiener_lcdf_unnorm(rt, a, t0, w, v, sv, INFTY, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lcdf_unnorm(rt, a, t0, w, v, sv, -INFTY, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lcdf_unnorm(rt, a, t0, w, v, sv, NAN, st0, 1e-4),
               std::domain_error);
}

// st0
TEST(mathPrimScalProbWienerFullLcdfPrecScal, invalid_st0) {
  using stan::math::INFTY;
  using stan::math::wiener_lcdf_unnorm;
  double rt = 1, a = 1, v = -1, w = 0.5, t0 = 0.1, sv = 0.2, sw = 0.2;
  EXPECT_THROW(wiener_lcdf_unnorm(rt, a, t0, w, v, sv, sw, -1, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lcdf_unnorm(rt, a, t0, w, v, sv, sw, INFTY, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lcdf_unnorm(rt, a, t0, w, v, sv, sw, -INFTY, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lcdf_unnorm(rt, a, t0, w, v, sv, sw, NAN, 1e-4),
               std::domain_error);
}

TEST(mathPrimCorrectValues, wiener_lcdf_unnorm) {
  /* Test concrete values. True values are computed in R using the R-package
   * WienR and the function WienerCDF() with its partial derivatives
   */
  std::vector<double> y_vec
      = {.2, .3, 2.2, .1, 2, 1.7, 3, 2.5, 3.7, 1.8, .9, 1.5, 2};
  std::vector<double> a_vec
      = {3, 1.7, 2.4, 2, 1.2, 1.5, 4, 1.5, 1.7, 6, 1, 2.5, 2};
  std::vector<double> t0_vec
      = {0, 0, 0, 0, 0.2, 0.3, 0.1, 0.7, 1, 0.6, .1, .2, .2};
  std::vector<double> w_vec
      = {.7, .92, .9, .4, 0.5, 0.3, 0.45, 0.6, 0.8, 0.35, .4, .5, .6};
  std::vector<double> v_vec
      = {-1, -7.3, -4.9, 2.5, 1, 0.5, -0.7, 0.9, -1.4, 2, .8, -1, -.9};
  std::vector<double> sv_vec = {0, 0, 0, 0, 0.5, 0, 0, 0.7, 0.8, 0, .4, .5, .2};
  std::vector<double> sw_vec = {0, 0, 0, 0, 0, 0.1, 0, 0.2, 0, 0.3, .2, 0, .2};
  std::vector<double> st0_vec
      = {0, 0, 0, 0, 0, 0, 0.1, 0, 0.1, 0.1, .2, .3, .25};

  std::vector<double> true_lcdf
      = {-4.09435375436108, -1.98560937871925,  -2.35200000058001,
         -6.09860875887167, -0.28490338793642,  -0.835235451149805,
         -3.60573984770431, -0.213132053911535, -0.887413009805869,
         -2.14898785081967, -0.574404607969948, -2.67530409425456,
         -1.55387294861195};
  std::vector<double> true_grad_y
      = {11.6955868009794,  0.000291754616022926, 2.81009792107004e-13,
         73.9485382296513,  0.00423167900126897,  0.170505395422671,
         0.301191971384078, 0.0216747718871863,   0.00332619998061515,
         3.66549265245825,  0.211550389634303,    0.732731869409894,
         0.117239643328064};
  std::vector<double> true_grad_a
      = {-1.90909592842206,  -1.1680106846017,   -0.979999994944343,
         -6.17343628534049,  0.189129525101159,  -0.0640548615428843,
         -1.01571721473935,  0.0439926894045548, -0.405445995424302,
         -0.875235676490553, 0.0169546689079121, -1.10780211022062,
         -0.724594716693696};
  std::vector<double> true_grad_t0
      = {-11.6955868009794,  -0.000291754616022926, -2.81009792107004e-13,
         -73.9485382296513,  -0.00423167900126897,  -0.170505395422671,
         -0.301191971384078, -0.0216747718871863,   -0.00332619998061515,
         -3.66549265245825,  -0.211550389634303,    -0.732731869409894,
         -0.117239643328064};
  std::vector<double> true_grad_w = {
      19.0909592842206, 24.8202270495257, 23.5200000150767, 20.578120951136,
      1.03356005252668, 2.72745006368943, 7.46860413650941, 0.749294067331125,
      3.99584607466979, 8.74288037387998, 1.83454127550223, 5.60378493291068,
      4.16252427817267};
  std::vector<double> true_grad_v
      = {1.04905306487445,  0.27197850836347,  0.479999997523509,
         0.977134030099708, 0.27971915432389,  0.762665629738307,
         3.35189640059343,  0.25195637817767,  0.558439280243674,
         1.59504189896887,  0.392934480883027, 1.62283843578951,
         1.22498956343152};
  std::vector<double> true_grad_sv = {0,
                                      0,
                                      0,
                                      0,
                                      -0.075059359595863,
                                      0,
                                      0,
                                      -0.118203610656295,
                                      0.137718761146619,
                                      0,
                                      -0.0215064672196581,
                                      1.05496167316742,
                                      0.220567350697758};
  std::vector<double> true_grad_sw = {0,
                                      0,
                                      0,
                                      0,
                                      0,
                                      -0.0276727265109265,
                                      0,
                                      -0.0201028131785313,
                                      0,
                                      1.43818889437335,
                                      -0.0415094057459971,
                                      0,
                                      0.25713151312856};
  std::vector<double> true_grad_st0 = {0,
                                       0,
                                       0,
                                       0,
                                       0,
                                       0,
                                       -0.151860277428223,
                                       0,
                                       -0.00171635548928614,
                                       -1.7784225450937,
                                       -0.12422577484349,
                                       -0.384925306866848,
                                       -0.0625008391233113};

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
    var lcdf = stan::math::wiener_lcdf_unnorm(y, a, t0, w, v, sv, sw, st0);
    lcdf.grad();
    EXPECT_NEAR(lcdf.val(), true_lcdf[i], err_tol_dens);
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
