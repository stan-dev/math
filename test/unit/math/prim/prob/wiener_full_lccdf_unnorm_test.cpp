#include <stan/math/prim/functor.hpp>
#include <stan/math/prim/fun.hpp>
#include <stan/math/prim/prob.hpp>
#include <stan/math/rev.hpp>

#include <gtest/gtest.h>
#include <vector>

TEST(mathPrimScalProbWienerFullLccdfScal, valid) {
  using stan::math::INFTY;
  using stan::math::wiener_lccdf_unnorm;
  double rt = 1, a = 1, v = -1, w = 0.5, t0 = 0.1, sv = 0.2, sw = 0.2,
         st0 = 0.1;
  EXPECT_NO_THROW(wiener_lccdf_unnorm(rt, a, t0, w, v, sv, sw, st0));
  rt = 5;
  a = 1;
  v = 1;
  w = 0.5;
  t0 = 0.0;
  sv = 0.0;
  sw = 0.0;
  st0 = 0.0;
  EXPECT_NO_THROW(wiener_lccdf_unnorm(rt, a, t0, w, v, sv, sw, st0));
}

// rt
TEST(mathPrimScalProbWienerFullLccdfScal, invalid_rt) {
  using stan::math::INFTY;
  using stan::math::wiener_lccdf_unnorm;
  double a = 1, v = -1, w = 0.5, t0 = 0.1, sv = 0.2, sw = 0.2, st0 = 0.1;
  EXPECT_THROW(wiener_lccdf_unnorm(0, a, t0, w, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(-1, a, t0, w, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(INFTY, a, t0, w, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(-INFTY, a, t0, w, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(NAN, a, t0, w, v, sv, sw, st0),
               std::domain_error);
}

// a
TEST(mathPrimScalProbWienerFullLccdfScal, invalid_a) {
  using stan::math::INFTY;
  using stan::math::wiener_lccdf_unnorm;
  double rt = 1, v = -1, w = 0.5, t0 = 0.1, sv = 0.2, sw = 0.2, st0 = 0.1;
  EXPECT_THROW(wiener_lccdf_unnorm(rt, 0, t0, w, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, -1, t0, w, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, INFTY, t0, w, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, -INFTY, t0, w, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, NAN, t0, w, v, sv, sw, st0),
               std::domain_error);
}

// v
TEST(mathPrimScalProbWienerFullLccdfScal, invalid_v) {
  using stan::math::INFTY;
  using stan::math::wiener_lccdf_unnorm;
  double rt = 1, a = 1, v = -1, w = 0.5, t0 = 0.1, sv = 0.2, sw = 0.2,
         st0 = 0.1;
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, w, INFTY, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, w, -INFTY, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, w, NAN, sv, sw, st0),
               std::domain_error);
}

// w
TEST(mathPrimScalProbWienerFullLccdfScal, invalid_w) {
  using stan::math::INFTY;
  using stan::math::wiener_lccdf_unnorm;
  double rt = 1, a = 1, v = -1, t0 = 0.1, sv = 0.2, sw = 0.2, st0 = 0.1;
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, -0.1, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, 0, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, 1, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, 1.1, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, INFTY, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, -INFTY, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, NAN, v, sv, sw, st0),
               std::domain_error);
}

// t0
TEST(mathPrimScalProbWienerFullLccdfScal, invalid_t0) {
  using stan::math::INFTY;
  using stan::math::wiener_lccdf_unnorm;
  double rt = 1, a = 1, v = -1, w = 0.5, sv = 0.2, sw = 0.2, st0 = 0.1;
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, 2, w, v, sv, sw, st0),
               std::domain_error);  // rt must be greater than t0
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, -1, w, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, INFTY, w, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, -INFTY, w, v, sv, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, NAN, w, v, sv, sw, st0),
               std::domain_error);
}

// sv
TEST(mathPrimScalProbWienerFullLccdfScal, invalid_sv) {
  using stan::math::INFTY;
  using stan::math::wiener_lccdf_unnorm;
  double rt = 1, a = 1, v = -1, w = 0.5, t0 = 0.1, sw = 0.2, st0 = 0.1;
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, w, v, -1, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, w, v, INFTY, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, w, v, -INFTY, sw, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, w, v, NAN, sw, st0),
               std::domain_error);
}

// sw
TEST(mathPrimScalProbWienerFullLccdfScal, invalid_sw) {
  using stan::math::INFTY;
  using stan::math::wiener_lccdf_unnorm;
  double rt = 1, a = 1, v = -1, w = 0.5, t0 = 0.1, sv = 0.2, st0 = 0.1;
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, w, v, sv, -1, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, 0.8, v, sv, 0.5, st0),
               std::domain_error);  // sw must be smaller than 2*(1-w)
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, 0.3, v, sv, 0.7, st0),
               std::domain_error);  // sw must be smaller than 2*w
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, w, v, sv, INFTY, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, w, v, sv, -INFTY, st0),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, w, v, sv, NAN, st0),
               std::domain_error);
}

// st0
TEST(mathPrimScalProbWienerFullLccdfScal, invalid_st0) {
  using stan::math::INFTY;
  using stan::math::wiener_lccdf_unnorm;
  double rt = 1, a = 1, v = -1, w = 0.5, t0 = 0.1, sv = 0.2, sw = 0.2;
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, w, v, sv, sw, -1),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, w, v, sv, sw, INFTY),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, w, v, sv, sw, -INFTY),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, w, v, sv, sw, NAN),
               std::domain_error);
}

TEST(mathPrimScalProbWienerFullLccdfPrecScal, valid) {
  using stan::math::INFTY;
  using stan::math::wiener_lccdf_unnorm;
  double rt = 1, a = 1, v = -1, w = 0.5, t0 = 0.1, sv = 0.2, sw = 0.2,
         st0 = 0.1;
  EXPECT_NO_THROW(wiener_lccdf_unnorm(rt, a, t0, w, v, sv, sw, st0, 1e-4));
  rt = 5;
  a = 1;
  v = 1;
  w = 0.5;
  t0 = 0.0;
  sv = 0.0;
  sw = 0.0;
  st0 = 0.0;
  EXPECT_NO_THROW(wiener_lccdf_unnorm(rt, a, t0, w, v, sv, sw, st0, 1e-4));
}

// rt
TEST(mathPrimScalProbWienerFullLccdfPrecScal, invalid_rt) {
  using stan::math::INFTY;
  using stan::math::wiener_lccdf_unnorm;
  double a = 1, v = -1, w = 0.5, t0 = 0.1, sv = 0.2, sw = 0.2, st0 = 0.1;
  EXPECT_THROW(wiener_lccdf_unnorm(0, a, t0, w, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(-1, a, t0, w, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(INFTY, a, t0, w, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(-INFTY, a, t0, w, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(NAN, a, t0, w, v, sv, sw, st0, 1e-4),
               std::domain_error);
}

// a
TEST(mathPrimScalProbWienerFullLccdfPrecScal, invalid_a) {
  using stan::math::INFTY;
  using stan::math::wiener_lccdf_unnorm;
  double rt = 1, v = -1, w = 0.5, t0 = 0.1, sv = 0.2, sw = 0.2, st0 = 0.1;
  EXPECT_THROW(wiener_lccdf_unnorm(rt, 0, t0, w, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, -1, t0, w, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, INFTY, t0, w, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, -INFTY, t0, w, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, NAN, t0, w, v, sv, sw, st0, 1e-4),
               std::domain_error);
}

// v
TEST(mathPrimScalProbWienerFullLccdfPrecScal, invalid_v) {
  using stan::math::INFTY;
  using stan::math::wiener_lccdf_unnorm;
  double rt = 1, a = 1, v = -1, w = 0.5, t0 = 0.1, sv = 0.2, sw = 0.2,
         st0 = 0.1;
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, w, INFTY, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, w, -INFTY, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, w, NAN, sv, sw, st0, 1e-4),
               std::domain_error);
}

// w
TEST(mathPrimScalProbWienerFullLccdfPrecScal, invalid_w) {
  using stan::math::INFTY;
  using stan::math::wiener_lccdf_unnorm;
  double rt = 1, a = 1, v = -1, t0 = 0.1, sv = 0.2, sw = 0.2, st0 = 0.1;
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, -0.1, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, 0, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, 1, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, 1.1, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, INFTY, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, -INFTY, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, NAN, v, sv, sw, st0, 1e-4),
               std::domain_error);
}

// t0
TEST(mathPrimScalProbWienerFullLccdfPrecScal, invalid_t0) {
  using stan::math::INFTY;
  using stan::math::wiener_lccdf_unnorm;
  double rt = 1, a = 1, v = -1, w = 0.5, sv = 0.2, sw = 0.2, st0 = 0.1;
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, 2, w, v, sv, sw, st0, 1e-4),
               std::domain_error);  // rt must be greater than t0
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, -1, w, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, INFTY, w, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, -INFTY, w, v, sv, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, NAN, w, v, sv, sw, st0, 1e-4),
               std::domain_error);
}

// sv
TEST(mathPrimScalProbWienerFullLccdfPrecScal, invalid_sv) {
  using stan::math::INFTY;
  using stan::math::wiener_lccdf_unnorm;
  double rt = 1, a = 1, v = -1, w = 0.5, t0 = 0.1, sw = 0.2, st0 = 0.1;
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, w, v, -1, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, w, v, INFTY, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, w, v, -INFTY, sw, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, w, v, NAN, sw, st0, 1e-4),
               std::domain_error);
}

// sw
TEST(mathPrimScalProbWienerFullLccdfPrecScal, invalid_sw) {
  using stan::math::INFTY;
  using stan::math::wiener_lccdf_unnorm;
  double rt = 1, a = 1, v = -1, w = 0.5, t0 = 0.1, sv = 0.2, st0 = 0.1;
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, w, v, sv, -1, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, 0.8, v, sv, 0.5, st0, 1e-4),
               std::domain_error);  // sw must be smaller than 2*(1-w)
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, 0.3, v, sv, 0.7, st0, 1e-4),
               std::domain_error);  // sw must be smaller than 2*w
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, w, v, sv, INFTY, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, w, v, sv, -INFTY, st0, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, w, v, sv, NAN, st0, 1e-4),
               std::domain_error);
}

// st0
TEST(mathPrimScalProbWienerFullLccdfPrecScal, invalid_st0) {
  using stan::math::INFTY;
  using stan::math::wiener_lccdf_unnorm;
  double rt = 1, a = 1, v = -1, w = 0.5, t0 = 0.1, sv = 0.2, sw = 0.2;
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, w, v, sv, sw, -1, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, w, v, sv, sw, INFTY, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, w, v, sv, sw, -INFTY, 1e-4),
               std::domain_error);
  EXPECT_THROW(wiener_lccdf_unnorm(rt, a, t0, w, v, sv, sw, NAN, 1e-4),
               std::domain_error);
}

TEST(mathPrimCorrectValues, wiener_lccdf_unnorm) {
  /* Test concrete values. True values are computed in R using the R-package
   * WienR and the function WienerCDF() with its partial derivatives:
   * ccdf = cdf(big_value) - cdf
   * lccdf = log(cdf(big_value) - cdf)
   * lccdf' = ccdf'/ccdf = (cdf(big_value)'-cdf')/ccdf
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

  std::vector<double> true_lccdf
      = {-1.92033254902358,   -13.5626831740558, -33.8122958054115,
         -0.0207304574971297, -7.08470781580343, -3.44486241813253,
         -4.18311596420185,   -4.95366166391584, -7.23109216780916,
         -0.125273337362733,  -3.79815495569855, -2.96111699037856,
         -4.14929265642466};
  std::vector<double> true_grad_y
      = {-1.33001519743633,  -31.1083518155787, -12.9346392526053,
         -0.169566706752058, -3.7986585460178,  -2.31784674820993,
         -0.536530793068282, -2.48164917576694, -1.89222689054599,
         -0.484444809563098, -5.31460285637481, -0.975151391419348,
         -1.57127046017878};
  std::vector<double> true_grad_a
      = {-0.432928212420575, -0.0286062052702194, 0.308724832214765,
         0.05132866977503,   10.5080654585197,    4.37077488076664,
         -0.193420148034756, 5.36345669388884,    5.14057131554871,
         0.11692460737807,   7.07166462218294,    0.365404340761721,
         1.53902298927608};
  std::vector<double> true_grad_t0 = {
      1.33001519743633, 31.1083518155787,  12.9346392526053,  0.169566706752058,
      3.7986585460178,  2.31784674820993,  0.536530793068282, 2.48164917576694,
      1.89222689054599, 0.484444809563098, 5.31460285637481,  0.975151391419348,
      1.57127046017878};
  std::vector<double> true_grad_w
      = {4.61303457108487,  0.611152450869066,  1.93288590604027,
         0.139815038667115, -0.911757809733564, 1.54221437873658,
         3.63420918438666,  -1.85952752467816,  -3.65499582974638,
         -1.12442042700282, 0.262989430764412,  1.3830938669116,
         0.518345884055846};
  std::vector<double> true_grad_v
      = {1.83079071623625,   2.56351255240926,    11.4026845637584,
         0.0274975992798594, -0.960807197217193,  0.132382852324283,
         5.47350811619098,   -0.652456751242452,  1.57258003231325,
         -0.207055671666776, -0.0823432112415456, 2.18370562291933,
         2.62134233757004};
  std::vector<double> true_grad_sv = {0,
                                      0,
                                      0,
                                      0,
                                      -0.209370036446517,
                                      0,
                                      0,
                                      -0.430772737121245,
                                      1.14730002231356,
                                      0,
                                      -0.293744554675593,
                                      1.73837609123659,
                                      0.997525744657476};
  std::vector<double> true_grad_sw = {0,
                                      0,
                                      0,
                                      0,
                                      0,
                                      -0.106375521997502,
                                      0,
                                      -0.114806591325169,
                                      0,
                                      -0.201319331904563,
                                      -0.180444908443688,
                                      0,
                                      -0.17718676934793};
  std::vector<double> true_grad_st0 = {0,
                                       0,
                                       0,
                                       0,
                                       0,
                                       0,
                                       0.270517553006863,
                                       0,
                                       0.976409725660212,
                                       0.235042776747307,
                                       3.12081986216078,
                                       0.512275314141889,
                                       0.837649445725816};

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
    var lccdf = stan::math::wiener_lccdf_unnorm(y, a, t0, w, v, sv, sw, st0);
    lccdf.grad();
    EXPECT_NEAR(lccdf.val(), true_lccdf[i], err_tol_dens);
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
