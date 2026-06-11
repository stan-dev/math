#include <stan/math/prim/functor.hpp>
#include <stan/math/prim/fun.hpp>
#include <stan/math/prim/prob.hpp>
#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>

#include <gtest/gtest.h>
#include <vector>

// CHECK THAT ALL VALID SCALAR TYPES ARE ACCEPTED
template <typename F>
inline void check_scalar_types(F& f, double value, double res, double deriv) {
  // - f: Function with a single parameter exposed, all others
  // have to be scalars
  // - value: value to be used for the parameter
  // - res: expected result of calling `f` with `value`
  // - deriv: expected result of partial of f with respect to
  // the parameter in `value`
  using stan::math::var;
  double err_tol = 2e-4;

  // type double
  EXPECT_NEAR(f(value), res, err_tol);

  // type var with derivative
  var value_var = value;
  var result_var = f(value_var);
  result_var.grad();
  EXPECT_NEAR(value_of(result_var), res, err_tol);
  EXPECT_NEAR(value_var.adj(), deriv, err_tol);
}

TEST_F(AgradRev, ProbWienerFull_wiener_full_all_scalar) {
  // tests all parameter types individually, with other
  // parameters set to double
  using stan::math::wiener_lpdf;

  std::vector<double> rt{1, 1, 1, 1, 1, 1, 1, 1, 1};
  std::vector<double> a{1, 1, 1, 1, 1, 1, 1, 1, 1};
  std::vector<double> v{1, 1, 1, 1, 1, 1, 1, 1, 1};
  std::vector<double> w{0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
  std::vector<double> t0{0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0};
  std::vector<double> sv{0.1, 0.1, 0.1, 0, 0.1, 0, 0, 0, 0};
  std::vector<double> sw{0.1, 0.1, 0, 0.1, 0, 0.1, 0, 0, 0};
  std::vector<double> st0{0.1, 0, 0.1, 0.1, 0, 0, 0.1, 0, 0};

  std::vector<double> result{-2.426204, -2.710348, -2.422505,
                             -2.422795, -2.70665,  -2.706812,
                             -2.419095, -2.703112, -3.790072};
  std::vector<double> drt{-5.437337, -5.436798, -5.437332, -5.434799, -5.436791,
                          -5.434801, -5.434802, -5.434802, -5.434802};
  std::vector<double> da{5.857318, 6.395031, 5.856484, 5.858549, 6.394195,
                         6.396512, 5.857722, 6.395684, 8.369604};
  std::vector<double> dv{-0.2428443, -0.2967977, -0.2436664,
                         -0.2446629, -0.297619,  -0.2991696,
                         -0.2454931, -0.3,       -0.5};
  std::vector<double> dw{-0.9891369, -0.9887305, -0.9973428,
                         -0.9915453, -0.9969335, -0.991674,
                         -0.9997794, -0.9999097, -0.9999953};
  std::vector<double> dt0{5.437337, 5.436798, 5.437332, 5.434799, 5.436791,
                          5.434801, 5.434802, 5.434802, 5.434802};
  std::vector<double> dsv{-0.06793703, -0.07047449, -0.06797882,
                          0.0,         -0.07050737, 0.0,
                          0.0,         0.0,         0.0};
  std::vector<double> dsw{-0.07406705, -0.07407386, 0.0, -0.07410901, 0.0,
                          -0.07410686, 0.0,         0.0, 0.0};
  std::vector<double> dst0{2.963915, 0.0,     2.963912, 2.962338, 0.0,
                           0.0,      2.96234, 0.0,      0.0};

  for (int i = 0; i < result.size(); i++) {
    // rt
    auto f_rt = [&](auto value) {
      return wiener_lpdf(value, a[i], t0[i], w[i], v[i], sv[i], sw[i], st0[i]);
    };
    check_scalar_types(f_rt, rt[i], result[i], drt[i]);
    // a
    auto f_a = [&](auto value) {
      return wiener_lpdf(rt[i], value, t0[i], w[i], v[i], sv[i], sw[i], st0[i]);
    };
    check_scalar_types(f_a, a[i], result[i], da[i]);
    // v
    auto f_v = [&](auto value) {
      return wiener_lpdf(rt[i], a[i], t0[i], w[i], value, sv[i], sw[i], st0[i]);
    };
    check_scalar_types(f_v, v[i], result[i], dv[i]);
    // w
    auto f_w = [&](auto value) {
      return wiener_lpdf(rt[i], a[i], t0[i], value, v[i], sv[i], sw[i], st0[i]);
    };
    check_scalar_types(f_w, w[i], result[i], dw[i]);
    // t0
    auto f_t0 = [&](auto value) {
      return wiener_lpdf(rt[i], a[i], value, w[i], v[i], sv[i], sw[i], st0[i]);
    };
    check_scalar_types(f_t0, t0[i], result[i], dt0[i]);
    // sv
    auto f_sv = [&](auto value) {
      return wiener_lpdf(rt[i], a[i], t0[i], w[i], v[i], value, sw[i], st0[i]);
    };
    check_scalar_types(f_sv, sv[i], result[i], dsv[i]);
    // sw
    auto f_sw = [&](auto value) {
      return wiener_lpdf(rt[i], a[i], t0[i], w[i], v[i], sv[i], value, st0[i]);
    };
    check_scalar_types(f_sw, sw[i], result[i], dsw[i]);
    // st0
    auto f_st0 = [&](auto value) {
      return wiener_lpdf(rt[i], a[i], t0[i], w[i], v[i], sv[i], sw[i], value);
    };
    check_scalar_types(f_st0, st0[i], result[i], dst0[i]);
  }
}

TEST_F(AgradRev, ProbWienerFullPrec_wiener_full_prec_all_scalar) {
  // tests all parameter types individually, with other parameters
  // set to double
  using stan::math::wiener_lpdf;
  double err_tol = 2e-6;

  std::vector<double> rt{1, 1, 1, 1, 1, 1, 1, 1, 1};
  std::vector<double> a{1, 1, 1, 1, 1, 1, 1, 1, 1};
  std::vector<double> v{1, 1, 1, 1, 1, 1, 1, 1, 1};
  std::vector<double> w{0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
  std::vector<double> t0{0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0};
  std::vector<double> sv{0.1, 0.1, 0.1, 0, 0.1, 0, 0, 0, 0};
  std::vector<double> sw{0.1, 0.1, 0, 0.1, 0, 0.1, 0, 0, 0};
  std::vector<double> st0{0.1, 0, 0.1, 0.1, 0, 0, 0.1, 0, 0};

  std::vector<double> result{-2.426204, -2.710348, -2.422505,
                             -2.422795, -2.70665,  -2.706812,
                             -2.419095, -2.703112, -3.790072};
  std::vector<double> drt{-5.437337, -5.436798, -5.437332, -5.434799, -5.436791,
                          -5.434801, -5.434802, -5.434802, -5.434802};
  std::vector<double> da{5.857318, 6.395031, 5.856484, 5.858549, 6.394195,
                         6.396512, 5.857722, 6.395684, 8.369604};
  std::vector<double> dv{-0.2428443, -0.2967977, -0.2436664,
                         -0.2446629, -0.297619,  -0.2991696,
                         -0.2454931, -0.3,       -0.5};
  std::vector<double> dw{-0.9891369, -0.9887305, -0.9973428,
                         -0.9915453, -0.9969335, -0.991674,
                         -0.9997794, -0.9999097, -0.9999953};
  std::vector<double> dt0{5.437337, 5.436798, 5.437332, 5.434799, 5.436791,
                          5.434801, 5.434802, 5.434802, 5.434802};
  std::vector<double> dsv{-0.06793703, -0.07047449, -0.06797882,
                          0.0,         -0.07050737, 0.0,
                          0.0,         0.0,         0.0};
  std::vector<double> dsw{-0.07406705, -0.07407386, 0.0, -0.07410901, 0.0,
                          -0.07410686, 0.0,         0.0, 0.0};
  std::vector<double> dst0{2.963915, 0.0,     2.963912, 2.962338, 0.0,
                           0.0,      2.96234, 0.0,      0.0};

  for (int i = 0; i < result.size(); i++) {
    // rt
    auto f_rt = [&](auto value) {
      return wiener_lpdf(value, a[i], t0[i], w[i], v[i], sv[i], sw[i], st0[i],
                         1e-6);
    };
    check_scalar_types(f_rt, rt[i], result[i], drt[i]);
    // a
    auto f_a = [&](auto value) {
      return wiener_lpdf(rt[i], value, t0[i], w[i], v[i], sv[i], sw[i], st0[i],
                         1e-6);
    };
    check_scalar_types(f_a, a[i], result[i], da[i]);
    // v
    auto f_v = [&](auto value) {
      return wiener_lpdf(rt[i], a[i], t0[i], w[i], value, sv[i], sw[i], st0[i],
                         1e-6);
    };
    check_scalar_types(f_v, v[i], result[i], dv[i]);
    // w
    auto f_w = [&](auto value) {
      return wiener_lpdf(rt[i], a[i], t0[i], value, v[i], sv[i], sw[i], st0[i],
                         1e-6);
    };
    check_scalar_types(f_w, w[i], result[i], dw[i]);
    // t0
    auto f_t0 = [&](auto value) {
      return wiener_lpdf(rt[i], a[i], value, w[i], v[i], sv[i], sw[i], st0[i],
                         1e-6);
    };
    check_scalar_types(f_t0, t0[i], result[i], dt0[i]);
    // sv
    auto f_sv = [&](auto value) {
      return wiener_lpdf(rt[i], a[i], t0[i], w[i], v[i], value, sw[i], st0[i],
                         1e-6);
    };
    check_scalar_types(f_sv, sv[i], result[i], dsv[i]);
    // sw
    auto f_sw = [&](auto value) {
      return wiener_lpdf(rt[i], a[i], t0[i], w[i], v[i], sv[i], value, st0[i],
                         1e-6);
    };
    check_scalar_types(f_sw, sw[i], result[i], dsw[i]);
    // st0
    auto f_st0 = [&](auto value) {
      return wiener_lpdf(rt[i], a[i], t0[i], w[i], v[i], sv[i], sw[i], value,
                         1e-6);
    };
    check_scalar_types(f_st0, st0[i], result[i], dst0[i]);
  }
}

// CHECK THAT ALL VALID Vector TYPES ARE ACCEPTED
template <typename F>
inline void check_vector_types(F& f, std::vector<double> value, double res) {
  // - f: Function where all inputs are vectors
  // - value: value to be used for the parameter
  // - res: expected result of calling `f` with `value`
  // - deriv: expected result of partial of f with respect to
  // the parameter in `value`
  using stan::math::var;
  double err_tol = 2e-4;

  // type double
  EXPECT_NEAR(f(value), res, err_tol);

  // type var with derivative
  var result_var = f(value);
  result_var.grad();
  EXPECT_NEAR(value_of(result_var), res, err_tol);
}

TEST_F(AgradRev, ProbWienerFull_wiener_full_all_vector) {
  // tests all parameter types individually, with other
  // parameters set to std::vector<double>
  using stan::math::wiener_lpdf;

  std::vector<double> rt{1, 1, 1, 1, 1, 1, 1, 1, 1};
  std::vector<double> a{1, 1, 1, 1, 1, 1, 1, 1, 1};
  std::vector<double> v{1, 1, 1, 1, 1, 1, 1, 1, 1};
  std::vector<double> w{0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
  std::vector<double> t0{0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0};
  std::vector<double> sv{0.1, 0.1, 0.1, 0, 0.1, 0, 0, 0, 0};
  std::vector<double> sw{0.1, 0.1, 0, 0.1, 0, 0.1, 0, 0, 0};
  std::vector<double> st0{0.1, 0, 0.1, 0.1, 0, 0, 0.1, 0, 0};

  double result{-24.307593};

  // rt
  auto f_rt = [&](auto value) {
    return wiener_lpdf(value, a, t0, w, v, sv, sw, st0);
  };
  check_vector_types(f_rt, rt, result);
  // a
  auto f_a = [&](auto value) {
    return wiener_lpdf(rt, value, t0, w, v, sv, sw, st0);
  };
  check_vector_types(f_a, a, result);
  // v
  auto f_v = [&](auto value) {
    return wiener_lpdf(rt, a, t0, w, value, sv, sw, st0);
  };
  check_vector_types(f_v, v, result);
  // w
  auto f_w = [&](auto value) {
    return wiener_lpdf(rt, a, t0, value, v, sv, sw, st0);
  };
  check_vector_types(f_w, w, result);
  // t0
  auto f_t0 = [&](auto value) {
    return wiener_lpdf(rt, a, value, w, v, sv, sw, st0);
  };
  check_vector_types(f_t0, t0, result);
  // sv
  auto f_sv = [&](auto value) {
    return wiener_lpdf(rt, a, t0, w, v, value, sw, st0);
  };
  check_vector_types(f_sv, sv, result);
  // sw
  auto f_sw = [&](auto value) {
    return wiener_lpdf(rt, a, t0, w, v, sv, value, st0);
  };
  check_vector_types(f_sw, sw, result);
  // st0
  auto f_st0
      = [&](auto value) { return wiener_lpdf(rt, a, t0, w, v, sv, sw, value); };
  check_vector_types(f_st0, st0, result);
}

TEST_F(AgradRev, mathRevCorrectValues_wiener_lccdf_unnorm) {
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
  double err_tol_dens = 1e-4;
  double err_tol = 1e-1;
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
