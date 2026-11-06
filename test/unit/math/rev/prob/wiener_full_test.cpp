#include <stan/math/prim/functor.hpp>
#include <stan/math/prim/fun.hpp>
#include <stan/math/prim/prob.hpp>
#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>

#include <gtest/gtest.h>
#include <cmath>
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
   * Except the third point: there the ccdf is 1.6e-14 of P and the
   * difference above loses digits (WienR value -33.8123, gradient w
   * 1.9329). The values for it are from the large-time series of P - F
   * (199 terms, 50 digits, mpmath), gradients by mpmath.diff.
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
      = {-1.92033254902358,   -13.5626831740558, -33.795600277824003,
         -0.0207304574971297, -7.08470781580343, -3.44486241813253,
         -4.18311596420185,   -4.95366166391584, -7.23109216780916,
         -0.125273337362733,  -3.79815495569855, -2.96111699037856,
         -4.14929265642466};
  std::vector<double> true_grad_y
      = {-1.33001519743633,  -31.1083518155787, -12.889970183576842,
         -0.169566706752058, -3.7986585460178,  -2.31784674820993,
         -0.536530793068282, -2.48164917576694, -1.89222689054599,
         -0.484444809563098, -5.31460285637481, -0.975151391419348,
         -1.57127046017878};
  std::vector<double> true_grad_a
      = {-0.432928212420575, -0.0286062052702194, 0.35604435682103767,
         0.05132866977503,   10.5080654585197,    4.37077488076664,
         -0.193420148034756, 5.36345669388884,    5.14057131554871,
         0.11692460737807,   7.07166462218294,    0.365404340761721,
         1.53902298927608};
  std::vector<double> true_grad_t0
      = {1.33001519743633,  31.1083518155787, 12.889970183576842,
         0.169566706752058, 3.7986585460178,  2.31784674820993,
         0.536530793068282, 2.48164917576694, 1.89222689054599,
         0.484444809563098, 5.31460285637481, 0.975151391419348,
         1.57127046017878};
  std::vector<double> true_grad_w
      = {4.61303457108487,  0.611152450869066,  2.1023852348638125,
         0.139815038667115, -0.911757809733564, 1.54221437873658,
         3.63420918438666,  -1.85952752467816,  -3.65499582974638,
         -1.12442042700282, 0.262989430764412,  1.3830938669116,
         0.518345884055846};
  std::vector<double> true_grad_v
      = {1.83079071623625,   2.56351255240926,    11.40027803089135,
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

TEST_F(AgradRev, ProbWienerFull_density_below_smallest_double) {
  /* Densities far below the smallest normal double at short times. With
   * st0 > y - t0 and sw -> 0, the density is F5(y - t0) / st0, where F5 is
   * the distribution function of the 5-parameter model. Hence
   * d/dst0 = -1 / st0 and d/dy = -d/dt0 = p5(y - t0) / F5(y - t0). The
   * references for log(F5 / st0) and p5 / F5 are from Boost Gauss-Kronrod
   * quadrature of the 5-parameter density, scaled by its largest value.
   * d/dsw (sw = 1e-8) is O(sw) as sw -> 0. Its reference is
   * int (w' - w) / sw dp5/dw' / int p5 over the box, by adaptive
   * Gauss-Kronrod quadrature (differentiation under the integral). The
   * end-point form (integral) / density - 1 / sw multiplies the integration
   * error by 1 / sw and is not accurate at sw = 1e-8.
   */
  using stan::math::var;
  const double a = 0.1, t0 = 0.0, w = 1e-4, v = -10.0, sv = 0.1, st0 = 1e-3;
  // y, log(F5(y) / st0), p5(y) / F5(y), d/dsw at sw = 1e-8
  const std::vector<std::vector<double>> refs = {
      {1e-6, -4998.0685384088774, 4999186937.1439915, 0.083341666361125277},
      {3e-6, -1665.4275870208894, 555376585.70015228, 0.0092735554080942972},
      {5.5435493922668729e-06, -901.03029980710346, 162610279.12970507,
       0.0027270180916692255},
      {7.5356592945287405e-06, -662.74912830365952, 87982443.012321636,
       0.0014831386423082942},
      {3e-05, -166.59848075414246, 5538827.1596829286, 0.00010898610669626657}};
  for (const auto& ref : refs) {
    for (double sw_value : {0.0, 1e-8}) {
      for (double precision : {1e-4, 1e-8}) {
        var y = ref[0];
        var t0_v = t0;
        var sw = sw_value;
        var st0_v = st0;
        var lp = stan::math::wiener_lpdf(y, a, t0_v, w, v, sv, sw, st0_v,
                                         precision);
        lp.grad();
        const double tol = precision < 1e-6 ? 1e-9 : 1e-5;
        EXPECT_NEAR(lp.val(), ref[1], tol * 10);
        EXPECT_NEAR(y.adj(), ref[2], tol * ref[2]);
        EXPECT_NEAR(t0_v.adj(), -ref[2], tol * ref[2]);
        EXPECT_FLOAT_EQ(st0_v.adj(), -1 / st0);
        if (sw_value == 0) {
          EXPECT_EQ(sw.adj(), 0);
        } else {
          EXPECT_NEAR(sw.adj(), ref[3], precision < 1e-6 ? 1e-8 : 1e-6);
        }
        stan::math::recover_memory();
      }
    }
  }
}

TEST_F(AgradRev, ProbWienerFull_sw_gradient_st0_beyond_y) {
  /* d/dsw when st0 > y - t0. The integral over the st0 variable must use the
   * limits of that variable, [0, (y - t0) / st0], not [0, 1]: with [0, 1]
   * the quadrature nodes miss the support if (y - t0) / st0 < 0.004, and
   * d/dsw = -1 / sw. The references are the 5-parameter density averaged
   * over the box by adaptive Gauss-Kronrod quadrature, and the gradients
   * from the weighted 5-parameter gradients and the end-point formulas
   * d/dsw = -1/sw + (G(w + sw/2) + G(w - sw/2)) / (2 J) and
   * d/dst0 = -1/st0 (st0 > y - t0).
   */
  using stan::math::var;
  // y, a, t0, w, v, sv, sw, st0
  const std::vector<std::vector<double>> points
      = {{0.501, 0.5, 0.5, 0.5, 1.0, 0.5, 0.2, 0.3},
         {0.501, 0.5, 0.5, 0.3, -1.0, 0.5, 0.1, 0.3}};
  // value, then d/dy, d/da, d/dt0, d/dw, d/dv, d/dsv, d/dsw, d/dst0
  const std::vector<std::vector<double>> refs
      = {{-23.72203172192728, 21439.651877785887, -85.330940582416858,
          -21439.651877785887, 104.16367572856434, 0.20369451548936232,
          0.020277897748002094, 47.08183786536835, -3.3333333333333335},
         {-57.294953191342536, 54282.626135080958, -217.73468386330191,
          -54282.626135080958, 165.94975734995256, 0.32888822649013677,
          0.053597130002326632, 72.974882133158005, -3.3333333333333335}};
  for (size_t i = 0; i < points.size(); ++i) {
    std::vector<var> x(points[i].begin(), points[i].end());
    var lp = stan::math::wiener_lpdf(x[0], x[1], x[2], x[3], x[4], x[5], x[6],
                                     x[7]);
    lp.grad();
    EXPECT_NEAR(lp.val(), refs[i][0], 1e-6);
    for (size_t k = 0; k < x.size(); ++k) {
      EXPECT_NEAR(x[k].adj(), refs[i][k + 1],
                  1e-6 * std::fmax(1.0, std::fabs(refs[i][k + 1])))
          << "point " << i << ", argument " << k;
    }
    stan::math::recover_memory();
  }
  // A density near the smallest normal double, sharply peaked in the st0
  // direction: the 2-D integration error is 6.4e-4 here, so the tolerance is
  // looser.
  var y = 0.200001;
  var sw = 0.25;
  var lp = stan::math::wiener_lpdf(y, 0.1, 0.2, 0.5, -30.0, 0.1, sw, 1e-3);
  lp.grad();
  EXPECT_NEAR(lp.val(), -708.03896211156052, 2e-3);
  EXPECT_NEAR(sw.adj(), 1875.1599297689354, 2e-3 * 1875.16);
  EXPECT_NEAR(y.adj(), 704621632.51963747, 2e-3 * 704621632.5);
}

namespace {
// value, then d/dy, d/da, d/dt0, d/dw, d/dv, d/dsv, d/dsw, d/dst0 of the
// 7-parameter function at each point (y, a, t0, w, v, sv, sw, st0)
template <typename F>
void expect_wiener7_refs(const F& f,
                         const std::vector<std::vector<double>>& points,
                         const std::vector<std::vector<double>>& refs,
                         double value_tol, double grad_tol) {
  using stan::math::var;
  ASSERT_EQ(refs.size(), points.size());
  for (size_t i = 0; i < points.size(); ++i) {
    ASSERT_EQ(refs[i].size(), points[i].size() + 1);
    std::vector<var> x(points[i].begin(), points[i].end());
    var lp = f(x);
    lp.grad();
    EXPECT_NEAR(lp.val(), refs[i][0], value_tol) << "point " << i;
    for (size_t k = 0; k < x.size(); ++k) {
      EXPECT_NEAR(x[k].adj(), refs[i][k + 1],
                  grad_tol * std::fmax(1.0, std::fabs(refs[i][k + 1])))
          << "point " << i << ", argument " << k;
    }
    stan::math::recover_memory();
  }
}
}  // namespace

TEST_F(AgradRev, ProbWienerFull_density_below_smallest_double_long_times) {
  /* Densities far below the smallest normal double at long times
   * (y - t0 = 10 and 50), with a = 0.1 and 8 and v = -30 and -10, at
   * precision 1e-8.
   * References: the 5-parameter density (drift variability in closed form)
   * averaged over the box by mpmath quadrature with 30 digits, gradients by
   * central differences of it; the derivative with respect to a zero
   * variability is 0. An adaptive Gauss-Kronrod quadrature of the same
   * density agrees to 1e-8.
   */
  using stan::math::var;
  // y, a, t0, w, v, sv, sw, st0
  const std::vector<std::vector<double>> points
      = {{10.2, 0.1, 0.2, 1e-4, -30.0, 0.1, 5e-5, 1e-3},
         {50.2, 0.1, 0.2, 1e-4, -10.0, 0.1, 5e-5, 0.0},
         {50.2, 8.0, 0.2, 1e-4, -30.0, 0.1, 5e-5, 1e-3}};
  const std::vector<std::vector<double>> refs
      = {{-9030.3379144485420, -865.364663872592, 98643.1370267705,
          865.364663872592, 10002.6700507408, 272.804010205697,
          7441.29375966739, 0.227254900812711, 494.322285309951},
         {-26343.862677660422, -515.70133184146, 493453.554734888,
          515.70133184146, 10000.6523835387, 333.399993194435, 11112.2222128716,
          0.0555106500343369, 0.0},
         {-15174.807064836688, -199.016665244638, -19.2309950333399,
          199.016665244638, 10156.2502489645, 1005.32591939905,
          101064.687112651, 13.3994390340451, 102.807012091621}};
  expect_wiener7_refs(
      [](const std::vector<var>& x) {
        return stan::math::wiener_lpdf(x[0], x[1], x[2], x[3], x[4], x[5], x[6],
                                       x[7], 1e-8);
      },
      points, refs, 1e-6, 1e-6);
}

TEST_F(AgradRev, ProbWienerFull_sw_st0_gradients_small_density) {
  /* Point 0: sw > 0 and st0 = 0 with a small density. The series tolerance
   * of d/dsw must be relative to the density; an absolute tolerance gives a
   * relative error of 2e3 here. Point 1: sw and st0 non-zero with a small
   * density. Precision 1e-8.
   * References: as in ProbWienerFull_density_below_smallest_double_long_times.
   */
  using stan::math::var;
  // y, a, t0, w, v, sv, sw, st0
  const std::vector<std::vector<double>> points
      = {{10.2, 8.0, 0.2, 1e-4, -10.0, 0.1, 5e-5, 0.0},
         {10.2, 8.0, 0.2, 0.5, -10.0, 0.1, 0.25, 1e-3}};
  const std::vector<std::vector<double>> refs
      = {{-539.34627674639635, -40.6169325975188, -9.3959721858679,
          40.6169325975188, 10070.6429283228, 98.1810756505509,
          963.043270681908, 6.03293734926645, 0.0},
         {-488.42625756357737, -41.1734783364799, -3.51612680187419,
          41.1734783364799, 72.4306217206936, 93.7323936174272,
          877.668111678012, 32.2153115886327, 20.7280690602107}};
  expect_wiener7_refs(
      [](const std::vector<var>& x) {
        return stan::math::wiener_lpdf(x[0], x[1], x[2], x[3], x[4], x[5], x[6],
                                       x[7], 1e-8);
      },
      points, refs, 1e-6, 1e-6);
}

TEST_F(AgradRev, ProbWienerFull_lcdf_lccdf_two_dimensions) {
  /* 7-parameter lcdf and lccdf with two non-zero variabilities (sv and st0,
   * sv and sw). The gradients with respect to a and w integrate the
   * small-time series of the 4-parameter cdf gradients. At point 0,
   * st0 > y - t0, so the lccdf also contains the mass of the trials whose
   * non-decision time ends after y.
   * References: the 5-parameter density integrated by mpmath quadrature
   * with 30 digits,
   *   cdf  = 1/(sw st0) int dw' int_0^d p5(s | w') min(d - s, st0) ds,
   *   ccdf = 1/(sw st0) int dw' int_0^inf p5(s | w')
   *          clamp(s - d + st0, 0, st0) ds,
   * with d = y - t0 and without the average over a zero variability;
   * gradients by central differences of it, and 0 for a zero variability.
   * An adaptive Gauss-Kronrod quadrature agrees to 1e-9.
   */
  using stan::math::var;
  // y, a, t0, w, v, sv, sw, st0
  const std::vector<std::vector<double>> points
      = {{0.20753147706690653, 0.96248012425856511, 0.11408699836935123,
          0.20587307426609611, 2.6537214716508419, 0.80835488843424108, 0.0,
          0.21064846346879404},
         {0.54520091427708517, 1.7495584716477295, 0.30505332833140775,
          0.29376778741039961, 2.2034724592953872, 0.56060550137195397,
          0.14611719349304031, 0.0}};
  const std::vector<std::vector<double>> lcdf_refs
      = {{-4.9647686969783789, 52.5586358384741, -6.37930101549765,
          -52.5586358384741, 7.95074990826158, 0.56917143940557,
          0.213094214615845, 0.0, -4.74724564106845},
         {-2.1004294695969092, 12.5103886780412, -2.41021208698703,
          -12.5103886780412, 6.05817203831149, 0.755694568011582,
          0.226041551017222, 0.323581506000001, 0.0}};
  const std::vector<std::vector<double>> lccdf_refs
      = {{-0.45507152221964691, -0.578239142061406, 0.60207029331835,
          0.578239142061406, 2.63491457462573, 0.201109932917652,
          -0.0498196829121461, 0.0, 0.052228205752569},
         {-0.28525919173166422, -2.03681867081548, 0.567349515332851,
          2.03681867081548, 0.144651269830577, 0.0364920550026578,
          -0.11785874504401, -0.149632279530551, 0.0}};
  expect_wiener7_refs(
      [](const std::vector<var>& x) {
        return stan::math::wiener_lcdf_unnorm(x[0], x[1], x[2], x[3], x[4],
                                              x[5], x[6], x[7]);
      },
      points, lcdf_refs, 1e-7, 1e-6);
  expect_wiener7_refs(
      [](const std::vector<var>& x) {
        return stan::math::wiener_lccdf_unnorm(x[0], x[1], x[2], x[3], x[4],
                                               x[5], x[6], x[7]);
      },
      points, lccdf_refs, 1e-7, 1e-6);
}

TEST_F(AgradRev, ProbWienerFullLccdf_late_start_and_small_ccdf) {
  /* Points 0 and 1: st0 > y - t0. A trial whose non-decision time ends
   * after y has not responded at y, so it contributes P (the probability
   * of the upper boundary) to the ccdf. Points 2 and 3: small ccdf, where
   * P - F loses digits (point 2 gave -inf).
   * References: with d = y - t0, the 5-parameter density p5 (drift
   * variability in closed form) integrated by adaptive Gauss-Kronrod
   * quadrature,
   *   ccdf = 1/(sw st0) int dw' int_{d - st0}^inf p5(s | w')
   *          clamp(s - d + st0, 0, st0) ds,
   * which has no cancellation; gradients by central finite differences of
   * this reference (step 1e-5 max(|x|, 0.05)). The s integral is split at
   * the kink of the weight, s = d.
   */
  using stan::math::var;
  // y, a, t0, w, v, sv, sw, st0
  const std::vector<std::vector<double>> points
      = {{0.6, 1.0, 0.5, 0.5, 1.0, 0.0, 0.0, 0.3},
         {0.6, 1.5, 0.5, 0.4, -1.0, 0.8, 0.2, 0.25},
         {3.0, 0.6, 0.1, 0.6, 5.0, 0.4, 0.0, 0.0},
         {2.5, 0.8, 0.2, 0.65, -3.0, 0.5, 0.1, 0.1}};
  // value, then d/dy, d/da, d/dt0, d/dw, d/dv, d/dsv, d/dsw, d/dst0
  const std::vector<std::vector<double>> refs
      = {{-0.34070836357439016, -0.850496615431, 0.3784844401, 0.850496615445,
          0.992216360718, 0.264038728681, 0, 0, 0.0927560195689},
         {-1.739587600152078, -0.071706798351, -0.379766643732, 0.0717067982769,
          3.13857732603, 1.04914003454, 0.584703035322, 0.0865195797228,
          0.00522743821563},
         {-64.741572165806787, -19.6190832883, 132.880948289, 19.6190832895,
          -3.07426977741, -9.85932868275, 38.0863232508, 0, 0},
         {-25.325740009943843, -9.49511893744, 42.1178255789, 9.49511893822,
          -0.138690615232, 4.61306469729, 9.90862435266, -0.103076095215,
          5.49202370514}};
  for (size_t i = 0; i < points.size(); ++i) {
    std::vector<var> x(points[i].begin(), points[i].end());
    var lp = stan::math::wiener_lccdf_unnorm(x[0], x[1], x[2], x[3], x[4], x[5],
                                             x[6], x[7]);
    lp.grad();
    EXPECT_NEAR(lp.val(), refs[i][0], 1e-6) << "point " << i;
    for (size_t k = 0; k < x.size(); ++k) {
      EXPECT_NEAR(x[k].adj(), refs[i][k + 1],
                  2e-5 * std::fmax(1.0, std::fabs(refs[i][k + 1])))
          << "point " << i << ", argument " << k;
    }
    stan::math::recover_memory();
  }
}

TEST_F(AgradRev, ProbWienerFull_lcdf_lccdf_small_values) {
  /* 7-parameter lccdf with drift variability only, where the ccdf is far
   * below 1e-16 times the probability P of the upper boundary, so that the
   * form P - F cancels; the last point has a ccdf of 6.6e-308.
   * References: the 5-parameter density (drift variability in closed form)
   * integrated from y - t0 to infinity by adaptive Gauss-Kronrod quadrature,
   * gradients by central finite differences (step 1e-5 max(|x|, 0.05));
   * the first and the last point by mpmath quadrature with 30 digits, which
   * agrees with the Gauss-Kronrod values to 3e-9.
   * Then a ccdf and a cdf whose integral underflows to 0: the value is -inf,
   * the gradient integrals, whose tolerances are relative to the value, are
   * skipped, and the gradients are set to 0, so they are finite.
   */
  using stan::math::var;
  // y, a, t0, w, v, sv
  const std::vector<std::vector<double>> points
      = {{3.1, 2.0, 0.1, 0.5, 5.0, 0.01},
         {9.0, 2.0, 0.1, 0.5, 5.0, 0.01},
         {22.1, 2.0, 0.1, 0.5, 5.0, 0.01},
         {52.0, 2.0, 0.1, 0.5, 5.0, 0.01}};
  // value, then d/dy, d/da, d/dt0, d/dw, d/dv, d/dsv
  const std::vector<std::vector<double>> refs = {
      {-39.052352532661321, -13.7265716431519, 5.29026166529988,
       13.7265716431519, -9.99697894788928, -14.3598466219783, 2.0326581579936},
      {-119.99567604004699, -13.7118474462, 12.5677181533, 13.7118474512,
       -9.99123493273, -43.8253261844, 19.1182682983},
      {-299.40720932606365, -13.679247658, 28.7261448207, 13.6792476439,
       -9.97817507482, -109.124647395, 118.86296619},
      {-707.30981035866349, -13.6053173797974, 65.6068611282751,
       13.6053173797974, -9.94849419259966, -257.529037001709,
       662.696324212265}};
  for (size_t i = 0; i < points.size(); ++i) {
    std::vector<var> x(points[i].begin(), points[i].end());
    var lp = stan::math::wiener_lccdf_unnorm(x[0], x[1], x[2], x[3], x[4], x[5],
                                             0.0, 0.0);
    lp.grad();
    EXPECT_NEAR(lp.val(), refs[i][0], 1e-8) << "point " << i;
    for (size_t k = 0; k < x.size(); ++k) {
      EXPECT_NEAR(x[k].adj(), refs[i][k + 1],
                  1e-6 * std::fmax(1.0, std::fabs(refs[i][k + 1])))
          << "point " << i << ", argument " << k;
    }
    stan::math::recover_memory();
  }

  // y, a, t0, w, v, sv
  std::vector<var> x_ccdf{60.0, 2.0, 0.1, 0.5, 5.0, 0.01};
  var lccdf = stan::math::wiener_lccdf_unnorm(x_ccdf[0], x_ccdf[1], x_ccdf[2],
                                              x_ccdf[3], x_ccdf[4], x_ccdf[5],
                                              0.0, 0.0);
  lccdf.grad();
  EXPECT_EQ(lccdf.val(), stan::math::NEGATIVE_INFTY);
  for (size_t k = 0; k < x_ccdf.size(); ++k) {
    EXPECT_TRUE(std::isfinite(x_ccdf[k].adj())) << "lccdf, argument " << k;
  }
  stan::math::recover_memory();

  std::vector<var> x_cdf{1e-6, 0.1, 0.0, 1e-4, -10.0, 0.1};
  var lcdf = stan::math::wiener_lcdf_unnorm(
      x_cdf[0], x_cdf[1], x_cdf[2], x_cdf[3], x_cdf[4], x_cdf[5], 0.0, 0.0);
  lcdf.grad();
  EXPECT_EQ(lcdf.val(), stan::math::NEGATIVE_INFTY);
  for (size_t k = 0; k < x_cdf.size(); ++k) {
    EXPECT_TRUE(std::isfinite(x_cdf[k].adj())) << "lcdf, argument " << k;
  }
  stan::math::recover_memory();
}
