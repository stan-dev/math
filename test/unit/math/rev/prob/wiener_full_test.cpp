#include <stan/math/prim/functor.hpp>
#include <stan/math/prim/fun.hpp>
#include <stan/math/prim/prob.hpp>
#include <stan/math/rev.hpp>

#include <gtest/gtest.h>
#include <vector>

// CHECK THAT ALL VALID SCALAR TYPES ARE ACCEPTED
template <typename F>
void check_scalar_types(F& f, double value, double res, double deriv) {
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

TEST(ProbWienerFull, wiener_full_all_scalar) {
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

TEST(ProbWienerFullPrec, wiener_full_prec_all_scalar) {
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
void check_vector_types(F& f, std::vector<double> value, double res) {
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

TEST(ProbWienerFull, wiener_full_all_vector) {
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
