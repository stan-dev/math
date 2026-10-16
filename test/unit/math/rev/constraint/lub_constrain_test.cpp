#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <type_traits>
#include <vector>

namespace lub_constrain_rev_test {
// Reference derivatives of y = lb + (ub - lb) * inv_logit(x) and of the
// log Jacobian log(ub - lb) + log(inv_logit(x)) + log(inv_logit(-x)) at
// lb = -2, ub = 3, computed with mpmath at 50 digits:
//   dy/dx = (ub - lb) * inv_logit(x) * inv_logit(-x), dy/dlb = inv_logit(-x),
//   dy/dub = inv_logit(x), dlogJ/dx = inv_logit(-x) - inv_logit(x),
//   dlogJ/dlb = -1 / (ub - lb), dlogJ/dub = 1 / (ub - lb).
// Entries below the smallest subnormal double are 0.
const std::vector<double> xs{-800, -40, -37, -30, -20, -1, 0,
                             1,    20,  30,  37,  40,  800};
const std::vector<double> dx{0.0,
                             2.1241771276457945e-17,
                             4.2665238128720322e-16,
                             4.6788114844192117e-13,
                             1.0305768069709247e-8,
                             9.8305966620740926e-1,
                             1.25,
                             9.8305966620740926e-1,
                             1.0305768069709247e-8,
                             4.6788114844192117e-13,
                             4.2665238128720322e-16,
                             2.1241771276457945e-17,
                             0.0};
const std::vector<double> dlb{1.0,
                              1.0,
                              9.9999999999999991e-1,
                              9.9999999999990642e-1,
                              9.9999999793884638e-1,
                              7.3105857863000488e-1,
                              0.5,
                              2.6894142136999512e-1,
                              2.0611536181902036e-9,
                              9.357622968839299e-14,
                              8.5330476257440651e-17,
                              4.248354255291589e-18,
                              0.0};
const std::vector<double> dub{0.0,
                              4.248354255291589e-18,
                              8.5330476257440651e-17,
                              9.357622968839299e-14,
                              2.0611536181902036e-9,
                              2.6894142136999512e-1,
                              0.5,
                              7.3105857863000488e-1,
                              9.9999999793884638e-1,
                              9.9999999999990642e-1,
                              9.9999999999999991e-1,
                              1.0,
                              1.0};
const std::vector<double> dlogj{1.0,
                                9.9999999999999999e-1,
                                9.9999999999999983e-1,
                                9.9999999999981285e-1,
                                9.9999999587769276e-1,
                                4.6211715726000976e-1,
                                0.0,
                                -4.6211715726000976e-1,
                                -9.9999999587769276e-1,
                                -9.9999999999981285e-1,
                                -9.9999999999999983e-1,
                                -9.9999999999999999e-1,
                                -1.0};
constexpr double lb = -2.0;
constexpr double ub = 3.0;

inline void expect_rel(double expected, double actual, const char* what,
                       size_t i) {
  EXPECT_NEAR(expected, actual, 1e-13 * std::fabs(expected))
      << what << " at x = " << xs[i];
}

// Checks the adjoint of the first element (or the scalar) of a.
template <typename T>
inline void expect_adj(const T& a, double expected, const char* what,
                       size_t i) {
  if constexpr (std::is_same<T, stan::math::var>::value) {
    expect_rel(expected, a.adj(), what, i);
  } else if constexpr (stan::is_autodiff_v<T>) {
    expect_rel(expected, a.adj()(0), what, i);
  }
}

template <typename TX, typename TL, typename TU>
inline void scalar_case(size_t i) {
  {
    TX x = xs[i];
    TL l = lb;
    TU u = ub;
    stan::math::var y = stan::math::lub_constrain(x, l, u);
    y.grad();
    expect_adj(x, dx[i], "dx", i);
    expect_adj(l, dlb[i], "dlb", i);
    expect_adj(u, dub[i], "dub", i);
    stan::math::recover_memory();
  }
  TX x = xs[i];
  TL l = lb;
  TU u = ub;
  stan::math::var lp = 0;
  stan::math::var y = stan::math::lub_constrain(x, l, u, lp);
  y.grad();
  expect_adj(x, dx[i], "dx, lp variant", i);
  expect_adj(l, dlb[i], "dlb, lp variant", i);
  expect_adj(u, dub[i], "dub, lp variant", i);
  stan::math::set_zero_all_adjoints();
  lp.grad();
  expect_adj(x, dlogj[i], "dlogJ/dx", i);
  expect_adj(l, -0.2, "dlogJ/dlb", i);
  expect_adj(u, 0.2, "dlogJ/dub", i);
  stan::math::recover_memory();
}

// x, lb, ub are vectors or scalars; the first element holds xs[i]. With
// tail, a second element with an infinite vector bound takes the
// overloads' mixed-infinity path.
template <typename T>
inline T make(const Eigen::VectorXd& v, double s) {
  if constexpr (stan::is_stan_scalar<T>::value) {
    return T(s);
  } else {
    return T(v);
  }
}

template <typename TX, typename TL, typename TU>
inline void vector_case(size_t i, bool tail) {
  using stan::math::var;
  constexpr bool l_vec = !stan::is_stan_scalar<TL>::value;
  constexpr bool u_vec = !stan::is_stan_scalar<TU>::value;
  const int n = tail ? 2 : 1;
  Eigen::VectorXd x_d(n), l_d(n), u_d(n);
  x_d(0) = xs[i];
  l_d(0) = lb;
  u_d(0) = ub;
  if (tail) {
    x_d(1) = 0.3;
    l_d(1) = l_vec ? stan::math::NEGATIVE_INFTY : lb;
    u_d(1) = u_vec && !l_vec ? stan::math::INFTY : ub;
  }
  for (bool use_lp : {false, true}) {
    TX x = x_d;
    TL l = make<TL>(l_d, lb);
    TU u = make<TU>(u_d, ub);
    var lp = 0;
    var s = use_lp ? stan::math::sum(stan::math::lub_constrain(x, l, u, lp))
                   : stan::math::sum(stan::math::lub_constrain(x, l, u));
    s.grad();
    expect_adj(x, dx[i], "dx", i);
    if (l_vec || !tail) {
      expect_adj(l, dlb[i], "dlb", i);
    }
    if (u_vec || !tail) {
      expect_adj(u, dub[i], "dub", i);
    }
    if (use_lp) {
      stan::math::set_zero_all_adjoints();
      lp.grad();
      expect_adj(x, dlogj[i], "dlogJ/dx", i);
      if (l_vec || !tail) {
        expect_adj(l, -0.2, "dlogJ/dlb", i);
      }
      if (u_vec || !tail) {
        expect_adj(u, 0.2, "dlogJ/dub", i);
      }
    }
    stan::math::recover_memory();
  }
}

template <typename V>
inline void vector_cases(size_t i) {
  using stan::math::var;
  for (bool tail : {false, true}) {
    vector_case<V, var, var>(i, tail);
    vector_case<V, V, var>(i, tail);
    vector_case<V, var, V>(i, tail);
    vector_case<V, V, V>(i, tail);
    vector_case<V, V, double>(i, tail);
    vector_case<V, double, V>(i, tail);
    vector_case<Eigen::VectorXd, V, V>(i, tail);
  }
  vector_case<V, var, double>(i, false);
  vector_case<V, double, var>(i, false);
  vector_case<Eigen::VectorXd, var, var>(i, false);
}
}  // namespace lub_constrain_rev_test

TEST_F(AgradRev, lub_constrain_scalar_derivatives_tails) {
  using lub_constrain_rev_test::scalar_case;
  using stan::math::var;
  for (size_t i = 0; i < lub_constrain_rev_test::xs.size(); ++i) {
    scalar_case<var, var, var>(i);
    scalar_case<var, var, double>(i);
    scalar_case<var, double, var>(i);
    scalar_case<var, double, double>(i);
    scalar_case<double, var, var>(i);
  }
}

TEST_F(AgradRev, lub_constrain_matrix_derivatives_tails) {
  using stan::math::var;
  for (size_t i = 0; i < lub_constrain_rev_test::xs.size(); ++i) {
    lub_constrain_rev_test::vector_cases<Eigen::Matrix<var, -1, 1>>(i);
  }
}

TEST_F(AgradRev, lub_constrain_var_matrix_derivatives_tails) {
  using stan::math::var_value;
  for (size_t i = 0; i < lub_constrain_rev_test::xs.size(); ++i) {
    lub_constrain_rev_test::vector_cases<var_value<Eigen::VectorXd>>(i);
  }
}

TEST_F(AgradRev, lub_constrain_std_vector_derivatives_tails) {
  using stan::math::var;
  for (size_t i = 0; i < lub_constrain_rev_test::xs.size(); ++i) {
    std::vector<var> x{lub_constrain_rev_test::xs[i]};
    var l = lub_constrain_rev_test::lb;
    var u = lub_constrain_rev_test::ub;
    std::vector<var> y = stan::math::lub_constrain(x, l, u);
    y[0].grad();
    lub_constrain_rev_test::expect_adj(x[0], lub_constrain_rev_test::dx[i],
                                       "dx", i);
    lub_constrain_rev_test::expect_adj(l, lub_constrain_rev_test::dlb[i], "dlb",
                                       i);
    stan::math::recover_memory();
  }
}
