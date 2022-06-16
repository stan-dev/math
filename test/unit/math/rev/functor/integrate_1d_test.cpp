#include <gtest/gtest.h>
#include <stan/math.hpp>
#include <test/unit/math/rev/fun/util.hpp>
#include <test/unit/util.hpp>
#include <vector>
#include <iostream>
#include <limits>
#include <sstream>

std::ostringstream *msgs = nullptr;

struct f1 {
  template <typename T1, typename T2, typename T3>
  inline stan::return_type_t<T1, T2, T3> operator()(
      const T1 &x, const T2 &xc, const std::vector<T3> &theta,
      const std::vector<double> &x_r, const std::vector<int> &x_i,
      std::ostream *msgs) const {
    return exp(x) + theta[0];
  }
};

struct f2 {
  template <typename T1, typename T2, typename T3>
  inline stan::return_type_t<T1, T2, T3> operator()(
      const T1 &x, const T2 &xc, const std::vector<T3> &theta,
      const std::vector<double> &x_r, const std::vector<int> &x_i,
      std::ostream *msgs) const {
    return exp(theta[0] * cos(2 * 3.141593 * x)) + theta[0];
  }
};

struct f3 {
  template <typename T1, typename T2, typename T3>
  inline stan::return_type_t<T1, T2, T3> operator()(
      const T1 &x, const T2 &xc, const std::vector<T3> &theta,
      const std::vector<double> &x_r, const std::vector<int> &x_i,
      std::ostream *msgs) const {
    return exp(x) + pow(theta[0], x_r[0]) + 2 * pow(theta[1], x_r[1])
           + 2 * theta[2];
  }
};

struct f4 {
  template <typename T1, typename T2, typename T3>
  inline stan::return_type_t<T1, T2, T3> operator()(
      const T1 &x, const T2 &xc, const std::vector<T3> &theta,
      const std::vector<double> &x_r, const std::vector<int> &x_i,
      std::ostream *msgs) const {
    return exp(-x) / sqrt(x);
  }
};

struct f5 {
  template <typename T1, typename T2, typename T3>
  inline stan::return_type_t<T1, T2, T3> operator()(
      const T1 &x, const T2 &xc, const std::vector<T3> &theta,
      const std::vector<double> &x_r, const std::vector<int> &x_i,
      std::ostream *msgs) const {
    return exp(-theta[0] * x) / sqrt(theta[1] * x);
  }
};

struct f6 {
  template <typename T1, typename T2, typename T3>
  inline stan::return_type_t<T1, T2, T3> operator()(
      const T1 &x, const T2 &xc, const std::vector<T3> &theta,
      const std::vector<double> &x_r, const std::vector<int> &x_i,
      std::ostream *msgs) const {
    return sqrt(x / (1 - theta[0] * x * x));
  }
};

struct f7 {
  template <typename T1, typename T2, typename T3>
  inline stan::return_type_t<T1, T2, T3> operator()(
      const T1 &x, const T2 &xc, const std::vector<T3> &theta,
      const std::vector<double> &x_r, const std::vector<int> &x_i,
      std::ostream *msgs) const {
    return exp(-theta[0] * x);
  }
};

struct f8 {
  template <typename T1, typename T2, typename T3>
  inline stan::return_type_t<T1, T2, T3> operator()(
      const T1 &x, const T2 &xc, const std::vector<T3> &theta,
      const std::vector<double> &x_r, const std::vector<int> &x_i,
      std::ostream *msgs) const {
    return exp(theta[0] * x);
  }
};

struct f10 {
  template <typename T1, typename T2, typename T3>
  inline stan::return_type_t<T1, T2, T3> operator()(
      const T1 &x, const T2 &xc, const std::vector<T3> &theta,
      const std::vector<double> &x_r, const std::vector<int> &x_i,
      std::ostream *msgs) const {
    return 1 / (1 + pow(x, x_i[0]) / x_r[0]);
  }
};

struct f11 {
  template <typename T1, typename T2, typename T3>
  inline stan::return_type_t<T1, T2, T3> operator()(
      const T1 &x, const T2 &xc, const std::vector<T3> &theta,
      const std::vector<double> &x_r, const std::vector<int> &x_i,
      std::ostream *msgs) const {
    return pow(x, theta[0] - 1.0)
           * pow((x > 0.5) ? xc : (1 - x), theta[1] - 1.0);
  }
};

struct f12 {
  template <typename T1, typename T2, typename T3>
  inline stan::return_type_t<T1, T2, T3> operator()(
      const T1 &x, const T2 &xc, const std::vector<T3> &theta,
      const std::vector<double> &x_r, const std::vector<int> &x_i,
      std::ostream *msgs) const {
    T3 mu = theta[0];
    T3 sigma = theta[1];
    return exp(-0.5 * stan::math::square((x - mu) / sigma))
           / (sigma * sqrt(2.0 * stan::math::pi()));
  }
};

struct f13 {
  template <typename T1, typename T2, typename T3>
  inline stan::return_type_t<T1, T2, T3> operator()(
      const T1 &x, const T2 &xc, const std::vector<T3> &theta,
      const std::vector<double> &x_r, const std::vector<int> &x_i,
      std::ostream *msgs) const {
    return x + theta[0] + theta[1];
  }
};

/*
 * Return the adjoint if the argument is a var
 *
 * @param v variable
 * @return adjoint of var
 */
double get_adjoint_if_var(stan::math::var v) { return v.adj(); }

/*
 * If the argument is not a var, return a NaN
 *
 * @param v variable
 * @return NaN
 */
double get_adjoint_if_var(double v) {
  return std::numeric_limits<double>::quiet_NaN();
}

/*
 * test_derivatives is a helper function to make it easy to test the
 * integrate_1d function.
 *
 * It takes in a callable function object, integration limits, parameters, real
 * and integer data. It integrates the provided function and compares the
 * computed integral and gradients against the provided integral value (val) and
 * gradients (grad, d_a, d_b).
 *
 * The first three template arguments specify how the left limit, right limit,
 * and parameters are tested. If the template arguments are vars, then the
 * gradients are checked against the provided gradients, otherwise the provided
 * gradients aren't used.
 *
 * The prototype for f is:
 *   struct f10 {
 *     template <typename T1, typename T2, typename T3>
 *     inline stan::return_type_t<T1, T2, T3> operator()(
 *         const T1& x, const T2& xc, const std::vector<T3>& theta, const
 * std::vector<double>& x_r, const std::vector<int>& x_i, std::ostream& msgs)
 * const {
 *     }
 *  };
 *
 * @tparam T_a Type of integration left limit
 * @tparam T_b Type of integration right limit
 * @tparam T_theta Type of parameters
 * @tparam F Type of f
 * @param f a functor with signature given above
 * @param a lower limit of integration
 * @param b upper limit of integration
 * @param param parameters to be passed to f (should be
 * std::vector<stan::math::var>)
 * @param x_r real data to be passed to f (should be std::vector<double>)
 * @param x_i integer data to be passed to f (should be std::vector<int>)
 * @param val correct value of integral
 * @param grad correct value of gradients (not used if T_theta is not var)
 * @param d_a correct value of gradient of integral with respect to left hand
 * limit (not used if T_a is not var)
 * @param d_b correct value of gradient of integral with respect to right hand
 * limit (not used if T_b is not var)
 */
template <typename T_a, typename T_b, typename T_theta, typename F>
void test_derivatives(const F &f, double a, double b,
                      std::vector<double> thetas,
                      const std::vector<double> &x_r,
                      const std::vector<int> &x_i, double val,
                      std::vector<double> grad, double d_a = 0.0,
                      double d_b = 0.0) {
  using stan::math::value_of;
  using stan::math::var;

  std::vector<double> tolerances = {1e-4, 1e-6, 1e-8};

  for (auto tolerance : tolerances) {
    // Reset autodiff stack
    stan::math::recover_memory();
    // Convert endpoints into given template types
    T_a a_(a);
    T_b b_(b);
    std::vector<T_theta> thetas_(thetas.size());
    for (size_t i = 0; i < thetas.size(); ++i)
      thetas_[i] = thetas[i];

    var integral = stan::math::integrate_1d(f, a_, b_, thetas_, x_r, x_i, msgs,
                                            tolerance);
    integral.grad();
    EXPECT_LE(std::abs(val - integral.val()), tolerance);
    if (stan::is_var<T_theta>::value) {
      for (size_t i = 0; i < grad.size(); ++i)
        EXPECT_LE(std::abs(grad[i] - get_adjoint_if_var(thetas_[i])),
                  tolerance);
    }
    if (stan::is_var<T_a>::value) {
      EXPECT_LE(std::abs(d_a - get_adjoint_if_var(a_)), tolerance);
    }
    if (stan::is_var<T_b>::value) {
      EXPECT_LE(std::abs(d_b - get_adjoint_if_var(b_)), tolerance);
    }
  }
}

TEST(StanMath_integrate_1d_rev, test_integer_arguments) {
  stan::math::var v;
  std::vector<stan::math::var> theta = {0.5};
  EXPECT_NO_THROW(
      v = stan::math::integrate_1d(f2{}, 0, 1, theta, {}, {}, msgs, 1e-6));
  EXPECT_NO_THROW(
      v = stan::math::integrate_1d(f2{}, 0.0, 1, theta, {}, {}, msgs, 1e-6));
  EXPECT_NO_THROW(
      v = stan::math::integrate_1d(f2{}, 0, 1.0, theta, {}, {}, msgs, 1e-6));
}

TEST(StanMath_integrate_1d_rev, TestDerivatives_easy) {
  // Easy integrals
  using stan::math::var;
  test_derivatives<double, double, var>(f1{}, 0.2, 0.7, {0.75}, {}, {},
                                        0.7923499493102901 + 0.5 * 0.75, {0.5});
  test_derivatives<var, var, var>(f2{}, 0.0, 1.0, {0.5}, {}, {},
                                  1.56348343527304, {1.25789445875152},
                                  -2.148721270700128, 2.14872127069993);
  test_derivatives<var, var, double>(f2{}, 0.0, 1.0, {0.5}, {}, {},
                                     1.56348343527304, {}, -2.148721270700128,
                                     2.14872127069993);
  test_derivatives<double, double, var>(f1{}, 0.0, 0.0, {0.75}, {}, {}, 0.0,
                                        {0.0});
  test_derivatives<double, double, var>(f2{}, 1.0, 1.0, {0.5}, {}, {}, 0.0,
                                        {0.0});
}

TEST(StanMath_integrate_1d_rev, TestDerivatives_zero_crossing) {
  // Zero crossing integral + test x_r + vars at endpoints
  using stan::math::var;
  test_derivatives<var, var, var>(f3{}, -1.0, 1.0, {0.5, 1.75, 3.9}, {2.5, 3.0},
                                  {},
                                  2.350402387287579 + 2.0 * pow(0.5, 2.5)
                                      + 4.0 * pow(1.75, 3.0) + 4.0 * 3.9,
                                  {5 * pow(0.5, 1.5), 12 * 1.75 * 1.75, 4.0},
                                  -19.06340613646808, 21.41380852375568);
}

TEST(StanMath_integrate_1d_rev, TestDerivatives_var_right_endpoint_var_params) {
  // Zero crossing integral + test x_r + vars at right endpoint
  using stan::math::var;
  test_derivatives<double, var, var>(
      f3{}, -1.0, 1.0, {0.5, 1.75, 3.9}, {2.5, 3.0}, {},
      2.350402387287579 + 2.0 * pow(0.5, 2.5) + 4.0 * pow(1.75, 3.0)
          + 4.0 * 3.9,
      {5 * pow(0.5, 1.5), 12 * 1.75 * 1.75, 4.0}, 0.0, 21.41380852375568);
}

TEST(StanMath_integrate_1d_rev, TestDerivatives_var_left_endpoint_var_params) {
  // Zero crossing integral + test x_r + var at left endpoint
  using stan::math::var;
  test_derivatives<var, double, var>(
      f3{}, -1.0, 1.0, {0.5, 1.75, 3.9}, {2.5, 3.0}, {},
      2.350402387287579 + 2.0 * pow(0.5, 2.5) + 4.0 * pow(1.75, 3.0)
          + 4.0 * 3.9,
      {5 * pow(0.5, 1.5), 12 * 1.75 * 1.75, 4.0}, -19.06340613646808, 0.0);
}

TEST(StanMath_integrate_1d_rev, TestDerivatives_no_param_vars) {
  // No param vars
  using stan::math::var;
  test_derivatives<var, var, double>(f3{}, -1.0, 1.0, {0.5, 1.75, 3.9},
                                     {2.5, 3.0}, {},
                                     2.350402387287579 + 2.0 * pow(0.5, 2.5)
                                         + 4.0 * pow(1.75, 3.0) + 4.0 * 3.9,
                                     {}, -19.06340613646808, 21.41380852375568);
}

TEST(StanMath_integrate_1d_rev, TestDerivatives_left_limit_var) {
  // No param vars, only left limit var
  using stan::math::var;
  test_derivatives<var, double, double>(f3{}, -1.0, 1.0, {0.5, 1.75, 3.9},
                                        {2.5, 3.0}, {},
                                        2.350402387287579 + 2.0 * pow(0.5, 2.5)
                                            + 4.0 * pow(1.75, 3.0) + 4.0 * 3.9,
                                        {}, -19.06340613646808, 0.0);
}

TEST(StanMath_integrate_1d_rev, TestDerivatives_right_limit_var) {
  // No param vars, only right limit var
  using stan::math::var;
  test_derivatives<double, var, double>(f3{}, -1.0, 1.0, {0.5, 1.75, 3.9},
                                        {2.5, 3.0}, {},
                                        2.350402387287579 + 2.0 * pow(0.5, 2.5)
                                            + 4.0 * pow(1.75, 3.0) + 4.0 * 3.9,
                                        {}, 0.0, 21.41380852375568);
}

TEST(StanMath_integrate_1d_rev, TestDerivatives_tricky1) {
  // Tricky integral from Boost docs + limit at infinity + no gradients
  using stan::math::var;
  test_derivatives<double, double, var>(f4{}, 0.0,
                                        std::numeric_limits<double>::infinity(),
                                        {}, {}, {}, 1.772453850905516, {});
}

TEST(StanMath_integrate_1d_rev, TestDerivatives_tricky2) {
  // Tricky integral from Boost docs + limit at infinity with gradients
  using stan::math::var;
  test_derivatives<double, double, var>(
      f5{}, 0.0, std::numeric_limits<double>::infinity(), {0.5, 3.0}, {}, {},
      1.772453850905516 / sqrt(0.5 * 3.0),
      {-1.772453850905516 * 3.0 / (2 * pow(0.5 * 3.0, 1.5)),
       -1.772453850905516 * 0.5 / (2 * pow(0.5 * 3.0, 1.5))});
}

TEST(StanMath_integrate_1d_rev, TestDerivatives_tricky3) {
  // Tricky integral from Boost docs
  using stan::math::var;
  test_derivatives<double, double, var>(
      f6{}, 0.0, 1.0, {0.75}, {}, {}, 0.851926727945904, {0.4814066053874294});
}

TEST(StanMath_integrate_1d_rev, TestDerivatives_zero_crossing2) {
  // Zero crossing integral + limit at infinity + var at left limit
  using stan::math::var;
  test_derivatives<var, double, var>(
      f7{}, -5.0, std::numeric_limits<double>::infinity(), {1.5}, {}, {},
      1205.361609637375, {5223.23364176196}, -1808.042414456063,
      std::numeric_limits<double>::quiet_NaN());
}

TEST(StanMath_integrate_1d_rev, TestDerivatives_zero_crossing3) {
  // Zero crossing integral + limit at negative infinity + var at right limit
  using stan::math::var;
  test_derivatives<double, var, var>(
      f8{}, -std::numeric_limits<double>::infinity(), 5.0, {1.5}, {}, {},
      1205.361609637375, {5223.23364176196},
      std::numeric_limits<double>::quiet_NaN(), 1808.042414456063);
}

TEST(StanMath_integrate_1d_rev, TestDerivatives_indefinite) {
  // Both limits at infinity + test x_r/x_i + no gradients
  using stan::math::var;
  test_derivatives<double, double, var>(
      f10{}, -std::numeric_limits<double>::infinity(),
      std::numeric_limits<double>::infinity(), {}, {1.7}, {4},
      2.536571480364399, {});
}

TEST(StanMath_integrate_1d_rev, TestDerivatives_endpoint_precision) {
  // Various integrals of beta function
  using stan::math::var;
  test_derivatives<double, double, var>(f11{}, 0.0, 1.0, {0.1, 0.1}, {}, {},
                                        19.71463948905016,
                                        {-101.229055967892, -101.229055967892});
  test_derivatives<double, double, var>(
      f11{}, 0.0, 1.0, {0.5, 0.51}, {}, {}, 3.098843783331868,
      {-4.346514423368675, -4.196150770134913});
  test_derivatives<double, double, var>(
      f11{}, 0.0, 1.0, {0.51, 0.5}, {}, {}, 3.098843783331868,
      {-4.196150770134913, -4.346514423368675});
  test_derivatives<double, double, var>(
      f11{}, 0.0, 1.0, {5.0, 3.0}, {}, {}, 0.00952380952380952,
      {-0.004852607709750566, -0.01040816326530613});
  test_derivatives<double, double, var>(
      f11{}, 0.0, 1.0, {3.0, 5.0}, {}, {}, 0.00952380952380952,
      {-0.01040816326530613, -0.004852607709750566});
}

TEST(StanMath_integrate_1d_rev, TestDerivatives_gaussian) {
  // Check Gaussian integrates to 1.0 always
  using stan::math::var;
  test_derivatives<double, double, var>(
      f12{}, -std::numeric_limits<double>::infinity(),
      std::numeric_limits<double>::infinity(), {5.7, 1}, {}, {}, 1.0,
      {0.0, 0.0});
}

TEST(StanMath_integrate_1d_rev, TestDerivativesSameVarAtEndpointAndInParams) {
  using stan::math::var;

  var a = 2.0;
  var b = 4.0;
  std::vector<var> thetas = {a, b};

  var integral
      = stan::math::integrate_1d(f13{}, a, b, thetas, {}, {}, msgs, 1e-8);
  integral.grad();

  EXPECT_LT(std::abs(18.0 - integral.val()), 1e-8);
  EXPECT_LT(std::abs(-6.0 - a.adj()), 1e-8);
  EXPECT_LT(std::abs(12.0 - b.adj()), 1e-8);
}

TEST(StanMath_integrate_1d_rev, TestBeta) {
  using stan::math::exp;
  using stan::math::integrate_1d;
  using stan::math::var;

  var alpha = 9.0 / 5;
  var beta = 13.0 / 7;
  std::vector<stan::math::var> theta = {alpha, beta};
  auto pdf = [](auto x, auto xc, auto theta, auto x_r, auto x_i,
                std::ostream *msgs) {
    return exp(stan::math::beta_lpdf(x, theta[0], theta[1]));
  };
  var I = integrate_1d(pdf, 0.0, 1.0, theta, {}, {}, msgs, 1e-8);
  EXPECT_FLOAT_EQ(1, I.val());

  std::vector<stan::math::var> x{alpha, beta};
  std::vector<double> g;
  I.grad(x, g);
  EXPECT_FLOAT_EQ(1, 1 + g[0]);
  EXPECT_FLOAT_EQ(1, 1 + g[1]);
}

TEST(StanMath_integrate_1d_rev, TestCauchy) {
  using stan::math::exp;
  using stan::math::integrate_1d;
  using stan::math::var;

  var mu = 9.0 / 5;
  var sigma = 13.0 / 7;
  std::vector<stan::math::var> theta = {mu, sigma};
  double b = std::numeric_limits<double>::infinity();
  double a = -b;
  auto pdf = [](auto x, auto xc, auto theta, auto x_r, auto x_i,
                std::ostream *msgs) {
    return exp(stan::math::cauchy_lpdf(x, theta[0], theta[1]));
  };
  var I = integrate_1d(pdf, a, b, theta, {}, {}, msgs, 1e-8);
  EXPECT_FLOAT_EQ(1, I.val());

  std::vector<stan::math::var> x{mu, sigma};
  std::vector<double> g;
  I.grad(x, g);
  EXPECT_FLOAT_EQ(1, 1 + g[0]);
  EXPECT_FLOAT_EQ(1, 1 + g[1]);
}

TEST(StanMath_integrate_1d_rev, TestChiSquare) {
  using stan::math::exp;
  using stan::math::integrate_1d;
  using stan::math::var;

  var nu = 9.0 / 5;
  std::vector<stan::math::var> theta = {nu};
  double b = std::numeric_limits<double>::infinity();
  double a = 0;
  auto pdf = [](auto x, auto xc, auto theta, auto x_r, auto x_i,
                std::ostream *msgs) {
    return exp(stan::math::chi_square_lpdf(x, theta[0]));
  };
  var I = integrate_1d(pdf, a, b, theta, {}, {}, msgs, 1e-8);
  EXPECT_FLOAT_EQ(1, I.val());

  std::vector<stan::math::var> x = {nu};
  std::vector<double> g;
  I.grad(x, g);
  EXPECT_FLOAT_EQ(1, 1 + g[0]);
}

TEST(StanMath_integrate_1d_rev, TestDoubleExponential) {
  using stan::math::exp;
  using stan::math::integrate_1d;
  using stan::math::var;

  var mu = 9.0 / 5;
  var sigma = 13.0 / 7;
  std::vector<stan::math::var> theta = {mu, sigma};
  double a = -std::numeric_limits<double>::infinity();
  double b = mu.val();
  auto pdf = [](auto x, auto xc, auto theta, auto x_r, auto x_i,
                std::ostream *msgs) {
    return exp(stan::math::double_exponential_lpdf(x, theta[0], theta[1]));
  };
  // requires two subintervals to achieve numerical accuracy
  var I = integrate_1d(pdf, a, b, theta, {}, {}, msgs, 1e-8)
          + integrate_1d(pdf, b, -a, theta, {}, {}, msgs, 1e-8);
  EXPECT_FLOAT_EQ(1, I.val());

  std::vector<stan::math::var> x{mu, sigma};
  std::vector<double> g;
  I.grad(x, g);
  EXPECT_FLOAT_EQ(1, 1 + g[0]);
  EXPECT_FLOAT_EQ(1, 1 + g[1]);
}

TEST(StanMath_integrate_1d_rev, TestExponential) {
  using stan::math::exp;
  using stan::math::integrate_1d;
  using stan::math::var;

  var beta = 9.0 / 5;
  std::vector<stan::math::var> theta = {beta};
  double b = std::numeric_limits<double>::infinity();
  double a = 0;
  auto pdf = [](auto x, auto xc, auto theta, auto x_r, auto x_i,
                std::ostream *msgs) {
    return exp(stan::math::exponential_lpdf(x, theta[0]));
  };
  var I = integrate_1d(pdf, a, b, theta, {}, {}, msgs, 1e-8);
  EXPECT_FLOAT_EQ(1, I.val());

  std::vector<stan::math::var> x = {beta};
  std::vector<double> g;
  I.grad(x, g);
  EXPECT_FLOAT_EQ(1, 1 + g[0]);
}

TEST(StanMath_integrate_1d_rev, TestFrechet) {
  using stan::math::exp;
  using stan::math::integrate_1d;
  using stan::math::var;

  var alpha = 9.0 / 5;
  var sigma = 13.0 / 7;
  std::vector<stan::math::var> theta = {alpha, sigma};
  double b = std::numeric_limits<double>::infinity();
  double a = 0;
  auto pdf = [](auto x, auto xc, auto theta, auto x_r, auto x_i,
                std::ostream *msgs) {
    return exp(stan::math::frechet_lpdf(x, theta[0], theta[1]));
  };
  var I = integrate_1d(pdf, a, b, theta, {}, {}, msgs, 1e-8);
  EXPECT_FLOAT_EQ(1, I.val());

  std::vector<stan::math::var> x = {alpha, sigma};
  std::vector<double> g;
  I.grad(x, g);
  EXPECT_FLOAT_EQ(1, 1 + g[0]);
  EXPECT_FLOAT_EQ(1, 1 + g[1]);
}

TEST(StanMath_integrate_1d_rev, TestGamma) {
  using stan::math::exp;
  using stan::math::integrate_1d;
  using stan::math::var;

  var alpha = 9.0 / 5;
  var beta = 13.0 / 7;
  std::vector<stan::math::var> theta = {alpha, beta};
  double b = std::numeric_limits<double>::infinity();
  double a = 0;
  auto pdf = [](auto x, auto xc, auto theta, auto x_r, auto x_i,
                std::ostream *msgs) {
    return exp(stan::math::gamma_lpdf(x, theta[0], theta[1]));
  };
  var I = integrate_1d(pdf, a, b, theta, {}, {}, msgs, 1e-8);
  EXPECT_FLOAT_EQ(1, I.val());

  std::vector<stan::math::var> x{alpha, beta};
  std::vector<double> g;
  I.grad(x, g);
  EXPECT_FLOAT_EQ(1, 1 + g[0]);
  EXPECT_FLOAT_EQ(1, 1 + g[1]);
}

TEST(StanMath_integrate_1d_rev, TestGumbel) {
  using stan::math::exp;
  using stan::math::integrate_1d;
  using stan::math::var;

  var mu = 9.0 / 5;
  var beta = 13.0 / 7;
  std::vector<stan::math::var> theta = {mu, beta};
  double b = std::numeric_limits<double>::infinity();
  double a = -b;
  auto pdf = [](auto x, auto xc, auto theta, auto x_r, auto x_i,
                std::ostream *msgs) {
    return exp(stan::math::gumbel_lpdf(x, theta[0], theta[1]));
  };
  var I = integrate_1d(pdf, a, b, theta, {}, {}, msgs, 1e-8);
  EXPECT_FLOAT_EQ(1, I.val());

  std::vector<stan::math::var> x{mu, beta};
  std::vector<double> g;
  I.grad(x, g);
  EXPECT_FLOAT_EQ(1, 1 + g[0]);
  EXPECT_FLOAT_EQ(1, 1 + g[1]);
}

TEST(StanMath_integrate_1d_rev, TestInvChiSquared) {
  using stan::math::exp;
  using stan::math::integrate_1d;
  using stan::math::var;

  var nu = 9.0 / 5;
  std::vector<stan::math::var> theta = {nu};
  double b = std::numeric_limits<double>::infinity();
  double a = 0;
  auto pdf = [](auto x, auto xc, auto theta, auto x_r, auto x_i,
                std::ostream *msgs) {
    return exp(stan::math::inv_chi_square_lpdf(x, theta[0]));
  };
  var I = integrate_1d(pdf, a, b, theta, {}, {}, msgs, 1e-8);
  EXPECT_FLOAT_EQ(1, I.val());

  std::vector<stan::math::var> x = {nu};
  std::vector<double> g;
  I.grad(x, g);
  EXPECT_FLOAT_EQ(1, 1 + g[0]);
}

TEST(StanMath_integrate_1d_rev, TestLogistic) {
  using stan::math::exp;
  using stan::math::integrate_1d;
  using stan::math::var;

  var mu = 9.0 / 5;
  var sigma = 13.0 / 7;
  std::vector<stan::math::var> theta = {mu, sigma};
  double b = std::numeric_limits<double>::infinity();
  double a = -b;
  auto pdf = [](auto x, auto xc, auto theta, auto x_r, auto x_i,
                std::ostream *msgs) {
    return exp(stan::math::logistic_lpdf(x, theta[0], theta[1]));
  };
  var I = integrate_1d(pdf, a, b, theta, {}, {}, msgs, 1e-8);
  EXPECT_FLOAT_EQ(1, I.val());

  std::vector<stan::math::var> x{mu, sigma};
  std::vector<double> g;
  I.grad(x, g);
  EXPECT_FLOAT_EQ(1, 1 + g[0]);
  EXPECT_FLOAT_EQ(1, 1 + g[1]);
}

TEST(StanMath_integrate_1d_rev, TestLogNormal) {
  using stan::math::exp;
  using stan::math::integrate_1d;
  using stan::math::var;

  var mu = 9.0 / 5;
  var sigma = 13.0 / 7;
  std::vector<stan::math::var> theta = {mu, sigma};
  double b = std::numeric_limits<double>::infinity();
  double a = 0;
  auto pdf = [](auto x, auto xc, auto theta, auto x_r, auto x_i,
                std::ostream *msgs) {
    return exp(stan::math::lognormal_lpdf(x, theta[0], theta[1]));
  };
  var I = integrate_1d(pdf, a, b, theta, {}, {}, msgs, 1e-8);
  EXPECT_FLOAT_EQ(1, I.val());

  std::vector<stan::math::var> x{mu, sigma};
  std::vector<double> g;
  I.grad(x, g);
  EXPECT_FLOAT_EQ(1, 1 + g[0]);
  EXPECT_FLOAT_EQ(1, 1 + g[1]);
}

TEST(StanMath_integrate_1d_rev, TestNormal) {
  using stan::math::exp;
  using stan::math::integrate_1d;
  using stan::math::var;

  var mu = 9.0 / 5;
  var sigma = 13.0 / 7;
  std::vector<stan::math::var> theta = {mu, sigma};
  double b = std::numeric_limits<double>::infinity();
  double a = -b;
  auto pdf = [](auto x, auto xc, auto theta, auto x_r, auto x_i,
                std::ostream *msgs) {
    return exp(stan::math::normal_lpdf(x, theta[0], theta[1]));
  };
  var I = integrate_1d(pdf, a, b, theta, {}, {}, msgs, 1e-8);
  EXPECT_FLOAT_EQ(1, I.val());

  std::vector<stan::math::var> x{mu, sigma};
  std::vector<double> g;
  I.grad(x, g);
  EXPECT_FLOAT_EQ(1, 1 + g[0]);
  EXPECT_FLOAT_EQ(1, 1 + g[1]);
}

TEST(StanMath_integrate_1d_rev, TestPareto) {
  using stan::math::exp;
  using stan::math::integrate_1d;
  using stan::math::var;

  var m = 9.0 / 5;
  var alpha = 13.0 / 7;
  std::vector<stan::math::var> theta = {m, alpha};
  double b = std::numeric_limits<double>::infinity();
  var a = m;
  auto pdf = [](auto x, auto xc, auto theta, auto x_r, auto x_i,
                std::ostream *msgs) {
    return exp(stan::math::pareto_lpdf(x, theta[0], theta[1]));
  };
  var I = integrate_1d(pdf, a, b, theta, {}, {}, msgs, 1e-8);
  EXPECT_FLOAT_EQ(1, I.val());

  std::vector<stan::math::var> x{m, alpha};
  std::vector<double> g;
  I.grad(x, g);
  EXPECT_FLOAT_EQ(1, 1 + g[0]);
  EXPECT_FLOAT_EQ(1, 1 + g[1]);
}

TEST(StanMath_integrate_1d_rev, TestPareto2) {
  using stan::math::exp;
  using stan::math::integrate_1d;
  using stan::math::var;

  var mu = 9.0 / 5;
  var lambda = 13.0 / 7;
  var alpha = 11.0 / 3;
  std::vector<stan::math::var> theta = {mu, lambda, alpha};
  double b = std::numeric_limits<double>::infinity();
  var a = mu;
  auto pdf = [](auto x, auto xc, auto theta, auto x_r, auto x_i,
                std::ostream *msgs) {
    return exp(stan::math::pareto_type_2_lpdf(x, theta[0], theta[1], theta[2]));
  };
  var I = integrate_1d(pdf, a, b, theta, {}, {}, msgs, 1e-8);
  EXPECT_FLOAT_EQ(1, I.val());

  std::vector<stan::math::var> x{mu, lambda, alpha};
  std::vector<double> g;
  I.grad(x, g);
  EXPECT_FLOAT_EQ(1, 1 + g[0]);
  EXPECT_FLOAT_EQ(1, 1 + g[1]);
  EXPECT_FLOAT_EQ(1, 1 + g[2]);
}

TEST(StanMath_integrate_1d_rev, TestRayleigh) {
  using stan::math::exp;
  using stan::math::integrate_1d;
  using stan::math::var;

  var sigma = 13.0 / 7;
  std::vector<stan::math::var> theta = {sigma};
  double b = std::numeric_limits<double>::infinity();
  double a = 0;
  auto pdf = [](auto x, auto xc, auto theta, auto x_r, auto x_i,
                std::ostream *msgs) {
    return exp(stan::math::rayleigh_lpdf(x, theta[0]));
  };
  var I = integrate_1d(pdf, a, b, theta, {}, {}, msgs, 1e-8);
  EXPECT_FLOAT_EQ(1, I.val());

  std::vector<stan::math::var> x{sigma};
  std::vector<double> g;
  I.grad(x, g);
  EXPECT_FLOAT_EQ(1, 1 + g[0]);
}

TEST(StanMath_integrate_1d_rev, TestScaledInvChiSquare) {
  using stan::math::exp;
  using stan::math::integrate_1d;
  using stan::math::var;

  var nu = 9.0 / 5;
  var s = 13.0 / 7;
  std::vector<stan::math::var> theta = {nu, s};
  double b = std::numeric_limits<double>::infinity();
  double a = 0;
  auto pdf = [](auto x, auto xc, auto theta, auto x_r, auto x_i,
                std::ostream *msgs) {
    return exp(stan::math::scaled_inv_chi_square_lpdf(x, theta[0], theta[1]));
  };
  var I = integrate_1d(pdf, a, b, theta, {}, {}, msgs, 1e-8);
  EXPECT_FLOAT_EQ(1, I.val());

  std::vector<stan::math::var> x{nu, s};
  std::vector<double> g;
  I.grad(x, g);
  EXPECT_FLOAT_EQ(1, 1 + g[0]);
  EXPECT_FLOAT_EQ(1, 1 + g[1]);
}

TEST(StanMath_integrate_1d_rev, TestStudentT) {
  using stan::math::exp;
  using stan::math::integrate_1d;
  using stan::math::var;

  var nu = 11.0 / 3;
  var mu = 9.0 / 5;
  var sigma = 13.0 / 7;
  std::vector<stan::math::var> theta = {nu, mu, sigma};
  double b = std::numeric_limits<double>::infinity();
  double a = -b;
  auto pdf = [](auto x, auto xc, auto theta, auto x_r, auto x_i,
                std::ostream *msgs) {
    return exp(stan::math::student_t_lpdf(x, theta[0], theta[1], theta[2]));
  };
  var I = integrate_1d(pdf, a, b, theta, {}, {}, msgs, 1e-8);
  EXPECT_FLOAT_EQ(1, I.val());

  std::vector<stan::math::var> x{nu, mu, sigma};
  std::vector<double> g;
  I.grad(x, g);
  EXPECT_FLOAT_EQ(1, 1 + g[0]);
  EXPECT_FLOAT_EQ(1, 1 + g[1]);
  EXPECT_FLOAT_EQ(1, 1 + g[2]);
}

TEST(StanMath_integrate_1d_rev, TestUniform) {
  using stan::math::exp;
  using stan::math::integrate_1d;
  using stan::math::var;

  var a = 9.0 / 5;
  var b = 13.0 / 7;
  std::vector<stan::math::var> theta = {a, b};
  auto pdf = [](auto x, auto xc, auto theta, auto x_r, auto x_i,
                std::ostream *msgs) {
    return exp(stan::math::uniform_lpdf(x, theta[0], theta[1]));
  };
  var I = integrate_1d(pdf, a, b, theta, {}, {}, msgs, 1e-8);
  EXPECT_FLOAT_EQ(1, I.val());

  std::vector<stan::math::var> x{a, b};
  std::vector<double> g;
  I.grad(x, g);
  EXPECT_FLOAT_EQ(1, 1 + g[0]);
  EXPECT_FLOAT_EQ(1, 1 + g[1]);
}

TEST(StanMath_integrate_1d_rev, TestVonMises) {
  using stan::math::exp;
  using stan::math::integrate_1d;
  using stan::math::var;

  var mu = 9.0 / 5;
  var kappa = 13.0 / 7;
  std::vector<stan::math::var> theta = {mu, kappa};
  double b = stan::math::pi() * 2;
  double a = 0;
  auto pdf = [](auto x, auto xc, auto theta, auto x_r, auto x_i,
                std::ostream *msgs) {
    return exp(stan::math::von_mises_lpdf(x, theta[0], theta[1]));
  };
  var I = integrate_1d(pdf, a, b, theta, {}, {}, msgs, 1e-8);
  EXPECT_FLOAT_EQ(1, I.val());

  std::vector<stan::math::var> x{mu, kappa};
  std::vector<double> g;
  I.grad(x, g);
  EXPECT_FLOAT_EQ(1, 1 + g[0]);
  EXPECT_FLOAT_EQ(1, 1 + g[1]);
}

TEST(StanMath_integrate_1d_rev, TestWeibull) {
  using stan::math::exp;
  using stan::math::integrate_1d;
  using stan::math::var;

  var alpha = 9.0 / 5;
  var sigma = 13.0 / 7;
  std::vector<stan::math::var> theta = {alpha, sigma};
  double b = std::numeric_limits<double>::infinity();
  double a = 0;
  auto pdf = [](auto x, auto xc, auto theta, auto x_r, auto x_i,
                std::ostream *msgs) {
    return exp(stan::math::weibull_lpdf(x, theta[0], theta[1]));
  };
  var I = integrate_1d(pdf, a, b, theta, {}, {}, msgs, 1e-8);
  EXPECT_FLOAT_EQ(1, I.val());

  std::vector<stan::math::var> x{alpha, sigma};
  std::vector<double> g;
  I.grad(x, g);
  EXPECT_FLOAT_EQ(1, 1 + g[0]);
  EXPECT_FLOAT_EQ(1, 1 + g[1]);
}
