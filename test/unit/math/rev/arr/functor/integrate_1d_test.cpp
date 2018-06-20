#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <iostream>
#include <sstream>
#include <vector>

std::ostringstream msgs;

struct f1 {
  template <typename T1, typename T2>
  inline typename stan::return_type<T1, T2>::type operator()(
      const T1& x, const std::vector<T2>& theta, const std::vector<double>& x_r, const std::vector<int>& x_i, std::ostream& msgs) const {
    return exp(x) + theta[0];
  }
};

struct f2 {
  template <typename T1, typename T2>
  inline typename stan::return_type<T1, T2>::type operator()(
      const T1& x, const std::vector<T2>& theta, const std::vector<double>& x_r, const std::vector<int>& x_i, std::ostream& msgs) const {
    return exp(theta[0] * cos(2 * 3.141593 * x)) + theta[0];
  }
};

struct f3 {
  template <typename T1, typename T2>
  inline typename stan::return_type<T1, T2>::type operator()(
      const T1& x, const std::vector<T2>& theta, const std::vector<double>& x_r, const std::vector<int>& x_i, std::ostream& msgs) const {
    return exp(x) + pow(theta[0], x_r[0]) + 2 * pow(theta[1], x_r[1]) + 2 * theta[2];
  }
};

struct f6 {
  template <typename T1, typename T2>
  inline typename stan::return_type<T1, T2>::type operator()(
      const T1& x, const std::vector<T2>& theta, const std::vector<double>& x_r, const std::vector<int>& x_i, std::ostream& msgs) const {
    return sqrt(x / (1 - theta[0] * x * x));
  }
};

/*
 * test_derivatives is a helper function to make it easy to test the integrate_1d
 * function.
 *
 * It takes in a callable function object, integration limits, parameters, real and
 * integer data. It integrates the provided function and compares the computed integral
 * and gradients against the provided integral value (val) and gradients (grad).
 *
 * The prototype for f is:
 *   struct f10 {
 *     template <typename T1, typename T2>
 *     inline typename stan::return_type<T1, T2>::type operator()(
 *         const T1& x, const std::vector<T2>& theta, const std::vector<double>& x_r, const std::vector<int>& x_i, std::ostream& msgs) const {
 *     }
 *  };
 *
 * @tparam F Type of f
 * @param f a functor with signature given above
 * @param a lower limit of integration
 * @param b upper limit of integration
 * @param param parameters to be passed to f (should be std::vector<stan::math::var>)
 * @param x_r real data to be passed to f (should be std::vector<double>)
 * @param x_i integer data to be passed to f (should be std::vector<int>)
 * @param val correct value of integral
 * @param grad correct value of gradients
 */
template <typename F>
void test_derivatives(const F& f, double a, double b, std::vector<double> thetas, const std::vector<double>& x_r, const std::vector<int>& x_i, double val, std::vector<double> grad) {
  using stan::math::var;

  std::vector<double> tolerances = { 1e-4, 1e-6, 1e-8 };

  for(auto tolerance : tolerances) {
    // Reset autodiff stack
    stan::math::recover_memory();
    // Convert parameters into vars for actual calls
    std::vector<var> theta_vars(thetas.size());
    for(size_t i = 0; i < thetas.size(); ++i) {
      theta_vars[i] = thetas[i];
    }
    
    var integral = stan::math::integrate_1d(f, a, b, theta_vars, x_r, x_i, msgs, tolerance);
    std::vector<double> computed_grad;
    integral.grad(theta_vars, computed_grad);
    EXPECT_LE(std::abs(val - integral.val()), tolerance);
    for(size_t i = 0; i < grad.size(); ++i) {
      EXPECT_LE(std::abs(grad[i] - computed_grad[i]), tolerance);
    }
  }
}

TEST(StanMath_integrate_1d, TestDerivatives) {
  // Easy integrals
  test_derivatives(f1{}, 0.2, 0.7, {0.75}, {}, {}, 0.7923499493102901 + 0.5 * 0.75, {0.5});
  test_derivatives(f2{}, 0.0, 1.0, {0.5}, {}, {}, 1.56348343527304, {1.25789445875152}) ;
  test_derivatives(f1{}, 0.0, 0.0, {0.75}, {}, {}, 0.0, {0.0});
  test_derivatives(f2{}, 1.0, 1.0, {0.5}, {}, {}, 0.0, {0.0}) ;
  // Zero crossing integral + test x_r
  test_derivatives(f3{}, -1.0, 1.0, {0.5, 1.75, 3.9}, {2.5, 3.0}, {}, 2.350402387287579 + 2.0 * pow(0.5, 2.5) + 4.0 * pow(1.75, 3.0) + 4.0 * 3.9, { 5 * pow(0.5, 1.5), 12 * 1.75 * 1.75, 4.0});
  // Tricky integral from Boost 1d integratioin docs
  test_derivatives(f6{}, 0.0, 1.0, {0.75}, {}, {}, 0.851926727945904, {0.4814066053874294});
}
