#include <stan/math/prim/functor.hpp>
#include <stan/math/prim/fun.hpp>

#include <gtest/gtest.h>
#include <vector>

namespace hcubature_test {

struct my_params {
  long double x;
  double y;
};

template <typename T_x, typename T_p>
stan::return_type_t<T_x, T_p> f1(const T_x& x, const T_p& p) {
  using T_x_ref = stan::ref_type_t<T_x>;
  T_x_ref x_ref = x;
  stan::scalar_seq_view<T_x_ref> x_vec(x_ref);
  return cos(x_vec[0]);
}

template <typename T_x, typename T_p>
stan::return_type_t<T_x, T_p> f2(const T_x& x, const T_p& p) {
  using T_x_ref = stan::ref_type_t<T_x>;
  T_x_ref x_ref = x;
  stan::scalar_seq_view<T_x_ref> x_vec(x_ref);
  return cos(x_vec[0]) * cos(x_vec[1]);
}

template <typename T_x, typename T_p>
stan::return_type_t<T_x, T_p> f3(const T_x& x, const T_p& p) {
  using T_x_ref = stan::ref_type_t<T_x>;
  T_x_ref x_ref = x;
  stan::scalar_seq_view<T_x_ref> x_vec(x_ref);
  my_params* pars = static_cast<my_params*>(p);
  long double radius = (pars->x);
  double result;
  if (std::pow(x_vec[0], 2) + std::pow(x_vec[1], 2) < std::pow(radius, 2)) {
    result = 1;
  } else {
    result = 0;
  }
  return result;
}

template <typename T_x, typename T_p>
stan::return_type_t<T_x, T_p> f4(const T_x& x, const T_p& p) {
  using T_x_ref = stan::ref_type_t<T_x>;
  using namespace stan::math;
  T_x_ref x_ref = x;
  stan::scalar_seq_view<T_x_ref> x_vec(x_ref);
  my_params* pars = static_cast<my_params*>(p);
  double sigma = (pars->y);
  double numerator = std::pow(x_vec[0] - 0.5, 2) + std::pow(x_vec[1] - 0.5, 2);
  return std::pow(TWO_OVER_SQRT_PI / (2.0 * sigma), 2)
         * exp(-numerator / std::pow(sigma, 2));
}

template <typename T_x, typename T_p>
stan::return_type_t<T_x, T_p> f5(const T_x& x, const T_p& p) {
  using T_x_ref = stan::ref_type_t<T_x>;
  using namespace stan::math;
  T_x_ref x_ref = x;
  stan::scalar_seq_view<T_x_ref> x_vec(x_ref);
  double val = std::pow((1 - x_vec[0]) / x_vec[0], 2)
               + std::pow((1 - x_vec[1]) / x_vec[1], 2)
               + std::pow((1 - x_vec[2]) / x_vec[2], 2);
  double scale = TWO_OVER_SQRT_PI / std::pow(x_vec[0], 2) * TWO_OVER_SQRT_PI
                 / std::pow(x_vec[1], 2) * TWO_OVER_SQRT_PI
                 / std::pow(x_vec[2], 2);
  return exp(-val) * scale;
}

template <typename T_x, typename T_p>
stan::return_type_t<T_x, T_p> f6(const T_x& x, const T_p& p) {
  using T_x_ref = stan::ref_type_t<T_x>;
  T_x_ref x_ref = x;
  stan::scalar_seq_view<T_x_ref> x_vec(x_ref);
  return 2 * x_vec[0] * 2 * x_vec[1] * 2 * x_vec[2];
}

template <typename T_x, typename T_p>
stan::return_type_t<T_x, T_p> f7(const T_x& x, const T_p& p) {
  using T_x_ref = stan::ref_type_t<T_x>;
  T_x_ref x_ref = x;
  stan::scalar_seq_view<T_x_ref> x_vec(x_ref);
  my_params* pars = static_cast<my_params*>(p);
  double a = (pars->y);
  double result = (a / (a + 1) * std::pow((a + 1) / (a + x_vec[0]), 2))
                  * (a / (a + 1) * std::pow((a + 1) / (a + x_vec[1]), 2))
                  * (a / (a + 1) * std::pow((a + 1) / (a + x_vec[2]), 2))
                  * (a / (a + 1) * std::pow((a + 1) / (a + x_vec[3]), 2));
  return result;
}
}  // namespace hcubature_test

/*
 * test_integration is a helper function to make it easy to test the
 * hcubature function.
 *
 * It takes in a callable function object, parameters, dimension,
 * integration limits, the maximal number of evaluations, and an absolute and a
 *relative error. It integrates the provided function and compares the computed
 *integral against the provided integral (val).
 *
 * The prototype for f is:
 * template <typename T_x, typename T_p>
 * stan::return_type_t<T_x, T_p> f1(const T_x& x, const T_p& p) {
 *	using T_x_ref = stan::ref_type_t<T_x>;
 *	T_x_ref x_ref = x;
 *	stan::scalar_seq_view<T_x_ref> x_vec(x_ref);
 *  my_params* pars = static_cast<my_params*>(p);
 *  type_1 var_1 = (pars->par_1);
 *	return ;
 * }
 * The parameters can be handed to f as a struct:
 * struct my_params {
 * type_1 par_1;
 * type_2 par_2;
 * };
 *
 * @tparam F Type of f
 * @param f a functor with signature given above
 * @param dim dimension of the integral
 * @param a lower limit of integration as vector
 * @param b upper limit of integration as vector
 * @param maxEval maximal number of evaluations
 * @param reqAbsError absolute error
 * @param reqRelError relative error as vector
 * @param pars parameters to be passed to f (should be
 * std::vector<stan::math::var>)
 * @param val correct value of integral
 */

template <typename F>
void test_integration(const F& f, hcubature_test::my_params* pars, int dim,
                      std::vector<double> a, std::vector<double> b, int maxEval,
                      double reqAbsError, std::vector<double> reqRelError,
                      double val) {
  using stan::math::hcubature;

  for (auto tolerance : reqRelError) {
    EXPECT_NEAR(hcubature(f, pars, dim, a, b, maxEval, reqAbsError, tolerance),
                val, tolerance);
  }
}

// Test values
TEST(StanMath_hcubature_prim, test1) {
  // Integrals from
  // https://www.quantargo.com/help/r/latest/packages/cubature/2.0.4.1/hcubature

  int dim = 1;
  std::vector<double> a = {0.0};
  std::vector<double> b = {1.0};
  std::vector<double> reqRelError = {1e-4, 1e-6, 1e-7};
  hcubature_test::my_params pars = {};
  test_integration(hcubature_test::f1<std::vector<double>, void*>, &pars, dim,
                   a, b, 6000, 0.0, reqRelError, 0.841471);

  dim = 2;
  a = {0.0, 0.0};
  b = {1.0, 1.0};
  reqRelError = {1e-4, 1e-6, 1e-7};
  test_integration(hcubature_test::f2<std::vector<double>, void*>, &pars, dim,
                   a, b, 6000, 0.0, reqRelError, 0.7080734);

  reqRelError = {1e-4};
  pars = {0.50124145262344534123412, 0.0};
  test_integration(hcubature_test::f3<std::vector<double>, void*>, &pars, dim,
                   a, b, 10000, 0.0, reqRelError, 0.1972807);

  // (Gaussian centered at 1/2)
  reqRelError = {1e-4, 1e-6, 1e-7};
  pars = {0.0, 0.1};
  test_integration(hcubature_test::f4<std::vector<double>, void*>, &pars, dim,
                   a, b, 6000, 0.0, reqRelError, 1);

  dim = 3;
  a = {0.0, 0.0, 0.0};
  b = {1.0, 1.0, 1.0};
  reqRelError = {1e-4, 1e-6};
  test_integration(hcubature_test::f5<std::vector<double>, void*>, &pars, dim,
                   a, b, 6000, 0.0, reqRelError, 1.00001);

  reqRelError = {1e-4, 1e-6, 1e-8};
  test_integration(hcubature_test::f6<std::vector<double>, void*>, &pars, dim,
                   a, b, 6000, 0.0, reqRelError, 1);

  // (Tsuda's example)
  dim = 4;
  a = {0.0, 0.0, 0.0, 0.0};
  b = {1.0, 1.0, 1.0, 1.0};
  reqRelError = {1e-4, 1e-6};
  pars = {0.0, (1 + sqrt(10.0)) / 9.0};
  test_integration(hcubature_test::f7<std::vector<double>, void*>, &pars, dim,
                   a, b, 20000, 0.0, reqRelError, 0.999998);
}
