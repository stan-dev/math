#include <stan/math/prim/functor.hpp>
#include <stan/math/prim/fun.hpp>

#include <gtest/gtest.h>
#include <vector>

namespace hcubature_test {

template <typename T_x>
inline auto f1(const T_x& x) {
  return stan::math::sum(stan::math::cos(x));
}

template <typename T_x>
inline auto f2(const T_x& x) {
  return stan::math::prod(stan::math::cos(x));
}

template <typename T_x>
inline auto f3(const T_x& x, double radius) {
  using stan::math::square;
  return stan::return_type_t<T_x>(stan::math::sum(square(x)) < square(radius));
}

template <typename T_x>
inline auto f4(const T_x& x, double sigma) {
  using stan::math::as_array_or_scalar;
  using stan::math::square;
  using stan::math::sum;
  using stan::math::TWO_OVER_SQRT_PI;

  auto numerator = sum(square(as_array_or_scalar(x) - 0.5));
  return square(TWO_OVER_SQRT_PI / (2.0 * sigma))
         * exp(-numerator / square(sigma));
}

template <typename T_x>
inline auto f5(const T_x& x) {
  using stan::math::prod;
  using stan::math::square;
  using stan::math::TWO_OVER_SQRT_PI;

  const auto& x_arr = stan::math::as_array_or_scalar(x);
  auto val = stan::math::sum(square((1 - x_arr) / x_arr));
  auto scale = prod(TWO_OVER_SQRT_PI / square(x_arr));
  return exp(-val) * scale;
}

template <typename T_x>
inline auto f6(const T_x& x) {
  const auto& x_arr = stan::math::as_array_or_scalar(x);
  return stan::math::prod(2 * x_arr);
}

template <typename T_x>
inline auto f7(const T_x& x, double a) {
  const auto& x_arr = stan::math::as_array_or_scalar(x);
  return stan::math::prod(a / (a + 1)
                          * stan::math::square((a + 1) / (a + x_arr)));
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

template <typename F, typename ArgsTupleT, typename T_a, typename T_b,
          typename T_relerr>
void test_integration(const F& f, const ArgsTupleT& pars, int dim, const T_a& a,
                      const T_b& b, int maxEval, double reqAbsError,
                      const T_relerr& reqRelError, double val) {
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
  const Eigen::VectorXd a{{0.0}};
  const Eigen::VectorXd b{{1.0}};
  const Eigen::VectorXd reqRelError{{1e-4, 1e-6, 1e-7}};
  test_integration([](auto&& x) { return hcubature_test::f1(x); },
                   std::make_tuple(), dim, a, b, 6000, 0.0, reqRelError,
                   0.841471);

  dim = 2;
  const Eigen::VectorXd a_2{{0.0, 0.0}};
  const Eigen::VectorXd b_2{{1.0, 1.0}};
  test_integration([](auto&& x) { return hcubature_test::f2(x); },
                   std::make_tuple(), dim, a_2, b_2, 6000, 0.0, reqRelError,
                   0.7080734);

  const Eigen::VectorXd reqRelError_2{{1e-4}};
  test_integration(
      [](auto&& x, auto&& radius) { return hcubature_test::f3(x, radius); },
      std::make_tuple(0.50124145262344534123412), dim, a_2, b_2, 10000, 0.0,
      reqRelError_2, 0.1972807);

  // (Gaussian centered at 1/2)
  test_integration(
      [](auto&& x, auto&& sigma) { return hcubature_test::f4(x, sigma); },
      std::make_tuple(0.1), dim, a_2, b_2, 6000, 0.0, reqRelError, 1);

  dim = 3;
  const Eigen::VectorXd a_3{{0.0, 0.0, 0.0}};
  const Eigen::VectorXd b_3{{1.0, 1.0, 1.0}};
  const Eigen::VectorXd reqRelError_3{{1e-4, 1e-6}};
  test_integration([](auto&& x) { return hcubature_test::f5(x); },
                   std::make_tuple(), dim, a_3, b_3, 6000, 0.0, reqRelError_3,
                   1.00001);

  const Eigen::VectorXd reqRelError_4{{1e-4, 1e-6, 1e-8}};
  test_integration([](auto&& x) { return hcubature_test::f6(x); },
                   std::make_tuple(), dim, a_3, b_3, 6000, 0.0, reqRelError_4,
                   1);

  // (Tsuda's example)
  dim = 4;
  const Eigen::VectorXd a_4{{0.0, 0.0, 0.0, 0.0}};
  const Eigen::VectorXd b_4{{1.0, 1.0, 1.0, 1.0}};
  test_integration([](auto&& x, auto&& a) { return hcubature_test::f7(x, a); },
                   std::make_tuple((1 + sqrt(10.0)) / 9.0), dim, a_4, b_4,
                   20000, 0.0, reqRelError_3, 0.999998);
}
