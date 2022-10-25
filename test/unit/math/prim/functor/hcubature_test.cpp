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
	return cos(x_vec[0])*cos(x_vec[1]);
}


template <typename T_x, typename T_p>
stan::return_type_t<T_x, T_p> f2(const T_x& x, const T_p& p) {
	using T_x_ref = stan::ref_type_t<T_x>;
	T_x_ref x_ref = x;
	stan::scalar_seq_view<T_x_ref> x_vec(x_ref);
	my_params *pars = static_cast<my_params*> (p);
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
stan::return_type_t<T_x, T_p> f3(const T_x& x, const T_p& p) {
	using T_x_ref = stan::ref_type_t<T_x>;
	using namespace stan::math;
	T_x_ref x_ref = x;
	stan::scalar_seq_view<T_x_ref> x_vec(x_ref);
	my_params *pars = static_cast<my_params*> (p);
	double a = (pars->y);
	double s = std::pow(x_vec[0] - 0.5, 2) + std::pow(x_vec[1] - 0.5, 2);
	return std::pow(TWO_OVER_SQRT_PI / (2.0 * a), 2) * exp(-s / std::pow(a, 2));
}


	

template <typename T_x, typename T_p>
stan::return_type_t<T_x, T_p> f4(const T_x& x, const T_p& p) {
	using T_x_ref = stan::ref_type_t<T_x>;
	using namespace stan::math;
	T_x_ref x_ref = x;
	stan::scalar_seq_view<T_x_ref> x_vec(x_ref);
	double val = std::pow((1-x_vec[0])/x_vec[0], 2) +
				 std::pow((1-x_vec[1])/x_vec[1], 2) +
				 std::pow((1-x_vec[2])/x_vec[2], 2);
	double scale = TWO_OVER_SQRT_PI/std::pow(x_vec[0], 2) *
				   TWO_OVER_SQRT_PI/std::pow(x_vec[1], 2) *
				   TWO_OVER_SQRT_PI/std::pow(x_vec[2], 2);
	return exp(-val) * scale;
}


template <typename T_x, typename T_p>
stan::return_type_t<T_x, T_p> f5(const T_x& x, const T_p& p) {
	using T_x_ref = stan::ref_type_t<T_x>;
	T_x_ref x_ref = x;
	stan::scalar_seq_view<T_x_ref> x_vec(x_ref);
	return 2*x_vec[0] * 2*x_vec[1] * 2*x_vec[2];
}


template <typename T_x, typename T_p>
stan::return_type_t<T_x, T_p> f6(const T_x& x, const T_p& p) {
	using T_x_ref = stan::ref_type_t<T_x>;
	T_x_ref x_ref = x;
	stan::scalar_seq_view<T_x_ref> x_vec(x_ref);
	my_params *pars = static_cast<my_params*> (p);
	double a = (pars->y);
	double result = (a / (a + 1) * std::pow((a + 1) / (a + x_vec[0]), 2)) *
					(a / (a + 1) * std::pow((a + 1) / (a + x_vec[1]), 2)) *
					(a / (a + 1) * std::pow((a + 1) / (a + x_vec[2]), 2)) *
					(a / (a + 1) * std::pow((a + 1) / (a + x_vec[3]), 2));
	return result;
}
}  // namespace hcubature_test


/*
 * test_integration is a helper function to make it easy to test the
 * hcubature function.
 *
 * It takes in a callable function object, parameters, dimension, 
 * integration limits, the maximal number of evaluations, and an absolute and a relative error. It integrates 
 * the provided function and compares the
 * computed integral against the provided integral (val).
 *
 * The prototype for f is:
 * template <typename T_x, typename T_p>
 * stan::return_type_t<T_x, T_p> f1(const T_x& x, const T_p& p) {
 *	using T_x_ref = stan::ref_type_t<T_x>;
 *	T_x_ref x_ref = x;
 *	stan::scalar_seq_view<T_x_ref> x_vec(x_ref);
 *	return ;
 * }
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
void test_integration(const F &f, 
					  hcubature_test::my_params* pars,
					  int dim,
					  std::vector<double> a, 
					  std::vector<double> b,
                      double maxEval,
                      double reqAbsError,
                      std::vector<double> reqRelError, 
					  double val) {
  using stan::math::hcubature;
  

  for (auto tolerance : reqRelError) {
    EXPECT_LE(std::abs(hcubature(f, pars, dim, a, b, maxEval, reqAbsError, 
								 tolerance)
                       - val),
              tolerance);		  
  }
}





/*
TEST(StanMath_hcubature_prim, TestThrows) {
  // Left limit of integration must be less than or equal to right limit
  EXPECT_THROW(stan::math::hcubature(hcubature_test::f2{}, 1.0, 0.0,
                                        std::vector<double>(), {}, {},
                                        hcubature_test::msgs, 1e-6),
               std::domain_error);
  // NaN limits not okay
  EXPECT_THROW(
      stan::math::hcubature(hcubature_test::f2{}, 0.0,
                               std::numeric_limits<double>::quiet_NaN(),
                               std::vector<double>(), {}, {},
                               hcubature_test::msgs, 1e-6),
      std::domain_error);
  EXPECT_THROW(
      stan::math::hcubature(
          hcubature_test::f2{}, std::numeric_limits<double>::quiet_NaN(),
          0.0, std::vector<double>(), {}, {}, hcubature_test::msgs, 1e-6),
      std::domain_error);
  EXPECT_THROW(
      stan::math::hcubature(
          hcubature_test::f2{}, std::numeric_limits<double>::quiet_NaN(),
          std::numeric_limits<double>::quiet_NaN(), std::vector<double>(), {},
          {}, hcubature_test::msgs, 1e-6),
      std::domain_error);
  // Two of the same inf limits not okay
  EXPECT_THROW(
      stan::math::hcubature(
          hcubature_test::f2{}, -std::numeric_limits<double>::infinity(),
          -std::numeric_limits<double>::infinity(), std::vector<double>(), {},
          {}, hcubature_test::msgs, 1e-6),
      std::domain_error);

  EXPECT_THROW(stan::math::hcubature(hcubature_test::f2{},
                                        std::numeric_limits<double>::infinity(),
                                        std::numeric_limits<double>::infinity(),
                                        std::vector<double>(), {}, {},
                                        hcubature_test::msgs, 1e-6),
               std::domain_error);
  // xc should be nan if there are infinite limits
  EXPECT_THROW(stan::math::hcubature(hcubature_test::f11{}, 0.0,
                                        std::numeric_limits<double>::infinity(),
                                        std::vector<double>(), {}, {},
                                        hcubature_test::msgs, 1e-6),
               std::runtime_error);
  EXPECT_THROW(stan::math::hcubature(hcubature_test::f11{},
                                        std::numeric_limits<double>::infinity(),
                                        0.0, std::vector<double>(), {}, {},
                                        hcubature_test::msgs, 1e-6),
               std::domain_error);
  EXPECT_THROW(stan::math::hcubature(hcubature_test::f11{},
                                        std::numeric_limits<double>::infinity(),
                                        std::numeric_limits<double>::infinity(),
                                        std::vector<double>(), {}, {},
                                        hcubature_test::msgs, 1e-6),
               std::domain_error);
  // But not otherwise
  EXPECT_NO_THROW(stan::math::hcubature(hcubature_test::f11{}, 0.0, 1.0,
                                           std::vector<double>(), {}, {},
                                           hcubature_test::msgs, 1e-6));
}

TEST(StanMath_hcubature_prim, test_integer_arguments) {
  double v;
  EXPECT_NO_THROW(v = stan::math::hcubature(hcubature_test::f2{}, 0, 1,
                                               std::vector<double>(), {}, {},
                                               hcubature_test::msgs, 1e-6));
  EXPECT_NO_THROW(v = stan::math::hcubature(hcubature_test::f2{}, 0.0, 1,
                                               std::vector<double>(), {}, {},
                                               hcubature_test::msgs, 1e-6));
  EXPECT_NO_THROW(v = stan::math::hcubature(hcubature_test::f2{}, 0, 1.0,
                                               std::vector<double>(), {}, {},
                                               hcubature_test::msgs, 1e-6));
}

*/

// Test values
TEST(StanMath_hcubature_prim, test1) {
  // Integrals from https://www.quantargo.com/help/r/latest/packages/cubature/2.0.4.1/hcubature

int dim = 2;
std::vector<double> a = {0.0, 0.0};
std::vector<double> b = {1.0, 1.0};
std::vector<double> reqRelError = {1e-4, 1e-6, 1e-7};
hcubature_test::my_params pars = {};

test_integration(hcubature_test::f1<std::vector<double>, void*>, 
				 &pars, dim, a, b, 6000, 0.0, reqRelError,
                 0.7080734);


reqRelError = {1e-4};
pars = {0.50124145262344534123412, 0.0};
test_integration(hcubature_test::f2<std::vector<double>, void*>, 
				 &pars, dim, a, b, 6000, 0.0, reqRelError,
                 0.1972807);


// (Gaussian centered at 1/2)
reqRelError = {1e-4, 1e-6, 1e-7};
pars = {0.0, 0.1};
test_integration(hcubature_test::f3<std::vector<double>, void*>, 
				 &pars, dim, a, b, 6000, 0.0, reqRelError,
                 1);


dim = 3;
a = {0.0, 0.0, 0.0};
b = {1.0, 1.0, 1.0};
reqRelError = {1e-4, 1e-6};
test_integration(hcubature_test::f4<std::vector<double>, void*>, 
				 &pars, dim, a, b, 6000, 0.0, reqRelError,
                 1.00001);


reqRelError = {1e-4, 1e-6, 1e-8};
test_integration(hcubature_test::f5<std::vector<double>, void*>, 
				 &pars, dim, a, b, 6000, 0.0, reqRelError,
                 1);


// (Tsuda's example)
dim = 4;
a = {0.0, 0.0, 0.0, 0.0};
b = {1.0, 1.0, 1.0, 1.0};
reqRelError = {1e-4, 1e-6};
pars = {0.0, (1 + sqrt(10.0)) / 9.0};
test_integration(hcubature_test::f6<std::vector<double>, void*>, 
				 &pars, dim, a, b, 19000, 0.0, reqRelError,
                 0.999998);
}






/*



template <typename T_y, typename T_alpha, typename T_delta, typename T_beta, typename T_t0, typename T_sv, typename T_sw, typename T_st0, typename T_err>
struct my_params {
	T_y y;
	T_alpha a;
	T_delta v;
	T_beta w;
	T_t0 t0;
	T_sv sv;
	T_sw sw;
	T_st0 st0;
	T_err lerr;
};


	using T_return_type = return_type_t<T_x, T_p>; 
	my_params<T_return_type, T_return_type, T_return_type, T_return_type, T_return_type, T_return_type, T_return_type, 
	T_return_type, T_return_type>  *params = static_cast<my_params<T_return_type, T_return_type, T_return_type, T_return_type, 
	T_return_type, T_return_type, T_return_type, T_return_type, T_return_type> *>(p);
	T_return_type y = (params->y);
	T_return_type a = (params->a);


internal::my_params<T_partials_return, T_partials_return, T_partials_return, T_partials_return, 
T_partials_return, T_partials_return, T_partials_return, T_partials_return, T_partials_return> params_new_error = {y_val, alpha_val, delta_val, beta_val, t0_val, sv_val, sw_val, st0_val, log(.1) + lerror_bound - log(2) + log(fabs(deriv))};
				





 1/dens*hcubature(internal::int_dsvddiff<std::vector<double>, void*>, &params_new_error, dim, xmin, xmax, Meval, abstol, reltol/2);
				



*/