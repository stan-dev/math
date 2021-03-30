#include <stan/math/rev.hpp>
#include <boost/numeric/odeint.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/prim/functor/harmonic_oscillator.hpp>
#include <test/unit/math/rev/functor/test_fixture_ode.hpp>
#include <test/unit/math/rev/functor/test_fixture_ode_sho.hpp>
#include <test/unit/math/rev/functor/ode_test_functors.hpp>

#include <boost/preprocessor/seq/for_each_product.hpp>
#include <boost/preprocessor/seq/to_tuple.hpp>
#include <boost/preprocessor/tuple/to_seq.hpp>
#include <boost/preprocessor/tuple/enum.hpp>
#include <boost/preprocessor/tuple/elem.hpp>
#include <boost/preprocessor/seq/for_each.hpp>
#include <boost/preprocessor/seq/enum.hpp>
#include <boost/preprocessor/seq/elem.hpp>

#define TEST_TYPE(R,SEQ_X) \
  std::tuple<BOOST_PP_TUPLE_ELEM(3,0,BOOST_PP_SEQ_TO_TUPLE(SEQ_X)),     \
             BOOST_PP_TUPLE_ELEM(3,0,BOOST_PP_SEQ_TO_TUPLE(SEQ_X)),     \
             double,                                                    \
             BOOST_PP_TUPLE_ELEM(3, 1, BOOST_PP_SEQ_TO_TUPLE(SEQ_X)),   \
             BOOST_PP_TUPLE_ELEM(3, 2, BOOST_PP_SEQ_TO_TUPLE(SEQ_X))>,
   
#define TEST_TYPE_PRODUCT(TUP_SOLVER,TUP_TY,TUP_TP)             \
  BOOST_PP_SEQ_FOR_EACH_PRODUCT (  \
    TEST_TYPE, \
     (BOOST_PP_TUPLE_TO_SEQ(TUP_SOLVER)) \
     (BOOST_PP_TUPLE_TO_SEQ(TUP_TY)) \
     (BOOST_PP_TUPLE_TO_SEQ(TUP_TP)) \
   )

#define FUNCTOR_TYPES (ode_ckrk_functor, ode_rk45_functor, ode_bdf_functor, ode_adams_functor)
#define TY_TYPES (double, stan::math::var, stan::math::var_value<double>)
#define TP_TYPES (double, stan::math::var, stan::math::var_value<double>)

using harmonic_oscillator_test_types = ::testing::Types<
  TEST_TYPE_PRODUCT(FUNCTOR_TYPES, TY_TYPES, TP_TYPES)
  std::tuple<ode_rk45_functor, ode_rk45_functor, double, stan::math::var_value<double>, stan::math::var_value<double>>>;

TYPED_TEST_SUITE_P(harmonic_oscillator_analytical_test);
TYPED_TEST_P(harmonic_oscillator_analytical_test, dv) {
  double omega = 0.5;
  double chi = 1.25;
  double t = 5;
  this->ode_sol_dv.ts = std::vector<double>{t};
  this->ode_sol_dv.y0[0] = chi;
  Eigen::VectorXd x(1);
  x << omega;
  this->test_analytical(this->ode_sol_dv, this->analy_sol_functor(),
                        this->analy_grad_omega_sol_functor(), 1e-5, x, t, omega,
                        chi);
}
TYPED_TEST_P(harmonic_oscillator_analytical_test, vd) {
  double omega = 0.5;
  double chi = 1.25;
  double t = 5;
  this->ode_sol_vd.ts = std::vector<double>{t};
  this->ode_sol_vd.theta[0] = omega;
  Eigen::VectorXd x(1);
  x << chi;
  this->test_analytical(this->ode_sol_vd, this->analy_sol_functor(),
                        this->analy_grad_chi_sol_functor(), 1e-5, x, t, omega,
                        chi);
}

TYPED_TEST_P(harmonic_oscillator_analytical_test, vv) {
  double omega = 0.5;
  double chi = 1.25;
  double t = 5;
  this->ode_sol_vv.ts = std::vector<double>{t};
  Eigen::VectorXd x(2);
  x << omega, chi;
  this->test_analytical(this->ode_sol_vv, this->analy_sol_functor(),
                        this->analy_grad_sol_functor(), 1e-5, x, t, omega, chi);
}

REGISTER_TYPED_TEST_SUITE_P(harmonic_oscillator_analytical_test, dv, vd, vv);
INSTANTIATE_TYPED_TEST_SUITE_P(StanOde, harmonic_oscillator_analytical_test,
                               harmonic_oscillator_test_types);
