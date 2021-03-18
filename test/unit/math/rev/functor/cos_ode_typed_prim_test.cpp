#include <stan/math/rev.hpp>
#include <boost/numeric/odeint.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/functor/test_fixture_ode.hpp>
#include <test/unit/math/rev/functor/test_fixture_ode_cos_scalar.hpp>
#include <test/unit/math/rev/functor/ode_test_functors.hpp>

using cos_arg_test_types = ::testing::Types<
  std::tuple<ode_rk45_functor, ode_rk45_functor>,
  std::tuple<ode_ckrk_functor, ode_ckrk_functor>,
  std::tuple<ode_bdf_functor, ode_bdf_functor>,
  std::tuple<ode_adams_functor, ode_adams_functor> >;

TYPED_TEST_SUITE_P(cos_arg_test);
TYPED_TEST_P(cos_arg_test, y0_error) {
  this->test_y0_error();
  this->test_y0_error_with_tol();
}
TYPED_TEST_P(cos_arg_test, t0_error) {
  this->test_t0_error();
  this->test_t0_error_with_tol();
}
TYPED_TEST_P(cos_arg_test, ts_error) {
  this->test_ts_error();
  this->test_ts_error_with_tol();
}
TYPED_TEST_P(cos_arg_test, one_arg_error) {
  this->test_one_arg_error();
  this->test_one_arg_error_with_tol();
}
TYPED_TEST_P(cos_arg_test, two_arg_error) {
  this->test_two_arg_error();
  this->test_two_arg_error_with_tol();
}
TYPED_TEST_P(cos_arg_test, rhs_wrong_size_error) {
  this->test_rhs_wrong_size_error();
  this->test_rhs_wrong_size_error_with_tol();
}
TYPED_TEST_P(cos_arg_test, error_name) {
  this->test_error_name();
  this->test_error_name_with_tol();
}
TYPED_TEST_P(cos_arg_test, tol_error) {
  this->test_rtol_error();
  this->test_atol_error();
  this->test_max_num_step_error();
  this->test_too_much_work();
}
REGISTER_TYPED_TEST_SUITE_P(cos_arg_test, y0_error, t0_error,
                            ts_error, one_arg_error, two_arg_error,
                            rhs_wrong_size_error, error_name, tol_error);
INSTANTIATE_TYPED_TEST_SUITE_P(StanOde, cos_arg_test,
                               cos_arg_test_types);
