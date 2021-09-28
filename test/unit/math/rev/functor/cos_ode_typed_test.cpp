#include <stan/math/rev.hpp>
#include <boost/mp11.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/functor/test_fixture_ode.hpp>
#include <test/unit/math/rev/functor/test_fixture_ode_cos_scalar.hpp>
#include <test/unit/math/rev/functor/ode_test_functors.hpp>

/**
 *
 * Use same solver functor type for both w & w/o tolerance control
 */
template <typename solve_type, typename... Ts>
using ode_test_tuple = std::tuple<solve_type, solve_type, Ts...>;

/**
 * Outer product of test types
 */
using cos_arg_test_types = boost::mp11::mp_product<
    ode_test_tuple,
    ::testing::Types<ode_adams_functor, ode_bdf_functor, ode_ckrk_functor,
                     ode_rk45_functor, ode_adjoint_functor> >;

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
TYPED_TEST_P(cos_arg_test, value) { this->test_value(); }
TYPED_TEST_P(cos_arg_test, grad) {
  this->test_grad_t0();
  this->test_grad_ts();
  this->test_grad_ts_repeat();
  this->test_scalar_arg();
  this->test_std_vector_arg();
  this->test_vector_arg();
  this->test_row_vector_arg();
  this->test_matrix_arg();
  this->test_scalar_std_vector_args();
  this->test_std_vector_std_vector_args();
  this->test_std_vector_vector_args();
  this->test_std_vector_row_vector_args();
  this->test_std_vector_matrix_args();
  this->test_arg_combos_test();
}
TYPED_TEST_P(cos_arg_test, tol_grad) {
  this->test_tol_int_t0();
  this->test_tol_int_ts();
  this->test_tol_t0();
  this->test_tol_ts();
  this->test_tol_ts_repeat();
  this->test_tol_scalar_arg();
  this->test_tol_scalar_arg_multi_time();
  this->test_tol_std_vector_arg();
  this->test_tol_vector_arg();
  this->test_tol_row_vector_arg();
  this->test_tol_matrix_arg();
  this->test_tol_scalar_std_vector_args();
  this->test_tol_std_vector_std_vector_args();
  this->test_tol_std_vector_vector_args();
  this->test_tol_std_vector_row_vector_args();
  this->test_tol_std_vector_matrix_args();
  this->test_tol_arg_combos_test();
}
REGISTER_TYPED_TEST_SUITE_P(cos_arg_test, y0_error, t0_error, ts_error,
                            one_arg_error, two_arg_error, rhs_wrong_size_error,
                            error_name, tol_error, value, grad, tol_grad);
INSTANTIATE_TYPED_TEST_SUITE_P(StanOde, cos_arg_test, cos_arg_test_types);
