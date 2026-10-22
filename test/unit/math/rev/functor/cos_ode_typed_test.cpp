#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <boost/mp11.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/functor/test_fixture_ode.hpp>
#include <test/unit/math/rev/functor/test_fixture_ode_cos_scalar.hpp>
#include <test/unit/math/rev/functor/ode_test_functors.hpp>
#include <cmath>

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
  EXPECT_NO_THROW(this->test_grad_t0());
  EXPECT_NO_THROW(this->test_grad_ts());
  EXPECT_NO_THROW(this->test_grad_ts_repeat());
  EXPECT_NO_THROW(this->test_scalar_arg());
  EXPECT_NO_THROW(this->test_std_vector_arg());
  EXPECT_NO_THROW(this->test_vector_arg());
  EXPECT_NO_THROW(this->test_row_vector_arg());
  EXPECT_NO_THROW(this->test_matrix_arg());
  EXPECT_NO_THROW(this->test_scalar_std_vector_args());
  EXPECT_NO_THROW(this->test_std_vector_std_vector_args());
  EXPECT_NO_THROW(this->test_std_vector_vector_args());
  EXPECT_NO_THROW(this->test_std_vector_row_vector_args());
  EXPECT_NO_THROW(this->test_std_vector_matrix_args());
  EXPECT_NO_THROW(this->test_arg_combos_test());
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

TEST(StanOdeAdjoint, quadrature_buffer_reuse) {
  // Repeated solves and reverse passes reuse CVODES' quadrature work buffers.
  // Every parameter has the same analytic derivative through their sum.
  const std::vector<double> times{0.45, 1.1};
  constexpr double rate = 1.5;
  for (int iteration = 0; iteration < 10; ++iteration) {
    SCOPED_TRACE(iteration);
    stan::math::nested_rev_autodiff nested;
    const Eigen::VectorXd initial = Eigen::VectorXd::Zero(1);
    std::vector<stan::math::var> parameters{0.25, 0.5, 0.75};
    ode_adjoint_functor solver;
    auto output = solver(stan::test::CosArg1(), initial, 0.0, times, nullptr,
                         parameters);
    for (int index : {0, 1, 0, 1}) {
      SCOPED_TRACE(index);
      nested.set_zero_all_adjoints();
      ASSERT_NO_THROW(output[index](0).grad());
      const double phase = rate * times[index];
      const double expected
          = (phase * std::cos(phase) - std::sin(phase)) / (rate * rate);
      EXPECT_NEAR(output[index](0).val(), std::sin(phase) / rate, 1e-8);
      for (const auto& parameter : parameters) {
        EXPECT_NEAR(parameter.adj(), expected, 1e-8);
      }
    }
  }
}
