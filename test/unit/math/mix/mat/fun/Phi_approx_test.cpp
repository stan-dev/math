#include <stan/math/mix/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/prim/mat/vectorize/prim_scalar_unary_test.hpp>
#include <test/unit/math/rev/mat/vectorize/rev_scalar_unary_test.hpp>
#include <test/unit/math/fwd/mat/vectorize/fwd_scalar_unary_test.hpp>
#include <test/unit/math/mix/mat/vectorize/mix_scalar_unary_test.hpp>
#include <test/unit/math/prim/mat/vectorize/vector_builder.hpp>
#include <vector>

/**
 * This is the structure for testing vectorized Phi_approx (defined in the
 * testing framework).
 */
struct Phi_approx_test {

  /**
   * Redefinition of function brought in from stan::math.  The reason
   * to do this is that it wraps it up in this static template class.
   *
   * This is the version that's being tested.  
   *
   * WARNING:  assumes that the scalar values for all instantiations
   * (prim, rev, fwd, mix) ***have already been tested***.
   *
   * @tparam R Return type.
   * @tparam T Argument type.
   */
  template <typename R, typename T>
  static R apply(const T& x) {
    using stan::math::Phi_approx;
    return Phi_approx(x);
  }

  /**
   * This defines the truth against which we're testing.
   *
   * Because this is *not an independent test*, this function just
   * delegates to the actual function defined in stan::math.
   *
   * Redundant definition of function from stan::math to apply to an
   * integer and return a double.
   *
   * This function delegates to apply(), defined above, directly.  
   *
   * WARNING:  this is *not an independent test*.
   */
  static double apply_base(int x) {
    return apply<double>(x);
  }

  /**
   * This is the generic version of the integer version defined
   * above.  For every other type, the return type is the same as the
   * reference type.
   *
   * WARNING:  this is *not an independent test of the underlying function*.
   */
  template <typename T>
  static T apply_base(const T& x) {
    return apply<T>(x);
  }

  /**
   * Return sequence of valid double-valued inputs.
   */
  static std::vector<double> valid_inputs() {
    return test::math::vector_builder<double>()
      .add(-15.2).add(0.0).add(1.3).build();
  }

  /**
   * Return sequence of invalid double-valued inputs.
   */
  static std::vector<double> invalid_inputs() {
    return test::math::vector_builder<double>()
      .build();
  }

  /**
   * Return sequence of valid integer inputs.
   */
  static std::vector<int> int_valid_inputs() {
    return test::math::vector_builder<int>()
      .add(-10).add(0).add(1).add(5).add(10).build();
  }

  /**
   * Return sequence of invalid integer inputs.
   */
  static std::vector<int> int_invalid_inputs() {
    return test::math::vector_builder<int>()
      .build();
  }
};

INSTANTIATE_TYPED_TEST_CASE_P(, prim_scalar_unary_test, Phi_approx_test);
INSTANTIATE_TYPED_TEST_CASE_P(, rev_scalar_unary_test, Phi_approx_test);
INSTANTIATE_TYPED_TEST_CASE_P(, fwd_scalar_unary_test, Phi_approx_test);
INSTANTIATE_TYPED_TEST_CASE_P(, mix_scalar_unary_test, Phi_approx_test);
