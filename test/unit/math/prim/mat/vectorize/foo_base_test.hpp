#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_FOO_BASE_TEST_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_FOO_BASE_TEST_HPP

#include <test/unit/math/prim/mat/vectorize/foo_fun.hpp>
#include <test/unit/math/prim/mat/vectorize/vector_builder.hpp>

/**
 * This is the structure for testing mock function foo (defined in the
 * testing framework).  See README.txt for more instructions.
 */
struct foo_base_test {

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
    using stan::math::foo;
    return foo(x);
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
    return stan::math::vector_builder<double>()
      .add(1.3).add(-2.6).add(0).add(-0.2).build();
  }

  /**
   * Return sequence of invalid double-valued inputs.
   */
  static std::vector<double> illegal_inputs() {
    return stan::math::vector_builder<double>()
      .add(10.6).add(10.6).add(25.7).add(100.25).build();
  }

  /**
   * Return sequence of valid integer inputs.
   */
  static std::vector<int> int_valid_inputs() {
    return stan::math::vector_builder<int>()
      .add(1).add(-2).add(0).add(3).build();
  }

  /**
   * Return sequence of invalid integer inputs.
   */
  static std::vector<int> int_illegal_inputs() {
    return stan::math::vector_builder<int>()
      .add(10).add(25).add(100).add(50).build();
  }
};
#endif
