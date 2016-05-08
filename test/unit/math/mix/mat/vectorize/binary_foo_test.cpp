#include <stan/math/mix/mat.hpp>
#include <gtest/gtest.h>
#include <boost/math/tools/promotion.hpp>
#include <test/unit/math/prim/mat/vectorize/prim_scalar_binary_test.hpp>
#include <test/unit/math/rev/mat/vectorize/rev_scalar_binary_test.hpp>
//#include <test/unit/math/fwd/mat/vectorize/fwd_scalar_unary_test.hpp>
//#include <test/unit/math/mix/mat/vectorize/mix_scalar_unary_test.hpp>
#include <test/unit/math/prim/mat/vectorize/binary_foo_fun.hpp>
#include <test/unit/math/prim/mat/vectorize/vector_builder.hpp>

/**
 * This is the structure for testing mock function foo (defined in the
 * testing framework).  See README.txt for more instructions.
 */
struct binary_foo_test {

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
  template <typename R, typename T1, typename T2>
  static R apply(const T1& x, const T2& y) {
    using stan::math::binary_foo;
    return binary_foo(x, y);
  }

  /**
   * This is the generic version of the integer version defined
   * above.  For every other type, the return type is the same as the
   * reference type.
   *
   * WARNING:  this is *not an independent test of the underlying function*.
   */
  static double apply_base(int x, int y) {
    return apply<double, int, int>(x, y);
  }

  /**
   * This is the generic version of the integer version defined
   * above.  For every other type, the return type is the same as the
   * reference type.
   *
   * WARNING:  this is *not an independent test of the underlying function*.
   */
  template <typename T1, typename T2>
  static typename boost::math::tools::promote_args<T1, T2>::type
  apply_base(const T1& x, const T2& y) {
    return apply<typename boost::math::tools::promote_args<T1, T2>::type, 
                 T1, T2>(x, y);
  }

  /**
   * Return sequence of valid double-valued inputs.
   */
  static std::vector<double> valid_inputs() {
    return test::math::vector_builder<double>()
      .add(1.3).add(-2.6).add(0).add(-0.2).build();
  }

  /**
   * Return sequence of valid double-valued inputs.
   */
  static std::vector<double> valid_inputs1() {
    return test::math::vector_builder<double>()
      .add(0.7).add(2.3).add(3.5).add(0).add(0).add(0)
      .add(-0.3).add(-5.3).add(-3.7).build();
  }

  /**
   * Return sequence of valid double-valued inputs.
   */
  static std::vector<double> valid_inputs2() {
    return test::math::vector_builder<double>()
      .add(1.3).add(-2.6).add(0).add(-1.2).add(2.3).add(0)
      .add(-3.7).add(0).add(2.3).build();
  }

  /**
   * Return sequence of invalid double-valued inputs.
   */
  static std::vector<double> invalid_inputs1() {
    return std::vector<double>();
  }

  /**
   * Return sequence of invalid double-valued inputs.
   */
  static std::vector<double> invalid_inputs2() {
    return std::vector<double>();
  }

  /**
   * Return sequence of valid integer inputs.
   */
  static std::vector<int> int_valid_inputs() {
    return test::math::vector_builder<int>()
      .add(1).add(-2).add(0).add(3).build();
  }

  /**
   * Return sequence of valid integer inputs.
   */
  static std::vector<int> int_valid_inputs1() {
    return test::math::vector_builder<int>()
      .add(6).add(3).add(7).add(0).add(0).add(0)
      .add(-2).add(-5).add(-4).build();
  }

  /**
   * Return sequence of valid integer inputs.
   */
  static std::vector<int> int_valid_inputs2() {
    return test::math::vector_builder<int>()
      .add(3).add(-2).add(0).add(5).add(0).add(-4)
      .add(0).add(6).add(-3).build();
  }

  /**
   * Return sequence of invalid integer inputs.
   */
  static std::vector<int> int_invalid_inputs1() {
    return std::vector<int>();
  }

  /**
   * Return sequence of invalid integer inputs.
   */
  static std::vector<int> int_invalid_inputs2() {
    return std::vector<int>();
  }
};

INSTANTIATE_TYPED_TEST_CASE_P(, prim_scalar_binary_test, binary_foo_test);
INSTANTIATE_TYPED_TEST_CASE_P(, rev_scalar_unary_test, binary_foo_test);
//INSTANTIATE_TYPED_TEST_CASE_P(, fwd_scalar_unary_test, foo_test);
//INSTANTIATE_TYPED_TEST_CASE_P(, mix_scalar_unary_test, foo_test);
