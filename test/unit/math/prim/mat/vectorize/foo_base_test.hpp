#ifndef TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_FOO_BASE_TEST_HPP
#define TEST_UNIT_MATH_PRIM_MAT_VECTORIZE_FOO_BASE_TEST_HPP

#include <stan/math/prim/mat/vectorize/apply_scalar_unary.hpp>
#include <stan/math/fwd/mat/vectorize/apply_scalar_unary.hpp>
#include <stan/math/rev/mat/vectorize/apply_scalar_unary.hpp>
#include <test/unit/math/prim/mat/vectorize/foo_fun.hpp>

struct foo_base_test {

  template <typename R, typename T>
  static R apply(const T& x) {
    using stan::math::foo;
    return foo(x);
  }

  static double apply_base(const int x) {
    return apply<double>(x);
  }

  template <typename T>
  static T apply_base(const T& x) {
    return apply<T>(x);
  }

  static std::vector<double> valid_inputs() {
    using std::vector;

    vector<double> valid_inputs;
    valid_inputs.push_back(1.3);
    valid_inputs.push_back(-2.6);
    valid_inputs.push_back(0);
    valid_inputs.push_back(-0.2);

    return valid_inputs;
  }

  static std::vector<double> illegal_inputs() {
    using std::vector;

    vector<double> illegal_inputs(2, 10.6);
    illegal_inputs.push_back(25.7);
    illegal_inputs.push_back(100.25);

    return illegal_inputs;
  }

  static std::vector<int> int_valid_inputs() {
    using std::vector;

    vector<int> valid_inputs;
    valid_inputs.push_back(1);
    valid_inputs.push_back(-2);
    valid_inputs.push_back(0);
    valid_inputs.push_back(3);

    return valid_inputs;
  }

  static std::vector<int> int_illegal_inputs() {
    using std::vector;

    vector<int> illegal_inputs;
    illegal_inputs.push_back(10);
    illegal_inputs.push_back(25);
    illegal_inputs.push_back(100);
    illegal_inputs.push_back(50);

    return illegal_inputs;
  }

};
#endif
