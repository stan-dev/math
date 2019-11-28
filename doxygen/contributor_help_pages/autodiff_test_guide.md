# Automatic Differentiation Testing Guide {#autodiff_test_guide}

The automatic differentiation testing framework tests the value and derivatives (up to third order) of a function against results derived from the primitive implementation with finite differences.  Multivariate functions are tested as unary functions by projecting to each output dimension.  All that needs to be provided is a working primitive implementation and the inputs to test.

Derivatives are tested

* with user supplied arguments plus all combinations of automatic differentiation variables up to third order:
    - first order: `var` (gradient), `fvar<double>` tested using `gradient()` functionals
    - second order: `fvar<var>`, `fvar<fvar<double>>` tested using `hessian()` functionals
    - third order: `fvar<fvar<var>>` tested using `grad_hessian()` functional

* for all possible output values for functions with multivariate outputs,

* for all possible combinations of argument instantiations with primitive (`int` or `double`) and automatic differentiation variables, and

* up to three levels of container nesting for elementwise vectorized functions.

Exceptions are tested insofar as if one version of the primitive or autodiff versions of the function throws an exception they must all throw an exception  That is, each autodiff call must throw exceptions for exactly the same argument values as the primitive implementation.

### Scope of tests

The framework can test unary and binary functions with any of the following input and output types.  

The basic types that can be handled by the testing framework are

* `int`: integer
* `double`: double-precision floating point
* `Eigen::Matrix<double, -1, 1>`: column vector
* `Eigen::Matrix<double, 1, -1>`: row vector
* `Eigen::Matrix<double, -1, -1>`: matrix

For any type `T` that can be handled, the framework can also handle the corresponding standard vector type (used for arrays in Stan)

* `std::vector<T>`: array of elements of type `T`

### Relative error comparisons

The tests are instrumented to test relative error between the values and derivatives produced by the function being tested and those produced by finite differences.  For non-zero inputs, the relative error is defined to be the error divided by the average absolute values of the arguments:

```cpp
error(u, v) = u - v
relative_error(u, v) = error(u, v) / (0.5 * (abs(u) + abs(v)))
```

If one of the inputs is zero, error is used rather than relative error (because relative error devolves if one or both inputs are zero)

### Tests and tolerances

The following table lists the complete set of tests, with their default tolerances, the type of automatic differentiation variable used, the name of the functional used to perform the test, and the order of the test.  The first column lists the names of the member variables in the simple struct `stan::test::ad_tolerances`.  Every test accepts an `ad_tolerances` object as a first argument, so the defaults may be overridden for particular tests.   The default tolerances are set relatively high because of the inaccuracy of simple finite difference calculations as used by these tests.  

| function                     | tolerance |  autodiff type       | functional     | (order) test            |
|-----------------------------:|:---------:|:--------------------:|:--------------:|:------------------------|
| `gradient_val_`              | 1e-8      | `var`                | `gradient`     | (0) value               |
| `gradient_grad_`             | 1e-4      | `var`                | `gradient`     | (1) gradient            |
| `gradient_fvar_val_`         | 1e-8      | `fvar<double>`       | `gradient`     | (0) value               |
| `gradient_fvar_grad_`        | 1e-4      | `fvar<double>`       | `gradient`     | (1) gradient            |
| `hessian_val_`               | 1e-8      | `fvar<var>`          | `hessian`      | (0) value               |
| `hessian_grad_`              | 1e-4      | `fvar<var>`          | `hessian`      | (1) gradient            |
| `hessian_hessian_`           | 1e-3      | `fvar<var>`          | `hessian`      | (2) Hessian             |
| `hessian_fvar_val_`          | 1e-8      | `fvar<fvar<double>>` | `hessian`      | (0) value               |
| `hessian_fvar_grad_`         | 1e-4      | `fvar<fvar<double>>` | `hessian`      | (1) gradient            |
| `hessian_fvar_hessian_`      | 1e-3      | `fvar<fvar<double>>` | `hessian`      | (2) Hessian             |
| `grad_hessian_val_`          | 1e-8      | `fvar<fvar<var>>`    | `grad_hessian` | (0) value               |
| `grad_hessian_hessian_`      | 1e-3      | `fvar<fvar<var>>`    | `grad_hessian` | (2) Hessian             |
| `grad_hessian_grad_hessian_` | 1e-2      | `fvar<fvar<var>>`    | `grad_hessian` | (3) gradient of Hessian |


### Examples

These examples are all included in the unit tests for the test framework in file `test/unit/math/test_ad_test.cpp`

#### Example: matrix inverse (unary function)

Let's start by seeing how to test matrix inverse.  In Stan, that's a unary function `inverse` that takes a matrix input and returns a matrix output.

```cpp
Matrix<T, -1, -1> inverse(const Matrix<T, -1, 1>& a);
```

We'll need to include the test framework (which includes both the Google test framework and all of the Stan math library (by including `<stan/math.hpp>).  The entire library is included to ensure that all functions work under the standard include system when all higher-order autodiff is included.

Each test is scoped within a Google test framework `TEST` macro.  Thus all of our test files will look like this:

```cpp
#include <test/unit/math/test_ad.hpp>

TEST(test_name, subtest1_name) {
  ... first test goes here...
}

TEST(test_name, subtest2_name) {
  ... second test goes here...
}

... more tests may follow ...
```

If there is more than one test, there will be more than one `TEST` macro invocation, each with a unique combination of test name and subtest name.

The body of a test will define a polymorphic functor to test, then a sequence of arguments to test with thest framework.  Although a functor is necessary in order for there to be a single type encapsulating all possible overloads, it can be a simple closure.

Here's what the first test looks like for `inverse`, repeating the macro wrapper and include

```cpp
#include <test/unit/math/test_ad.hpp>

TEST(inverse, two_by_two) {
  auto inverse_functor = [](const auto& u) {
    return stan::math::inverse(u);
  };

  Eigen::MatrixXd x(2, 2);
  x << 1.9, 0.3,
       0.3, 1.7;

  stan::test::expect_ad(inverse_functor, x);
}
```

That's it.  The `inverse_functor` is defined using a lambda (with no closure --- the `[ ]` bit) taking a `const auto&` argument and returning the result of applying `stan::math::inverse`.  Generic lambdas are a C++14 feature (in `-std=c++1y` for the compiler).  They're convenient here because they define a polymorphic functor defined for any argument for which the body is defined, which here means any argument `u` for which `stan::math::inverse(u)` is defined.  Bare functions aren't enough for the testing framework functions like `expect_ad` to perform type inference as they do not have a single type.  The functor defined by the lambda closure has a single type.

The second bit of the code defines the Matrix to test, here an invertible positive-definite matrix with first row `[1.9  0.3]` and second row `[0.3  1.9]`.  

The final bit of the code simply invokes the test framework for the functor and argument(s) being tested.  The `expect_ad` function is both variadic and overloaded to deal with arbitrary Stan-relevant input and output types.  

#### Example 1a: Adding tolerances

If we wanted to lower the tolerance for the gradient test using reverse-mode autodiff, we could use

```cpp
ad_tolerances tols;  // constructed with default tolerances
tols.gradient_grad_ = 1e-3;
expect_ad(tols, inverse_functor, x);
```

#### Example 2: Binary Function (multiply matrix/scalar)

Binary functions are handled the same way as unary functions.  The framework tests all combinations of input autodiff types, e.g., `f(double, var)`, `f(var, double)`, and `f(var, var)` for testing reverse-mode instantiations.  When there is more than one autodiff variable argument, all autodiff arguments must be of the same type.  

```cpp
#include <test/unit/math/test_ad.hpp>

TEST(multiply_matrix_scalar, test1) {
  auto g = [ ] (const auto& u, const auto& v) {
    return stan::math::multiply(u, v);
  };

  Eigen::MatrixXd x(2, 3);
  x << 1, 2, 3, 4, 5, 6;
  double y = -3;
  stan::test::expect_ad(g, x, y);

  double y2 = std::numeric_limits<double>::quiet_NaN();
  stan::test::expect_ad(g, x, y2);
}
```

The lambda here abstracts two arguments to which it applies the `multiply` function.  There are two tests.  Here The first takes `y = -3` and the second `y = NaN`.  

#### Example 3: vectorized unary scalar function

Most of the unary scalar functions in Stan are vectorized in that they can be applied elementwise to any container type.

The only difference in calls is that the `expect_ad_vectorized` function is used in place of `expect_ad`.  The `log10` function is vectorized to operate elementwise, and the following code will test it.

```cpp
#include <test/unit/math/test_ad.hpp>

TEST(test_unit_math_test_ad, expect_ad_vectorized) {
  auto g = [ ] (const auto& u) {
    return stan::math::log10(u);
  };

  stan::test::expect_ad_vectorized(g, 3.2);  // finite return
  stan::test::expect_ad_vectorized(g, 0);  // infinite return
  stan::test::expect_ad_vectorized(g, -1.7);  // NaN return
}
```

#### Example 4: Mismatched autodiff and double functions

The following function has an autodiff implementation that differs by a negation from the primitive version, throwing off value and derivative tests.

```cpp
#include <test/unit/math/test_ad.hpp>

TEST(test_ad, mismatch) {
  double x = 3.2;
  auto g = [](const auto& u) { return f_mismatch(u); };
  // include following line to show exception error behavior
  stan::test::expect_ad(g, x);
}
```

```cpp
~/cmdstan/stan/lib/stan_math$ ./runTests.py -j4 test/unit/math/test_ad_test.cpp

...

[ RUN      ] test_ad.mismatch
./test/unit/math/mix/mat/util/autodiff_tester.hpp:75: Failure
The difference between x1 and x2 is 12.800000000000001, which exceeds tol, where
x1 evaluates to -6.4000000000000004,
x2 evaluates to 6.4000000000000004, and
tol evaluates to 1.0000000000000001e-09.
expect_near(-6.4000000000000004, 6.4000000000000004)
test_gradient fx = fx_ad

./test/unit/math/mix/mat/util/autodiff_tester.hpp:75: Failure
The difference between x1 and x2 is 4.0000000000000764, which exceeds tol, where
x1 evaluates to -2.0000000000000759,
x2 evaluates to 2, and
tol evaluates to 9.9999999999999995e-08.
expect_near(-2.0000000000000759, 2)
expect_near elt x1(i) = x2(i)
test gradient grad_fd == grad_ad

[  FAILED  ] test_ad.mismatch (1 ms)
```

#### Example 5: Mismatched Exceptions

```cpp
#include <test/unit/math/test_ad.hpp>

template <typename T>
T f_misthrow(const T& x) {
  return -2 * x;
}

double f_misthrow(const double& x) {
  throw std::runtime_error("f_misthrow(double) called");
  return -2 * x;
}

TEST(test_ad, misthrow) {
  double x = 1.73;
  auto h = [](const auto& u) { return f_misthrow(u); };
  stan::test::expect_ad(h, x);
}
```

Uncommenting the test above in the test file and running produces the exception diagnostic:

```cpp
~/cmdstan/stan/lib/stan_math$ ./runTests.py -j4 test/unit/math/test_ad_test.cpp

[ RUN      ] test_ad.misthrow
./test/unit/math/mix/mat/util/autodiff_tester.hpp:234: Failure
Failed
double throws, expect autodiff version to throw
./test/unit/math/mix/mat/util/autodiff_tester.hpp:234: Failure
Failed
double throws, expect autodiff version to throw
./test/unit/math/mix/mat/util/autodiff_tester.hpp:234: Failure
Failed
double throws, expect autodiff version to throw
./test/unit/math/mix/mat/util/autodiff_tester.hpp:234: Failure
Failed
double throws, expect autodiff version to throw
./test/unit/math/mix/mat/util/autodiff_tester.hpp:234: Failure
Failed
double throws, expect autodiff version to throw
[  FAILED  ] test_ad.misthrow (0 ms)
```
