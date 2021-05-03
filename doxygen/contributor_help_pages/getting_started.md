# Getting Started Guide {#getting_started}

This is meant as a basic guide for writing and maintaining functions in the
[stan-dev/math](https://github.com/stan-dev/math) repository. Stan Math is the
automatic differentiation (autodiff) library behind the Stan language, and
the vast majority of the functions exposed at the Stan language level are
implemented here in C++.

In the course of the Math library's existence, C++ has changed substantially.
Math was originally written before C++11. It currently targets C++14. In the near
future it will transition to C++17. With this in mind, there are many different
ways to write Math functions. This guide tries to document best practices,
conventions which not all functions in Math follow, but should be followed for
new code to keep the code from getting unwieldy (the old patterns will be
updated eventually).

The title contains "Current State" to emphasize that if any information here is
out of date or any advice does not work, it should be reported as a bug (the
[example-models](https://github.com/stan-dev/example-models)).

This document builds on the information in the [Stan Math Library](https://arxiv.org/abs/1509.07164) which should be skimmed over before .

## Preliminary Resources:

Before contributing to Stan Math it's recommended you become familiar with C++. While this guide attempts to cover the edges and odd things we do in the math library, we recommend the resources below for getting started on autodiff in C++.

Books:
- [Effective Modern C++](https://www.amazon.com/Effective-Modern-Specific-Ways-Improve/dp/1491903996)
- [Bob's autodiff book](https://github.com/bob-carpenter/ad-handbook/blob/master/ad-handbook-draft.pdf)
- [C++ Templates: The Complete Guide](https://www.amazon.com/C-Templates-Complete-Guide-2nd/dp/0321714121)
- [Numerical Recipes](https://www.amazon.com/Numerical-Recipes-3rd-Scientific-Computing/dp/0521880688)

Talks:
- [Give me 15 minutes & I'll change your view of GDB](https://www.youtube.com/watch?v=PorfLSr3DDI)
- [What Has My Compiler Done for Me Lately?](https://www.youtube.com/watch?v=bSkpMdDe4g4)
- [Going Nowwhere Faster](https://www.youtube.com/watch?v=2EWejmkKlxs)

Blog Posts:
- [Thinking About Automatic Differentiation in Fun New Ways](https://blog.mc-stan.org/2020/11/23/thinking-about-automatic-differentiation-in-fun-new-ways/)

A generally good resource for making minimal examples for bugs is [godbolt.org](godbolt.org)

## Code Structure

The Stan Math library is spit into 4 main source folders that hold functions, type traits, and classes.

- prim: General `Scalar`, `Matrix`, and `std::vector<T>` types
- rev: Specializations for reverse mode automatic differentiation
- fwd: Specializations for forward mode automatic differentiation.
- mix: Sources to allow mixes of forward and reverse mode.
- opencl: Sources for doing reverse mode automatic differentiation on GPUs.

Within each of those folders you will find any one of the following folders

- core: Base implementations of custom scalars or backend setup.
  - Ex: in `prim` this is operators for complex numbers and the setup for threading, `rev`'s core is the scalar and its base operators for reverse mode and the stack allocator, and `fwd` has the scalar for forward mode autodiff and its operators.
- err: Functions that perform a check and if true throw an error.
- fun: The math functions exposed to the Stan language.
- functor: Functions that take in other functions and data as input such as `reduce_sum()`
- meta: Type traits for compile time deduction on types.

Any function callable from the math library, excluding those in the `internal` namespace, should be compatible with Stan's reverse and forward mode autodiff types. Any new function introduced to the Math library is expected to support higher order automatic differentiation. Though exceptions to supporting higher order autodiff have been made in the past such as the differential equations solvers only supporting reverse mode autodiff.

The structure for contributing a function:
- Functions which operate on arithmetic and autodiff types go in `stan/math/prim`.
- Functions specializations for reverse mode go in `stan/math/rev`.
- Functions specializations for forward mode autodiff go in `stan/math/fwd`
- Functions specializations for a combination of both forward and reverse mode go in `stan/mat/mix/fun`.

Adding a function to Stan can be as simple as adding a generic templated function to `prim/fun`. The `rev`, `fwd`, `mix`, and `opencl` folders are used for adding specialization for:

1. Performance
2. Better numerical behavior
4. Utilize particular hardware

## Adding a Simple Example function

A generic function that performs a dot product on itself could be written in `prim/fun` as

```cpp
// stan/math/prim/fun/dot_self.hpp
template <typename EigVec>
inline double dot_self(const EigVec& x) {
  const auto& x_ref = to_ref(x); // (1)
  value_type_t<EigVec> sum_x = 0.0; // (2)
  for (int i = 0; i < x.size(); ++i) {
    sum_x += x_ref.coeff(i) * x_ref.coeff(i); // (3)
  }
  return sum_x;
}
```

This is pretty standard C++ besides (1), (2), and (3) which I'll go over here.

#### (1) Safe elementwise access for Eigen expressions

TL;DR: Before accessing individual coefficients of an Eigen type, use `to_ref()` to make sure it's a type that's safe to access by coefficient.


In the Stan math library we allow functions to accept Eigen expressions. This is rather nice as for instance the code

```cpp
Eigen::MatrixXd x = multiply(add(A, multiply(B, C)), add(D, E));
```
Actually delays evaluation when making `x` and produces an object of type

```cpp
Eigen::Product<Eigen::Add<Matrix, Eigen::Product<Matrix, Matrix>>, Eigen::Add<Matrix, Matrix>>
```

Using lazily evaluated expressions allows Eigen to avoid redundant copies, reads, and writes to our data. However, this comes at a cost.

In (2), when we access the coefficients of `x`, if its type is similar to the wacky expression above we can get incorrect results as Eigen does not guarantee any safety of results when performing coefficient level access on a expression type that transforms its inputs. So `to_ref()` looks at its input type, and if the input type is an Eigen expression that it evaluates the expression so that all the computations are performed and the return object is then safe to access.

But! Suppose our input is something like

```cpp
Eigen::Matrix<var, Dynamic, 1> A = Eigen::VectorXd::Random(20);
var b = dot_self(A.segment(1, 5));
```

Here, we will have an expression `Eigen::Block<Eigen::MatrixXd>` as the input type. But an expression that is only a view into an object is safe to read coefficent-wise. In this case `to_ref()` knows this and, as long as you used `const auto&` as the object type returned from `to_ref(x)`, then `to_ref(x)` will deduce the correct type `Block` instead of casting to an evaluated vector. We use `const auto&` here even when the output will store an object as C++'s [lifetime](https://en.cppreference.com/w/cpp/language/lifetime) rules allow for this.

> The lifetime of a temporary object may be extended by binding to a const lvalue reference or to an rvalue reference (since C++11), see reference initialization for details.

If we used `auto` here then if the return type of `to_ref()` was an `Eigen::MatrixXd` then it would make a copy. But with `const auto&` we do not make a copy the lifetime of the object referenced by the `const auto&` will be the same as if we had done `auto`


#### (2) Using `value_type_t<>` and Friends

The type trait `value_type_t` is one of several type traits we use in the library to query information about types. `value_type_t` will return the inner type of a container,
so `value_type_t<Eigen::Matrix<double, -1, -1>>` will return `double`, `value_type_t<std::vector<std::vector<double>>>` will return a `std::vector<double>` and `value_type_t<double>` will simply return a double.
See the module on type traits for a guide on Stan's specific type traits.

#### (3) Accessing Eigen matrices though `.coeff()` and `.coeffRef()`

This is a small bit, but Eigen performs bounds checks by default when using `[]` or `()` to access elements. However `.coeff()` and `.coeffRef()` do not. Because Stan model's perform bounds checking at a higher level its safe to remove the bounds checks here.

#### Testing

From there you can use the unit [test suite for automatic differentiation](@ref autodiff_test_guide) to verify your code produces the correct gradients and higher order derivatives and your done!

```cpp
// Example tests in test/unit/math/mix/fun/dot_self_test.cpp
TEST(MathMixMatFun, dotSelf) {
  auto f = [](const auto& y) { return stan::math::dot_self(y); };

  Eigen::VectorXd x0(0);
  Eigen::VectorXd x1(1);
  x1 << 2;
  Eigen::VectorXd x2(2);
  x2 << 2, 3;
  Eigen::VectorXd x3(3);
  x3 << 2, 3, 4;
  for (const auto& a : std::vector<Eigen::VectorXd>{x0, x1, x2, x3}) {
    stan::test::expect_ad(f, a);
  }
}
```

That's it, your done! Make a PR, add some docs, have a cup of tea all is well.

But, while your function will work for all of our types, it probably won't be as fast as it could be. Indeed, since the above is using Stan's basic scalar reverse mode, for a vector of size N there will be N individual callback functions on the reverse pass callback stack. This is not the worst, Stan's base scalar operations are as efficient as they can be, but by writing a specialization for reverse mode we can just have one callback with much more efficient memory access. The next section may seem a bit advanced, but those who enjoy adventure, let's talk about how to make this function go faster.

## Specializing for Reverse Mode

At this point I'll assume you have looked over the preliminary resources listed above describing Stan math's automatic differentiation.

For reverse mode we can do the math by hand or use a site like [matrixcalculus.org](http://www.matrixcalculus.org/) to find the adjoint method for a self dot product is

```cpp
x.adj() += sum_x.adj() * 2.0 * x.val()
```

and so we want to add a specialization in `stan/math/rev/fun` that calculates reverse mode's adjoint directly. The methods on `x` above use custom Eigen methods for matrices, vectors, and arrays, `.val()` and `.adj()`. For an `Eigen::Matrix<var, Rows, Cols>` the internal Eigen representation comes down to a struct holding an array of `var` types and information on the row and column sizes.

```
struct DenseStorage {
  var* m_data; // array of Stan's var type
  Eigen::Index m_rows;
  Eigen::Index m_cols;
}
```

It's very common for us to want to access just the `val_` or `adj_` of the `vari`'s inside of the `var`'s and so we have written custom methods `.adj()` and `.val()` using Eigen's plugin system which returns an expression of an `Eigen::Matrix<double, Rows, Cols>` representing the `var`'s values and adjoints, respectively.

#### Using type traits to Expose Our New Function

First we need to stop our current function from compiling for reverse mod autodiff. This can be done using the [requires](@ref require_meta) type traits in stan. The only thing that has changed in the function above and below is the new non-type template parameter ([NTTP](https://en.cppreference.com/w/cpp/language/template_parameters)) that is now in the template definition.

```cpp
template <typename EigVec, require_not_st_var_t<EigVec>* = nullptr> // (1)
inline double dot_self(const EigVec& x) {
  auto x_ref = to_ref(x);
  var sum = 0.0;
  for (int i = 0; i < x.size(); ++i) {
    sum += x_ref.coeff(i) * x_ref.coeff(i);
  }
  return sum;
}
```

TL;DR: Use Stan's `require` type trait library for limiting a function's scope to certain types.

Reading from left to right, `require_not_st_var_t` requires that a signature's template does _not have a scalar type of a `var`_. This line is exactly the same as

```cpp
template <typename EigVec, std::enable_if_t<is_var<scalar_type_t<EigVec>>::value>* = nullptr>
```

but we tend to use these type traits so frequently in the math library it made sense to write specializations for them.

But how does that line work? Note that, on success, `std::enable_if_t` returns a `void` type, and on failure the functions is removed from the set of possible signatures though Substitution Failure is not an Error [SFINAE](https://en.cppreference.com/w/cpp/language/sfinae). So on success the template signature is

```cpp
template <typename EigVec, void* = nullptr>
```

where `void* = nullptr` follows a nice little rule in the C++ standard that says any void pointer non-type template parameters with a default value of a `nullptr` will be ignored by template instantiation. So good! We get exactly what we want at the expense of a some annoying legalize and a bit of an odd `nullptr` style syntax.

So we stopped our original function for compiling for `var` types and now we can add a specialization that does work for `var`.

### Working in Reverse Mode

So now let's write our reverse mode specialization in `rev/fun`.

```cpp
/**
 * Returns the dot product of a vector of var with itself.
 *
 * @tparam EigVec An Eigen type with compile time rows or columns equal to 1 and a scalar var type.
 * @param[in] v Vector.
 * @return Dot product of the vector with itself.
 */
template <typename EigVec, require_eigen_vt<is_var, EigVec>* = nullptr>
inline var dot_self(const EigVec& v) {
  // (1) put v into our arena
  arena_t<EigVec> arena_v(v);
   // calculate forward pass
  var res = arena_v.val().dot(arena_v.val());
  // (2) Place a callback for the reverse pass on the callback stack.
  reverse_pass_callback([res, arena_v]() mutable {
    arena_v.adj().array() += (2.0 * res.adj().array()) * arena_v.val().array();
  });
  return res;
}
```

There's a lot going on here particular to Stan math and so we will talk about each special bit.

### (1) Arena memory

The type `arena_t<T>` and the function `to_arena()` takes an object and places it on Stan's arena allocator.

Remember from [Thinking About Automatic Differentiation in Fun New Ways](https://blog.mc-stan.org/2020/11/23/thinking-about-automatic-differentiation-in-fun-new-ways/) reverse mode automatic differentiation comes down to.

1. Calculating a function z = f(x, y)‘s output (forward pass)
2. Storing the input and output into a memory arena
3. Adding a callback to a callback stack that calculates the adjoint of the function

This process allows us to call `grad()` and from the last in, call each of
the callbacks in the callback stack and accumulates the gradients through
all the outputs and inputs. (reverse pass)

Reverse mode autodiff in Math requires a huge number of temporaries and
intermediates to be computed and saved for the reverse pass. There are
so many allocations that the overhead of `malloc()` becomes noticeable. To
avoid this, Math provides its own memory arena. The assumptions of
the Math arena are:

1. Objects allocated in the arena are [TriviallyCopyable](https://en.cppreference.com/w/cpp/named_req/TriviallyCopyable) and [TriviallyDestructable](https://en.cppreference.com/w/cpp/language/destructor#Trivial_destructor). These objects do not need their destructors
called (although there are ways to destruct something allocated in the arena,
it just isn't automatic).
2. All objects allocated on the arena have the same lifetime (that is the
lifetime of the stack - until `recover_memory()` is called).

The arena memory pattern fits nicely into automatic differentiation for MCMC
as we will most commonly be calling gradients with data of the same size over
and over again.

Functions that use reverse mode autodiff are allowed to save what they need to
the arena and then the objects saved here will be available during the reverse
pass.

There is a utility for saving objects into the arena that do need their
destructors called. This is discussed under `make_callback_ptr()`. Arena types
are helper functions and types for easily saving variables to the arena.

`arena_t<T>` defines, for a given type `T`, the type of a lightweight
object pointing to memory in the arena that can store a `T`. Memory for
`arena_t<T>` objects is only recovered when `recover_memory` is called.
No destructors are called for objects stored with `arena_t<T>()`.

Copying a variable `x` to the arena can be done with either `to_arena()` or
`arena_t<T>`:

```cpp
template <typename T>
auto myfunc(const T& x) {
  auto arena_x = to_arena(x);
  //...
}
```

```cpp
template <typename T>
auto myfunc(const T& x) {
  arena_t<T> arena_x = x;
  // ...
}
```

## (2) Setting up the Reverse Pass

Once we have stored the data we need for the reverse pass we need to actually write that reverse pass! We need to take our adjoint calculation and put it onto the callback stack so that when the users call `grad()` the adjoints are propagated upwards properly.

For this we have a function called `reverse_pass_callback()`. Calling `reverse_pass_callback()` with a functor `f` creates an object on the callback stack that will call `f`.

```cpp
reverse_pass_callback([res, arena_v]() mutable { // (2)
  v.adj() += (2.0 * res.adj()) * v.val();
});
```

Notice that the lambda we create *copies* the objects `res` and `arena_v`. The trick here is that, because anything allocated with Stan's arena using `arena_t` or `to_arena()` are [TriviallyCopyable](https://en.cppreference.com/w/cpp/named_req/TriviallyCopyable) (as are arithmetic scalars). So the memory is safe and sound within our arena and so these are safe and fast to copy.

An alternative approach to the above function could also have been

```cpp
template <typename EigVec, require_eigen_vt<is_var, EigVec>* = nullptr>
inline var dot_self(const EigVec& v) {
  arena_t<EigVec> arena_v(v); // (1) put v into our arena
  return make_callback_var(arena_v.val().dot(arena_v.val()), [arena_v](auto& vi) {
    arena_v.adj() += (2.0 * vi.adj()) * arena_v.val();
  });
}
```

`make_callback_var(x, f)` is similar to `reverse_pass_callback()`, but it constructs
an autodiff type from the argument `x`. The function `f` is called on the
reverse pass similarly to `reverse_pass_callback()`. `make_callback_var()` should
only be used when returning back a `var<double>` or `var<Matrix>`.

## Argument Types

The functions in Stan Math are used in the C++ code generated by the Stan compiler. This
means that all types in the Stan language map to `C++` types. The map between Stan types and types in `C++` is a one to many
relationship because of the different blocks in a Stan program. Objects created in `data`, `transformed data`, and `generated quantities` will hold arithmetic types while objects created in `parameters`, `transformed parameters`, and the `model` block will hold autodiff types. The basic Stan types are listed in the table below as well as arrays of any of these types.

| Stan Type           | Data                                                              | Model                                                                 |
|---------------------|-------------------------------------------------------------------|-----------------------------------------------------------------------|
| `int`               | `int`                                                             | `int`                                                                 |
| `int[N]`            | `std::vector<int>`                                                | `std::vector<int>`                                                    |
| `real`              | `double`                                                          | `var/fvar<T>`                                                         |
| `real[N]`           | `std::vector<double>`                                             | `std::vector<var/fvar<T>>`                                            |
| `vector[N]`         | `Eigen::Matrix<double, Dynamic, 1>(N)`                            | `Eigen::Matrix<var/fvar<T> Dynamic, 1>(N)`                            |
| `row_vector[N]`     | `Eigen::Matrix<double, 1, Dynamic>(N)`                            | `Eigen::Matrix<var/fvar<T> 1, Dynamic>(N)`                            |
| `matrix[N,M]`       | `Eigen::Matrix<double, Dynamic, Dynamic>(N, M)`                   | `Eigen::Matrix<var/fvar<T> Dynamic, Dynamic>(N, M)`                   |
| `matrix[N, M] x[K]` | `std::vector<Eigen::Matrix<double, Dynamic, Dynamic>>(K, {M, N})` | `std::vector<Eigen::Matrix<var/fvar<T> Dynamic, Dynamic>>(K, {M, N})` |


Any array type (`int[]`, `real[]`, or `T[]` for `T` any of the above types)
map to `std::vector<C>` types where `C` is the C++ equivalent of `T`.

### Unary Functions Working Over Arrays

The Stan function `sin()` expects the math library has the following signatures available.

```cpp
real sin(real);
vector sin(vector);
array vector sin(array vector)
```

where `array` can actually be multiple arrays within one another. The above translates to the following C++ signatures

```cpp
// reals
double sin(double);
var sin(var);
// vector
template <typename T>
Eigen::Vector<T, Dynamic, 1> sin(Eigen::Vector<T, Dynamic, 1>);
// array vector
template <typename T>
std::vector<Eigen::Vector<T, Dynamic, 1>> sin(std::vector<Eigen::Vector<T, Dynamic, 1>>);
// array array vector
template <typename T>
std::vector<std::vector<Eigen::Vector<T, Dynamic, 1>>> sin(std::vector<std::vector<Eigen::Vector<T, Dynamic, 1>>>);
```

A function that accepts an autodiff type should return an autodiff type. Notice that because we have to support arithmetic and autodiff types the signatures expand out for both arithmetic and autodiff types. It would be tedious and grueling to write out this method for all of the array types, so Stan math provides helper classes `apply_scalar_unary` and `apply_vector_unary` to apply a function over each array.

```cpp
/**
 * Structure to wrap sin() so it can be vectorized.
 *
 * @tparam T type of argument
 * @param x angle in radians
 * @return Sine of x.
 */
struct sin_fun {
  template <typename T>
  static inline T fun(const T& x) {
    using std::sin;
    return sin(x);
  }
};

/**
 * Vectorized version of sin().
 *
 * @tparam Container type of container
 * @param x angles in radians
 * @return Sine of each value in x.
 */
template <typename T, require_not_container_st<std::is_arithmetic, T>* = nullptr>
inline auto sin(const T& x) {
  return apply_scalar_unary<sin_fun, T>::apply(x);
}

/**
 * Version of sin() that accepts std::vectors containing Scalars or Eigen Matrix/Array types.
 *
 * @tparam Container Type of x
 * @param x Container
 * @return Sine of each value in x.
 */
template <typename Container, require_container_st<std::is_arithmetic, Container>* = nullptr>
inline auto sin(const Container& x) {
  return apply_vector_unary<Container>::apply(
      x, [&](const auto& v) { return v.array().sin(); });
}
```

In Stan math, a `container` type is an `std::vector` that holds either other `std::vector`s or `Eigen::Matrix` types. So in the first function we use `apply_scalar_unary` to apply `sin()` to either `Scalar`s, `std::vector`s holding scalars, or  `Eigen::Matrix` types. The second function which uses `apply_vector_unary()` will apply its lambda to the container's whose elements are also containers or Eigen matrices.


### Binary Functions (and least-upper-bound return type semantics)

Multiple argument functions are more involved. Consider
`real atan2(real, real)` which requires four C++ overloads:

```cpp
double atan2(double, double);
var atan2(double, var);
var atan2(var, double);
var atan2(var, var);
```

In this case, the return types are a function of both input types.
`return_type_t` is actually written to take any number of input types and
compute from them a single scalar return type. The return type is the simplest
type that can be constructed from all of the input types. For reverse mode
autodiff, this means that if the scalar type of any input argument is a
`var`, then the return type is a var.

A generic signature for `real atan2(real, real)` is:

```cpp
template <typename T, typename S>
return_type_t<T, S> atan2(T, S);
```

`return_type_t` can be used to construct non-scalar return types as well.
For the function `vector add(vector, vector)`, the following generic
C++ signature could be used:

```cpp
template <typename T, typename S>
Eigen::Matrix<return_type_t<T, S>, Eigen::Dynamic, 1> atan2(T, S);
```

Alternativly, with C++14 one can also delay deducing the return type with `auto` and the compiler can deduce the return type from the code (which would then internally use `return_type_t<T, S>`)

### Higher Order Autodiff

Adding the forward mode autodiff type simply adds more overloads. The forward mode
autodiff type itself takes a template argument to represent the types of the
value and gradient and allows recursive templating. This technically defines
an infinite number of new types, but the ones of interest in Math (and the
ones that should be tested are) `fvar<double>`, `fvar<fvar<double>>`,
`fvar<var>` and `fvar<fvar<var>>`. These are useful for computing various
high order derivatives. Going back to the `real sin(real)` function, Math
is now expected to implement the following, rather expanded list of functions:

```cpp
double sin(double);
var sin(var);
fvar<double> sin(fvar<double>);
fvar<fvar<double>> sin(fvar<fvar<double>>);
fvar<var> sin(fvar<var>);
fvar<fvar<var>> sin(fvar<fvar<var>>);
```

For the sake of performance, it may be desirable to also define `var sin(var)`,
or similarly `fvar<double> sin(fvar<double>)` etc.

`return_type_t` is defined similarly as for `var`. In general, autodiff types
should not be mixed, and so the `return_type_t` does not need to account
for various combinations of `var`, `fvar<double>`, etc. Sometimes it is
useful to mix autodiff types, but it is somewhat uncommon.

### Exceptions

Math functions are expected to validate inputs and throw exceptions.
It should be impossible to cause a segfault by calling a
user-facing Math function or trigger an exception that is not handled
by Stan itself. This means that argument sizes should be checked to be
compatible. For instance, a matrix multiply function would need to check
that the first matrix can in fact be multiplied by the second, if this
is not so it should throw.

If a scalar argument should be positive, an exception should be thrown
if it is not. If a matrix argument should be positive definite, the
function should check this before operating on the matrix. This can
lead to significant overhead in the checks, but the current library
policy is to get informative errors back to the user, which means doing
these checks.

Error checks should be done with the functions in
[stan/math/prim/err](https://github.com/stan-dev/math/tree/develop/stan/math/prim/err)
which throw consistently formatted errors.

### Tests

Tests in Math are written with the Google test framework. At this
point there are enough functions in Math that the way to get started
is by just copying the tests of a similar Function and editing them to
suit the new function. There are two basic sets of functions that every
Stan function should have. First, for a function named `foo`, there
should be a file in `test/unit/math/prim/fun` named `foo_test.cpp`.

This should include tests that:

1. Check that every exception is throwing in the right way at least once
and that all edge cases are handled correctly.

  For instance, if there is a `check_size_match()` function, then there
  should be a test that makes sure that if the sizes don't match, this
  error gets thrown. It is not necessary to check the error message in
  the tests, just check that an exception is thrown and the exception
  type is correct.

2. Check that the function returns the expect values for a couple
suitable inputs.

  For some functions these checks should be more extensive that others,
  but there should at least be a two or three even in simple cases.

The `prim` tests should only have arithmetic types. No autodiff
is tested here. The `prim` tests check that the function handles edge
cases correctly and produces the correct output for different function
values. The second set of tests should be in `stan/math/mix/fun` with
the same name `foo_test.cpp`. These tests will be used to test all
forms of the autodiff.

The tests in `mix/fun` should use the `expect_ad()` function to test
all forms of autodiff. `expect_ad()` takes as its first argument a
function to test and then up to three other arguments, all with
arithmetic scalar types. For a given function `f`, and two arguments
`a` and `b`, `expect_ad()` tests that:

1. For every combination of `a` and `b`, the returned value of `f`
matches the returned value with non-autodiff types.

2. The finite difference gradients of `f` match those computed with
autodiff.

3. Any time a `f` throws with arithmetic arguments, it also throws
with autodiff arguments.

## Reverse Mode

### Values

`value_of(x)` returns the values of the variable `x`. If `x` is not an
autodiff type, it simply forwards `x` along. The values of `x` have
the same shape as `x`. For reverse mode autodiff, values are `double`
variables. For higher order autodiff (`fvar<var>`, etc.) this is not
always the case. See `value_of_rec()` for a different behavior.

### Values and adjoint extensions to Eigen

The AOS matrix and vector autodiff types come with an extra `.val()`,
`.val_op()`, `.adj()`, and `.adj_op()` member functions.

`.val()` and `.val_op()` return expressions that evaluate to the values
of the autodiff matrix.

`.adj()` and `.adj_op()` return expressions that evaluate to the adjoints
of the autodiff matrix.

The non-`_op` versions of the functions should be preferred. If the
code does not compile (usually an error about const-ness) try the `_op`
versions. This can happen when multiplying the values/adjoints of an AOS
autodiff type against other matrices.

It is not safe to use `.noalias()` with the AOS value/adjoint expressions. In Eigen, `.noalias()` is similar to the keyworld [`restrict`](https://en.cppreference.com/w/c/language/restrict) which assumes that a pointer is unique within a given function. But when performing operations on Eigen matrices of `var` while using `.val()` and `.adj()` on the left and right hand side of an operation, the pointer holding the underlying `var` types will be used on both sides. Since the pointer is used on both sides it will not be unique in the operation and so using `.noalias()` will lead to undefined behavior. See the [Aliasing](https://eigen.tuxfamily.org/dox/group__TopicAliasing.html) and [Writing Efficient Matrix Product Expressions](https://eigen.tuxfamily.org/dox/TopicWritingEfficientProductExpression.html) Eigen docs for more info.
