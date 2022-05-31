# Adding New Functions Guide {#getting_started}

This is meant as a basic guide for writing and maintaining functions in the
[stan-dev/math](https://github.com/stan-dev/math) repository. Stan Math is the
automatic differentiation (autodiff) library behind the Stan language, and the
math functions exposed at the Stan language level are implemented here in C++.

In the course of the Math library's existence, C++ has changed substantially.
Math was originally written before C++11. It currently targets C++14. In the near
future it will transition to C++17. With this in mind, there are many different
ways to write Math functions. This guide tries to document best practices,
conventions which not all functions in Math follow, but should be followed for
new code to keep the code from getting unwieldy (the old patterns will be
updated eventually).

The title contains "Current State" to emphasize that if any information here is
out of date or any advice does not work, it should be reported as a bug (
[issue tracker link](https://github.com/stan-dev/math/issues/new)).

This document builds on the information in the [Stan Math Library](https://arxiv.org/abs/1509.07164). While some parts of the paper are out of date now it is very useful for understanding the patterns used in the math library.

## Preliminary Resources:

Before contributing to Stan Math it's recommended you become familiar with C++. While this guide attempts to cover the edges and odd things we do in the math library, we recommend the resources below for getting started on autodiff in C++.

Books and Papers:
- [Effective Modern C++](https://www.amazon.com/Effective-Modern-Specific-Ways-Improve/dp/1491903996)
- [C++ Templates: The Complete Guide](https://www.amazon.com/C-Templates-Complete-Guide-2nd/dp/0321714121)

Talks:
- [Give me 15 minutes & I'll change your view of GDB](https://www.youtube.com/watch?v=PorfLSr3DDI)
- [What Has My Compiler Done for Me Lately?](https://www.youtube.com/watch?v=bSkpMdDe4g4)
- [Going Nowwhere Faster](https://www.youtube.com/watch?v=2EWejmkKlxs)

Blog Posts and Articles:
- [Thinking About Automatic Differentiation in Fun New Ways](https://blog.mc-stan.org/2020/11/23/thinking-about-automatic-differentiation-in-fun-new-ways/)
- [AD HandBook](https://github.com/bob-carpenter/ad-handbook/blob/master/ad-handbook-draft.pdf)
- [Collected matrix derivative results for forward
and reverse mode AD](https://people.maths.ox.ac.uk/gilesm/files/AD2008.pdf)

A generally good resource for making minimal examples for bugs is [godbolt.org](https://godbolt.org/). The godbolt link [here](https://godbolt.org/z/GEhdcz8eo) ([backup gist](https://gist.github.com/SteveBronder/803227cccc66035f7b1988d0156e562a)) contains a function for printing C++ types and has Eigen and boost headers included.


The example function section is an amalgamation of the @ref common_pitfalls section. It's recommended you skim over the common pitfalls before going over the examples below.

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
- meta: Type traits and helper functions used internally.

Any function callable from the math library, excluding those in the `internal` namespace, should be compatible with Stan's reverse and forward mode autodiff types. Any new function introduced to the Math library is expected to support higher order automatic differentiation. Though exceptions to supporting higher order autodiff have been made in the past such as the differential equations solvers only supporting reverse mode autodiff.

You can contribute a new function by writing the function in the `prim/fun` folder. Though the function will need to accept arithmetic, reverse mode (@ref stan::math::var  ), and forward mode (`fvar<T>`) scalar types and should be tested using the `expect_ad()` testing framework. The `rev`, `fwd`, `mix`, and `opencl` folders are used for adding specialization for:

1. Performance
2. Better numerical behavior of gradient calculations
4. Utilize particular hardware

## Adding a Simple Example function

A generic function that takes in an Eigen vector and performs a dot product on itself could be written in `prim/fun` as

```cpp
// stan/math/prim/fun/dot_self.hpp
template <typename EigVec, require_eigen_vector_t<EigVec>* = nullptr> // (1)
inline value_type_t<EigVec> dot_self(const EigVec& x) {
  const auto& x_ref = to_ref(x); // (2)
  value_type_t<EigVec> sum_x = 0.0; // (3)
  for (int i = 0; i < x.size(); ++i) {
    sum_x += x_ref.coeff(i) * x_ref.coeff(i); // (4)
  }
  return sum_x;
}
```

Do keep in mind that in reality we would just want to return `x.squaredNorm()` on the input vector `x`. But unrolling this `dot_self()` example gives us a nice example with a few gotchas we need to watch for.

There are several odd parts to this function labeled (1), (2), (3), and (4) which I'll go over below.

#### (1) Using Stan's require type traits to specialize overloads.

The Stan math library requires that all functions are able to accept Eigen expression templates, hence the general `EigVec` template above. Though we want general templates parameters to the functions, we also want to ensure that only a restricted subset of types are able to be accepted for a particular function. For example, without the require in the above function it would be perfectly valid to pass in an `Eigen::MatrixXd`, but that would not be good! In order to restrict the types that can go into a function the math library uses the `require` type traits to detecting which function a particular type should use. This is very similar conceptually to [C++20 requires](https://en.cppreference.com/w/cpp/language/constraints#Constraints).

To avoid document duplication I recommend you skim over @ref require_meta_doc The guide there explains the different styles and patterns of requires that we use throughout the math library as well as in this example. At the bottom of that page you will see subsections which will show you the `requires` we have already setup for types we commonly see.

We want to restrict `self_dot()` to only accept Eigen types which have one row or one column at compile time, which for us is `require_eigen_vector_t`. A very odd piece of the template syntax is how we set our require to be a pointer template parameter with a default value of nullptr.

```cpp
template <typename EigVec, require_eigen_vector_t<EigVec>* = nullptr>
```

One C++ trick is that any non-type template parameter of type `void*` with a default of `nullptr` will be ignored by the compiler. Under the hood, if any of the `requires` are successful they will return back a type `void`. So in the case we are passing a type that our requires accepts we get back a `void* = nullptr` which is safely ignored by the compiler. In the case that the type does not satisfy the `require` then the function is removed from the set of possible functions the caller could use via SFINAE. So with this scheme we end up having a very nice pattern for writing generic templates for functions while also being able to restrict the set of types that a function can be used for

#### (2) Ensure Eigen Matrices are Only Evaluated Once

TL;DR: Before accessing individual coefficients of an Eigen type, use @ref stan::math::to_ref to make sure accessing the coefficients will not cause a compiler error or forcing the coefficient to be computed multiple times.

In the Stan math library we allow functions to accept Eigen expressions. This can lead to some nice outcomes such as in the below example. Here the additions that go into `subtract()` can all be combined together to form one large computation instead of many smaller ones.

```cpp
Eigen::MatrixXd x = subtract(add(A, multiply(B, C)), add(D, E));
```

The type created on the RHS of the above is actually an eigen expression that looks like the below.

```cpp
Eigen::Subtract<Eigen::Add<Matrix, Eigen::Product<Matrix, Matrix>>, Eigen::Add<Matrix, Matrix>>
```

Using lazily evaluated expressions allows Eigen to avoid redundant copies, reads, and writes to our data. However, this comes at the cost of code complication time and being careful around handling inputs to functions.

In (3), when we access the coefficients of `x`, if its type is similar to the wacky expression above we can get incorrect results as Eigen does not guarantee any safety of results when performing coefficient level access on a expression type that transforms its inputs. So @ref stan::math::to_ref  looks at its input type, and if the input type is an Eigen expression that performs a calculation then it returns as if we called `.eval()` on our input, otherwise it returns the same type as `x`.

Another example to show the need for @ref stan::math::to_ref  is in the below (which you can copy and paste into godbolt.org to try yourself)

```cpp
#include <Eigen/Dense>
#include <string_view>
#include <iostream>

template <typename T>
double dot_self(const T& A) {
  double sum_a = 0;
  for (int i = 0; i < A.size(); ++i) {
    sum_a += A.coeff(i) * A.coeff(i);
  }
  return sum_a;
}

int main() {
  Eigen::MatrixXd A_mat = Eigen::MatrixXd::Random(5, 5);
  Eigen::VectorXd A_vec = Eigen::VectorXd::Random(5);
  // This works fine
  //double b2 = dot_self((A_mat * A_vec).eval());
  // This fails
  double b2 = dot_self(A_mat * A_vec);
  std::cout << b2;
}
```

Pasting the above into godbolt will show a compilation error! That's because Eigen does now allow coefficient levels access to product functions. Instead that expression needs to be evaluated before it is accessed elementwise. Alternatively, say we have something such as the below where we pass in an expression for elementwise multiplication.

```cpp
int main() {
  Eigen::VectorXd A_vec1 = Eigen::VectorXd::Random(5);
  Eigen::VectorXd A_vec2 = Eigen::VectorXd::Random(5);
  double b2 = dot_self(A_vec1.array() * A_vec2.array());
  std::cout << b2;
}
```


This will appear fine, but in our call to `dot_self`, when we call `A.coeff(i)` twice in the for loop we will end up evaluating the expression twice.

Now Suppose our input is something like

```cpp
Eigen::Matrix<var, Dynamic, 1> A = Eigen::VectorXd::Random(20);
var b = dot_self(A.segment(1, 5));
```

Here, we will have an expression `Eigen::Block<Eigen::MatrixXd>` as the input type. But an expression that is only a view into an object is safe to read coefficent-wise. In this case @ref stan::math::to_ref  knows this and, as long as you used `const auto&` as the object type returned from `to_ref(x)`, then `to_ref(x)` will deduce the correct type `Block` instead of casting to an evaluated vector. We use `const auto&` here even when the output will store an object as C++'s [lifetime](https://en.cppreference.com/w/cpp/language/lifetime) rules allow for this.

> The lifetime of a temporary object may be extended by binding to a const lvalue reference or to an rvalue reference (since C++11), see reference initialization for details.

If we used `auto` here then if the return type of @ref stan::math::to_ref  was an `Eigen::MatrixXd` then it would make a copy. But with `const auto&` we do not make a copy the lifetime of the object referenced by the `const auto&` will be the same as if we had done `auto`


#### (3) Using value_type_t and Friends

The type trait @ref stan::value_type_t  is one of several type traits we use in the library to query information about types. @ref stan::value_type_t  will return the inner type of a container,
so `value_type_t<Eigen::Matrix<double, -1, -1>>` will return `double`, `value_type_t<std::vector<std::vector<double>>>` will return a `std::vector<double>` and `value_type_t<double>` will simply return a double.

See the docs for @ref stan::value_type, @ref stan::scalar_type and @ref stan::base_type for more information on these.

#### (4) Accessing Eigen matrices though `.coeff()` and `.coeffRef()`

This is a small bit, but Eigen performs bounds checks by default when using `[]` or `()` to access elements. However `.coeff()` and `.coeffRef()` do not. Because Stan model's perform bounds checking at a higher level its safe to remove the bounds checks here. Note that `.coeff()` should be used when accessing values of an Eigen types and `.coeffRef()` should be used when assigning values of an Eigen type.

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

#### Using type traits to Expose Our New Function

First we need to stop our current function from compiling for reverse mod autodiff. This can be done using the [requires](@ref require_meta) type traits in stan. The only thing that has changed in the function above and below is the require on the function.

```cpp
template <typename EigVec, require_eigen_vector_t<EigVec>* = nullptr,
 require_not_st_var<EigVec>* = nullptr>
inline double dot_self(const EigVec& x) {
  auto x_ref = to_ref(x);
  var sum = 0.0;
  for (int i = 0; i < x.size(); ++i) {
    sum += x_ref.coeff(i) * x_ref.coeff(i);
  }
  return sum;
}
```

The require in the above stops the function from working specifically with Eigen vectors whose inner `Scalar` type is a @ref stan::math::var type. Now we can write the signature of our function that will go into the `stan/math/rev/fun/` folder

```cpp
template <typename EigVec, require_eigen_vector_vt<is_var, EigVec>* = nullptr>
inline double dot_self(const EigVec& x) {
  // ... Will fill this out in a moment
}

```

The above signature is set to only work with Eigen vector types whose @ref stan::value_type is @ref stan::math::var .


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
  // (2) calculate forward pass using (3) the .val() method for matrices of var types  
  var res = dot_self(arena_v.val());
  // (4) Place a callback for the reverse pass on the callback stack.
  reverse_pass_callback([res, arena_v]() mutable {
    arena_v.adj().array() += (2.0 * res.adj().array()) * arena_v.val().array();
  });
  return res;
}
```

There's a lot going on here particular to Stan math and so we will talk about each special bit.

### (1) Arena memory

The type `arena_t<T>` and the function @ref stan::math::to_arena takes an object and places it on Stan's arena allocator.

Remember from [Thinking About Automatic Differentiation in Fun New Ways](https://blog.mc-stan.org/2020/11/23/thinking-about-automatic-differentiation-in-fun-new-ways/) reverse mode automatic differentiation comes down to.

1. Calculating a function z = f(x, y)‘s output (forward pass)
2. Storing the input and output into a memory arena
3. Adding a callback to a callback stack that calculates the adjoint of the function

This process allows us to call @ref stan::math::grad , call each of
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
lifetime of the stack - until @ref stan::math::recover_memory is called).

The arena memory pattern fits nicely into automatic differentiation for MCMC
as we will most commonly be calling gradients with data of the same size over
and over again.

Functions that use reverse mode autodiff are allowed to save what they need to
the arena and then the objects saved here will be available during the reverse
pass.

There is a utility for saving objects into the arena that do need their
destructors called. This is discussed under @ref stan::math::make_chainable_ptr . Arena types
are helper functions and types for easily saving variables to the arena.

@ref stan::arena_t defines, for a given type `T`, the type of a lightweight
object pointing to memory in the arena that can store a `T`. Memory for
@ref stan::arena_t objects is only recovered when @ref stan::math::recover_memory is called.
No destructors are called for objects stored with @ref stan::arena_t.

Copying a variable `x` to the arena can be done with either @ref stan::math::to_arena or
@ref stan::arena_t :

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

## (2) Calculating the forward pass

This one does not require too much explanation. The only thing I would like to call out here is that because of our previous `requires` we can now call `dot_self()` using the values of the matrix which will then know to call the previous definition of `dot_self` we wrote that can work with `double` scalar types.

## (3) Accessing the Values and Adjoints of Matrices.

The calculation for forward and reverse pass in the above uses custom Eigen methods for matrices, vectors, and arrays. In the above, `.val()` will return a view of all of the values of the @ref stan::math::var scalars in the matrix while `.adj()` will return a view for all of the adjoints of the matrices @ref stan::math::var scalars. You can think of `.val()` and `.adj()` as returning an Eigen matrix of doubles such as `Eigen::Matrix<double, Rows, Cols>` where both functions's scalars are references to the values or adjoints of the original matrix with @ref stan::math::var scalars.

## (4) Setting up the Reverse Pass

Once we have stored the data we need for the reverse pass we need to actually write that reverse pass! We need to take our adjoint calculation and put it onto the callback stack so that when the users call @ref stan::math::grad the adjoints are propagated upwards properly.

For reverse mode we can do the math by hand or use a site like [matrixcalculus.org](http://www.matrixcalculus.org/) to find the adjoint method for a self dot product is

```cpp
x.adj() += sum_x.adj() * 2.0 * x.val()
```

Once we know the adjoint calculation we can use @ref stan::math::reverse_pass_callback to actually place the adjoint calculation onto the callback stack for the reverse pass. Calling @ref stan::math::reverse_pass_callback with a lambda will place the lambda on Stan's callback stack so that it can be calculated in the reverse pass.

```cpp
reverse_pass_callback([res, arena_v]() mutable { // (2)
  v.adj() += (2.0 * res.adj()) * v.val();
});
```

Notice that for the lambda we create *copies* of the objects `res` and `arena_v`. The trick here is that anything allocated with Stan's arena using @ref stan::arena_t or @ref stan::math::to_arena are [TriviallyCopyable](https://en.cppreference.com/w/cpp/named_req/TriviallyCopyable) (as are arithmetic scalars). So the memory is safe and sound within our arena and so these are safe and fast to copy.

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

`make_callback_var(x, f)` is similar to @ref stan::math::reverse_pass_callback , but it constructs
an autodiff type from the argument `x`. The function `f` is called on the
reverse pass similarly to @ref stan::math::reverse_pass_callback . @ref stan::math::make_callback_var should
only be used when returning back a `var_value<double>` or `var_value<Matrix>`.

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

A function that accepts an autodiff type should return an autodiff type. Notice that because we have to support arithmetic and autodiff types the signatures expand out for both arithmetic and autodiff types. It would be tedious and grueling to write out this method for all of the array types, so Stan math provides helper classes @ref stan::math::apply_scalar_unary and @ref stan::math::apply_vector_unary to apply a function over each array.

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

In Stan math, a `container` type is an `std::vector` that holds either other `std::vector`s or `Eigen::Matrix` types. So in the first function we use @ref stan::math::apply_scalar_unary to apply `sin()` to either `Scalar`s, `std::vector`s holding scalars, or  `Eigen::Matrix` types. The second function which uses @ref stan::math::apply_vector_unary will apply its lambda to the container's whose elements are also containers or Eigen matrices.


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
@ref stan::return_type_t is actually written to take any number of input types and
compute from them a single scalar return type. The return type is the simplest
type that can be constructed from all of the input types. For reverse mode
autodiff, this means that if the scalar type of any input argument is a
@ref stan::math::var , then the return type is a var.

A generic signature for `real atan2(real, real)` is:

```cpp
template <typename T, typename S>
return_type_t<T, S> atan2(T, S);
```

@ref stan::return_type_t can be used to construct non-scalar return types as well.
For the function `vector add(vector, vector)`, the following generic
C++ signature could be used:

```cpp
template <typename T, typename S>
Eigen::Matrix<return_type_t<T, S>, Eigen::Dynamic, 1> atan2(T, S);
```

Alternatively, with C++14 one can also delay deducing the return type with `auto` and the compiler can deduce the return type from the code (which would then internally use `return_type_t<T, S>`)

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

`return_type_t` is defined similarly as for @ref stan::math::var . In general, autodiff types
should not be mixed, and so the @ref stan::return_type_t does not need to account
for various combinations of @ref stan::math::var , `fvar<double>`, etc. Sometimes it is
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

  For instance, if there is a @ref stan::math::check_size_match function, then there
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

@ref stan::math::value_of returns the values of the variable `x`. If `x` is not an
autodiff type, the values are simply the input and the input is forwards along. For reverse mode autodiff, values are `double`
variables. For higher order autodiff (`fvar<var>`, etc.) this is not
always the case. See `value_of_rec()` for a different behavior. The values of `x` have
the same shape as `x`.

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

It is not safe to use `.noalias()` with the AOS value/adjoint expressions. In Eigen, `.noalias()` is similar to the keyworld [`restrict`](https://en.cppreference.com/w/c/language/restrict) which assumes that a pointer is unique within a given function. But when performing operations on Eigen matrices of @ref stan::math::var   while using `.val()` and `.adj()` on the left and right hand side of an operation, the pointer holding the underlying @ref stan::math::var   types will be used on both sides. Since the pointer is used on both sides it will not be unique in the operation and so using `.noalias()` will lead to undefined behavior. See the [Aliasing](https://eigen.tuxfamily.org/dox/group__TopicAliasing.html) and [Writing Efficient Matrix Product Expressions](https://eigen.tuxfamily.org/dox/TopicWritingEfficientProductExpression.html) Eigen docs for more info.
