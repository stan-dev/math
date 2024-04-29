## Common pitfalls {#common_pitfalls}



### Eigen expressions (and handling argument types with SFINAE)

The additional complexity of the Math library is that `vector`, `row_vector`, and `matrix` do not map to just templated `Eigen::Matrix` types, they can also map to a variety of Eigen expressions.

For instance, in the following code, the result `c` is not a vector but a vector-like type representing the not-yet-calculated sum of `a` and `b`:

```cpp
Eigen::VectorXd a(5);
Eigen::VectorXd b(5);
// auto is actually a
// Eigen::CwiseBinaryOp<Eigen::internal::scalar_sum_op<double, double>, const Eigen::Matrix<double, -1, 1>, const Eigen::Matrix<double, -1, 1>>
auto c = a + b;
```

Similarly, there are other representations of vectors and matrices yet discussed, in particular matrices stored in Math's special arena memory or matrices stored on a GPU.
In either case, it is expedient to not write explicit overloads for all the different types that a function might accept but limit them with template metaprograms.

For instance, if only base Eigen types are allowed, then a function that takes column vector types could be defined as follows:

```cpp
template <typename T>
T norm(const Eigen::Matrix<T, Eigen::Dynamic, 1>&);
```

A typical problem with a function like this is that `norm` can similarly be defined for a row vector (and the implementation is the same).
In this case there are two signatures:

```cpp
template <typename T>
T norm(const Eigen::Matrix<T, Eigen::Dynamic, 1>&);
template <typename T>
T norm(const Eigen::Matrix<T, 1, Eigen::Dynamic>&);
```

The simple solution to this problem is to move to a fully templated signature:

```cpp
template <typename T>
return_type_t<T> norm(const T&);
```

But what if we want different overloads to handle different types?
The most common situation is when there is one template overload that works with any autodiff type, and then a second faster overload that only works with a specific autodiff type.

There can easily be ambiguities between the two function signatures.
The previous examples took advantage of simple overloads.
The more general solution is template metaprograms that additionally make use of C++ substitution failure
is not an error (SFINAE) semantics.
For the norm function above, SFINAE could be used to limit one signature to work with reverse
mode autodiff types and one to work with anything else:

```cpp
// For Eigen vector types with var Scalar types
template <typename T,
          require_eigen_vector_vt<is_var, T>* = nullptr>
return_type_t<T> norm(const T&);
// For Eigen vector types with a non var scalar type.
template <typename T,
          require_eigen_vector_t<T>* = nullptr,
          require_not_vt_var<T>* = nullptr>
return_type_t<T> norm(const T&);
```

SFINAE should be thought of as a filter on what functions are visible to the compiler when it does a name lookup for a specific function.
The type trait `require_st_var` should be read "require the @ref stan::scalar_type of the argument to be a @ref stan::math::var".
The metaprogram `require_st_var<T>` will evaluate to a valid type if the scalar type of `T` is a @ref stan::math::var.
If the scalar type of `T` is not a @ref stan::math::var, then `require_st_var<T>` does not evaluate to a valid type and the compiler treats it as if the signature does not exist.
This is how SFINAE (substitution failure is not an error) works.
Because the substitution does not work, the signature is ignored.
This is all built based on C++'s [std::enable_if](https://en.cppreference.com/w/cpp/types/enable_if) template metaprogram.

Again, there are many ways to solve a problem in C++.
In particular there are cases where simple overloads or template specializations can achieve the same thing as the SFINAE template metaprograms.
Unless there is a specific reason though, new functions in Math should use the SFINAE metaprograms to handle different implementations.
These are tested to work with the broad set of C++ types that the relatively compact set of Stan types map to.

All of the SFINAE template metaprograms in Math are special to Stan or Math.
Documentation can be found in the @ref require_meta docs.

### Reference types

Any time a function takes a vector or matrix type, it needs to be able to also handle an expression that evaluates to that type.
Some expressions can be expensive to evaluate, so each expression should be evaluated only once.
If the results of an expression is needed in multiple places, they should be saved.
Eigen expressions can result from many places, including arithmetic operations,
a matrix multiply, a sum, or a variety of other things.
Some expressions are cheap to evaluate, any of the Eigen views qualify here (transpose, block access, segment, etc).
In this case, evaluating the expression is not necessary -- it would only lead to another allocation and copy.

The Math function @ref stan::math::to_ref is a solution for this.
For a vector or matrix variable `x`, `to_ref(x)` returns an Eigen reference to the input variable `x`.
If `x` was an expensive expression, it will be evaluated.
If it was a cheap expression, the reference type won't evaluate.
If it isn't an Eigen type, the `to_ref` just forwards.

```cpp
template <typename T,
          require_eigen_t<T>* = nullptr>
auto myfunc(const T& x) {
  const auto& x_ref = to_ref(x);
}
```

### Avoiding hanging pointers in return expressions with `make_holder`

If a function returns an Eigen expression, it may be necessary to use the @ref stan::math::make_holder utility.

Returning expressions is tricky because the expression may have a longer lifetime than the objects it operates on.
For instance, returning an expression on a local variable leads to this situation:

```cpp
template <typename T>
auto add_zero(const T& a) {
  Eigen::VectorXd b = Eigen::VectorXd(a.size());
  return a + b;
}
```


Because the return type of this function is `auto`, it is will return an Eigen expression.
The following code will segfault on the construction of `c` because by this time the temporary `b` will have been destructed.

```
Eigen::VectorXd a(5);
Eigen::VectorXd c = add_zero(a);
std::cout << c(0);
```

[godbolt example](https://godbolt.org/z/6Te856Mv7)

When we return back an **Eigen expression** containing objects created in the local scope we need a way to store those locally created objects till the expression is evaluated.
Stan's @ref stan::math::make_holder can be used to extend the lifetime of objects created in a local scope so that it matches the expression it is used in.

```cpp
template <typename T>
auto add_zero(T&& a) {
  Eigen::VectorXd b = Eigen::VectorXd(a.size());
  return make_holder([](auto&& l, auto&& r) { return l + r; }, a, std::move(b));
}
```

Returning expressions is an advanced feature, and it is easy to make mistakes.
In this regard, it is simplest to start development not returning expressions (in this case holders are unnecessary) and only add expression return types later.

It is always possible to return a non-expression type by evaluating the Eigen expression.
For convenience, there is an `eval` function in Math that will evaluate Eigen expressions and forward anything else.
This is convenient because it works on non-Eigen types as well (as compared to the built-in `.eval()` member function on Eigen types).

The implementation of @ref stan::math::make_holder is [here](https://github.com/stan-dev/math/blob/develop/stan/math/prim/meta/holder.hpp).

### Move Semantics

Move semantics generally work as

```cpp

class my_big_type {
  double* m_data;

  explicit my_big_type(size_t N) m_data(static_cast<double*>(malloc(N))) {}
  // Standard copy constructor
  my_big_type(const my_big_type& other) : m_data(static_cast<double*>(malloc(N))) {
    std::copy(other.m_data, this->m_data);
  }
  // Standard move constructor
  my_big_type(my_big_type&& other) : m_data(other.m_data) {
    other.m_data = nullptr;
  }
};
```

We can see in the above that the standard style of a move (the constructor taking an rvalue reference) is to copy the pointer and then null out the original pointer.
But in Stan, particularly for reverse mode, we need to keep memory around even if it's only temporary for when we call the gradient calculations in the reverse pass.
And since memory for reverse mode is stored in our arena allocator no copying happens in the first place.

Functions for Stan Math's reverse mode autodiff should use [_perfect forwarding_](https://drewcampbell92.medium.com/understanding-move-semantics-and-perfect-forwarding-part-3-65575d523ff8) arguments. Perfect forwarding arguments use a template parameter wit no attributes such as `const` and `volatile` and have a double ampersand `&&` next to them.

```c++
  template <typename T>
  auto my_function(T&& x) {
    return my_other_function(std::forward<T>(x));
  }
```

The `std::forward<T>` in the in the code above tells the compiler that if `T` is deduced to be an rvalue reference (such as `Eigen::MatrixXd&&`), then it should be moved to `my_other_function`, where there it can possibly use another objects move constructor to reuse memory.
A perfect forwarding argument of a function accepts any reference type as its input argument. 
The above signature is equivalent to writing out several functions with different reference types

```c++
  // Accepts a plain lvalue reference
  auto my_function(Eigen::MatrixXd& x) {
    return my_other_function(x);
  }
  // Accepts a const lvalue reference
  auto my_function(const Eigen::MatrixXd& x) {
    return my_other_function(x);
  }
  // Accepts an rvalue reference
  auto my_function(Eigen::MatrixXd&& x) {
    return my_other_function(std::move(x));
  }
  // Accepts a const rvalue reference
  auto my_function(const Eigen::MatrixXd&& x) {
    return my_other_function(std::move(x));
  }
```

In Stan, perfect forwarding is used in reverse mode functions which can accept an Eigen matrix type.

```c++
template <typename T, require_eigen_vt<is_var, T>* = nullptr>
inline auto sin(T&& x) {
  // Store `x` on the arena
  arena_t<T> x_arena(std::forward<T>(x));
  arena_t<T> ret(x_arena.val().array().sin().matrix());
  reverse_pass_callback([x_arena, ret] mutable {
    x_arena.adj() += ret.adj().cwiseProduct(x_arena.val().array().cos().matrix());
  });
  return ret;
}
```

Let's go through the above line by line.

```c++
template <typename T, require_eigen_vt<is_var, T>* = nullptr>
inline auto sin(T&& x) {
```

The signature for this function has a template `T` that is required to be an Eigen type with a `value_type` that is a `var` type. 
The template parameter `T` is then used in the signature as an perfect forwarding argument. 

```c++
  // Store `x` on the arena
  arena_t<T> x_arena(std::forward<T>(x));
```

The input is stored in the arena, which is where the perfect forwarding magic actually occurs. 
If `T` is an lvalue type such as `Eigen::MatrixXd&` then `arena_matrix` will use it's copy constructor, creating new memory in Stan's arena allocator and then copying the values of `x` into that memory. 
But if `T` was a temporary rvalue type such as `Eigen::MatrixXd&&`, then the `arena_matrix` class will use it's move constructor to place the temporary matrix in Stan's `var_alloc_stack_`. 
The `var_alloc_stick_` is used to hold objects that were created outside of the arena allocator but need to be deleted when the arena allocator is cleared. 
This allows the `arena_matrix` to reuse the memory from the temporary matrix. Then the matrix will be deleted once arena allocator requests memory to be cleared.

```c++
  arena_t<T> ret(x_arena.val().array().sin().matrix());
```

This construction of an `arena_matrix` will *not* use the move constructor for `arena_matrix`. 
Here, `x_arena` is an `arena_matrix<T>`, which is then wrapped in an expression to compute the elementwise `sin`. 
That expression will be evaluated into new memory allocated in the arena allocator and then a pointer to it will be stored in the `arena_matrix.`

```c++
  reverse_pass_callback([x_arena, ret] mutable {
    x_arena.adj() += ret.adj().cwiseProduct(x_arena.val().array().cos().matrix());
  });
  return ret;
```

The rest of this code follows the standard format for the rest of Stan Math's reverse mode that accepts Eigen types as input. 
The `reverse_pass_callback` function accepts a lambda as input and places the lambda in Stan's callback stack to be called later when `grad()` is called by the user.
Since `arena_matrix` types only store a pointer to memory allocated elsewhere they are copied into the lambda. 
The body of the lambda holds the gradient calculation needed for the reverse mode pass.

Then finally `ret`, the `arena_matrix` type is returned by the function.

When working with arithmetic types, keep in mind that moving Scalars is often less optimal than simply taking their copy.
For instance, Stan's `var` type uses the pointer to implementation (PIMPL) pattern, so it simply holds a pointer of size 8 bytes.
A `double` is also 8 bytes which just so happens to fit exactly in a [word](https://en.wikipedia.org/wiki/Word_(computer_architecture)) of most modern CPUs with at least 64-byte cache lines.
While a reference to a double is also 8 bytes, unless the function is inlined by the compiler, the computer will have to place the reference into a cache, then go fetch the value that is being referenced which now takes up two words instead of one!

The general rules to follow for passing values to a function are:

1. If a type is smaller than two words (16 bytes), pass it by value
2. If you are writing a function for reverse mode, pass values by `const&`
3. In prim, if you are confident and working with larger types, use perfect forwarding to pass values that can be moved from. Otherwise simply pass values by `const&`.

### Using auto is Dangerous With Eigen Matrix Functions in Reverse Mode

The use of auto with the Stan Math library should be used with care, like in [Eigen](https://eigen.tuxfamily.org/dox/TopicPitfalls.html). 
Along with the cautions mentioned in the Eigen docs, there are also memory considerations when using reverse mode automatic differentiation. 
When returning from a function in the Stan Math library with an Eigen matrix output with a scalar `var` type, the actual returned type will often be an `arena_matrix<Eigen::Matrix<...>>`. 
The `arena_matrix` class is an Eigen matrix where the underlying array of memory is located in Stan's memory arena. 
The `arena_matrix` that is returned by Math functions is normally the same one resting in the callback used to calculate gradients in the reverse pass.
Directly changing the elements of this matrix would also change the memory the reverse pass callback sees which would result in incorrect calculations.

The simple solution to this is that when you use a math library function that returns a matrix and then want to assign to any of the individual elements of the matrix, assign to an actual Eigen matrix type instead of using auto. 
In the below example, we see the first case which uses auto and will change the memory of the `arena_matrix` returned in the callback for multiply's reverse mode.
Directly below it is the safe version, which just directly assigns to an Eigen matrix type and is safe to do element insertion into.

```c++
Eigen::Matrix<var, -1, 1> y;
Eigen::Matrix<var, -1, -1> X;
// Bad!! Will change memory used by reverse pass callback within multiply!
auto mu = multiply(X, y);
mu(4) = 1.0;
// Good! Will not change memory used by reverse pass callback within multiply
Eigen::Matrix<var, -1, 1> mu_good = multiply(X, y);
mu_good(4) = 1.0;
```

The reason we do this is for cases where function returns are passed to other functions. 
An `arena_matrix` will always make a shallow copy when being constructed from another `arena_matrix`, which lets the functions avoid unnecessary copies.

```c++
Eigen::Matrix<var, -1, 1> y1;
Eigen::Matrix<var, -1, -1> X1;
Eigen::Matrix<var, -1, 1> y2;
Eigen::Matrix<var, -1, -1> X2;
auto mu1 = multiply(X1, y1);
auto mu2 = multiply(X2, y2);
// Inputs not copied in this case!
auto z = add(mu1, mu2);
```


### Passing variables that need destructors called after the reverse pass (`make_chainable_ptr`)

When possible, non-arena variables should be copied to the arena to be used in the reverse pass.
The two tools for that are @ref stan::arena_t and @ref stan::math::to_arena .

When these tools do not work, there is @ref stan::math::make_chainable_ptr .
@ref stan::math::make_chainable_ptr constructs a copy of its argument onto stan's arena allocator and returns
a pointer to that copy.
The copy of the argument will only be destructed when @ref stan::math::recover_memory is called.

The pointer is cheap to copy around and is safe to copy into lambdas for @ref stan::math::reverse_pass_callback and @ref stan::math::make_callback_var .

As an example, see the implementation of @ref stan::math::mdivide_left [here](https://github.com/stan-dev/math/blob/develop/stan/math/rev/fun/mdivide_left.hpp) where @ref stan::math::make_chainable_ptr is used to save the result of an Eigen Householder QR decomposition for use in the reverse pass.

The implementation is in [stan/math/rev/core/chainable_object.hpp](https://github.com/stan-dev/math/blob/develop/stan/math/rev/core/chainable_object.hpp)

### Returning arena types

Suppose we have a function such as

```cpp
/**
 * Returns the dot product of a vector of var with itself.
 *
 * @tparam EigVec An Eigen type with compile time rows or columns equal to 1 and a scalar var type.
 * @param[in] v Vector.
 * @return Dot product of the vector with itself.
 */
template <typename EigVec, require_eigen_vt<is_var, EigVec>* = nullptr>
inline var cool_fun(const EigVec& v) {
  arena_t<EigVec> arena_v(v);
  arena_t<EigVec> res = arena_v.val().array() * arena_v.val().array();
  reverse_pass_callback([res, arena_v]() mutable {
    arena_v.adj().array() += (2.0 * res.adj().array()) * arena_v.val().array();
  });
  return res;
}
```

There's a deceptive problem in this return!
We are returning back a @ref stan::math::arena_matrix, which is an Eigen matrix whose dynamic memory sits in our arena allocator.
The problem here is that we've also passed `res` to our callback, and now suppose a user does something like the following.

```cpp
Eigen::Matrix<var, -1, 1> x = Eigen::Matrix<double, -1, 1>::Random(5);
auto innocent_return = cool_fun(x);
innocent_return.coeffRef(3) = var(3.0);
auto other_return = cool_fun2(innocent_return);
grad();
```

Now `res` is `innocent_return` and we've changed one of the elements of `innocent_return`, but that is also going to change the element of `res` which is being used in our reverse pass callback!

Care must be taken by end users of Stan Math by using `auto` with caution. 
When a user wishes to manipulate the coefficients of a matrix that is a return from a function in Stan Math, they should assign the matrix to a plain Eigen type.

```c++
Eigen::Matrix<var, -1, 1> x = Eigen::Matrix<double, -1, 1>::Random(5);
Eigen::MatrixXd actually_innocent_return = cool_fun(x);
actually_innocent_return.coeffRef(3) = var(3.0);
auto still_unsafe_return = cool_fun2(actually_innocent_return);
grad();
```

### Const correctness, reverse mode autodiff, and arena types

In general, it's good to have arithmetic and integral types as `const`, for instance pulling out the size of a vector to reuse later such as `const size_t x_size = x.size();`.
However, vars and anything that can contain a var should not be `const`.
This is because in the reverse mode we want to update the value of the `adj_` inside of the var, which requires that it is non-const.

## Handy tricks

### @ref stan::math::forward_as

In functions such as [Stan's distributions](https://github.com/stan-dev/math/blob/1bf96579de5ca3d06eafbc2eccffb228565b4607/stan/math/prim/prob/exponential_cdf.hpp#L64) you will see code which uses a little function called @ref stan::math::forward_as inside of if statements whose values are known at compile time.
In the following code, `one_m_exp` can be either an Eigen vector type or a scalar.

```cpp
T_partials_return cdf(1.0);
if (is_vector<T_y>::value || is_vector<T_inv_scale>::value) {
  cdf = forward_as<T_partials_array>(one_m_exp).prod();
} else {
  cdf = forward_as<T_partials_return>(one_m_exp);
}
```

Since the if statements values are known at compile time, the compiler will always remove the unused side of the `if` during the dead code elimination pass.
But the dead code elimination pass does not happen until all the code is instantiated and verified as compilable.
So @ref stan::math::forward_as exists to trick the compiler into believing both sides of the `if` will compile.
If we used C++17, the above would become

```cpp
T_partials_return cdf(1.0);
if constexpr (is_vector<T_y>::value || is_vector<T_inv_scale>::value) {
  cdf = one_m_exp.prod();
} else {
  cdf = one_m_exp;
}
```


Where [`if constexpr`](https://en.cppreference.com/w/cpp/language/if) is run before any tests are done to verify the code can be compiled.

Using `forward_as<TheTypeIWant>(the_obj)` will, when `the_obj` matches the type the user passes, simply pass back a reference to `the_obj`.
But when `TheTypeIWant` and `the_obj` have different types it will throw a runtime error.
This function should only be used inside of `if` statements like the above where the conditionals of the `if` are known at compile time.
