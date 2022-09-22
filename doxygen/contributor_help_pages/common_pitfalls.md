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

In general, Stan Math does not use move semantics very often.
This is because of our arena allocator.
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

When working with arithmetic types, keep in mind that moving Scalars is often less optimal than simply taking their copy.
For instance, Stan's `var` type uses the pointer to implementation (PIMPL) pattern, so it simply holds a pointer of size 8 bytes.
A `double` is also 8 bytes which just so happens to fit exactly in a [word](https://en.wikipedia.org/wiki/Word_(computer_architecture)) of most modern CPUs with at least 64-byte cache lines.
While a reference to a double is also 8 bytes, unless the function is inlined by the compiler, the computer will have to place the reference into a cache, then go fetch the value that is being referenced which now takes up two words instead of one!

The general rules to follow for passing values to a function are:

1. If a type is smaller than two words (16 bytes), pass it by value
2. If you are writing a function for reverse mode, pass values by `const&`
3. In prim, if you are confident and working with larger types, use perfect forwarding to pass values that can be moved from. Otherwise simply pass values by `const&`.

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
The answer for this is simple but sadly requires a copy.

```cpp
template <typename EigVec, require_eigen_vt<is_var, EigVec>* = nullptr>
inline var cool_fun(const EigVec& v) {
  arena_t<EigVec> arena_v(v);
  arena_t<EigVec> res = arena_v.val().array() * arena_v.val().array();
  reverse_pass_callback([res, arena_v]() mutable {
    arena_v.adj().array() += (2.0 * res.adj().array()) * arena_v.val().array();
  });
  return plain_type_t<EigVec>(res);
}
```

we make a deep copy of the return whose inner `vari` will not be the same, but the `var` will produce a new copy of the pointer to the `vari`.
Now the user code above will be protected, and it is safe for them to assign to individual elements of the `auto` returned matrix.

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
