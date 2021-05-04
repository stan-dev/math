## Common Pitfalls {#common_pitfalls}



### Eigen Expressions (and handling argument types with SFINAE)

An additional complexity of the Math library is that `vector`, `row_vector`,
and `matrix` do not map to just templated `Eigen::Matrix` types, they can
also map to a variety of Eigen expressions.

For instance, in the following code, the result `c` is not a vector, but
a vector-like type representing the not-yet-calculated sum of `a` and `b`:

```cpp
Eigen::VectorXd a(5);
Eigen::VectorXd b(5);
// auto is actually a
// Eigen::CwiseBinaryOp<Eigen::internal::scalar_sum_op<double, double>, const Eigen::Matrix<double, -1, 1>, const Eigen::Matrix<double, -1, 1>>
auto c = a + b;
```

Similarly, there are other representations of vectors and matrices
yet discussed, in particular matrices stored in Math's special arena
memory or matrices stored on a GPU. In either case, it is expedient
to not write explicit overloads for all the different types that a
function might accept, but limit them with template meta-programs.

For instance, if only base Eigen types are allowed, then a function
that takes vector types could be defined as follows:

```cpp
template <typename T>
T norm(const Eigen::Matrix<T, Eigen::Dynamic, 1>&);
```

A typical problem with a function like this is that `norm` can
similarly be defined for a row vector (and the implementation is
the same). In this case there are two signatures:

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

The immediate problem here is what if there actually are differences in
the implementations? The most common situation is that encountered when
there is one generic templated implementation that works with any autodiff
type, and then a second faster overload that only works with a specific
autodiff type.

There can easily be ambiguities between the two function signatures. The
previous examples took advantage of simple overloads. The more general solution
are template metaprograms that additionally make use of C++ substitution failure
is not an error (SFINAE) semantics. For instance, for the norm function above,
SFINAE could be used to limit one signature to work with reverse
mode autodiff types and one to work with anything else:

```cpp
template <typename T,
          require_st_var<T>* = nullptr>
return_type_t<T> norm(const T&);
template <typename T,
          require_not_st_var<T>* = nullptr>
return_type_t<T> norm(const T&);
```

SFINAE should be thought of as filters on what functions are visible to the
compiler when it does name lookup for a specific function. `require_st_var`
should be read "require the scalar type of the argument to be a `var`". The
meta-program `require_st_var<T>` will evaluate to a valid type if the scalar
type of `T` is a `var`. If the scalar type of `T` is not a `var`, then
`require_st_var<T>` does not evaluate to a valid type, and it is not possible
to do template substitution on this signature, and the compiler treats it as
if it does not exist. This is how SFINAE (substitution failure is not an error)
works. Because the substitution does not work, the signature is ignored.
This is all built based on C++'s
[std::enable_if](https://en.cppreference.com/w/cpp/types/enable_if)
template metaprogram.

Again, there are many ways to solve a problem in C++. In particular there
are cases where simple overloads or template specializations can achieve the
same thing as the SFINAE template meta-programs. Unless there is a specific
reason though, new functions in Math should use the SFINAE meta-programs to
handle different implementations. These are tested to work with the broad
set of C++ types that the relatively compact set of Stan types map to.

The vast majority of the SFINAE template meta-programs in Math are special
to Stan or Math. Documentation can be found in the
[Doxygen](https://mc-stan.org/math/) docs. If functionality is not available,
it new meta-programs should be added.



### Reference types

Any time a function takes a vector or matrix type, it needs to be able
to also handle an expression that evaluates to that type. Some expressions
can be expensive to evaluate, so each expression should be evaluated only once.
If the results of an expression is needed in multiple places, they should be saved.
Eigen expressions can result from many places, including arithmetic operations,
a matrix multiply, a sum, or a variety of other things. Some expressions are cheap
to evaluate, any of the Eigen views qualify here (transpose, block access, segment, etc).
In this case, evaluating the expression is not necessary -- it would only lead to another
allocate and copy.

The Math function `to_ref` is a solution for this. For a vector
or matrix variable `x`, `to_ref(x)` returns an Eigen reference to the
input variable `x`. If `x` was an expensive expression, it will be evaluated.
If it was a cheap expression, the reference type won't evaluate. If it
isn't an Eigen type, the `to_ref` just forwards.

```cpp
template <typename T,
          require_eigen_t<T>* = nullptr>
auto myfunc(const T& x) {
  const auto& x_ref = to_ref(x);
}
```

### Holders

If a function returns an Eigen expression, it may be necessary to use
the `make_holder` utility.

Returning expressions is tricky because the expression may have a longer
lifetime than the objects it operates on. For instance, returning an
expression on a local variable leads to this situation:

```cpp
template <typename T>
auto add_zero(const T& a) {
  Eigen::VectorXd b = Eigen::VectorXd(a.size());
  return a + b;
}
```

Because the return type of this function is `auto`, it is possible
that it returns an Eigen expression.

The following code will segfault on the construction of `c` because
by this time the temporary `b` will have been destructed.

```
Eigen::VectorXd a(5);
Eigen::VectorXd c = add_zero(a);
```

`make_holder` can be used to extend the lifetime of `b` so that it matches
the expression it is used in.

```cpp
template <typename T>
auto add_zero(const T& a) {
  Eigen::VectorXd b = Eigen::VectorXd(a.size());
  return make_holder([](const auto& l, const auto& r) { return l + r; }, a, std::move(b));
}
```

Returning expressions is an advanced feature and it is easy to make
mistakes. In this regard, it is simplest to start development not
returning expressions (in this case holders are unnecessary), and only
add expression return types later.

It is always possible to return a non-expression type by evaluating
the Eigen expression. For convenience, there is an `eval` function
in Math that will evaluate Eigen expressions and forward anything
else. This is convenient because it works on non-Eigen types as well
(as compared to the built in `.eval()` member function on Eigen
types).

The implementation of `make_holder` is
[here](https://github.com/stan-dev/math/blob/develop/stan/math/prim/meta/holder.hpp).

### Move Semantics

In general, Stan math does not use move semantics very often. This is because of our arena allocator. Move semantics generally work as

```cpp

class my_big_type {
  double* m_data;

  explicit my_big_type(size_t N) m_data(static_cast<double*>(malloc(N))) {}
  // Standard copy constructor
  my_big_type(const my_big_type& other) m_data(static_cast<double*>(malloc(N))) {
    std::copy(other.m_data, this->m_data);
  }
  // Standard move constructor
  my_big_type(my_big_type&& other) m_data(other.m_data) {
    other.m_data = nullptr;
  }
};
```

We can see in the above that the standard style of a move (the constructor taking an rvalue reference) is to copy the pointer and then null out the original pointer. But in Stan, particularly for reverse mode, we need to keep memory around even if  it's only a temporary for when we call the gradient calculations in the reverse pass. And since memory for reverse mode is stored in our arena allocator no copying happens in the first place.

When working with arithmetic types, keep in mind that moving Scalars is often less optimal than simply taking their copy. For instance, Stan's `var` type is a PIMPL implimentation, so it simply holds a pointer of size 8 bytes. A `double` is also 8 bytes which just so happens to fit exactly in a [word](https://en.wikipedia.org/wiki/Word_(computer_architecture)) of most modern CPUs. While a reference to a double is also 8 bytes, unless the function is inlined by the compiler, the computer will have to place the reference into cache, then go fetch the value that is being referenced which now takes up two words instead of one!

The general rules to follow for passing values to a function are:

1. If a type is smaller than two words (16 bytes), pass it by value
2. If you are writing a function for reverse mode, pass values by `const&`
3. In prim, if you are confident and working with larger types, use perfect forwarding to pass values that can be moved from. Otherwise simply pass values by `const&`.

### Copying non-arena variables to lambdas used in the reverse pass (`make_callback_ptr`)

When possible, non-arena variables should be copied to the arena to be
used in the reverse pass. The two tools for that are `arena_t<T>` and
`to_arena`.

When these tools do not work, there is `make_callback_ptr(x)`.
`make_callback_ptr(x)` constructs a copy of the argument `x` and returns
a pointer to that copy. The copy of `x` will only be destructed when the
enclosing `recover_memory` is called.

The pointer is cheap to copy around and is safe to copy into lambdas for
`reverse_pass_callback` and `make_callback_var`.

As an example, see the implementation of `mdivide_left`
[here](https://github.com/stan-dev/math/blob/develop/stan/math/rev/fun/mdivide_left.hpp)
where `make_callback_ptr` is used to save the result of an Eigen
Householder QR decomposition for use in the reverse pass.

The implementation is in
[stan/math/rev/core/chainable_object.hpp](https://github.com/stan-dev/math/blob/develop/stan/math/rev/core/chainable_object.hpp)

### Returning arena types

### `.val()` vs. `.val_op()`

### Returning expressions

### Const correctness, reverse mode autodiff, and arena types

## Handy tricks

### `forward_as`

### Nested Stacks

### Scoped Stacks

As an aside, in order to get the actual vector result, Eigen defines the
`.eval()` operator, so `c.eval()` would return an actual `Eigen::VectorXd`.
Similarly, Math defines the `eval` function which will evaluate an Eigen
expression and forward anything that is not an Eigen expression.

### `value_of_rec`
