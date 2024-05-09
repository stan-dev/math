## Using requires<> template metaprograms for partial template specialization {#require_meta_doc}

Most of the functions in the Stan Math library are implemented as
[templated
functions](https://en.cppreference.com/w/cpp/language/function_template)
which allow for arguments that can be C++ primitive types
(e.g. `double`, `int`), Stan's reverse-mode or forward-mode automatic
differentiation (autodiff) types, or containers and expressions of
either primitive or autodiff types. We use templated functions rather
than overloaded functions for a number of reasons including the sheer
number of implementations we would need to write to handle the
combinations of arguments that are allowed in the Stan language.

Many functions in the Stan Math library have a single function
template defined in the `stan/math/prim` folder that is flexible
enough to accept both primitive and autodiff types. Some of these
functions also have template specializations (usually [partial
template
specializations](https://en.cppreference.com/w/cpp/language/partial_specialization))
that define an implementation where at least one template parameters
is restricted to a specific type. The source code for the function
template specializations are found in either `stan/math/rev/` for
reverse-mode implementations, `stan/math/fwd/` for forward-mode
implementations, or `stan/math/mix/` for implementations that are for
nested autodiff. This pattern of a primary template function
definition with specialization is commonly used in templated C++.

In the Stan Math library, we have also adopted another technique which
allows multiple template function definitions, while restricting the
definition to apply only to when certain criteria are met. Since this
technique is used repeatedly through the Math library and this is not a
common use of template metaprogramming, we'll describe it in the next
subsection.

### Using SFINAE to allow multiple template function definitions

In the Stan Math library, each function exposed through to the Stan
language must have definitions that allow for both primitive and
autodiff types. As the language has grown, we have also added
broadcasting, which allows users to mix scalars arguments with vector
or array arguments. 

The typical C++ method of defining a primary function template and a
set of function template specialization is untenable for many
functions in the Math library. For example, a single argument of the
function may take 7 distinct C++ types: `double`, `stan::math::var`,
`std::vector<double>`, `std::vector<stan::math::var>`,
`Eigen::Matrix<double, -1, 1>`, `Eigen::Matrix<stan::math::var, -1,
1>`, or `stan::math::var_value<Eigen::Matrix<double, -1, -1>>`.
For a 3-argument function, we would need to define 343 (7^3) different
function template specializations to handle all the autodiff types.


In the Math library, we use a technique similar to C++20's `require`
keyword that allows the definition of multiple template functions
where each handles a subset of allowable types.

When the compiler attempts to resolve which function should be called
from a set of templated function signatures there must be only one
possibly valid function signature available. This is called the [One
Definition
Rule](https://en.cppreference.com/w/cpp/language/definition). For
example, the following code would fail to compile because the compiler
is unable to differentiate between the two function signatures.

```c++
template <typename T>
T foo(T x) {
  return x;
}

template <typename K>
K foo(K x) {
  return x;
}
```

The compiler needs a way to differentiate between the two signatures
to select one and satisfy the One Definition Rule. One trick to have a
single valid definition is to utilize [Substitution Failure Is Not An
Error (SFNIAE)](https://en.cppreference.com/w/cpp/language/sfinae) to
purposefully create conditions where only one signature is valid
because all of the other conditions fail to compile. The simplest way
to do this is to start with a type trait like the below
`enable_if`. The `enable_if` is only defined for the case where `B` is
`true` and so if `B` is ever false the compiler would throw an error
saying that `enable_if` is not well defined.

```
// Define enable_if. Note: `type` is not a member typedef.
template<bool B, class T = void>
struct enable_if {};

// Only define member typedef `type` when B is true
template<class T>
struct enable_if<true, T> { typedef T type; };

template <bool B, typename T>
using enable_if_t = typename enable_if<B, T>::type;
```

Attempting to construct this `enable_if` with `B` being `false`
anywhere else in the program would cause the compiler to crash.
Using it in the template of a function signature allows SFINAE to
deduce which signature we would like to use.

```c++
// foo only works with floating point types 
template <typename T,  enable_if_t<std::is_floating_point<T>::value>>* = nullptr>
T foo(T x) {
  return x;
}

// foo only works with integer types
template <typename K,  enable_if_t<std::is_intergral<K>::value>>* = nullptr>
K foo(K x) {
  return x;
}

// Calls the first signature
double x_dbl = 1.0;
double y_dbl = foo(x_dbl); 

// Calls the second signature
int x = 1;
int y = foo(x);

```

The second template argument is referred to as a [non-type template
parameter](https://en.cppreference.com/w/cpp/language
template_parameters#Non-type_template_parameter) and has a default
value of `void`.  When the templated signature has the correct type
the `enable_if_t ` produces a `void` type which is then made into a
pointer and assigned a default value of `nullptr`.  When the templated
signature does not have the correct type, the compiler utilizes
[Substitution Failure Is Not An Error
(SFNIAE)](https://en.cppreference.com/w/cpp/language/sfinae), to
remove the offending signature from the list of possible matches while
continuing to search for the correct signature.


For convenience in using this technique, the Math library has
implemented a set of [`requires` type traits](@ref require_meta).
When we pass a type that satisfies the `requires` type trait, the
trait evaluates to `void`. When a type that does not satisfy the
`requires` template parameter is passed, there is a substitution
failure. These traits are used in the template functions by adding a
call to a `requires` type trait to the parameter list.

Below is an example to illustrate how this technique is used. After
the example, the rest of this page describes what the requires type
traits are, how to use them, and how to add new ones if necessary.

#### Example

Here's a function that would have two different template functions for
`stan::math::var` and `double`:

```c++
template <typename T, requires_not_var_t<T>* = nullptr>
T foo(const T& t) {
  // handles primitives
  return t;
}

template <typename T, requires_var_t<T>* = nullptr>
T foo(const T& t) {
  // handles stan::math::var
  return t;
}
```

When `foo()` is called with a `stan::math::var`, the first template
function matches but not the second. This works because
`requires_var_t<stan::math::var>` evaluates to `void` whereas
`requires_not_var_t<stan::math::var>` is a subsitution error causing
the compiler to omit the second definition for `stan::math::var`.

When `foo()` is called with `double` or `int`, the second template
function matches, but not the first.


### Requires<> type traits

The Stan Math library defines boolean type traits--template
metaprograms that operate on types at compile time--in the
`stan/math{prim, rev, fwd}/meta` folders. Each of these type traits
are named `is_{condition}` and the struct contains a `value` that is
`true` or `false` at compile time. For example, `is_var<T>::value` is
`true` if and only if the type `T` is `stan::math::var_value`.

We provide [`requires<>` type traits](@ref require_meta) based on the boolean
`is_{condition}` type traits. When types satisfy the condition, the
`requires<>` will evaluate to `void`. When the types do not
satisfy the condition, `requires<>` is an invalid subsitution
and is not used. (See @ref requires_impl for more details.)

Note: every possible requires<> type trait is not implemented in the
Stan Math library. If one of the missing requires<> type trait is
missing, we can implement it and include it. Please see @ref
requires_dev_guide for more information.


#### All requires<> type traits

For any boolean type trait, below is the list of possible requires<>
type traits.  Any `*` should be thought of as a wildcard where a type
traits name is put in its place. For example, for `is_var`, we can
substitute `var` for `*`.

1. `require_*_t`: A template parameter `T` must satisfy the `is_*`
type trait. This means `require_var_t<stan::math::var>` is
`void`, but `require_var_t<double>` is an invalid subsitution.

2. `require_not_*_t`: A template parameter `T` must not satisfy the
`is_*` type trait.

    *NOTE:* The `not` version of the `requires` template parameters
should be used sparingly.  Often a `requires` template parameter is
used to specify what types a function should accept.  Defining a
function by the types it cannot accept can make understanding what
goes into a function more difficult and error prone.

3. `require_all_*_t`: All types in the parameter pack of types must
satisfy the `is_*` type trait.

4. `require_any_*_t`: Any type in the parameter pack of types must
satisfy the `is_*` type trait.

5. `require_any_not_*_t`: Any type in the parameter pack must not
satisfy the `is_*` type trait.

6. `require_all_not_*_t`: All types in the parameter pack must not
satisfy the `is_*` type trait.

    `std::vector` and `Eigen` types have additional `requires`
    template parameters to detect if the @ref stan::value_type (the
    type of the elements of either `std::vector` or the `Eigen` type)
    or the @ref stan::scalar_type (the underlying scalar type after
    recursively walking through the container types) satisfy a
    condition to enable a class or function.

    The container `requires` template parameters have an ending at
    their signature of `_vt` and `_st` to symbolize whether you want
    to inspect the @ref stan::value_type or @ref stan::scalar_type.

    In the next requires traits, `is_type` is used to represent any
    boolean type trait.
    
7. `require_*_vt<is_type, T>`: A template parameter `T` must satisfy
the `is_*` type trait and `is_type<value_type<T>>::value` must
evaluate to true.

8. `require_not_*_vt<is_type, T>`: A template parameter `T` must not
satisfy the `is_*` type trait or `is_type<value_type<T>>::value` must
not evaluate to true.

9. `require_all_*_vt<is_type, T>`: All types in the parameter pack of
types must satisfy the `is_*` type trait and all
`is_type<value_type<T>>::value` must evaluate to true.

10. `require_any_*_vt<is_type, T>`: Any type in the parameter pack of
types must satisfy the `is_*` type trait and any
`is_type<value_type<T>>::value` must evaluate to true.

11. `require_any_not_*_vt<is_type, T>`: At least one type in the
parameter pack must not satisfy the `is_*` type trait and one of
`is_type<value_type<T>>::value` must evaluate to false.

12. `require_all_not_*_vt<is_type, T>`: None of the types in the
parameter pack must satisfy the `is_*` type trait and none of
`is_type<value_type<T>>::value` must evaluate to true.
    
13. `require_*_st<is_type, T>`: A template parameter `T` must satisfy
the `is_*` type trait and `is_type<scalar_type<T>>::value` must
evaluate to true.

14. `require_not_*_st<is_type, T>`: A template parameter `T` must not
satisfy the `is_*` type trait or `is_type<scalar_type<T>>::value` must
not evaluate to true.

15. `require_all_*_st<is_type, T>`: All types in the parameter pack of
types must satisfy the `is_*` type trait and all
`is_type<scalar_type<T>>::value` must evaluate to true.

16. `require_any_*_st<is_type, T>`: Any type in the parameter pack of
types must satisfy the `is_*` type trait and any
`is_type<scalar_type<T>>::value` must evaluate to true.

17. `require_any_not_*_st<is_type, T>`: At least one type in the
parameter pack must not satisfy the `is_*` type trait and one of
`is_type<scalar_type<T>>::value` must evaluate to false.

18. `require_all_not_*_st<is_type, T>`: None of the types in the
parameter pack must satisfy the `is_*` type trait and none of
`is_type<scalar_type<T>>::value` must evaluate to true.



### Implementation details of requires<> type traits {#requires_impl}

The [`requires` template parameters](@ref require_meta) type traits
are aliases for `std::enable_if_t` that have premade conditions for
turning on and off function definitions during compilation. These are
useful for having generalized templates while still overloading a
function or class.  You can think of these as "legacy concepts."
These are used in a very similar fashion to C++20's `requires`
keyword.

`requires` template parameters are
[`std::enable_if_t`](https://en.cppreference.com/w/cpp/types/enable_if)
aliases such as the following example definition of @ref
stan::require_t.

```cpp
template <typename T>
using require_t = std::enable_if_t<T::value>;
```

This differs from `std::enable_if_t` in that `std::enable_if_t`'s
argument must be boolean, but the alias @ref stan::require_t 's
template type `T` must have a valid boolean member named `value`.
This allows us to directly call @ref stan::require_t with type traits
instead of having to do the extra step of accessing the type traits
boolean member struct value explicity with calls such as
`a_type_trait::value`.

The most common use case for a `requires` template parameters is to
overload a function or declare specializations of a class.  For
example, the function below will only work on types derived from
[`Eigen::DenseBase`](https://eigen.tuxfamily.org/dox/classEigen_1_1DenseBase.html)
with only 1 row or column at compile time such as
`Eigen::Matrix<double, -1, 1>` or `Eigen::Matrix<double, 1, -1>`.

```cpp
template <typename EigVec, require_eigen_vector_t<EigVec>* = nullptr>
auto my_func(const EigVec& x) {
  // ...
}
```

For overloading classes and structs with this scheme we create an
initial forward definition with a `void` non-type template parameter.
Then the class overloads use the `requires` template parameter in
place of the non-type template parameter.

```cpp
template <typename T, typename = void>
class a_class;

/**
 * Overload for standard vectors with \ref stan::scalar_type
 *  of \ref stan::math::var
 */
template <typename T>
class a_class<T, require_std_vector_st<is_var, T>> {
  // ...
};
```

In the above example, `a_class` has an overload specifically for
standard vectors with a @ref stan::scalar_type of @ref stan::math::var.

There are also `requires` template parameters for generically checking
if a type's @ref stan::value_type or @ref stan::scalar_type is
correct.  To differentiate them from the Eigen and standard library
vector checks the `vt` and `st` come *before* the type such as
`require_vt_var<T>` which checks if a type `T`'s @ref stan::value_type
satisfies @ref stan::is_var.

The `requires` template parameters type traits allow Stan to have more
generic types so that the library can forward Eigen expression and
have better move semantics.  For instance, the code below will accept
any arbitrary Eigen expression that, if it's an rvalue, can be
forwarded to another function.

```cpp
 template <typename Mat1, typename Mat2,
  require_all_eigen_vt<is_arithmetic, Mat1, Mat2>* = nullptr>
 inline auto a_func(Mat1&& m1, Mat2&& m2) {
   check_not_nan(m1);
   check_not_nan(m2);
   return another_func(std::forward<Mat1>(m1), std::forward<Mat2>(m2));
```

## Developing new requires type traits {#requires_dev_guide}

Every requires type trait is not implemented for every boolean type
trait available. This was done intentionally to allow us to identify
which requires type traits are currently in use. If you
need a requires type trait and it is not currently available, please
feel free to implement the one you need and add a pull request.


### Adding a new boolean type trait

If you are adding a new boolean type trait, please add the primary
function template to `stan/math/prim/meta/`, then add any autodiff
specialization to the appropriate `stan/math/{rev, fwd, mix}/meta/`
folder.

### Adding a new requires


The Stan Math library requires a strict API to ensure consistency for
the `requires`. The below go over all of the possible API
configurations a developer should use when writing a new `requires`.

For the API docs below, let `T` represent the type parameter we want
to check, `is_type` is a generic type trait which will be replaced by
the developer, and `InnerCheck` is a type trait used to check either
the @ref value_type or @ref scalar_type of `T`.

Each requires ends in `_t`, `_vt`, or `_st`. They differ in the
following ways

* `_t` uses `Check` to test the type `T` passed in

    Ex:

    ```
    // Always decay types coming into the requires
    template <typename T>
    require_autodiff_t = require_t<is_autodiff<std::decay_t<T>>>;
    ```

* `_vt` uses `Check` to test the type `T` passed in and uses
  `InnerCheck` to test the @ref value_type of `T`

    ```
    template <template <class...> class TypeCheck, class... Check>   
    require_std_vector_vt = require_vt<is_std_vector, TypeCheck, std::decay_t<Check>...>;

    // Ex: Used to define a signature for `std::vectors` with a value type that is @ref stan::math::var
    template <typename StdVec, require_std_vector_vt<is_var, StdVec>* = nullptr>
    auto my_func(StdVec&& vec);
    ```

* `_st` uses `Check` to test the type `T` passed in and uses
  `InnerCheck` to test the @ref scalar_type of `T`

```cpp
template <template <class...> class TypeCheck, class... Check>   
require_std_vector_st = require_st<is_std_vector, TypeCheck, std::decay_t<Check>...>;

// Ex: Used to define a signature for `std::vectors` with a scalar type that is autodiffable
template <typename StdVec, require_std_vector_st<is_var, StdVec>* = nullptr>
auto my_func(StdVec&& vec);

The variant of the requires that places the `vt` or `st` before the
type trait name only checks the @ref value_type or @ref `scalar_type`
of `T` without testing `T`.

```
// Require the scalar type is an std::vector
template <typename T>   
require_st_std_vector = require_t<is_std_vector<scalar_type_t<std::decay_t<T>>>>;

// Ex: Used to define a signature for `std::vectors` with a scalar type that is autodiffable
template <typename StdVec, require_std_vector_st<is_var, StdVec>* = nullptr>
auto my_func(StdVec&& vec);
```

In the below, `{TYPE_TRAIT}` represents the name of the trait the
requires checks. Each new require _must_ follow this standard API.

1. `require_{TYPE_TRAIT}_t`: The template parameter must return `true`
when passed to the type trait

2. `require_not_{TYPE_TRAIT}_t`: The template parameter must return
`false` when passed to the type trait

3. `require_all_{TYPE_TRAIT}_t`: The template parameters must all return
`true` when passed to the type trait

4. `require_all_not_{TYPE_TRAIT}_t`: The template parameters must all
return `false` when passed to the type trait

5. `require_any_{TYPE_TRAIT}_t`: At least one of the template parameters
must return `true` when passed to the type trait

6. `require_any_not_{TYPE_TRAIT}_t`: At least one of the template
parameters must return `false` when passed to the type trait

In addition to all the requires with an `_t` at the end, the requires
also have `_st`, `_vt` variants where in addition to the logic above,
the @ref value_type or @ref scalar_type must follow the same logic as
the type for `T`. The `_st_`, and `_vt_` variants must also follow the
same logic but for checking only the inner @ref value_type or
@scalar_type.



