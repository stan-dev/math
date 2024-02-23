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

In the Stan Math library, all the primary function templates are
declared in `stan/math/prim/`. The role of the primary function
template is to define a generic implementation for all valid argument
types. The definition of the function often restricts the valid
argument types through [substitution failure is not an error
(SFINAE)](https://en.cppreference.com/w/cpp/language/sfinae).

For many functions, we write specializations of the function
templates, usually for readability of code and computationaly
efficiency. To take advantage of generic template types while
limiting the implementation to a subset of valid types, we use
[partial template
specialization](https://en.cppreference.com/w/cpp/language/partial_specialization)
and specifically [SFINAE in partial
specialization](https://en.cppreference.com/w/cpp/language/sfinae#SFINAE_in_partial_specializations)
to enable this behavior. The source code for function specializations are
found in either `stan/math/rev/` for reverse-mode implementations,
`stan/math/fwd/` for forward-mode implementations, or
`stan/math/mix/` for implementations that are for nested
autodiff.

In the partial template specializations, we include a `requires`
template parameter as a pointer [non-type template
parameter](https://en.cppreference.com/w/cpp/language/template_parameters#Non-type_template_parameter)
with a default value of `nullptr`; any `void*` non-template parameter
with a default of `nullptr` is ignored by the compiler. When we pass a
type that satisfies the `requires`, the `requires` trait evaluates to
`void`, the template parameter evaluates to `void*=nullptr`, and the
template argument is safely ignored by the compiler. When a type that
does not satisfy the `requires` template parameter is passed, the
definition is removed from the set of possible functions via
SFINAE. With this scheme we end up having a very nice pattern for
writing generic templates for functions while also being able to
restrict the set of types that a function can be used for.

For convenience in writing partial template specializations, we have
implemented a set of [`requires` type traits](@ref require_meta).

The rest of this page describes what the requires type traits are, how
to use them, and how to add new ones if necessary.

### Requires<> type traits

The Stan Math library defines boolean type traits--template
metaprograms that operate on types at compile time--in the
`stan/math{prim, rev, fwd}/meta` folders. Each of these type traits
are named `is_{condition}` and the struct contains a `value` that is
`true` or `false` at compile time. For example, `is_var<T>::value` is
`true` if and only if the type `T` is `stan::math::var_value`.

For ease of use in partial template specialization, we provide
`requires<>` type traits based on the boolean `is_{condition}` type
traits. When types satisfy the condition, the `requires<>::value` will
evaluate to `void`. When the types do not satisfy the condition,
`requires<>::value` is an invalid subsitution and is not used. (See
@ref requires_impl for more details.)

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
type trait. This means `require_var_t<stan::math::var>::value` is
`void`, but `require_var_t<double>::value` is an invalid subsitution.

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
which requires type traits are currently in use (as of 2024). If you
need a requires type trait and it is not currently available, please
feel free to implement the one you need and add a pull request.

We expect this to happen very infrequently. If it is happening often,
we can go back to having all possible combinations of requires traits
available.

### Adding a new boolean type trait

If you are adding a new boolean type trait, please add the primary
template function to `stan/math/prim/meta/`, then add any autodiff
specialization to the appropriate `stan/math/{rev, fwd, mix}/meta/`
folder.

### Adding a new requires

If you need to add a new `requires`, please submit a pull request! 

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

    // Ex: Used to define a signature for `std::vectors` with a value type that is autodiffable
    template <typename StdVec, require_std_vector_vt<is_var, StdVec>* = nullptr>
    auto my_func(StdVec&& vec);
    ```

* `_st` uses `Check` to test the type `T` passed in and uses
  `InnerCheck` to test the @ref scalar_type of `T`

    ```
    template <template <class...> class TypeCheck, class... Check>   
    require_std_vector_st = require_st<is_std_vector, TypeCheck, std::decay_t<Check>...>;

    // Ex: Used to define a signature for `std::vectors` with a scalar type that is autodiffable
    template <typename StdVec, require_std_vector_st<is_var, StdVec>* = nullptr>
    auto my_func(StdVec&& vec);
    ```

There variant of the requires that places the `vt` or `st` before the
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



