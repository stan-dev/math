## Using requires<> for general overloads {#require_meta_doc}

The [`requires` template parameters](@ref require_meta) type traits are aliases for `std::enable_if_t` that have premade conditions for turning on and off function definitions during compilation.
These are useful for having generalized templates while still overloading a function or class.
You can think of these as "legacy concepts".
These are used in a very similar fashion to C++20's `requires` keyword.

`requires` template parameters are [`std::enable_if_t`](https://en.cppreference.com/w/cpp/types/enable_if) aliases such as the following example definition of @ref stan::require_t.

```cpp
template <typename T>
using require_t = std::enable_if_t<T::value>;
```

This differes from `std::enable_if_t` in that `std::enable_if_t`'s argument must be boolean, but the alias @ref stan::require_t 's template type `T` must have a valid boolean member named `value`.
This allows us to directly call @ref stan::require_t with type traits instead of having to do the extra step of accessing the type traits boolean member struct value explicity with calls such as `a_type_trait::value`.

The most common use case for a `requires` template parameters is to overload a function or declare specializations of a class.
For example, the function below will only work on types derived from [`Eigen::DenseBase`](https://eigen.tuxfamily.org/dox/classEigen_1_1DenseBase.html) with only 1 row or column at compile time such as `Eigen::Matrix<double, -1, 1>` or `Eigen::Matrix<double, 1, -1>`.

```cpp
template <typename EigVec, require_eigen_vector_t<EigVec>* = nullptr>
auto my_func(const EigVec& x) {
  // ...
}
```

The `requires` template parameter is included in the template as a pointer [non-type template parameter](https://en.cppreference.com/w/cpp/language/template_parameters#Non-type_template_parameter) with a default value of `nullptr`.
This might look a bit odd, but it uses the fact that any non-type template parameter of type `void*` with a default of `nullptr` will be ignored by the compiler.
Under the hood, if any of the `requires` template parameters are successful they will return back a type `void`.
So when we pass a type that the `requires` template parameter accepts we get back a `void* = nullptr`, which is safely ignored by the compiler.
In the case that the type does not satisfy the `requires` template parameter then the function is removed from the set of possible functions the caller could use via SFINAE.
With this scheme we end up having a very nice pattern for writing generic templates for functions while also being able to restrict the set of types that a function can be used for.

For overloading classes and structs with this scheme we create an initial forward definition with a `void` non-type template parameter.
Then the class overloads use the `requires` template parameter in place of the non-type template parameter.

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

In the above example, `a_class` has an overload specifically for standard vectors with a @ref stan::scalar_type of @ref stan::math::var .

The examples below cover the general themes for all of the [`requires` template parameters](@ref require_meta) found in the Stan math library.
Any `*` should be thought of as a wildcard where a type traits name is put in its place.

- `requires_*_t`: A template parameter `T` must satisfy the `requires` template parameter in order for
the overload to be available.

```cpp
 // Works for stan::math::var types
 template <typename T, require_var_t<T>* = nullptr>
 auto add(T&& x, T&& y) { return x + y; }
```

- `require_not_*_t` : A template parameter `T` must *not* satisfy the `requires` template parameter in order for the overload to be availabe.

*NOTE:* The `not` version of the `requires` template parameters should be used sparingly.
Often a `requires` template parameter is used to specify what types a function should accept.
Defining a function by the types it cannot accept can make understanding what goes into a function more difficult and error prone.

```cpp
 // Works for anything that is not a std::vector
 template <typename T, require_not_std_vector_t<T>* = nullptr>
 auto add(T&& x, T&& y) { return x + y; }
```

- `require_all_*_t` : Takes a parameter pack of types to enable if all types satisfy the check.

```cpp
 // Works if T1 and T2 are complex
 template <typename T1, typename T2,
   require_all_complex_t<T1, T2>* = nullptr>
 auto add(T1&& x, T2&& y) { return x + y; }
```

- `require_any_*_t` : Takes a parameter pack of types to enable if any of the types satisfy the check.

```cpp
 // Works if either T1 or T2 enherit from EigenBase
 template <typename T1, typename T2, require_any_eigen_t<T1, T2>* = nullptr>
 auto add(T1&& x, T2&& y) { return x + y; }
```

- `require_not_any_*_t` : Takes a parameter pack of types to enable if any one of the types are not satisfied.

```cpp
 // Works if either neither T1 or T2 are arithmetic
 template <typename T1, typename T2,
   require_not_any_eigen_row_vector_t<T1, T2>* = nullptr>
 auto add(T1 x, T2 y) { return x + y; }
```

- `require_not_all_*_t` : Takes a parameter pack of types to enable if all of the types are not satisfied.

```cpp
 // Works if neither T1 and T2 are arithmetic
 template <typename T1, typename T2,
   require_not_all_arithmetic_t<T1, T2>* = nullptr>
 auto add(T1 x, T2 y) { return x + y; }
```

`std::vector` and `Eigen` types have additional `requires` template parameters to detect if the @ref stan::value_type (the first underlying type) or the  @ref stan::scalar_type (the containers underlying scalar type) satisfy a condition to enable a class or function.

The container `requires` template parameters have an ending at their signature of _vt and _st to symbolize whether you want to inspect the @ref stan::value_type or @ref stan::scalar_type.
A function that accepts eigen matrices with floating point value types can be defined as

```cpp
 template <typename Mat1, typename Mat2,
   require_all_eigen_vt<std::is_floating_point, Mat1, Mat2>* = nullptr>
 auto add(Mat1&& A, Mat2&& B) { return A + B;}
```

A function that accepts standard vectors of Eigen vectors whose scalar type is @ref stan::math::var types can be defined as

```cpp
 template <typename Vec1, typename Vec2,
   require_all_std_vector_vt<is_eigen_vector, Vec1, Vec2>* = nullptr,
   require_all_std_vector_st<is_var, Vec1, Vec2>* = nullptr>
 auto add(Vec1&& A, Vec2&& B) {
   std::vector<decltype<A[0] + B[0]>> return_vec;
   std::transform(A.begin(), A.end(), B.begin(), return_vec.begin(),
     [](auto&& x, auto&& y) {
         return x + y;
     });
   return return_vec;
 }
```

There are also `requires` template parameters for generically checking if a type's @ref stan::value_type or @ref stan::scalar_type is correct.
To differentiate them from the Eigen and standard library vector checks the `vt` and `st` come *before* the type such as `require_vt_var<T>` which checks if a type `T`'s @ref stan::value_type satisfies @ref stan::is_var.

The `requires` template parameters type traits allow Stan to have more generic types so that the library can forward Eigen expression and have better move semantics.
For instance, the code below will accept any arbitrary Eigen expression that, if it's an rvalue, can be forwarded to another function.

```cpp
 template <typename Mat1, typename Mat2,
  require_all_eigen_vt<is_arithmetic, Mat1, Mat2>* = nullptr>
 inline auto a_func(Mat1&& m1, Mat2&& m2) {
   check_not_nan(m1);
   check_not_nan(m2);
   return another_func(std::forward<Mat1>(m1), std::forward<Mat2>(m2));
```

#### Adding a new requires

For a full list of predefined `requires` template parameters please see @ref require_meta.

If you are adding a new type trait that contains a `bool` member named value, you can add it to the set of known `requires` template parameters by using the macros defined in @ref stan/math/prim/meta/require_helpers.hpp .

For an example, we will use the `is_double_only` type trait below that only has a `bool value` member equal to `true` if the type `T` is `double`.

```cpp
template <typename T>
struct is_double_only
    : bool_constant<std::is_same<double, std::decay_t<T>>::value> {};

```

We can add the associated `requires` template parameters for the standard requires of `require_double_only_t` using the @ref STAN_ADD_REQUIRE_UNARY macro.
We supply the name to use in the require (`double_only`), The type trait name (`is_double_only`) and which group this type trait is in `require_stan_scalar_real`.

```cpp
STAN_ADD_REQUIRE_UNARY(double_only, is_double_only,
                       require_stan_scalar_real);
```

Similary, we can use the @ref STAN_ADD_REQUIRE_UNARY_INNER macro to get the `requires` template parameters such as `require_vt_double_only`.

```cpp
STAN_ADD_REQUIRE_UNARY_INNER(double_only, is_double_only,
                             require_stan_scalar_real);
```
