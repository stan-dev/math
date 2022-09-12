#ifndef STAN_MATH_PRIM_META_HPP
#define STAN_MATH_PRIM_META_HPP

/**
 * \defgroup type_trait Type Traits
 * The type traits in Stan math are a mix of custom traits for detecting
 * value and container types of Eigen matrices, standard vectors, standard
 * complex numbers, and backports of C++17 traits.
 */

/**
 * \ingroup type_traits
 * \defgroup require_meta Pseudo-Concepts with requires<>
 * The requires type traits are wrappers for turning on and off definitions
 *  during compile time. You can think of these as "legacy concepts" aka I think
 *  the closest thing you can get to C++20 concepts without having the requires
 *  keyword.
 *
 * Legacy concepts are done through `enable_if_t` aliases which are called
 *  `require_*`. These are either defined in templates as an anonymous template
 * typename parameter or as an anonymous void* parameter with a default value
 * of `nullptr` There are 5 basic types of requires.
 *
 *    `requires_*_t`: This means the function requires a certain type to turn
 *on. For instance
 *~~~~~{.cpp}
 * // Works for arithmetic types
 * template <typename T, require_arithmetic_t<T>* = nullptr>
 * auto add(T&& x, T&& y) { return x + y; }
 *~~~~~
 *    `require_not_*_t` : When we want to disable for a type
 *~~~~~{.cpp}
 * // Works for non-arithmetic types
 * template <typename T, require_not_arithmetic_t<T>* = nullptr>
 * auto add(T&& x, T&& y) { return x + y; }
 *~~~~~
 *    `require_all_*_t` : Takes a parameter pack of types to enable if all types
 *satisfy the check.
 *~~~~~{.cpp}
 * // Works if T1 and T2 are arithmetic
 * template <typename T1, typename T2,
 *   require_all_arithmetic_t<T1, T2>* = nullptr>
 * auto add(T1&& x, T2&& y) { return x + y; }
 *~~~~~
 *    `require_any_*_t` : Takes a parameter pack of types to enable if any of
 *the types satisfy the check.
 *~~~~~{.cpp}
 * // Works if either T1 or T2 enherit from EigenBase
 * template <typename T1, typename T2, require_any_eigen_t<T1, T2>* = nullptr>
 * auto add(T1&& x, T2&& y) { return x + y; }
 *~~~~~
 *    `require_not_any_*_t` : Takes a parameter pack of types to enable if any
 *one of the types are not satisfied.
 *~~~~~{.cpp}
 * // Works if either neither T1 or T2 are arithmetic
 * template <typename T1, typename T2,
 *   require_not_any_arithmetic_t<T1, T2>* = nullptr>
 * auto add(T1 x, T2 y) { return x + y; }
 *~~~~~
 *    `require_not_all_*_t` : Takes a parameter pack of types to enable if all
 *of the types are not satisfied.
 *~~~~~{.cpp}
 * // Works if neither T1 and T2 are arithmetic
 * template <typename T1, typename T2,
 *   require_not_all_arithmetic_t<T1, T2>* = nullptr>
 * auto add(T1 x, T2 y) { return x + y; }
 *~~~~~
 *
 * `std::vector` and `Eigen` types have additional requires to
 *detect if the `value_type` (the first underlying type) or the `scalar_type`
 *(the containers underlying scalar type) satisfy a condition to enable a class
 *or function.
 *
 * The container requires have a _vt and _st to symbolize the above. A function
 *that accepts eigen matrices of whose whose value type is floating point types
 *can be defined as
 *
 *~~~~~{.cpp}
 * template <typename Mat1, typename Mat2,
 *   require_all_eigen_vt<std::is_floating_point, Mat1, Mat2>* = nullptr>
 * auto add(Mat1&& A, Mat2&& B) { return A + B;}
 *~~~~~
 *  A function that accepts standard vectors of Eigen vectors whose scalar type
 *is floating point types can be defined as
 *
 *~~~~~{.cpp}
 * template <typename Vec1, typename Vec2,
 *   require_all_std_vector_vt<is_eigen_vector, Vec1, Vec2>* = nullptr,
 *   require_all_std_vector_st<is_floating_point, Vec1, Vec2>* = nullptr>
 * auto add(Vec1&& A, Vec2&& B) {
 *   std::vector<decltype<A[0] + B[0]>> return_vec;
 *   std::transform(A.begin(), A.end(), B.begin(), return_vec.begin(),
 *     [](auto&& x, auto&& y) {
 *         return x + y;
 *     });
 *   return return_vec;
 * }
 *~~~~~
 *
 * There are also requires for generically checking if a type's `value_type` or
 *  `scalar_type` is correct. To differentiate them from the Eigen and Standard
 * vector checks the `vt` and `st` comes *before* the type such as
 * `require_vt_var<T>` which checks if a type T's `value_type` satisfies
 * `is_var`.
 *
 * The `requires` type traits allow Stan to have more generic types so that the
 * library can forward along Eigen expression and have better move semantics.
 * For instance, the code below will accept any arbitrary Eigen expression
 * that, if it's an rvalue, can be forwarded to another function.
 *
 *~~~~~{.cpp}
 * template <typename Mat1, typename Mat2,
 * require_all_eigen_vt<is_arithmetic, Mat1, Mat2>* = nullptr>
 * inline auto a_func(Mat1&& m1, Mat2&& m2) {
 *   check_nan(m1);
 *   check_nan(m2);
 *   // If m1 and/or m2 is an rvalue it will be moved over to this function
 *   // instead of copied to an lvalue
 *   auto B = another_func(std::forward<Mat1>(m1), std::forward<Mat2>(m2)); //
 *(3) return B;
 *~~~~~
 */

/**
 * \ingroup require_meta
 * \defgroup require_stan_scalar Scalar types
 * `require` type traits for types that are either `arithmetic`, `var`, `fvar`,
 * or `Complex`.
 */

/**
 * \ingroup require_stan_scalar
 * \defgroup require_stan_scalar_real Real types
 * `require` type traits for types that are either `arithmetic`, `var`, or
 * `fvar`.
 */

/**
 * \ingroup require_stan_scalar
 * \defgroup require_stan_scalar_complex Complex types
 * `require` type traits for types that are `Complex<T>`.
 */

/**
 * \ingroup require_meta
 * \defgroup require_std Standard library types and traits
 * `require` type traits that come from the standard library.
 */

/**
 * \ingroup require_meta
 * \defgroup require_eigens_types Eigen
 * `require` type traits to detect Eigen types.
 */

/**
 * \ingroup require_meta
 * \defgroup general_types General Types
 * `require` type traits for general types.
 */

/**
 * \ingroup require_meta
 * \defgroup macro_helpers Require Macro Generators
 * These macros are used on type traits to define the set of `requires`
 */

/**
 * \ingroup require_meta
 * \defgroup require_opencl_types OpenCL
 * `require` type traits to detect types used with OpenCL.
 */

#include <stan/math/prim/meta/compiler_attributes.hpp>
#include <stan/math/prim/meta/ad_promotable.hpp>
#include <stan/math/prim/meta/append_return_type.hpp>
#include <stan/math/prim/meta/base_type.hpp>
#include <stan/math/prim/meta/child_type.hpp>
#include <stan/math/prim/meta/contains_fvar.hpp>
#include <stan/math/prim/meta/contains_std_vector.hpp>
#include <stan/math/prim/meta/error_index.hpp>
#include <stan/math/prim/meta/forward_as.hpp>
#include <stan/math/prim/meta/holder.hpp>
#include <stan/math/prim/meta/include_summand.hpp>
#include <stan/math/prim/meta/index_type.hpp>
#include <stan/math/prim/meta/index_apply.hpp>
#include <stan/math/prim/meta/is_autodiff.hpp>
#include <stan/math/prim/meta/is_arena_matrix.hpp>
#include <stan/math/prim/meta/is_base_pointer_convertible.hpp>
#include <stan/math/prim/meta/is_dense_dynamic.hpp>
#include <stan/math/prim/meta/is_double_or_int.hpp>
#include <stan/math/prim/meta/is_complex.hpp>
#include <stan/math/prim/meta/is_constant.hpp>
#include <stan/math/prim/meta/is_container.hpp>
#include <stan/math/prim/meta/is_container_or_var_matrix.hpp>
#include <stan/math/prim/meta/is_eigen.hpp>
#include <stan/math/prim/meta/is_eigen_dense_base.hpp>
#include <stan/math/prim/meta/is_eigen_dense_dynamic.hpp>
#include <stan/math/prim/meta/is_eigen_matrix.hpp>
#include <stan/math/prim/meta/is_eigen_matrix_base.hpp>
#include <stan/math/prim/meta/is_eigen_sparse_base.hpp>
#include <stan/math/prim/meta/is_fvar.hpp>
#include <stan/math/prim/meta/is_kernel_expression.hpp>
#include <stan/math/prim/meta/is_matrix_cl.hpp>
#include <stan/math/prim/meta/is_matrix.hpp>
#include <stan/math/prim/meta/is_plain_type.hpp>
#include <stan/math/prim/meta/is_string_convertible.hpp>
#include <stan/math/prim/meta/is_tuple.hpp>
#include <stan/math/prim/meta/is_var.hpp>
#include <stan/math/prim/meta/is_var_and_matrix_types.hpp>
#include <stan/math/prim/meta/is_var_matrix.hpp>
#include <stan/math/prim/meta/is_var_dense_dynamic.hpp>
#include <stan/math/prim/meta/is_var_eigen.hpp>
#include <stan/math/prim/meta/is_rev_matrix.hpp>
#include <stan/math/prim/meta/is_vari.hpp>
#include <stan/math/prim/meta/is_var_or_arithmetic.hpp>
#include <stan/math/prim/meta/is_vector.hpp>
#include <stan/math/prim/meta/is_vector_like.hpp>
#include <stan/math/prim/meta/is_stan_scalar.hpp>
#include <stan/math/prim/meta/is_stan_scalar_or_eigen.hpp>
#include <stan/math/prim/meta/partials_return_type.hpp>
#include <stan/math/prim/meta/partials_type.hpp>
#include <stan/math/prim/meta/plain_type.hpp>
#include <stan/math/prim/meta/possibly_sum.hpp>
#include <stan/math/prim/meta/promote_args.hpp>
#include <stan/math/prim/meta/promote_scalar_type.hpp>
#include <stan/math/prim/meta/ref_type.hpp>
#include <stan/math/prim/meta/require_generics.hpp>
#include <stan/math/prim/meta/return_type.hpp>
#include <stan/math/prim/meta/scalar_type.hpp>
#include <stan/math/prim/meta/scalar_type_pre.hpp>
#include <stan/math/prim/meta/seq_view.hpp>
#include <stan/math/prim/meta/static_select.hpp>
#include <stan/math/prim/meta/value_type.hpp>
#include <stan/math/prim/meta/void_t.hpp>
#include <stan/math/prim/meta/StdVectorBuilder.hpp>
#include <stan/math/prim/meta/VectorBuilder.hpp>

#endif
