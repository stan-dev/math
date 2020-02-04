#ifndef STAN_MATH_HPP
#define STAN_MATH_HPP

/**
 * \defgroup prob_dists Probability Distributions
 */

/**
 * \ingroup prob_dists
 * \defgroup multivar_dists Multivariate Distributions
 * Distributions with Matrix inputs
 */
/**
 * \ingroup prob_dists
 * \defgroup univar_dists Univariate Distributions
 * Distributions with scalar, vector, or array input.
 */

/**
 * \defgroup type_trait Type Traits
 * The type traits in Stan math are a mix of custom traits for detecting
 *  value and container types of Eigen matrices and standard vectors as well
 *  as backports of C++17 traits.
 */

/**
 * \ingroup type_traits
 * \defgroup require_meta Pseudo-Concepts with requires<>
 * The requires type traits are wrappers for turning on and off definitions
 *  during compile time. You can think of these as "legacy concepts" aka I think
 *  the closest thing you can get to C++20 concepts without having the requires
 *  keyword.
 *
 * Legacy concepts are done through enable_if_t aliases which are called
 *  require_*. There are 5 basic types of requires.
 *
 *    requires_*_t: This means the function requires a certain type to turn on.
 *For instance
 *~~~~~{.cpp}
 * // Works for arithmetic types
 * template <typename T, require_arithmetic_t<T>...>
 * auto add(T&& x, T&& y) { return x + y; }
 *~~~~~
 *    require_not_*_t : When we want to disable for a type
 *~~~~~{.cpp}
 * // Works for non-arithmetic types
 * template <typename T, require_not_arithmetic_t<T>...>
 * auto add(T&& x, T&& y) { return x + y; }
 *~~~~~
 *    require_all_*_t : Takes a parameter pack of types to enable if all types
 *satisfy the check.
 *~~~~~{.cpp}
 * // Works if T1 and T2 are arithmetic
 * template <typename T1, typename T2, require_all_arithmetic_t<T1, T2>...>
 * auto add(T1&& x, T2&& y) { return x + y; }
 *~~~~~
 *    require_any_*_t : Takes a parameter pack of types to enable if any of the
 *types satisfy the check.
 *~~~~~{.cpp}
 * // Works if either T1 or T2 enherit from EigenBase
 * template <typename T1, typename T2, require_any_eigen_t<T1, T2>...>
 * auto add(T1&& x, T2&& y) { return x + y; }
 *~~~~~
 *    require_not_any_*_t : Takes a parameter pack of types to enable if any one
 *of the types are not satisfied.
 *~~~~~{.cpp}
 * // Works if either neither T1 or T2 are arithmetic
 * template <typename T1, typename T2, require_not_any_arithmetic_t<T1, T2>...>
 * auto add(T1 x, T2 y) { return x + y; }
 *~~~~~
 *    require_not_all_*_t : Takes a parameter pack of types to enable if all of
 *the types are not satisfied.
 *~~~~~{.cpp}
 * // Works if neither T1 and T2 are arithmetic
 * template <typename T1, typename T2, require_not_all_arithmetic_t<T1, T2>...>
 * auto add(T1 x, T2 y) { return x + y; }
 *~~~~~
 *
 * In addition, std::vector and Eigen types have additional requires to detect
 *if the value_type (the first underlying type) or the scalar_type (the
 *containers underlying scalar type) satisfy a condition to enable a class or
 *function.
 *
 * The container requires have a _vt and _st to symbolize the above. A function
 *that accepts eigen matrices of whose whose value type is floating point types
 *can be defined as
 *
 *~~~~~{.cpp}
 * template <typename Mat1, typename Mat2,
 * require_all_eigen_vt<std::is_floating_point, Mat1, Mat2>...>
 * auto add(Mat1&& A, Mat2&& B) { return A + B;}
 *~~~~~
 *  A function that accepts standard vectors of Eigen vectors whose scalar type
 *is floating point types can be defined as
 *
 *~~~~~{.cpp}
 * template <typename Vec1, typename Vec2,
 *   require_all_std_vector_vt<is_eigen_vector, Vec1, Vec2>...,
 *   require_all_std_vector_st<is_floating_point, Vec1, Vec2>...>
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
 * In general, these methods allow Stan to have more generic types so that the
 * library can forward along Eigen expression and have better move semantics.
 * For instance, the code below will accept any arbitrary Eigen expression
 * that, if it's an rvalue, can be forwarded to another function.
 *
 *~~~~~{.cpp}
 * template <typename Mat1, typename Mat2,
 * require_all_eigen_vt<is_arithmetic, Mat1, Mat2>...>
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
 * \defgroup require_base_types Basic Type Pseudo-Concepts
 * These type traits check the type properties of simple objects and are denoted
 * with a _t at the end. They check only the direct type and do not inspect
 * anything about type properties composed within a type.
 */

/**
 * \ingroup require_meta
 * \defgroup require_container_types Container Pseudo-Concepts
 * These type traits check the type properties of objects that act as
 * containers. The @c _vt methods check that a containers @c value_type fulfills
 * certain conditions while @c _st methods check that the objects @c scalar_type
 * fulfills certain conditions. @c value_type and @c scalar_type differ in that
 * @c value_type is the first level of a container while @c scalar_type
 * recursively goes through containers of containers till it comes to a simple
 * type.
 */

#include <stan/math/rev.hpp>

#endif
