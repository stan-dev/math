#ifndef STAN_MATH_PRIM_META_RETURN_TYPE_HPP
#define STAN_MATH_PRIM_META_RETURN_TYPE_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/promote_args.hpp>
#include <stan/math/prim/meta/base_type.hpp>
#include <stan/math/prim/meta/scalar_type.hpp>
#include <complex>
#include <vector>

namespace stan {

/**
 * Provides a member type alias named `type`, the value of which is
 * the least type under Stan's assignability relation that can be
 * assigned a `double` and all of the base types of the specified
 * arguments after removing qualifiers (`const` and `volatile`) and
 * decaying (lvalue to rvalue by removing references) and array to
 * pointer).   This type trait is used to calculate the return type
 * of real-valued functions involving heterogeneous arguments.
 *
 * @tparam Ts sequence of template parameter types
 * @see return_type for a definition of the type ordering
 * @see base_type for definition of base type extraction
 * @ingroup type_trait
 */
template <typename... Ts>
struct real_return {
  using type = double;
};

template <typename T, typename... Ts>
struct real_return<T, Ts...> {
  using type
      = promote_args_t<base_type_t<T>, typename real_return<Ts...>::type>;
};

/**
 * Convenience type to calculate the real return type.
 *
 * @tparam Ts sequence of argument types
 * @see real_return
 * @ingroup type_trait
 */
template <typename... Ts>
using real_return_t = typename real_return<Ts...>::type;

/**
 * Convenience type to calculate the complex return type, which wraps
 * `std::complex` around the return type of the specified template parameters.
 *
 * @tparam Ts sequence of argument types
 * @see real_return
 * @ingroup type_trait
 */
template <typename... Ts>
using complex_return_t = std::complex<real_return_t<Ts...>>;

/**
 * Convenience type to calculate the complex return type, which wraps
 * `std::vector` around the return type of the specified template parameters.
 *
 * @tparam Ts sequence of argument types
 * @see real_return
 * @ingroup type_trait
 */
template <typename... Ts>
using std_vector_return_t = std::vector<real_return_t<Ts...>>;

/**
 * Convenience type to calculate the complex return type, which wraps
 * `Eigen::Matrix< , -1, -1>` around the return type of the specified
 * template parameters.
 *
 * @tparam Ts sequence of argument types
 * @see real_return
 * @ingroup type_trait
 */
template <typename... Ts>
using matrix_return_t = Eigen::Matrix<real_return_t<Ts...>, -1, -1>;

/**
 * Convenience type to calculate the complex return type, which wraps
 * `Eigen::Matrix< , -1, 1>` around the return type of the specified
 * template parameters.
 *
 * @tparam Ts sequence of argument types
 * @see real_return
 * @ingroup type_trait
 */
template <typename... Ts>
using vector_return_t = Eigen::Matrix<real_return_t<Ts...>, -1, 1>;

/**
 * Convenience type to calculate the complex return type, which wraps
 * `Eigen::Matrix< , 1, -1>` around the return type of the specified
 * template parameters.
 *
 * @tparam Ts sequence of argument types
 * @see real_return
 * @ingroup type_trait
 */
template <typename... Ts>
using row_vector_return_t = Eigen::Matrix<real_return_t<Ts...>, 1, -1>;

/**
 * Defines a member type named `type` that is the least scalar type to
 * which both template parameter scalar types are assignable in Stan.
 * The ordering of types for which this is the least-upper-bound
 * operation is defined in the class documentation for `return_type`.
 *
 * @tparam T1 first scalar type
 * @tparam T2 second scalar type
 * @see return_type
 * @ingroup type_trait
 */
template <typename T1, typename T2>
struct scalar_lub {
  using type = promote_args_t<T1, T2>;
};

template <typename T1, typename T2>
struct scalar_lub<std::complex<T1>, T2> {
  using type = std::complex<promote_args_t<T1, T2>>;
};

template <typename T1, typename T2>
struct scalar_lub<T1, std::complex<T2>> {
  using type = std::complex<promote_args_t<T1, T2>>;
};

template <typename T1, typename T2>
struct scalar_lub<std::complex<T1>, std::complex<T2>> {
  using type = std::complex<promote_args_t<T1, T2>>;
};

/**
 * Convenience type for the least upper bound of the specified
 * template parameters in Stan's assignment ordering.
 *
 * @param T1 first type
 * @param T2 second type
 * @see scalar_lub
 * @ingroup type_trait
 */
template <typename T1, typename T2>
using scalar_lub_t = typename scalar_lub<T1, T2>::type;

/**
 * Template metaprogram to calculate the base scalar return type
 * resulting from promoting all the scalar types of the template
 * parameters to the least type to which all the base types of the
 * arguments are assignable. The metaprogram can take an arbitrary
 * number of template parameters.
 *
 * Complex numbers (instances of `std::complex<T>`) are considered
 * scalars.  All C++ primitive types (except `long double`) are
 * automatically promoted to `double`.
 *
 * The set of autodiff types is defined to be the smallest such that
 * - `var` is an autodiff type,
 * - `fvar<double>` is an autodiff type, and
 * - `fvar<T>` is an autodiff type if `T` is an autodiff type.

 * The set of scalar types is defined to be the smallest such that
 * - `double` and `long double` are scalar types,
 * - `T` is a scalar type if `T` is an autodiff type,
 * - `complex<double>` is a scalar type, and
 * - `complex<T>` is a scalar type if `T` is an autodiff type.
 *
 * The assignability relation among scalar types is defined to be
 * the smallest such that
 * - `double` is assignable to any type `T`,
 * - `T1` is assignable to `std::complex<T2>` if `T1` is assignable to
 * - `T2`, and
 * - `std::complex<T1>` is assignable to `std::complex<T2>` if `T1` is
 * assignable to `T2`.
 *
 * Example usage:
 *  - `return_type_t<int, double, float>` is `double`
 *  - `return_type_t<double, var>` is `var`
 *
 * @tparam Ts sequence of input types
 * @ingroup type_trait
 */
template <typename... Ts>
struct return_type {
  using type = double;
};

template <typename T, typename... Ts>
struct return_type<T, Ts...> {
  using type
      = scalar_lub_t<scalar_type_t<T>, typename return_type<Ts...>::type>;
};

/**
 * Convenience type for the return type of the specified template
 * parameters.
 *
 * @tparam Ts sequence of types
 * @see return_type
 * @ingroup type_trait
 */
template <typename... Ts>
using return_type_t = typename return_type<Ts...>::type;

/**
 * \ingroup require_stan_scalar_real
 * \addtogroup autodiff_types
 * \brief Require return type from parameter pack satisfies `Check`
 */
template <template <class...> class Check, typename... Ts>
using require_return_type_t = require_t<Check<return_type_t<Ts...>>>;

/**
 * \ingroup require_stan_scalar_real
 * \addtogroup autodiff_types
 * \brief Require return type from parameter pack does not satisfy `Check`
 */
template <template <class...> class Check, typename... Ts>
using require_not_return_type_t = require_not_t<Check<return_type_t<Ts...>>>;

}  // namespace stan
#endif
