#ifndef STAN_MATH_PRIM_META_IS_CONSTANT_HPP
#define STAN_MATH_PRIM_META_IS_CONSTANT_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/bool_constant.hpp>
#include <stan/math/prim/meta/conjunction.hpp>
#include <stan/math/prim/meta/is_vector.hpp>
#include <stan/math/prim/meta/is_eigen.hpp>
#include <stan/math/prim/meta/scalar_type.hpp>
#include <stan/math/prim/meta/value_type.hpp>
#include <stan/math/prim/meta/require_generics.hpp>
#include <type_traits>
#include <vector>
namespace stan {

/** \ingroup type_trait
 * Metaprogramming struct to detect whether a given type is constant
 * in the mathematical sense (not the C++ <code>const</code>
 * sense). If the parameter type is constant, <code>value</code>
 * will be equal to <code>true</code>.
 *
 * The baseline implementation in this abstract base class is to
 * classify a type <code>T</code> as constant if it can be converted
 * (i.e., assigned) to a <code>double</code>.  This baseline should
 * be overridden for any type that should be treated as a variable.
 *
 * @tparam T Type being tested.
 */
template <typename T, typename = void>
struct is_constant : bool_constant<std::is_convertible<T, double>::value> {};

/** \ingroup type_trait
 * Metaprogram defining an enum <code>value</code> which
 * is <code>true</code> if all of the type parameters
 * are constant (i.e., primitive types) and
 * <code>false</code> otherwise.
 */
template <typename... T>
using is_constant_all = math::conjunction<is_constant<T>...>;

/** \ingroup type_trait
 * Defines a static member named value and sets it to true
 * if the type of the elements in the provided std::vector
 * is constant, false otherwise. This is used in
 * the is_constant_all metaprogram.
 * @tparam type of the elements in the std::vector
 */
template <typename T>
struct is_constant<T, require_std_vector_t<T>>
    : bool_constant<is_constant<typename std::decay_t<T>::value_type>::value> {
};

/** \ingroup type_trait
 * Defines a public enum named value and sets it to true
 * if the type of the elements in the provided Eigen Matrix
 * is constant, false otherwise. This is used in
 * the is_constant_all metaprogram.
 *
 * @tparam T type of the Eigen Matrix
 */
template <typename T>
struct is_constant<T, require_eigen_t<T>>
    : bool_constant<is_constant<typename std::decay_t<T>::Scalar>::value> {};

//STAN_ADD_REQUIRE_UNARY(constant, is_constant, require_stan_scalar_real);
template <typename T>
using require_constant_t = require_t<is_constant<std::decay_t<T>>>;

template <typename T>
using require_not_constant_t
    = require_not_t<is_constant<std::decay_t<T>>>;

template <typename... Types>
using require_all_constant_t
    = require_all_t<is_constant<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_constant_t
    = require_any_t<is_constant<std::decay_t<Types>>...>;

template <typename... Types>
using require_all_not_constant_t
    = require_all_not_t<is_constant<std::decay_t<Types>>...>;

template <typename... Types>
using require_any_not_constant_t
    = require_any_not_t<is_constant<std::decay_t<Types>>...>;

  
//STAN_ADD_REQUIRE_UNARY_INNER(constant, is_constant, require_stan_scalar_real);
template <typename T>
using require_vt_constant
    = require_t<is_constant<value_type_t<std::decay_t<T>>>>;

template <typename T>
using require_not_vt_constant
    = require_not_t<is_constant<value_type_t<std::decay_t<T>>>>;

template <typename... Types>
using require_all_vt_constant
    = require_all_t<is_constant<value_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_vt_constant
    = require_any_t<is_constant<value_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_all_not_vt_constant
    = require_all_not_t<is_constant<value_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_not_vt_constant
    = require_any_not_t<is_constant<value_type_t<std::decay_t<Types>>>...>;

template <typename T>
using require_st_constant
    = require_t<is_constant<scalar_type_t<std::decay_t<T>>>>;

template <typename T>
using require_not_st_constant
    = require_not_t<is_constant<scalar_type_t<std::decay_t<T>>>>;

template <typename... Types>
using require_all_st_constant
    = require_all_t<is_constant<scalar_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_st_constant
    = require_any_t<is_constant<scalar_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_all_not_st_constant
    = require_all_not_t<is_constant<scalar_type_t<std::decay_t<Types>>>...>;

template <typename... Types>
using require_any_not_st_constant
    = require_any_not_t<is_constant<scalar_type_t<std::decay_t<Types>>>...>;


}  // namespace stan
#endif
