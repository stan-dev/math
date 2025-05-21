#ifndef STAN_MATH_PRIM_META_IS_EIGEN_HPP
#define STAN_MATH_PRIM_META_IS_EIGEN_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/is_eigen_matrix_base.hpp>
#include <stan/math/prim/meta/disjunction.hpp>
#include <stan/math/prim/meta/scalar_type.hpp>
#include <stan/math/prim/meta/value_type.hpp>
#include <type_traits>

namespace stan {

/**
 * Check if type derives from `EigenBase`
 * @tparam T Type to check if it is derived from `EigenBase`
 * @tparam Enable used for SFINAE deduction.
 * @ingroup type_trait
 **/
template <typename T>
struct is_eigen
    : bool_constant<is_base_pointer_convertible<Eigen::EigenBase, T>::value> {};

namespace internal {
// primary template handles types that have no nested ::type member:
template <class, class = void>
struct has_internal_trait : std::false_type {};

// specialization recognizes types that do have a nested ::type member:
template <class T>
struct has_internal_trait<T,
                          std::void_t<Eigen::internal::traits<std::decay_t<T>>>>
    : std::true_type {};

// primary template handles types that have no nested ::type member:
template <class, class = void>
struct has_scalar_trait : std::false_type {};

// specialization recognizes types that do have a nested ::type member:
template <class T>
struct has_scalar_trait<T, std::void_t<typename std::decay_t<T>::Scalar>>
    : std::true_type {};

}  // namespace internal

/**
 * Template metaprogram defining the base scalar type of
 * values stored in an Eigen matrix.
 *
 * @tparam T type to check.
 * @ingroup type_trait
 */
template <typename T>
struct scalar_type<T,
                   std::enable_if_t<is_eigen<T>::value
                                    && internal::has_scalar_trait<T>::value>> {
  using type = scalar_type_t<typename std::decay_t<T>::Scalar>;
};

/**
 * Template metaprogram defining the type of values stored in an
 * Eigen matrix, vector, or row vector.
 *
 * @tparam T type to check
 * @ingroup type_trait
 */
template <typename T>
struct value_type<T,
                  std::enable_if_t<is_eigen<T>::value
                                   && internal::has_scalar_trait<T>::value>> {
  using type = typename std::decay_t<T>::Scalar;
};

/**
 * Template metaprogram defining the base scalar type of
 * values stored in an Eigen matrix.
 *
 * @tparam T type to check.
 * @ingroup type_trait
 */
template <typename T>
struct scalar_type<T,
                   std::enable_if_t<is_eigen<T>::value
                                    && !internal::has_scalar_trait<T>::value>> {
  using type = scalar_type_t<
      typename Eigen::internal::traits<std::decay_t<T>>::Scalar>;
};

/**
 * Template metaprogram defining the type of values stored in an
 * Eigen matrix, vector, or row vector.
 *
 * @tparam T type to check
 * @ingroup type_trait
 */
template <typename T>
struct value_type<T,
                  std::enable_if_t<is_eigen<T>::value
                                   && !internal::has_scalar_trait<T>::value>> {
  using type = typename Eigen::internal::traits<std::decay_t<T>>::Scalar;
};

/*! \ingroup require_eigens_types */
/*! \defgroup eigen_types eigen  */
/*! \addtogroup eigen_types */
/*! @{ */

/*! \brief Require type satisfies @ref is_eigen */
/*! @tparam T the type to check */
template <typename T>
using require_eigen_t = require_t<is_eigen<std::decay_t<T>>>;

/*! \brief Require type does not satisfy @ref is_eigen */
/*! @tparam T the type to check */
template <typename T>
using require_not_eigen_t = require_not_t<is_eigen<std::decay_t<T>>>;

/*! \brief Require all of the types satisfy @ref is_eigen */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_all_eigen_t = require_all_t<is_eigen<std::decay_t<Types>>...>;

/*! \brief Require any of the types satisfy @ref is_eigen */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_any_eigen_t = require_any_t<is_eigen<std::decay_t<Types>>...>;

/*! \brief Require none of the types satisfy @ref is_eigen */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_all_not_eigen_t
    = require_all_not_t<is_eigen<std::decay_t<Types>>...>;

/*! \brief Require at least one of the types do not satisfy @ref is_eigen */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_any_not_eigen_t
    = require_any_not_t<is_eigen<std::decay_t<Types>>...>;
/*! @} */

/*! \ingroup require_eigens_types */
/*! \defgroup eigen_types eigen  */
/*! \addtogroup eigen_types */
/*! @{ */

/*! \brief Require type satisfies @ref is_eigen */
/*! and value type satisfies `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the value type against */
/*! @tparam Check The type to test @ref is_eigen for and whose @ref value_type
 * is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_eigen_vt = require_t<
    container_type_check_base<is_eigen, value_type_t, TypeCheck, Check...>>;

/*! \brief Require type does not satisfy @ref is_eigen or */
/*! value type does not satisfy `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the value type against */
/*! @tparam Check The type to test @ref is_eigen for and whose @ref value_type
 * is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_not_eigen_vt = require_not_t<
    container_type_check_base<is_eigen, value_type_t, TypeCheck, Check...>>;

/*! \brief Require any of the types satisfy @ref is_eigen */
/*! and any of the value types satisfy `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the value type against */
/*! @tparam Check The type to test @ref is_eigen for and whose @ref value_type
 * is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_any_eigen_vt = require_any_t<
    container_type_check_base<is_eigen, value_type_t, TypeCheck, Check>...>;

/*! \brief Require at least one of the types does not satisfy @ref is_eigen */
/*! and none of the value types satisfy `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the value type against */
/*! @tparam Check The type to test @ref is_eigen for and whose @ref value_type
 * is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_any_not_eigen_vt = require_any_not_t<
    container_type_check_base<is_eigen, value_type_t, TypeCheck, Check>...>;

/*! \brief Require all of the types satisfy @ref is_eigen */
/*! and all of the value types satisfy `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the value type against */
/*! @tparam Check The type to test @ref is_eigen for and whose @ref value_type
 * is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_all_eigen_vt = require_all_t<
    container_type_check_base<is_eigen, value_type_t, TypeCheck, Check>...>;

/*! \brief Require none of the types satisfy @ref is_eigen */
/*! and none of the value types satisfy `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the value type against */
/*! @tparam Check The type to test @ref is_eigen for and whose @ref value_type
 * is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_all_not_eigen_vt = require_all_not_t<
    container_type_check_base<is_eigen, value_type_t, TypeCheck, Check>...>;

/*! \brief Require type satisfies @ref is_eigen */
/*! and scalar type satisfies `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the scalar type against */
/*! @tparam Check The type to test @ref is_eigen for and whose @ref scalar_type
 * is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_eigen_st = require_t<
    container_type_check_base<is_eigen, scalar_type_t, TypeCheck, Check...>>;

/*! \brief Require type does not satisfy @ref is_eigen */
/*! or scalar type does not satisfy `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the scalar type against */
/*! @tparam Check The type to test @ref is_eigen for and whose @ref scalar_type
 * is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_not_eigen_st = require_not_t<
    container_type_check_base<is_eigen, scalar_type_t, TypeCheck, Check...>>;

/*! \brief Require any of the types satisfy @ref is_eigen */
/*! and any scalar type satisfies `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the scalar type against */
/*! @tparam Check The type to test @ref is_eigen for and whose @ref scalar_type
 * is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_any_eigen_st = require_any_t<
    container_type_check_base<is_eigen, scalar_type_t, TypeCheck, Check>...>;

/*! \brief Require at least one of the types does not satisfy @ref is_eigen */
/*! and any scalar type does not satisfy `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the scalar type against */
/*! @tparam Check The type to test @ref is_eigen for and whose @ref scalar_type
 * is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_any_not_eigen_st = require_any_not_t<
    container_type_check_base<is_eigen, scalar_type_t, TypeCheck, Check>...>;

/*! \brief Require all of the types does not satisfy @ref is_eigen */
/*! and all scalar types satisfy `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the scalar type against */
/*! @tparam Check The type to test @ref is_eigen for and whose @ref scalar_type
 * is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_all_eigen_st = require_all_t<
    container_type_check_base<is_eigen, scalar_type_t, TypeCheck, Check>...>;

/*! \brief Require none of the types satisfy @ref is_eigen */
/*! and none of the scalar types satisfy `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the scalar type against */
/*! @tparam Check The type to test @ref is_eigen for and whose @ref scalar_type
 * is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_all_not_eigen_st = require_all_not_t<
    container_type_check_base<is_eigen, scalar_type_t, TypeCheck, Check>...>;
/*! @} */

/**
 * Check if a type is derived from `Eigen::ArrayBase`
 * @tparam T type to check
 * @ingroup type_trait
 */
template <typename T>
struct is_eigen_array
    : bool_constant<is_base_pointer_convertible<Eigen::ArrayBase, T>::value> {};

/*! \ingroup require_eigens_types */
/*! \defgroup eigen_array_types eigen_array  */
/*! \addtogroup eigen_array_types */
/*! @{ */

/*! \brief Require type satisfies @ref is_eigen_array */
/*! @tparam T the type to check */
template <typename T>
using require_eigen_array_t = require_t<is_eigen_array<std::decay_t<T>>>;

/*! \brief Require type does not satisfy @ref is_eigen_array */
/*! @tparam T the type to check */
template <typename T>
using require_not_eigen_array_t
    = require_not_t<is_eigen_array<std::decay_t<T>>>;

/*! \brief Require any of the types satisfy @ref is_eigen_array */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_any_eigen_array_t
    = require_any_t<is_eigen_array<std::decay_t<Types>>...>;
/*! @} */

/**
 * Check if a type is derived from `Eigen::MatrixBase` or `Eigen::ArrayBase`
 * @tparam T type to check.
 * @ingroup type_trait
 */
template <typename T>
using is_eigen_matrix_or_array
    = math::disjunction<is_eigen_matrix_base<T>, is_eigen_array<T>>;

/*! \ingroup require_eigens_types */
/*! \defgroup eigen_array_types eigen_array  */
/*! \addtogroup eigen_array_types */
/*! @{ */

/*! \brief Require type satisfies @ref is_eigen_array */
/*! and value type satisfies `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the value type against */
/*! @tparam Check The type to test @ref is_eigen_array for and whose @ref
 * value_type is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_eigen_array_vt
    = require_t<container_type_check_base<is_eigen_array, value_type_t,
                                          TypeCheck, Check...>>;
/*! @} */

namespace internal {
template <typename T>
struct is_eigen_contiguous_map_impl : std::false_type {};
template <typename T, int Opts>
struct is_eigen_contiguous_map_impl<Eigen::Map<T, Opts, Eigen::Stride<0, 0>>>
    : std::true_type {};

}  // namespace internal

/**
 * Check if a type is an `Eigen::Map` with contiguous stride
 * @ingroup type_trait
 */
template <typename T>
struct is_eigen_contiguous_map
    : internal::is_eigen_contiguous_map_impl<std::decay_t<T>> {};

}  // namespace stan
#endif
