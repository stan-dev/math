#ifndef STAN_MATH_PRIM_META_IS_VECTOR_HPP
#define STAN_MATH_PRIM_META_IS_VECTOR_HPP

#include <stan/math/prim/meta/bool_constant.hpp>
#include <stan/math/prim/meta/is_eigen.hpp>
#include <stan/math/prim/meta/is_var.hpp>
#include <stan/math/prim/meta/scalar_type.hpp>
#include <stan/math/prim/meta/value_type.hpp>
#include <stan/math/prim/meta/require_helpers.hpp>
#include <type_traits>
#include <vector>

namespace stan {

/** \ingroup type_trait
 * Base implementation for checking if type is std vector
 */
template <typename T, typename = void>
struct is_std_vector : std::false_type {};

namespace internal {
/** \ingroup type_trait
 * Underlying implementation for detecting if an Eigen Matrix is a column
 * vector.
 */
template <typename T, bool = is_eigen<T>::value>
struct is_eigen_col_vector_impl
    : bool_constant<std::decay_t<T>::ColsAtCompileTime == 1> {};

/** \ingroup type_trait
 * Specialization for when type is not an eigen vector.
 */
template <typename T>
struct is_eigen_col_vector_impl<T, false> : std::false_type {};

/** \ingroup type_trait
 * Underlying implementation for detecting if an Eigen Matrix is a row vector.
 */
template <typename T, bool = is_eigen<T>::value>
struct is_eigen_row_vector_impl
    : std::integral_constant<bool, std::decay_t<T>::RowsAtCompileTime == 1> {};

/** \ingroup type_trait
 * Specialization for when type is not an eigen vector.
 */
template <typename T>
struct is_eigen_row_vector_impl<T, false> : std::false_type {};

/** \ingroup type_trait
 * Underlying implementation for detecting if a Matrix is a column vector.
 */
template <typename T, bool = bool_constant<is_eigen<T>::value
                                           || (is_eigen<value_type_t<T>>::value
                                               && is_var<T>::value)>::value>
struct is_col_vector_impl
    : bool_constant<std::decay_t<T>::ColsAtCompileTime == 1> {};

/** \ingroup type_trait
 * Specialization for when type is not a vector.
 */
template <typename T>
struct is_col_vector_impl<T, false> : std::false_type {};

/** \ingroup type_trait
 * Underlying implementation for detecting if a Matrix is a row vector.
 */
template <typename T, bool = bool_constant<is_eigen<T>::value
                                           || (is_eigen<value_type_t<T>>::value
                                               && is_var<T>::value)>::value>
struct is_row_vector_impl
    : std::integral_constant<bool, std::decay_t<T>::RowsAtCompileTime == 1> {};

/** \ingroup type_trait
 * Specialization for when type is not an vector.
 */
template <typename T>
struct is_row_vector_impl<T, false> : std::false_type {};

}  // namespace internal

/** \ingroup type_trait
 * If the input type T is an eigen matrix with 1 row at compile time this
 * has a static member with a value of true. Else this has a static
 * member with a value of false.
 */
template <typename T>
struct is_eigen_col_vector : internal::is_eigen_col_vector_impl<T> {};

/*! \ingroup require_eigens_types */
/*! \defgroup eigen_col_vector_types eigen_col_vector  */
/*! \addtogroup eigen_col_vector_types */
/*! @{ */

/*! \brief Require type satisfies @ref is_eigen_col_vector */
/*! @tparam T the type to check */
template <typename T>
using require_eigen_col_vector_t
    = require_t<is_eigen_col_vector<std::decay_t<T>>>;

/*! \brief Require type does not satisfy @ref is_eigen_col_vector */
/*! @tparam T the type to check */
template <typename T>
using require_not_eigen_col_vector_t
    = require_not_t<is_eigen_col_vector<std::decay_t<T>>>;

/*! \brief Require all of the types satisfy @ref is_eigen_col_vector */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_all_eigen_col_vector_t
    = require_all_t<is_eigen_col_vector<std::decay_t<Types>>...>;
/*! @} */

/*! \ingroup require_eigens_types */
/*! \defgroup eigen_col_vector_types eigen_col_vector  */
/*! \addtogroup eigen_col_vector_types */
/*! @{ */

/*! \brief Require type satisfies @ref is_eigen_col_vector */
/*! and value type satisfies `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the value type against */
/*! @tparam Check The type to test @ref is_eigen_col_vector for and whose @ref
 * value_type is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_eigen_col_vector_vt
    = require_t<container_type_check_base<is_eigen_col_vector, value_type_t,
                                          TypeCheck, Check...>>;
/*! @} */

/** \ingroup type_trait
 * If the input type T has a static comple time constant type
 * `ColsAtCompileTime` equal to 1 this
 * has a static member with a value of true. Else this has a static
 * member with a value of false.
 */
template <typename T>
struct is_col_vector : internal::is_col_vector_impl<T> {};

/*! \ingroup require_eigens_types */
/*! \defgroup col_vector_types col_vector  */
/*! \addtogroup col_vector_types */
/*! @{ */

/*! \brief Require type satisfies @ref is_col_vector */
/*! @tparam T the type to check */
template <typename T>
using require_col_vector_t = require_t<is_col_vector<std::decay_t<T>>>;

/*! \brief Require type does not satisfy @ref is_col_vector */
/*! @tparam T the type to check */
template <typename T>
using require_not_col_vector_t = require_not_t<is_col_vector<std::decay_t<T>>>;

/*! \brief Require all of the types satisfy @ref is_col_vector */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_all_col_vector_t
    = require_all_t<is_col_vector<std::decay_t<Types>>...>;
/*! @} */

/** \ingroup type_trait
 * If the input type T is an eigen matrix with 1 column at compile time this
 * has a static member with a value of true. Else this has a static
 * member with a value of false.
 */
template <typename T>
struct is_eigen_row_vector : internal::is_eigen_row_vector_impl<T> {};

/*! \ingroup require_eigens_types */
/*! \defgroup eigen_row_vector_types eigen_row_vector  */
/*! \addtogroup eigen_row_vector_types */
/*! @{ */

/*! \brief Require type satisfies @ref is_eigen_row_vector */
/*! @tparam T the type to check */
template <typename T>
using require_eigen_row_vector_t
    = require_t<is_eigen_row_vector<std::decay_t<T>>>;
/*! @} */

/** \ingroup type_trait
 * If the input type T has a static comple time constant type
 * `RowsAtCompileTime` equal to 1 this
 * has a static member with a value of true. Else this has a static
 * member with a value of false.
 */
template <typename T>
struct is_row_vector : internal::is_row_vector_impl<T> {};

/*! \ingroup require_eigens_types */
/*! \defgroup row_vector_types row_vector  */
/*! \addtogroup row_vector_types */
/*! @{ */

/*! \brief Require type satisfies @ref is_row_vector */
/*! @tparam T the type to check */
template <typename T>
using require_row_vector_t = require_t<is_row_vector<std::decay_t<T>>>;
/*! @} */

/** \ingroup type_trait
 * If the input type T is an eigen matrix with 1 column or 1 row at compile time
 * this has a static member with a value of true. Else this has a static
 * member with a value of false.
 */
template <typename T>
struct is_eigen_vector : bool_constant<is_eigen_col_vector<T>::value
                                       || is_eigen_row_vector<T>::value> {};

/*! \ingroup require_eigens_types */
/*! \defgroup eigen_vector_types eigen_vector  */
/*! \addtogroup eigen_vector_types */
/*! @{ */

/*! \brief Require type satisfies @ref is_eigen_vector */
/*! @tparam T the type to check */
template <typename T>
using require_eigen_vector_t = require_t<is_eigen_vector<std::decay_t<T>>>;

/*! \brief Require type does not satisfy @ref is_eigen_vector */
/*! @tparam T the type to check */
template <typename T>
using require_not_eigen_vector_t
    = require_not_t<is_eigen_vector<std::decay_t<T>>>;

/*! \brief Require all of the types satisfy @ref is_eigen_vector */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_all_eigen_vector_t
    = require_all_t<is_eigen_vector<std::decay_t<Types>>...>;

/*! \brief Require any of the types satisfy @ref is_eigen_vector */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_any_eigen_vector_t
    = require_any_t<is_eigen_vector<std::decay_t<Types>>...>;

/*! \brief Require none of the types satisfy @ref is_eigen_vector */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_all_not_eigen_vector_t
    = require_all_not_t<is_eigen_vector<std::decay_t<Types>>...>;

/*! \brief Require at least one of the types do not satisfy @ref is_eigen_vector
 */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_any_not_eigen_vector_t
    = require_any_not_t<is_eigen_vector<std::decay_t<Types>>...>;
/*! @} */

/*! \ingroup require_eigens_types */
/*! \defgroup eigen_vector_types eigen_vector  */
/*! \addtogroup eigen_vector_types */
/*! @{ */

/*! \brief Require type satisfies @ref is_eigen_vector */
/*! and value type satisfies `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the value type against */
/*! @tparam Check The type to test @ref is_eigen_vector for and whose @ref
 * value_type is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_eigen_vector_vt
    = require_t<container_type_check_base<is_eigen_vector, value_type_t,
                                          TypeCheck, Check...>>;

/*! \brief Require type does not satisfy @ref is_eigen_vector or */
/*! value type does not satisfy `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the value type against */
/*! @tparam Check The type to test @ref is_eigen_vector for and whose @ref
 * value_type is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_not_eigen_vector_vt
    = require_not_t<container_type_check_base<is_eigen_vector, value_type_t,
                                              TypeCheck, Check...>>;

/*! \brief Require any of the types satisfy @ref is_eigen_vector */
/*! and any of the value types satisfy `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the value type against */
/*! @tparam Check The type to test @ref is_eigen_vector for and whose @ref
 * value_type is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_any_eigen_vector_vt
    = require_any_t<container_type_check_base<is_eigen_vector, value_type_t,
                                              TypeCheck, Check>...>;

/*! \brief Require at least one of the types does not satisfy @ref
 * is_eigen_vector */
/*! and none of the value types satisfy `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the value type against */
/*! @tparam Check The type to test @ref is_eigen_vector for and whose @ref
 * value_type is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_any_not_eigen_vector_vt
    = require_any_not_t<container_type_check_base<is_eigen_vector, value_type_t,
                                                  TypeCheck, Check>...>;

/*! \brief Require all of the types satisfy @ref is_eigen_vector */
/*! and all of the value types satisfy `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the value type against */
/*! @tparam Check The type to test @ref is_eigen_vector for and whose @ref
 * value_type is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_all_eigen_vector_vt
    = require_all_t<container_type_check_base<is_eigen_vector, value_type_t,
                                              TypeCheck, Check>...>;

/*! \brief Require none of the types satisfy @ref is_eigen_vector */
/*! and none of the value types satisfy `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the value type against */
/*! @tparam Check The type to test @ref is_eigen_vector for and whose @ref
 * value_type is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_all_not_eigen_vector_vt
    = require_all_not_t<container_type_check_base<is_eigen_vector, value_type_t,
                                                  TypeCheck, Check>...>;

/*! \brief Require type satisfies @ref is_eigen_vector */
/*! and scalar type satisfies `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the scalar type against */
/*! @tparam Check The type to test @ref is_eigen_vector for and whose @ref
 * scalar_type is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_eigen_vector_st
    = require_t<container_type_check_base<is_eigen_vector, scalar_type_t,
                                          TypeCheck, Check...>>;

/*! \brief Require type does not satisfy @ref is_eigen_vector */
/*! or scalar type does not satisfy `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the scalar type against */
/*! @tparam Check The type to test @ref is_eigen_vector for and whose @ref
 * scalar_type is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_not_eigen_vector_st
    = require_not_t<container_type_check_base<is_eigen_vector, scalar_type_t,
                                              TypeCheck, Check...>>;

/*! \brief Require any of the types satisfy @ref is_eigen_vector */
/*! and any scalar type satisfies `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the scalar type against */
/*! @tparam Check The type to test @ref is_eigen_vector for and whose @ref
 * scalar_type is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_any_eigen_vector_st
    = require_any_t<container_type_check_base<is_eigen_vector, scalar_type_t,
                                              TypeCheck, Check>...>;

/*! \brief Require at least one of the types does not satisfy @ref
 * is_eigen_vector */
/*! and any scalar type does not satisfy `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the scalar type against */
/*! @tparam Check The type to test @ref is_eigen_vector for and whose @ref
 * scalar_type is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_any_not_eigen_vector_st
    = require_any_not_t<container_type_check_base<
        is_eigen_vector, scalar_type_t, TypeCheck, Check>...>;

/*! \brief Require all of the types does not satisfy @ref is_eigen_vector */
/*! and all scalar types satisfy `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the scalar type against */
/*! @tparam Check The type to test @ref is_eigen_vector for and whose @ref
 * scalar_type is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_all_eigen_vector_st
    = require_all_t<container_type_check_base<is_eigen_vector, scalar_type_t,
                                              TypeCheck, Check>...>;

/*! \brief Require none of the types satisfy @ref is_eigen_vector */
/*! and none of the scalar types satisfy `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the scalar type against */
/*! @tparam Check The type to test @ref is_eigen_vector for and whose @ref
 * scalar_type is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_all_not_eigen_vector_st
    = require_all_not_t<container_type_check_base<
        is_eigen_vector, scalar_type_t, TypeCheck, Check>...>;
/*! @} */

/**
 * Require `Row` is a row vector and `Col` is a column vector.
 * @ingroup require_eigen_types
 */
template <typename Row, typename Col>
using require_eigen_row_and_col_t = require_t<
    math::conjunction<is_eigen_row_vector<Row>, is_eigen_col_vector<Col>>>;

/**
 * Require `Row` is not a row vector and `Col` is not a column vector.
 * @ingroup require_eigen_types
 */
template <typename Row, typename Col>
using require_not_eigen_row_and_col_t = require_not_t<
    math::conjunction<is_eigen_row_vector<Row>, is_eigen_col_vector<Col>>>;

/**
 * Require `Row` is a row vector and `Col` is a column vector.
 * @ingroup require_eigen_types
 */
template <typename Row, typename Col>
using require_row_and_col_vector_t
    = require_t<math::conjunction<is_row_vector<Row>, is_col_vector<Col>>>;

/**
 * Require `Row` is not a row vector and `Col` is not a column vector.
 * @ingroup require_eigen_types
 */
template <typename Row, typename Col>
using require_not_row_and_col_vector_t
    = require_not_t<math::conjunction<is_row_vector<Row>, is_col_vector<Col>>>;

/** \ingroup type_trait
 * If the input type T is either an eigen matrix with 1 column or 1 row at
 * compile time or a standard vector, this has a static member with a value
 * of true. Else this has a static member with a value of false.
 *
 */
template <typename T>
struct is_vector
    : bool_constant<is_eigen_vector<T>::value || is_std_vector<T>::value
                    || (is_var<T>::value
                        && is_eigen_vector<value_type_t<T>>::value)> {};

/*! \ingroup require_std */
/*! \defgroup vector_types vector  */
/*! \addtogroup vector_types */
/*! @{ */

/*! \brief Require type satisfies @ref is_vector */
/*! @tparam T the type to check */
template <typename T>
using require_vector_t = require_t<is_vector<std::decay_t<T>>>;

/*! \brief Require type does not satisfy @ref is_vector */
/*! @tparam T the type to check */
template <typename T>
using require_not_vector_t = require_not_t<is_vector<std::decay_t<T>>>;

/*! \brief Require all of the types satisfy @ref is_vector */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_all_vector_t = require_all_t<is_vector<std::decay_t<Types>>...>;

/*! \brief Require any of the types satisfy @ref is_vector */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_any_vector_t = require_any_t<is_vector<std::decay_t<Types>>...>;

/*! \brief Require none of the types satisfy @ref is_vector */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_all_not_vector_t
    = require_all_not_t<is_vector<std::decay_t<Types>>...>;

/*! \brief Require at least one of the types do not satisfy @ref is_vector */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_any_not_vector_t
    = require_any_not_t<is_vector<std::decay_t<Types>>...>;
/*! @} */

/*! \ingroup require_std */
/*! \defgroup vector_types vector  */
/*! \addtogroup vector_types */
/*! @{ */

/*! \brief Require type satisfies @ref is_vector */
/*! and value type satisfies `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the value type against */
/*! @tparam Check The type to test @ref is_vector for and whose @ref value_type
 * is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_vector_vt = require_t<
    container_type_check_base<is_vector, value_type_t, TypeCheck, Check...>>;

/*! \brief Require type does not satisfy @ref is_vector or */
/*! value type does not satisfy `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the value type against */
/*! @tparam Check The type to test @ref is_vector for and whose @ref value_type
 * is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_not_vector_vt = require_not_t<
    container_type_check_base<is_vector, value_type_t, TypeCheck, Check...>>;

/*! \brief Require any of the types satisfy @ref is_vector */
/*! and any of the value types satisfy `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the value type against */
/*! @tparam Check The type to test @ref is_vector for and whose @ref value_type
 * is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_any_vector_vt = require_any_t<
    container_type_check_base<is_vector, value_type_t, TypeCheck, Check>...>;

/*! \brief Require at least one of the types does not satisfy @ref is_vector */
/*! and none of the value types satisfy `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the value type against */
/*! @tparam Check The type to test @ref is_vector for and whose @ref value_type
 * is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_any_not_vector_vt = require_any_not_t<
    container_type_check_base<is_vector, value_type_t, TypeCheck, Check>...>;

/*! \brief Require all of the types satisfy @ref is_vector */
/*! and all of the value types satisfy `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the value type against */
/*! @tparam Check The type to test @ref is_vector for and whose @ref value_type
 * is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_all_vector_vt = require_all_t<
    container_type_check_base<is_vector, value_type_t, TypeCheck, Check>...>;

/*! \brief Require none of the types satisfy @ref is_vector */
/*! and none of the value types satisfy `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the value type against */
/*! @tparam Check The type to test @ref is_vector for and whose @ref value_type
 * is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_all_not_vector_vt = require_all_not_t<
    container_type_check_base<is_vector, value_type_t, TypeCheck, Check>...>;

/*! \brief Require type satisfies @ref is_vector */
/*! and scalar type satisfies `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the scalar type against */
/*! @tparam Check The type to test @ref is_vector for and whose @ref scalar_type
 * is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_vector_st = require_t<
    container_type_check_base<is_vector, scalar_type_t, TypeCheck, Check...>>;

/*! \brief Require type does not satisfy @ref is_vector */
/*! or scalar type does not satisfy `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the scalar type against */
/*! @tparam Check The type to test @ref is_vector for and whose @ref scalar_type
 * is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_not_vector_st = require_not_t<
    container_type_check_base<is_vector, scalar_type_t, TypeCheck, Check...>>;

/*! \brief Require any of the types satisfy @ref is_vector */
/*! and any scalar type satisfies `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the scalar type against */
/*! @tparam Check The type to test @ref is_vector for and whose @ref scalar_type
 * is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_any_vector_st = require_any_t<
    container_type_check_base<is_vector, scalar_type_t, TypeCheck, Check>...>;

/*! \brief Require at least one of the types does not satisfy @ref is_vector */
/*! and any scalar type does not satisfy `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the scalar type against */
/*! @tparam Check The type to test @ref is_vector for and whose @ref scalar_type
 * is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_any_not_vector_st = require_any_not_t<
    container_type_check_base<is_vector, scalar_type_t, TypeCheck, Check>...>;

/*! \brief Require all of the types does not satisfy @ref is_vector */
/*! and all scalar types satisfy `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the scalar type against */
/*! @tparam Check The type to test @ref is_vector for and whose @ref scalar_type
 * is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_all_vector_st = require_all_t<
    container_type_check_base<is_vector, scalar_type_t, TypeCheck, Check>...>;

/*! \brief Require none of the types satisfy @ref is_vector */
/*! and none of the scalar types satisfy `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the scalar type against */
/*! @tparam Check The type to test @ref is_vector for and whose @ref scalar_type
 * is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_all_not_vector_st = require_all_not_t<
    container_type_check_base<is_vector, scalar_type_t, TypeCheck, Check>...>;
/*! @} */

namespace internal {

/** \ingroup type_trait
 * This underlying implementation is used when the type is not an std vector.
 */
template <typename T>
struct is_std_vector_impl : std::false_type {};

/** \ingroup type_trait
 * This specialization implementation has a static member named value when the
 * template type is an std vector.
 */
template <typename... Args>
struct is_std_vector_impl<std::vector<Args...>> : std::true_type {};

}  // namespace internal

/** \ingroup type_trait
 * Checks if the decayed type of T is a standard vector.
 */
template <typename T>
struct is_std_vector<
    T, std::enable_if_t<internal::is_std_vector_impl<std::decay_t<T>>::value>>
    : std::true_type {};

/** \ingroup type_trait
 * Specialization of scalar_type for vector to recursively return the inner
 * scalar type.
 *
 * @tparam T type of standard vector
 */
template <typename T>
struct scalar_type<T, std::enable_if_t<is_std_vector<T>::value>> {
  using type = scalar_type_t<typename std::decay_t<T>::value_type>;
};

/** \ingroup type_trait
 * Template metaprogram class to compute the type of values stored
 * in a standard vector.
 *
 * @tparam T type of elements in standard vector.
 */
template <typename T>
struct value_type<T, std::enable_if_t<is_std_vector<T>::value>> {
  using type = typename std::decay_t<T>::value_type;
};

/*! \ingroup require_std */
/*! \defgroup std_vector_types std_vector  */
/*! \addtogroup std_vector_types */
/*! @{ */

/*! \brief Require type satisfies @ref is_std_vector */
/*! @tparam T the type to check */
template <typename T>
using require_std_vector_t = require_t<is_std_vector<std::decay_t<T>>>;

/*! \brief Require type does not satisfy @ref is_std_vector */
/*! @tparam T the type to check */
template <typename T>
using require_not_std_vector_t = require_not_t<is_std_vector<std::decay_t<T>>>;

/*! \brief Require all of the types satisfy @ref is_std_vector */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_all_std_vector_t
    = require_all_t<is_std_vector<std::decay_t<Types>>...>;

/*! \brief Require any of the types satisfy @ref is_std_vector */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_any_std_vector_t
    = require_any_t<is_std_vector<std::decay_t<Types>>...>;

/*! \brief Require none of the types satisfy @ref is_std_vector */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_all_not_std_vector_t
    = require_all_not_t<is_std_vector<std::decay_t<Types>>...>;

/*! \brief Require at least one of the types do not satisfy @ref is_std_vector
 */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_any_not_std_vector_t
    = require_any_not_t<is_std_vector<std::decay_t<Types>>...>;
/*! @} */

/*! \ingroup require_std */
/*! \defgroup std_vector_types std_vector  */
/*! \addtogroup std_vector_types */
/*! @{ */

/*! \brief Require type satisfies @ref is_std_vector */
/*! and value type satisfies `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the value type against */
/*! @tparam Check The type to test @ref is_std_vector for and whose @ref
 * value_type is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_std_vector_vt
    = require_t<container_type_check_base<is_std_vector, value_type_t,
                                          TypeCheck, Check...>>;

/*! \brief Require type does not satisfy @ref is_std_vector or */
/*! value type does not satisfy `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the value type against */
/*! @tparam Check The type to test @ref is_std_vector for and whose @ref
 * value_type is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_not_std_vector_vt
    = require_not_t<container_type_check_base<is_std_vector, value_type_t,
                                              TypeCheck, Check...>>;

/*! \brief Require any of the types satisfy @ref is_std_vector */
/*! and any of the value types satisfy `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the value type against */
/*! @tparam Check The type to test @ref is_std_vector for and whose @ref
 * value_type is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_any_std_vector_vt
    = require_any_t<container_type_check_base<is_std_vector, value_type_t,
                                              TypeCheck, Check>...>;

/*! \brief Require at least one of the types does not satisfy @ref is_std_vector
 */
/*! and none of the value types satisfy `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the value type against */
/*! @tparam Check The type to test @ref is_std_vector for and whose @ref
 * value_type is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_any_not_std_vector_vt
    = require_any_not_t<container_type_check_base<is_std_vector, value_type_t,
                                                  TypeCheck, Check>...>;

/*! \brief Require all of the types satisfy @ref is_std_vector */
/*! and all of the value types satisfy `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the value type against */
/*! @tparam Check The type to test @ref is_std_vector for and whose @ref
 * value_type is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_all_std_vector_vt
    = require_all_t<container_type_check_base<is_std_vector, value_type_t,
                                              TypeCheck, Check>...>;

/*! \brief Require none of the types satisfy @ref is_std_vector */
/*! and none of the value types satisfy `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the value type against */
/*! @tparam Check The type to test @ref is_std_vector for and whose @ref
 * value_type is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_all_not_std_vector_vt
    = require_all_not_t<container_type_check_base<is_std_vector, value_type_t,
                                                  TypeCheck, Check>...>;

/*! \brief Require type satisfies @ref is_std_vector */
/*! and scalar type satisfies `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the scalar type against */
/*! @tparam Check The type to test @ref is_std_vector for and whose @ref
 * scalar_type is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_std_vector_st
    = require_t<container_type_check_base<is_std_vector, scalar_type_t,
                                          TypeCheck, Check...>>;

/*! \brief Require type does not satisfy @ref is_std_vector */
/*! or scalar type does not satisfy `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the scalar type against */
/*! @tparam Check The type to test @ref is_std_vector for and whose @ref
 * scalar_type is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_not_std_vector_st
    = require_not_t<container_type_check_base<is_std_vector, scalar_type_t,
                                              TypeCheck, Check...>>;

/*! \brief Require any of the types satisfy @ref is_std_vector */
/*! and any scalar type satisfies `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the scalar type against */
/*! @tparam Check The type to test @ref is_std_vector for and whose @ref
 * scalar_type is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_any_std_vector_st
    = require_any_t<container_type_check_base<is_std_vector, scalar_type_t,
                                              TypeCheck, Check>...>;

/*! \brief Require at least one of the types does not satisfy @ref is_std_vector
 */
/*! and any scalar type does not satisfy `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the scalar type against */
/*! @tparam Check The type to test @ref is_std_vector for and whose @ref
 * scalar_type is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_any_not_std_vector_st
    = require_any_not_t<container_type_check_base<is_std_vector, scalar_type_t,
                                                  TypeCheck, Check>...>;

/*! \brief Require all of the types does not satisfy @ref is_std_vector */
/*! and all scalar types satisfy `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the scalar type against */
/*! @tparam Check The type to test @ref is_std_vector for and whose @ref
 * scalar_type is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_all_std_vector_st
    = require_all_t<container_type_check_base<is_std_vector, scalar_type_t,
                                              TypeCheck, Check>...>;

/*! \brief Require none of the types satisfy @ref is_std_vector */
/*! and none of the scalar types satisfy `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the scalar type against */
/*! @tparam Check The type to test @ref is_std_vector for and whose @ref
 * scalar_type is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_all_not_std_vector_st
    = require_all_not_t<container_type_check_base<is_std_vector, scalar_type_t,
                                                  TypeCheck, Check>...>;
/*! @} */

}  // namespace stan
#endif
