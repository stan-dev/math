#ifndef STAN_MATH_PRIM_META_IS_EIGEN_DENSE_DYNAMIC_HPP
#define STAN_MATH_PRIM_META_IS_EIGEN_DENSE_DYNAMIC_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/bool_constant.hpp>
#include <stan/math/prim/meta/is_eigen_matrix.hpp>
#include <stan/math/prim/meta/is_eigen_dense_base.hpp>
#include <stan/math/prim/meta/is_base_pointer_convertible.hpp>
#include <stan/math/prim/meta/require_helpers.hpp>
#include <type_traits>

namespace stan {

/**
 * Checks whether type T is derived from Eigen::DenseBase and has dynamic rows
 * and columns. If true this will have a static member function named value with
 * a type of true, else value is false.
 * @tparam T Type to check
 * @ingroup type_trait
 */
template <typename T>
using is_eigen_dense_dynamic = stan::internal::is_eigen_matrix_dynamic_impl<
    std::decay_t<T>, stan::is_eigen_dense_base<std::decay_t<T>>::value>;

/*! \ingroup require_eigens_types */
/*! \defgroup eigen_dense_dynamic_types eigen_dense_dynamic  */
/*! \addtogroup eigen_dense_dynamic_types */
/*! @{ */

/*! \brief Require type satisfies @ref is_eigen_dense_dynamic */
/*! @tparam T the type to check */
template <typename T>
using require_eigen_dense_dynamic_t
    = require_t<is_eigen_dense_dynamic<std::decay_t<T>>>;

/*! \brief Require all of the types satisfy @ref is_eigen_dense_dynamic */
/*! @tparam Types The types that are checked */
template <typename... Types>
using require_all_eigen_dense_dynamic_t
    = require_all_t<is_eigen_dense_dynamic<std::decay_t<Types>>...>;
/*! @} */

/*! \ingroup require_eigens_types */
/*! \defgroup eigen_dense_dynamic_types eigen_dense_dynamic  */
/*! \addtogroup eigen_dense_dynamic_types */
/*! @{ */

/*! \brief Require type satisfies @ref is_eigen_dense_dynamic */
/*! and value type satisfies `TypeCheck` */
/*! @tparam TypeCheck The type trait to check the value type against */
/*! @tparam Check The type to test @ref is_eigen_dense_dynamic for and whose
 * @ref value_type is checked with `TypeCheck` */
template <template <class...> class TypeCheck, class... Check>
using require_eigen_dense_dynamic_vt
    = require_t<container_type_check_base<is_eigen_dense_dynamic, value_type_t,
                                          TypeCheck, Check...>>;
/*! @} */

}  // namespace stan

#endif
