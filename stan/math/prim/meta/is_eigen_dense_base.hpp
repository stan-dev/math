#ifndef STAN_MATH_PRIM_META_IS_EIGEN_DENSE_BASE_HPP
#define STAN_MATH_PRIM_META_IS_EIGEN_DENSE_BASE_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/bool_constant.hpp>
#include <stan/math/prim/meta/is_base_pointer_convertible.hpp>
#include <stan/math/prim/meta/require_helpers.hpp>
#include <type_traits>

namespace stan {

/**
 * Checks whether type T is derived from Eigen::DenseBase.
 * If true this will have a static member function named value with a type
 * of true, else value is false.
 * @tparam T Type to check if it is derived from `DenseBase`
 * @tparam Enable used for SFINAE deduction.
 * @ingroup type_trait
 */
template <typename T>
struct is_eigen_dense_base
    : bool_constant<is_base_pointer_convertible<Eigen::DenseBase, T>::value> {};

/*! \ingroup require_eigens_types */
/*! \defgroup eigen_dense_base_types eigen_dense_base  */
/*! \addtogroup eigen_dense_base_types */
/*! @{ */

/*! \brief Require type satisfies @ref is_eigen_dense_base */
/*! @tparam T the type to check */
template <typename T>
using require_eigen_dense_base_t
    = require_t<is_eigen_dense_base<std::decay_t<T>>>;
/*! @} */

}  // namespace stan

#endif
