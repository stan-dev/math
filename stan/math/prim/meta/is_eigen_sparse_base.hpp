#ifndef STAN_MATH_PRIM_META_IS_EIGEN_SPARSE_BASE_HPP
#define STAN_MATH_PRIM_META_IS_EIGEN_SPARSE_BASE_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/bool_constant.hpp>
#include <stan/math/prim/meta/is_base_pointer_convertible.hpp>
#include <stan/math/prim/meta/require_helpers.hpp>
#include <type_traits>

namespace stan {

/**
 * Checks whether type T is derived from Eigen::SparseMatrixBase.
 * If true this will have a static member function named value with a type
 * of true, else value is false.
 * @tparam T Type to check if it is derived from `SparseMatrixBase`
 * @tparam Enable used for SFINAE deduction.
 * @ingroup type_trait
 */
template <typename T>
struct is_eigen_sparse_base
    : bool_constant<
          is_base_pointer_convertible<Eigen::SparseMatrixBase, T>::value> {};

/*! \ingroup require_eigens_types */
/*! \defgroup eigen_sparse_base_types eigen_sparse_base  */
/*! \addtogroup eigen_sparse_base_types */
/*! @{ */

/*! \brief Require type satisfies @ref is_eigen_sparse_base */
/*! @tparam T the type to check */
template <typename T>
using require_eigen_sparse_base_t
    = require_t<is_eigen_sparse_base<std::decay_t<T>>>;
/*! @} */

}  // namespace stan

#endif
