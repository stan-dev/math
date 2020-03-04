#ifndef STAN_MATH_PRIM_META_IS_EIGEN_SPARSE_MATRIX_HPP
#define STAN_MATH_PRIM_META_IS_EIGEN_SPARSE_MATRIX_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/bool_constant.hpp>
#include <stan/math/prim/meta/disjunction.hpp>
#include <type_traits>

namespace stan {
/** \addtogroup type_trait
 *  @{
 */

/*
 * Checks whether type T is derived from Eigen::SparseMatrixBase .If true this
 * will have a static member function named value with a type of true, else
 * value is false.
 * @tparam T Type to check if it is derived from `EigenBase`
 * @tparam Enable used for SFINAE deduction.
 */
template <typename T, typename Enable = void>
struct is_eigen_sparse_matrix : std::false_type {};

template <typename T>
struct is_eigen_sparse_matrix<
    T, std::enable_if_t<std::is_base_of<
           Eigen::SparseMatrixBase<typename std::decay_t<T>::PlainObject>,
           typename std::decay_t<T>::PlainObject>::value>> : std::true_type {};

template <typename T>
struct is_eigen_sparse_matrix<
    T, std::enable_if_t<std::is_base_of<
           Eigen::SparseMatrixBase<typename std::decay_t<T>::MatrixType>,
           typename std::decay_t<T>::MatrixType>::value>> : std::true_type {};
/** @}*/

}  // namespace stan
#endif
