#ifndef STAN_MATH_PRIM_MAT_META_ENABLE_IF_EIGEN_VECTOR_HPP
#define STAN_MATH_PRIM_MAT_META_ENABLE_IF_EIGEN_VECTOR_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/scal/meta/conjunction.hpp>
#include <stan/math/prim/scal/meta/disjunction.hpp>

#include <type_traits>

namespace stan {

template <typename T>
using enable_if_eigen_vector = std::enable_if_t<
    std::is_base_of<Eigen::EigenBase<std::decay_t<T>>, std::decay_t<T>>::value
    && T::RowsAtCompileTime + T::ColsAtCompileTime == 0>;

template <typename T>
using enable_if_eigen_row_vector = std::enable_if_t<
    std::is_base_of<Eigen::EigenBase<std::decay_t<T>>, std::decay_t<T>>::value
    && T::ColsAtCompileTime == 1>;

template <typename T>
using enable_if_not_eigen_row_vector = std::enable_if_t<
    !std::is_base_of<Eigen::EigenBase<std::decay_t<T>>, std::decay_t<T>>::value
    && !T::ColsAtCompileTime == 1>;

template <typename T>
using enable_if_eigen_col_vector = std::enable_if_t<
    std::is_base_of<Eigen::EigenBase<std::decay_t<T>>, std::decay_t<T>>::value
    && T::RowsAtCompileTime == 1>;

template <typename T>
using enable_if_not_eigen_col_vector = std::enable_if_t<
    !std::is_base_of<Eigen::EigenBase<std::decay_t<T>>, std::decay_t<T>>::value
    && !T::RowsAtCompileTime == 1>;

template <typename T1, typename T2>
using enable_if_dot_product = std::enable_if_t<(
    T1::RowsAtCompileTime == 1 && T1::ColsAtCompileTime == -1
    && T2::RowsAtCompileTime == -1 && T2::ColsAtCompileTime == 1)>;

template <typename T1, typename T2>
using enable_if_not_dot_product = std::enable_if_t<!(
    T1::RowsAtCompileTime == 1 && T1::ColsAtCompileTime == -1
    && T2::RowsAtCompileTime == -1 && T2::ColsAtCompileTime == 1)>;

}  // namespace stan
#endif
