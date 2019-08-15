#ifndef STAN_MATH_PRIM_MAT_META_ENABLE_IF_EIGEN_SIZES_MATH_HPP
#define STAN_MATH_PRIM_MAT_META_ENABLE_IF_EIGEN_SIZES_MATH_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/scal/meta/conjunction.hpp>
#include <stan/math/prim/scal/meta/disjunction.hpp>

#include <type_traits>

namespace stan {

template <typename T1, typename T2>
using enable_if_eigen_rows_match
    = std::enable_if_t<T1::RowsAtCompileTime == T2::RowsAtCompileTime>;

template <typename T1, typename T2>
using enable_if_not_eigen_rows_match
    = std::enable_if_t<!(T1::RowsAtCompileTime == T2::RowsAtCompileTime)>;

template <typename T1, typename T2>
using enable_if_eigen_cols_match
    = std::enable_if_t<T1::ColsAtCompileTime == T2::ColsAtCompileTime>;

template <typename T1, typename T2>
using enable_if_not_eigen_cols_match
    = std::enable_if_t<!(T1::ColsAtCompileTime == T2::ColsAtCompileTime)>;

template <typename T1, typename T2>
using enable_if_eigen_size_match
    = std::enable_if_t<T1::ColsAtCompileTime == T2::ColsAtCompileTime
                       || T1::RowsAtCompileTime == T2::RowsAtCompileTime>;

template <typename T1, typename T2>
using enable_if_not_eigen_size_match
    = std::enable_if_t<!(T1::ColsAtCompileTime == T2::ColsAtCompileTime
                         || T1::RowsAtCompileTime == T2::RowsAtCompileTime)>;

}  // namespace stan
#endif
