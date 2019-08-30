#ifndef STAN_MATH_FWD_MAT_FUN_VALUE_OF_HPP
#define STAN_MATH_FWD_MAT_FUN_VALUE_OF_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/fwd/meta.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/mat/fun/value_of_rec.hpp>

namespace stan {
namespace math {

/**
 * Convert a matrix of type T to a matrix of doubles.
 *
 * T must implement value_of. See
 * test/math/fwd/mat/fun/value_of.cpp for fvar and var usage.
 *
 * @tparam T Scalar type in matrix
 * @tparam R Rows of matrix
 * @tparam C Columns of matrix
 * @param[in] x Matrix to be converted
 * @return Matrix of values
 **/
template <typename T, enable_if_eigen<T>* = nullptr,
          enable_if_fvar<scalar_type_t<std::decay_t<T>>>* = nullptr>
inline auto value_of(T&& x) {
  return value_of(x.val().eval()).eval();
}

}  // namespace math
}  // namespace stan

#endif
