#ifndef STAN_MATH_PRIM_FUN_DIAG_MATRIX_HPP
#define STAN_MATH_PRIM_FUN_DIAG_MATRIX_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Return a square diagonal matrix with the specified vector of
 * coefficients as the diagonal values.
 *
 * @tparam EigVec type of the vector (must be derived from \c Eigen::MatrixBase
 * and have one compile time dimension equal to 1)
 * @param[in] v Specified vector.
 * @return Diagonal matrix with vector as diagonal values.
 */
template <typename EigVec, require_eigen_vector_t<EigVec>* = nullptr>
inline auto diag_matrix(EigVec&& v) {
  if constexpr (std::is_lvalue_reference_v<EigVec&&>) {
    return v.asDiagonal();
  } else {
    using diagonal_t = std::decay_t<decltype(v.asDiagonal())>;
    return typename diagonal_t::PlainObject(std::forward<EigVec>(v));
  }
}

}  // namespace math
}  // namespace stan

#endif
