#ifndef STAN_MATH_REV_FUN_LOG_DETERMINANT_SPD_HPP
#define STAN_MATH_REV_FUN_LOG_DETERMINANT_SPD_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/core/typedefs.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/fun/typedefs.hpp>

namespace stan {
namespace math {

/**
 * Returns the log det of a symmetric, positive-definite matrix
 *
 * @tparam T Type is an Eigen Matrix or `var_value` with inner Eigen matrix type
 * @param M a symmetric, positive-definite matrix
 * @return The log determinant of the specified matrix
 */
template <typename T, require_rev_matrix_t<T>* = nullptr>
inline var log_determinant_spd(const T& M) {
  if (M.size() == 0) {
    return var(0.0);
  }
  check_symmetric("log_determinant", "M", M);

  arena_t<T> arena_M = M;
  matrix_d M_d = arena_M.val();
  auto M_ldlt = M_d.ldlt();
  if (M_ldlt.info() != Eigen::Success) {
    constexpr double y = 0;
    throw_domain_error("log_determinant_spd", "matrix argument", y,
                       "failed LDLT factorization");
  }
  // compute the inverse of A (needed for the derivative)
  M_d.setIdentity(M.rows(), M.cols());
  M_ldlt.solveInPlace(M_d);
  auto arena_M_inv_transpose = to_arena(M_d.transpose());

  if (M_ldlt.isNegative() || (M_ldlt.vectorD().array() <= 1e-16).any()) {
    constexpr double y = 0;
    throw_domain_error("log_determinant_spd", "matrix argument", y,
                       "matrix is negative definite");
  }

  var log_det = sum(log(M_ldlt.vectorD()));

  reverse_pass_callback([arena_M, log_det, arena_M_inv_transpose]() mutable {
    arena_M.adj() += log_det.adj() * arena_M_inv_transpose;
  });
  return log_det;
}

}  // namespace math
}  // namespace stan
#endif
