#ifndef STAN_MATH_PRIM_ERR_CHECK_NONEMPTY_HPP
#define STAN_MATH_PRIM_ERR_CHECK_NONEMPTY_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/invalid_argument.hpp>
#include <sstream>

namespace stan {
namespace math {

/**
 * Check if the specified matrix is not empty. This check does not allow 0x0
 * matrices.
 * @tparam EigMat A type derived from `EigenBase` 
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Matrix to test
 * @throw <code>std::invalid_argument</code> if the matrix is empty
 */
template <typename EigMat, require_eigen_t<EigMat>* = nullptr>
inline void check_nonempty(const char* function, const char* name, const EigMat& y) {
  if (y.rows() > 0 && y.cols() > 0) {
    return;
  }
  invalid_argument(function, "Expecting a non empty matrix", name, "", ".");
}

}  // namespace math
}  // namespace stan
#endif
