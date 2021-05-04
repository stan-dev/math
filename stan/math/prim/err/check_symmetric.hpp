#ifndef STAN_MATH_PRIM_ERR_CHECK_SYMMETRIC_HPP
#define STAN_MATH_PRIM_ERR_CHECK_SYMMETRIC_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/throw_domain_error.hpp>
#include <stan/math/prim/err/check_square.hpp>
#include <stan/math/prim/err/constraint_tolerance.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/abs.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <sstream>
#include <string>
#include <cmath>

namespace stan {
namespace math {

/**
 * Check if the specified matrix is symmetric.
 * The error message is either 0 or 1 indexed, specified by
 * <code>stan::error_index::value</code>.
 * @tparam EigMat Type of matrix
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Matrix to test
 * @throw <code>std::invalid_argument</code> if the matrix is not square.
 * @throw <code>std::domain_error</code> if any element not on the
 *   main diagonal is <code>NaN</code>
 */
template <typename EigMat, require_matrix_t<EigMat>* = nullptr>
inline void check_symmetric(const char* function, const char* name,
                            const EigMat& y) {
  check_square(function, name, y);
  using std::fabs;

  Eigen::Index k = y.rows();
  if (k <= 1) {
    return;
  }
  const auto& y_ref = to_ref(y);
  for (Eigen::Index m = 0; m < k; ++m) {
    for (Eigen::Index n = m + 1; n < k; ++n) {
      if (!(fabs(value_of(y_ref(m, n)) - value_of(y_ref(n, m)))
            <= CONSTRAINT_TOLERANCE)) {
        [&]() STAN_COLD_PATH {
          std::ostringstream msg1;
          msg1 << "is not symmetric. " << name << "["
               << stan::error_index::value + m << ","
               << stan::error_index::value + n << "] = ";
          std::string msg1_str(msg1.str());
          std::ostringstream msg2;
          msg2 << ", but " << name << "[" << stan::error_index::value + n << ","
               << stan::error_index::value + m << "] = " << y_ref(n, m);
          std::string msg2_str(msg2.str());
          throw_domain_error(function, name, y_ref(m, n), msg1_str.c_str(),
                             msg2_str.c_str());
        }();
      }
    }
  }
}

}  // namespace math
}  // namespace stan
#endif
