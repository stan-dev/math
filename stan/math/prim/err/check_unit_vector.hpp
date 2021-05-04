#ifndef STAN_MATH_PRIM_ERR_CHECK_UNIT_VECTOR_HPP
#define STAN_MATH_PRIM_ERR_CHECK_UNIT_VECTOR_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/check_nonzero_size.hpp>
#include <stan/math/prim/err/throw_domain_error.hpp>
#include <stan/math/prim/err/constraint_tolerance.hpp>
#include <stan/math/prim/fun/abs.hpp>
#include <sstream>
#include <string>
#include <cmath>

namespace stan {
namespace math {

/**
 * Check if the specified vector is unit vector.
 * A valid unit vector is one where the square of the elements
 * summed is equal to 1. This function tests that the sum is within the
 * tolerance specified by <code>CONSTRAINT_TOLERANCE</code>. This
 * function only accepts Eigen vectors, statically typed vectors,
 * not general matrices with 1 column.
 * @tparam EigVec A type derived from `EigenBase` with either dynamic rows or
 * columns
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param theta Vector to test
 * @throw <code>std::invalid_argument</code> if <code>theta</code>
 *   is a 0-vector
 * @throw <code>std::domain_error</code> if the vector is not a unit
 *   vector or if any element is <code>NaN</code>
 */
template <typename EigVec, require_eigen_vector_t<EigVec>* = nullptr>
void check_unit_vector(const char* function, const char* name,
                       const EigVec& theta) {
  check_nonzero_size(function, name, theta);
  using std::fabs;
  value_type_t<EigVec> ssq = theta.squaredNorm();
  if (!(fabs(1.0 - ssq) <= CONSTRAINT_TOLERANCE)) {
    [&]() STAN_COLD_PATH {
      std::stringstream msg;
      msg << "is not a valid unit vector."
          << " The sum of the squares of the elements should be 1, but is ";
      std::string msg_str(msg.str());
      throw_domain_error(function, name, ssq, msg_str.c_str());
    }();
  }
}

}  // namespace math
}  // namespace stan
#endif
