#ifndef STAN_MATH_PRIM_ERR_CHECK_LDLT_FACTOR_HPP
#define STAN_MATH_PRIM_ERR_CHECK_LDLT_FACTOR_HPP

#include <stan/math/prim/err/throw_domain_error.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/LDLT_factor.hpp>
#include <sstream>
#include <string>

namespace stan {
namespace math {

/**
 * Raise domain error if the specified LDLT factor is invalid.  An
 * <code>LDLT_factor</code> is invalid if it was constructed from
 * a matrix that is not positive definite.
 *
 * @tparam T Type matrix of LDLT
 * @param[in] function name of function for error messages
 * @param[in] name variable name for error messages
 * @param[in] A the LDLT factor to check for validity
 * @throws <code>std::domain_error</code> if the LDLT factor is invalid
 */
template <typename T>
inline void check_ldlt_factor(const char* function, const char* name,
                              LDLT_factor<T>& A) {
  if (!(A.ldlt().info() == Eigen::Success && A.ldlt().isPositive()
        && (A.ldlt().vectorD().array() > 0).all())) {
    std::ostringstream msg;
    msg << "is not positive definite.  last conditional variance is ";
    std::string msg_str(msg.str());
    auto too_small = A.ldlt().vectorD().tail(1)(0);
    throw_domain_error(function, name, too_small, msg_str.c_str(), ".");
  }
}

}  // namespace math
}  // namespace stan
#endif
