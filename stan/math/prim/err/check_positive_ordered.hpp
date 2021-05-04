#ifndef STAN_MATH_PRIM_ERR_CHECK_POSITIVE_ORDERED_HPP
#define STAN_MATH_PRIM_ERR_CHECK_POSITIVE_ORDERED_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/throw_domain_error.hpp>
#include <stan/math/prim/err/check_ordered.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <sstream>
#include <string>

namespace stan {
namespace math {

/**
 * Check if the specified vector contains non-negative values and is sorted into
 * strictly increasing order.
 * @tparam EigVec A type derived from `EigenBase` with 1 compile time row or
 * column
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Vector to test
 * @throw <code>std::domain_error</code> if the vector contains non-positive
 *   values, if the values are not ordered, if there are duplicated
 *   values, or if any element is <code>NaN</code>.
 */
template <typename EigVec, require_eigen_vector_t<EigVec>* = nullptr>
void check_positive_ordered(const char* function, const char* name,
                            const EigVec& y) {
  if (y.size() == 0) {
    return;
  }

  const auto& y_ref = to_ref(y);
  if (y_ref[0] < 0) {
    [&]() STAN_COLD_PATH {
      std::ostringstream msg;
      msg << "is not a valid positive_ordered vector."
          << " The element at " << stan::error_index::value << " is ";
      std::string msg_str(msg.str());
      throw_domain_error(function, name, y_ref[0], msg_str.c_str(),
                         ", but should be postive.");
    }();
  }
  check_ordered(function, name, y_ref);
}
}  // namespace math
}  // namespace stan
#endif
