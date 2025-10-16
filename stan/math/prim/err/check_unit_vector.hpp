#ifndef STAN_MATH_PRIM_ERR_CHECK_UNIT_VECTOR_HPP
#define STAN_MATH_PRIM_ERR_CHECK_UNIT_VECTOR_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/check_nonzero_size.hpp>
#include <stan/math/prim/err/constraint_tolerance.hpp>
#include <stan/math/prim/err/make_iter_name.hpp>
#include <stan/math/prim/err/throw_domain_error.hpp>
#include <stan/math/prim/fun/abs.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <sstream>
#include <string>
#include <cmath>

namespace stan {
namespace math {

/**
 * Throw an exception if the specified vector does not have unit Euclidiean
 * length. A valid unit vector is one where the square of the elements summed is
 * equal to 1. This function tests that the sum is within the tolerance
 * specified by `CONSTRAINT_TOLERANCE`. This function only accepts Eigen
 * vectors, statically typed vectors, not general matrices with 1 column.
 * @tparam Vec A type derived from `Eigen::EigenBase` with either dynamic rows
 * or columns but not both
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param theta Vector to test
 * @throw `std::invalid_argument` if `theta` is a 0-vector
 * @throw `std::domain_error` if the vector is not a unit vector or if any
 * element is `NaN`
 */
template <typename Vec, require_vector_t<Vec>* = nullptr,
          require_not_std_vector_t<Vec>* = nullptr>
void check_unit_vector(const char* function, const char* name,
                       const Vec& theta) {
  check_nonzero_size(function, name, theta);
  using std::fabs;
  scalar_type_t<Vec> ssq = value_of_rec(theta).squaredNorm();
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

/**
 * Throw an exception if the each element in a standard vector does not have
 * unit Euclidiean length.
 * @tparam StdVec A standard vector with inner type inheriting from
 * `Eigen::EigenBase`
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param theta Vector to test
 * @throw `std::invalid_argument` if `theta`
 *   is a 0-vector
 * @throw `std::domain_error` if the vector is not a unit vector or if any
 * element is `NaN`
 */
template <typename StdVec, require_std_vector_t<StdVec>* = nullptr>
void check_unit_vector(const char* function, const char* name,
                       const StdVec& theta) {
  for (size_t i = 0; i < theta.size(); ++i) {
    check_unit_vector(function, internal::make_iter_name(name, i).c_str(),
                      theta[i]);
  }
}

}  // namespace math
}  // namespace stan
#endif
