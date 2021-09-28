#ifndef STAN_MATH_PRIM_ERR_CHECK_POSITIVE_ORDERED_HPP
#define STAN_MATH_PRIM_ERR_CHECK_POSITIVE_ORDERED_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/check_ordered.hpp>
#include <stan/math/prim/err/make_iter_name.hpp>
#include <stan/math/prim/err/throw_domain_error.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <sstream>
#include <string>

namespace stan {
namespace math {

/**
 * Throw an exception if the specified the vector contains negative values or
 * is not sorted into strictly increasing order.
 * @tparam Vec A type derived from `Eigen::EigenBase` with 1 compile time row or
 * column
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Vector to test
 * @throw `std::domain_error` if the vector contains non-positive values, if the
 * values are not ordered, if there are duplicated values, or if any element is
 * `NaN`
 */
template <typename Vec, require_vector_t<Vec>* = nullptr,
          require_not_std_vector_t<Vec>* = nullptr>
void check_positive_ordered(const char* function, const char* name,
                            const Vec& y) {
  if (y.size() == 0) {
    return;
  }
  auto&& y_ref = to_ref(value_of_rec(y));
  if (y_ref.coeff(0) < 0) {
    [&]() STAN_COLD_PATH {
      std::ostringstream msg;
      msg << "is not a valid positive_ordered vector."
          << " The element at " << stan::error_index::value << " is ";
      std::string msg_str(msg.str());
      throw_domain_error(function, name, value_of_rec(to_ref(y).coeff(0)),
                         msg_str.c_str(), ", but should be postive.");
    }();
  }
  check_ordered(function, name, y_ref);
}

/**
 * Throw an exception if any of the vectors in a standard vector contains
 * negative values or is not sorted into strictly increasing order.
 * @tparam StdVec A standard vector type with an `value_type` inheriting from
 * `Eigen::EigenBase` with 1 compile time row or column
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Vector to test
 * @throw `std::domain_error` if the vector contains non-positive values, if the
 * values are not ordered, if there are duplicated values, or if any element is
 * `NaN`
 */
template <typename StdVec, require_std_vector_t<StdVec>* = nullptr>
void check_positive_ordered(const char* function, const char* name,
                            const StdVec& y) {
  for (size_t i = 0; i < y.size(); ++i) {
    check_positive_ordered(function, internal::make_iter_name(name, i).c_str(),
                           y[i]);
  }
}
}  // namespace math
}  // namespace stan
#endif
