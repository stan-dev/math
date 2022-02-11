#ifndef STAN_MATH_PRIM_ERR_CHECK_ORDERED_HPP
#define STAN_MATH_PRIM_ERR_CHECK_ORDERED_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/make_iter_name.hpp>
#include <stan/math/prim/err/throw_domain_error.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <sstream>
#include <string>
#include <vector>

namespace stan {
namespace math {

/**
 * Throw an exception if the specified vector is not sorted into strictly
 * increasing order.
 * @tparam T_y A type inheriting from EigenBase with either 1 compile time row
 * or column
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Vector to test
 * @throw `std::domain_error` if the vector elements are not ordered, if there
 * are duplicated values, or if any element is `NaN`
 */
template <typename T_y, require_vector_t<T_y>* = nullptr,
          require_not_std_vector_t<T_y>* = nullptr>
void check_ordered(const char* function, const char* name, const T_y& y) {
  const auto& y_ref = to_ref(value_of_rec(y));
  for (Eigen::Index n = 1; n < y_ref.size(); n++) {
    if (!(y_ref[n] > y_ref[n - 1])) {
      [&]() STAN_COLD_PATH {
        std::ostringstream msg1;
        msg1 << "is not a valid ordered vector."
             << " The element at " << stan::error_index::value + n << " is ";
        std::string msg1_str(msg1.str());
        std::ostringstream msg2;
        msg2 << ", but should be greater than the previous element, "
             << y_ref[n - 1];
        std::string msg2_str(msg2.str());
        throw_domain_error(function, name, y_ref[n], msg1_str.c_str(),
                           msg2_str.c_str());
      }();
    }
  }
}

/**
 * Throw an exception if the specified vector is not sorted into strictly
 * increasing order.
 * @tparam T_y A standard vector with inner scalar type
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y `std::vector` to test
 * @throw `std::domain_error` if the vector elements are not ordered, if there
 * are duplicated values, or if any element is `NaN`
 */
template <typename T_y, require_std_vector_vt<is_stan_scalar, T_y>* = nullptr>
void check_ordered(const char* function, const char* name, const T_y& y) {
  for (size_t n = 1; n < y.size(); n++) {
    if (!(y[n] > y[n - 1])) {
      [&]() STAN_COLD_PATH {
        std::ostringstream msg1;
        msg1 << "is not a valid ordered vector."
             << " The element at " << stan::error_index::value + n << " is ";
        std::string msg1_str(msg1.str());
        std::ostringstream msg2;
        msg2 << ", but should be greater than the previous element, "
             << y[n - 1];
        std::string msg2_str(msg2.str());
        throw_domain_error(function, name, y[n], msg1_str.c_str(),
                           msg2_str.c_str());
      }();
    }
  }
}

/**
 * Throw an exception if each vector in a standard vector is not sorted into
 * strictly increasing order.
 * @tparam T_y A standard vector with an inner vector type
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y `std::vector` to test
 * @throw `std::domain_error` if the vector elements are not ordered, if there
 * are duplicated values, or if any element is `NaN`
 */
template <typename T_y, require_std_vector_t<T_y>* = nullptr,
          require_not_vt_stan_scalar<T_y>* = nullptr>
void check_ordered(const char* function, const char* name, const T_y& y) {
  for (size_t i = 0; i < y.size(); ++i) {
    check_ordered(function, internal::make_iter_name(name, i).c_str(), y[i]);
  }
}

}  // namespace math
}  // namespace stan
#endif
