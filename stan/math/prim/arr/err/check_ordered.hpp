#ifndef STAN_MATH_PRIM_ARR_ERR_CHECK_ORDERED_HPP
#define STAN_MATH_PRIM_ARR_ERR_CHECK_ORDERED_HPP

#include <stan/math/prim/arr/meta/index_type.hpp>
#include <stan/math/prim/scal/err/domain_error.hpp>
#include <stan/math/prim/scal/meta/error_index.hpp>
#include <sstream>
#include <string>
#include <vector>

namespace stan {
  namespace math {

    /**
     * Check if the specified vector is sorted into
     * strictly increasing order.
     *
     * @tparam T_y Type of scalar
     *
     * @param function Function name (for error messages)
     * @param name Variable name (for error messages)
     * @param y <code>std::vector</code> to test
     *
     * @throw <code>std::domain_error</code> if the vector elements are
     *   not ordered, if there are duplicated
     *   values, or if any element is <code>NaN</code>.
     */
    template <typename T_y>
    void check_ordered(const char* function,
                       const char* name,
                       const std::vector<T_y>& y) {
      for (size_t n = 1; n < y.size(); n++) {
        if (!(y[n] > y[n - 1])) {
          std::ostringstream msg1;
          msg1 << "is not a valid ordered vector."
               << " The element at " << stan::error_index::value + n
               << " is ";
          std::string msg1_str(msg1.str());
          std::ostringstream msg2;
          msg2 << ", but should be greater than the previous element, "
               << y[n - 1];
          std::string msg2_str(msg2.str());
          domain_error(function, name, y[n],
                       msg1_str.c_str(), msg2_str.c_str());
        }
      }
    }
  }
}
#endif
