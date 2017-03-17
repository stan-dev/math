#ifndef STAN_MATH_PRIM_SCAL_ERR_CHECK_2F1_CONVERGES_HPP
#define STAN_MATH_PRIM_SCAL_ERR_CHECK_2F1_CONVERGES_HPP

#include <stdexcept>
#include <stringstream>
#include <cmath>

namespace stan {
  namespace math {

    /**
     * Check if the hypergeometric function (2F1) called with 
     * supplied arguments will converge. 
     *
     * @tparam T_a1 Type of a1
     * @tparam T_a2 Type of a2
     * @tparam T_b1 Type of b1
     * @tparam T_z Type of z
     *
     * @param function Name of function ultimately relying on 2F1 (for error 
     &   messages)
     * @param name Variable name (for error messages)
     * @param a1 Variable to check
     * @param a2 Variable to check
     * @param b1 Variable to check
     * @param z Variable to check
     *
     * @throw <code>domain_error</code> if 2F1(a1, a2, b1, z)
     *   does not meet convergence conditions, or if any coefficient is NaN. 
     */
    template <typename T_a1, typename T_a2, typename T_b1, typename T_z>
    inline void check_2F1_converges(const char* function, const char* name,
      const T_a1& a1, const T_a2& a2, const T_b1& b1, const T_z& z
    ) {
      using std::fabs;
      if (fabs(z) < 1)
        return;
      if (fabs(z) == 1) {
        if (b1 > a1 + a2)
          return;
      }     
      std::stringstream msg;
      msg << "hypergeometric function 2F1 does not meet convergence " 
          << "conditions with given arguments. " 
          << "a1: " << a1 << ", a2: " << a2  
          << "b1: " << b1 << ", z: " << z;
      throw std::domain_error(msg.str());
    }

  }
}
#endif
