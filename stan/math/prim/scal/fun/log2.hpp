#ifndef STAN_MATH_PRIM_SCAL_FUN_LOG2_HPP
#define STAN_MATH_PRIM_SCAL_FUN_LOG2_HPP

#include <stan/math/prim/scal/fun/constants.hpp>
#include <boost/math/tools/promotion.hpp>
#include <cmath>

namespace stan {
  namespace math {

    /**
     * Returns the base 2 logarithm of the argument (C99).
     *
     * The function is defined by:
     *
     * <code>log2(a) = log(a) / std::log(2.0)</code>.
     *
     * @tparam T type of scalar
     * @param a Value.
     * @return Base 2 logarithm of the value.
     */
    template <typename T>
    inline typename boost::math::tools::promote_args<T>::type
    log2(const T& a) {
      using std::log;
      return log(a) / LOG_2;
    }

    /**
     * Return the base two logarithm of the specified
     * argument.  This version is required to disambiguate
     * <code>log2(int)</code>.
     *
     * @param[in] x Argument.
     * @return Base two logarithm of the argument.
     */
    inline double log2(int x) {
      return log2(static_cast<double>(x));
    }

    /**
     * Return natural logarithm of two.
     *
     * @return Natural logarithm of two.
     */
    inline double log2() {
      return LOG_2;
    }

  }
}

#endif
