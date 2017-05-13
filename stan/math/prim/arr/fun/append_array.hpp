#ifndef STAN_MATH_PRIM_ARR_FUN_APPEND_ARRAY_HPP
#define STAN_MATH_PRIM_ARR_FUN_APPEND_ARRAY_HPP

#include <boost/math/tools/promotion.hpp>
#include <vector>

namespace stan {
  namespace math {

    /**
     * Return the concatenation of two specified integer vectors in the
     *   order of the arguments
     *
     * @param x First vector.
     * @param y Second vector.
     * @return A vector of x and y concatenated together (in that order).
     */
    inline std::vector<int>
    append_array(const std::vector<int>& x, const std::vector<int>& y) {
      std::vector<int> z;
      z.reserve(x.size() + y.size());
      z.insert(z.end(), x.begin(), x.end());
      z.insert(z.end(), y.begin(), y.end());
      return z;
    }

    /**
     * Return the concatenation of two specified vectors in the order of
     *   the arguments
     *
     * @tparam T1 Scalar type of first vector.
     * @tparam T2 Scalar Type of second vector.
     * @param x First vector.
     * @param y Second vector.
     * @return A vector of x and y concatenated together (in that order).
     */
    template <typename T1, typename T2>
    inline std::vector<typename boost::math::tools::promote_args<T1, T2>::type>
    append_array(const std::vector<T1>& x, const std::vector<T2>& y) {
      std::vector<typename boost::math::tools::promote_args<T1, T2>::type> z;
      z.reserve(x.size() + y.size());
      z.insert(z.end(), x.begin(), x.end());
      z.insert(z.end(), y.begin(), y.end());
      return z;
    }
  }
}
#endif
