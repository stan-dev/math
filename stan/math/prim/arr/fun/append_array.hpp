#ifndef STAN_MATH_PRIM_ARR_FUN_APPEND_ARRAY_HPP
#define STAN_MATH_PRIM_ARR_FUN_APPEND_ARRAY_HPP

#include <vector>

namespace stan {
  namespace math {

    /**
     * Concatenate two vectors into one.
     *
     * @tparam T Type of elements contained in vector.
     * @param x First vector.
     * @param y Second vector.
     * @return A vector of x and y concatenated together (in that order).
     */
    template <typename T>
    inline std::vector<T> append_array(std::vector<T> x, std::vector<T> y) {
      std::vector<T> z;
      z.reserve(x.size() + y.size());
      z.insert(z.end(), x.begin(), x.end());
      z.insert(z.end(), y.begin(), y.end());
      return z;
    }
  }
}
#endif
