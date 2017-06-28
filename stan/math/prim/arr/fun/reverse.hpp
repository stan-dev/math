#ifndef STAN_MATH_PRIM_ARR_FUN_REVERSE_HPP
#define STAN_MATH_PRIM_ARR_FUN_REVERSE_HPP

#include <algorithm>
#include <vector>

namespace stan {
  namespace math {

    /**
     * Return a copy of the specified vector in reversed order.
     *
     * @tparam T type of elements in vector
     * @param xs vector to order
     * @return vector in reversed order.
     */
    template <typename T>
    inline std::vector<T> reverse(const std::vector<T>& xs) {
      std::vector<T> reversed(xs.size());
      std::reverse_copy(xs.begin(), xs.end(), reversed.begin());
      return reversed;
    }
  }
}
#endif
