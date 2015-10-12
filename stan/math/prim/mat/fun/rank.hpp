#ifndef STAN_MATH_PRIM_MAT_FUN_RANK_HPP
#define STAN_MATH_PRIM_MAT_FUN_RANK_HPP

#include <stan/math/prim/mat/err/check_range.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <vector>

namespace stan {
  namespace math {

    /**
     * Return the number of components of v less than v[s].
     *
     * @tparam T Type of elements.
     * @param[in] v Input vector.
     * @param[in] s Position in vector.
     * @return Number of components of v less than v[s].
     */
    template <typename T>
    inline int rank(const std::vector<T> & v, int s) {
      using stan::math::check_range;
      int size = static_cast<int>(v.size());
      check_range("rank", "v", size, s);
      --s;
      int count(0);
      T compare(v[s]);
      for (int i = 0; i < size; ++i)
        if (v[i] < compare)
          ++count;
      return count;
    }

    /**
     * Return the number of components of v less than v[s].
     *
     * @tparam T Type of elements of the vector.
     * @param[in] v Input vector.
     * @param s Index for input vector.
     * @return Number of components of v less than v[s].
     */
    template <typename T, int R, int C>
    inline int rank(const Eigen::Matrix<T, R, C> & v, int s) {
      using stan::math::check_range;
      int size = v.size();
      check_range("rank", "v", size, s);
      --s;
      const T * vv = v.data();
      int count(0);
      T compare(vv[s]);
      for (int i = 0; i < size; ++i)
        if (vv[i] < compare)
          ++count;
      return count;
    }

  }
}
#endif
