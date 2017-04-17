#ifndef STAN_MATH_PRIM_ARR_META_CONTAINER_VIEW_HPP
#define STAN_MATH_PRIM_ARR_META_CONTAINER_VIEW_HPP

#include <stan/math/prim/scal/meta/container_view.hpp>
#include <vector>

namespace stan {
  namespace math {

    /**
     * Template specialization for scalar view of
     * array y with scalar type T2 with proper indexing
     * inferred from input vector x of scalar type T1     
     *
     * @tparam T1 scalar type of input vector
     * @tparam T2 scalar type returned by view.
     */
    template <typename T1, typename T2>
    class container_view<std::vector<T1>, T2> {
      public:
        /**
         * Constructor
         *
         * @param x input vector
         * @param y underlying array 
         */
        container_view(const std::vector<T1>& x, T2* y)
         : y_(y) { }

        /**
         * operator[](int i) returns reference to scalar
         * view indexed at i 
         *
         * @param i index of scalar element
         */
        T2& operator[](int i) {
          return y_[i];
        }
      private:
        T2* y_;
    };
  }
}

#endif
