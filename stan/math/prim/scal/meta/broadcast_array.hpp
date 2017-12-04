#ifndef STAN_MATH_PRIM_SCAL_META_BROADCAST_ARRAY_HPP
#define STAN_MATH_PRIM_SCAL_META_BROADCAST_ARRAY_HPP

#include <stdexcept>

namespace stan {
  namespace math {
    namespace internal {
      template <typename T>
      class broadcast_array {
      private:
        T& prim_;

      public:
        explicit broadcast_array(T& prim): prim_(prim) {}

        T& operator[](int /*i*/) {
          return prim_;
        }

        void operator= (Eigen::Matrix<T, -1, -1> A) {
          this[0][0] = A(0,0);
        }
      };

      template <typename T, typename S>
      class empty_broadcast_array {
      public:
        empty_broadcast_array() {}
       /**
        * Not implemented so cannot be called.
        */
        T& operator[] (int /*i*/);

        void operator= (Eigen::Matrix<T, -1, -1> /*A*/);
      };
    }
  }
}

#endif
