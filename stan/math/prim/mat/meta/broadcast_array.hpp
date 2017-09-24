#include <stan/math/prim/scal/meta/broadcast_array.hpp>
#include <Eigen/Dense>
#include <stdexcept>

#ifndef STAN_MATH_PRIM_MAT_META_BROADCAST_ARRAY_HPP
#define STAN_MATH_PRIM_MAT_META_BROADCAST_ARRAY_HPP

namespace stan {
  namespace math {
    namespace internal {
      template <typename ViewElt, typename OpElt, int R, int C>
      class empty_broadcast_array<ViewElt, Eigen::Matrix<OpElt, R, C>> {
      public:
        empty_broadcast_array() {}
        // We provide stub methods for the empty_broadcast_array which should
        // never be called.
        ViewElt& operator[](int /*i*/);
        Eigen::Matrix<ViewElt, 1, C>& row(int /*i*/);
        Eigen::Matrix<ViewElt, R, 1>& col(int /*i*/);
      };
    }
  }
}

#endif
