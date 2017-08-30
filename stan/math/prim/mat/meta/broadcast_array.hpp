#ifndef STAN_MATH_PRIM_MAT_META_BROADCAST_ARRAY_HPP
#define STAN_MATH_PRIM_MAT_META_BROADCAST_ARRAY_HPP

#include <stan/math/prim/scal/meta/broadcast_array.hpp>
#include <stdexcept>
#include <Eigen/Dense>

namespace stan {
  namespace math {
    namespace internal {
      template <typename ViewElt, typename OpElt, int R, int C>
      class empty_broadcast_array<ViewElt, Eigen::Matrix<OpElt, R, C>> {
      public:
        empty_broadcast_array() {}
        ViewElt& operator[] (int /*i*/) {
          throw std::logic_error("Don't do this");
        }
        Eigen::Matrix<ViewElt, 1, C>& row(int /*i*/) {
          throw std::logic_error("Don't do this");
        }
      };
    }
  }
}

#endif