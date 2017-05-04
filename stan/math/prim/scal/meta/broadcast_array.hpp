#ifndef STAN_MATH_PRIM_SCAL_META_BROADCAST_ARRAY_HPP
#define STAN_MATH_PRIM_SCAL_META_BROADCAST_ARRAY_HPP

#include <stdexcept>

namespace stan {
  namespace math {
    namespace detail {
      template <typename T>
      class broadcast_array {
      private:
        T& prim_;

      public:
        explicit broadcast_array(T& prim): prim_(prim) {}

        T& operator[] (int /*i*/) {
          return prim_;
        }
        int size() { return 1; }
      };

      template <>
      class broadcast_array<void> {
      public:
        broadcast_array() {}
        double& operator[] (int /*i*/) {
          throw std::logic_error("Don't do this");
        }
        int size() { return 0; }
      };
    }  // end namespace detail
  }  // end namespace math
}  // end namespace stan

#endif
