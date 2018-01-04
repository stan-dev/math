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
  explicit broadcast_array(T& prim) : prim_(prim) {}

  T& operator[](int /*i*/) { return prim_; }
};

template <typename T, typename S>
class empty_broadcast_array {
 public:
  empty_broadcast_array() {}
  /**
   * Not implemented so cannot be called.
   */
  T& operator[](int /*i*/);
};

template <typename T, typename S>
class empty_nested_broadcast_array {
  typedef empty_broadcast_array<T, S> empty_broadcast_array_t;

 public:
  empty_nested_broadcast_array() {}
  /**
   * Not implemented so cannot be called.
   */
  empty_broadcast_array_t operator[](int /*i*/) {
    return empty_broadcast_array_t();
  }
};
}  // namespace internal
}  // namespace math
}  // namespace stan

#endif
