#ifndef STAN_MATH_PRIM_SCAL_META_BROADCAST_ARRAY_HPP
#define STAN_MATH_PRIM_SCAL_META_BROADCAST_ARRAY_HPP

#include <stdexcept>

namespace stan {
namespace math {
namespace internal {
template <typename T, typename = void>
class broadcast_array {
 private:
  T prim_;

 public:
   template <typename T1, typename T2>
   using is_same_op = std::is_same<std::decay_t<T1>, std::decay_t<T2>>;

  template <typename T1, typename = std::enable_if_t<is_same_op<T, T1>::value>>
  explicit broadcast_array(T1&& prim) : prim_(std::forward<T1>(prim)) {}

  auto&& operator[](int /*i*/) { return prim_; }

  /**
   * We can assign any right hand side which allows for indexing to a
   * broadcast_array. The idea is that the entry for the first index is what
   * gets assigned. The most common use-case should be where the rhs is some
   * container of length 1.
   */
  template <typename Y, typename = std::enable_if_t<is_same_op<T, Y>::value>>
  void operator=(Y&& m) {
    prim_ = std::forward<decltype(m[0])>(m[0]);
  }
};

template <typename ViewElt, typename OpElt, typename = void>
class empty_broadcast_array {
 public:
  empty_broadcast_array() {}
  /**
   * Not implemented so cannot be called.
   */
  ViewElt& operator[](int /*i*/);

  /**
   * Not implemented so cannot be called.
   */
  template <typename Y>
  void operator=(const Y& /*A*/);
};
}  // namespace internal
}  // namespace math
}  // namespace stan

#endif
