#ifndef STAN_MATH_PRIM_META_BROADCAST_ARRAY_HPP
#define STAN_MATH_PRIM_META_BROADCAST_ARRAY_HPP

#include <stan/math/prim/meta/require_generics.hpp>
#include <stan/math/prim/meta/promote_scalar_type.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
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

  /** \ingroup type_trait
   * We can assign any right hand side which allows for indexing to a
   * broadcast_array. The idea is that the entry for the first index is what
   * gets assigned. The most common use-case should be where the rhs is some
   * container of length 1.
   */
  template <typename Y>
  void operator=(const Y& m) {
    prim_ = m[0];
  }
};

template <typename T, typename S, typename Enable = void>
class empty_broadcast_array {
 public:
  empty_broadcast_array() {}
  /** \ingroup type_trait
   * Not implemented so cannot be called.
   */
  T& operator[](int /*i*/);

  /** \ingroup type_trait
   * Not implemented so cannot be called.
   */
  template <typename Y>
  void operator=(const Y& /*A*/);
};

template <typename ViewElt, typename T>
class empty_broadcast_array<ViewElt, T, require_eigen_t<T>> {
  enum { R = T::RowsAtCompileTime, C = T::ColsAtCompileTime };
  using T_arg = promote_scalar_t<ViewElt, T>;

 public:
  empty_broadcast_array() {}
  /** \ingroup type_trait
   * Not implemented so cannot be called.
   */
  ViewElt& operator[](int /*i*/);
  /** \ingroup type_trait
   * Not implemented so cannot be called.
   */
  ViewElt& operator()(int /*i*/);
  /** \ingroup type_trait
   * Not implemented so cannot be called.
   */
  void operator=(const T_arg& /*A*/);
  /** \ingroup type_trait
   * Not implemented so cannot be called.
   */
  void operator+=(T_arg /*A*/);
  /** \ingroup type_trait
   * Not implemented so cannot be called.
   */
  void operator-=(T_arg /*A*/);
  /** \ingroup type_trait
   * Not implemented so cannot be called.
   */
  T& row(int /*i*/);
  /** \ingroup type_trait
   * Not implemented so cannot be called.
   */
  T& col(int /*i*/);
};
}  // namespace internal
}  // namespace math
}  // namespace stan
#endif
