#ifndef STAN_MATH_PRIM_META_BROADCAST_ARRAY_HPP
#define STAN_MATH_PRIM_META_BROADCAST_ARRAY_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
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

template <typename T, typename S>
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

template <typename ViewElt, typename OpElt, int R, int C>
class empty_broadcast_array<ViewElt, Eigen::Matrix<OpElt, R, C> > {
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
  void operator=(const Eigen::Matrix<ViewElt, R, C>& /*A*/);
  /** \ingroup type_trait
   * Not implemented so cannot be called.
   */
  void operator+=(Eigen::Matrix<ViewElt, R, C> /*A*/);
  /** \ingroup type_trait
   * Not implemented so cannot be called.
   */
  void operator-=(Eigen::Matrix<ViewElt, R, C> /*A*/);
  /** \ingroup type_trait
   * Not implemented so cannot be called.
   */
  Eigen::Matrix<ViewElt, 1, C>& row(int /*i*/);
  /** \ingroup type_trait
   * Not implemented so cannot be called.
   */
  Eigen::Matrix<ViewElt, R, 1>& col(int /*i*/);
};
}  // namespace internal
}  // namespace math
}  // namespace stan
#endif
