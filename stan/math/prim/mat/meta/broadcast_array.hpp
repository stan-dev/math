#ifndef STAN_MATH_PRIM_MAT_META_BROADCAST_ARRAY_HPP
#define STAN_MATH_PRIM_MAT_META_BROADCAST_ARRAY_HPP

#include <stan/math/prim/scal/meta/broadcast_array.hpp>
#include <stan/math/prim/scal/meta/require_generics.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stdexcept>

namespace stan {
namespace math {
namespace internal {
template <typename ViewElt, typename T>
class empty_broadcast_array<ViewElt, T, require_eigen_t<T>> {
  enum { R = T::RowsAtCompileTime, C = T::ColsAtCompileTime };
  using T_arg = std::conditional_t<
      std::is_same<typename Eigen::internal::traits<T>::XprKind,
                   Eigen::MatrixXpr>::value,
      Eigen::Matrix<ViewElt, R, C>, Eigen::Array<ViewElt, R, C>>;

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
