#ifndef STAN_MATH_PRIM_MAT_META_BROADCAST_ARRAY_HPP
#define STAN_MATH_PRIM_MAT_META_BROADCAST_ARRAY_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/meta/is_eigen.hpp>
#include <stan/math/prim/scal/meta/broadcast_array.hpp>
#include <stdexcept>
#include <type_traits>

namespace stan {
namespace math {
namespace internal {
template <typename ViewElt, typename OpElt>
class empty_broadcast_array<ViewElt, OpElt,
                            std::enable_if_t<is_eigen<OpElt>::value>> {
 public:
  empty_broadcast_array() {}
  /**
   * Not implemented so cannot be called.
   */
  ViewElt& operator[](int /*i*/);
  /**
   * Not implemented so cannot be called.
   */
  ViewElt& operator()(int /*i*/);
  /**
   * Not implemented so cannot be called.
   */
   template <typename Eig>
  void operator=(const Eig& /*A*/);
  /**
   * Not implemented so cannot be called.
   */
  template <typename Eig>
  void operator+=(Eig /*A*/);
  /**
   * Not implemented so cannot be called.
   */
   template <typename Eig>
   void operator-=(Eig /*A*/);
  /**
   * Not implemented so cannot be called.
   */
  Eigen::Matrix<ViewElt, OpElt::RowsAtCompileTime, OpElt::ColsAtCompileTime>& row(int /*i*/);
  /**
   * Not implemented so cannot be called.
   */
  Eigen::Matrix<ViewElt, OpElt::RowsAtCompileTime, OpElt::ColsAtCompileTime>& col(int /*i*/);
};
}  // namespace internal
}  // namespace math
}  // namespace stan
#endif
