#ifndef STAN_MATH_PRIM_SCAL_META_RM_ZEROING_HPP
#define STAN_MATH_PRIM_SCAL_META_RM_ZEROING_HPP

namespace stan {
namespace math {

template <class T>
struct zeroing;  // forward declare

}  // namespace math

/**
 * rm_zeroing replaces <code>zeroing<T></code> with
 * T. This is the general version of the trait.
 */
template <class T>
struct rm_zeroing {
  typedef T type;
};

/**
 * rm_zeroing replaces <code>zeroing<t></code> with
 * T. This is the removal version of the trait.
 */
template <class T>
struct rm_zeroing<math::zeroing<T>> {
  typedef T type;
};

template <class T>
using rm_zeroing_t = typename rm_zeroing<T>::type;

}  // namespace stan
#endif
