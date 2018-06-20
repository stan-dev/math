#ifndef STAN_MATH_PRIM_SCAL_META_RM_ZEROING_HPP
#define STAN_MATH_PRIM_SCAL_META_RM_ZEROING_HPP

namespace stan {
namespace math {

template <class T>
struct zeroing;  // forward declare

}  // namespace math

/**
 * trait to change z_fvar to fvar and z<var> to var
 * 
 * this is the end condition
 */
template <class T>
struct rm_zeroing {
  typedef T type;
};

/**
 * trait to change z_fvar to fvar and z<var> to var
 */
template <class T>
struct rm_zeroing<math::zeroing<T>> {
  typedef T type;
};

template <class T>
using rm_zeroing_t = typename rm_zeroing<T>::type;

}  // namespace stan
#endif
