#ifndef STAN_MATH_PRIM_MAT_META_ENABLE_IF_GENERICS_HPP
#define STAN_MATH_PRIM_MAT_META_ENABLE_IF_GENERICS_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>

namespace stan {
template <typename T>
struct eigen_traits {
  const int T::RowsAtCompileTime rows;
  const int T::ColsAtCompileTime cols;
  typedef typename T::PlainObject plain_object;
  typedef typename T::PlainObjectBase plain_object_base;
}
}  // namespace stan

#endif
