#ifndef STAN_MATH_REV_CORE_BUILD_VARI_ARRAY_HPP
#define STAN_MATH_REV_CORE_BUILD_VARI_ARRAY_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/core/vari.hpp>

namespace stan {
namespace math {

template <int R, int C>
vari** build_vari_array(const Eigen::Matrix<var, R, C>& x) {
  vari** x_vi_ = ChainableStack::instance().memalloc_.alloc_array<vari*>(x.size());
  for (int i = 0; i < x.size(); ++i) {
    x_vi_[i] = x.data()[i].vi_;
  }
  return x_vi_;
}

}  // namespace math
}  // namespace stan
#endif
