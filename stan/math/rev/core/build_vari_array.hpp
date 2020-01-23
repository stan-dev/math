#ifndef STAN_MATH_REV_CORE_BUILD_VARI_ARRAY_HPP
#define STAN_MATH_REV_CORE_BUILD_VARI_ARRAY_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/core/vari.hpp>

namespace stan {
namespace math {

/**
 * Allocates and populates a flat array of vari pointers in the autodiff arena
 * with the varis pointed to by the vars in the input Eigen matrix
 *
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 * @param x input matrix
 * @return Flat array of vari pointers
 */
template <int R, int C>
vari** build_vari_array(const Eigen::Matrix<var, R, C>& x) {
  vari** x_vi_
      = ChainableStack::instance_->memalloc_.alloc_array<vari*>(x.size());
  for (int i = 0; i < x.size(); ++i) {
    x_vi_[i] = x(i).vi_;
  }
  return x_vi_;
}

}  // namespace math
}  // namespace stan
#endif
