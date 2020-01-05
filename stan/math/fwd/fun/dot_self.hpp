#ifndef STAN_MATH_FWD_FUN_DOT_SELF_HPP
#define STAN_MATH_FWD_FUN_DOT_SELF_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/fun/dot_product.hpp>

namespace stan {
namespace math {

template <typename T, int R, int C>
inline fvar<T> dot_self(const Eigen::Matrix<fvar<T>, R, C>& v) {
  check_vector("dot_self", "v", v);
  return dot_product(v, v);
}

}  // namespace math
}  // namespace stan
#endif
