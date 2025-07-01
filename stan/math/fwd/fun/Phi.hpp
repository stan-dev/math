#ifndef STAN_MATH_FWD_FUN_PHI_HPP
#define STAN_MATH_FWD_FUN_PHI_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/fun/pow.hpp>
#include <stan/math/fwd/fun/exp.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/Phi.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T, require_fvar_t<T>* = nullptr>
inline auto Phi(T&& x) {
  return std::decay_t<T>(Phi(x.val_), x.d_ * exp(x.val_ * x.val_ / -2.0) * INV_SQRT_TWO_PI);
}

}  // namespace math
}  // namespace stan
#endif
