#ifndef STAN_MATH_FWD_PROB_STD_NORMAL_LOG_QF_HPP
#define STAN_MATH_FWD_PROB_STD_NORMAL_LOG_QF_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/prob/std_normal_log_qf.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T>
inline fvar<T> std_normal_log_qf(const fvar<T>& p) {
  const T xv = std_normal_log_qf(p.val_);
  return fvar<T>(xv, p.d_ * exp(p.val_ - std_normal_lpdf(xv)));
}
}  // namespace math
}  // namespace stan
#endif
