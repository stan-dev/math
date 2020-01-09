#ifndef STAN_MATH_PRIM_FUN_GET_LP_HPP
#define STAN_MATH_PRIM_FUN_GET_LP_HPP

#include <stan/math/prim/fun/accumulator.hpp>

namespace stan {
namespace math {

template <typename T_lp, typename T_lp_accum>
inline return_type_t<T_lp, T_lp_accum> get_lp(
    const T_lp& lp, const accumulator<T_lp_accum>& lp_accum) {
  return lp + lp_accum.sum();
}

}  // namespace math
}  // namespace stan
#endif
