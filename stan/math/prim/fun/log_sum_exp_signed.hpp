#ifndef STAN_MATH_PRIM_FUN_LOG_SUM_EXP_SIGNED_HPP
#define STAN_MATH_PRIM_FUN_LOG_SUM_EXP_SIGNED_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/log1p_exp.hpp>
#include <cmath>
#include <vector>

namespace stan {
namespace math {

template <typename T1, typename T2,
          require_all_stan_scalar_t<T1, T2>* = nullptr>
inline std::tuple<return_type_t<T1, T2>, int>
  log_sum_exp_signed(const T1& a, int a_sign, const T2& b, int b_sign) {
  if (a_sign == b_sign) {
    return std::make_tuple(log_sum_exp(a, b), a_sign);
  }
  bool a_larger = (value_of_rec(a) > value_of_rec(b));
  return std::make_tuple(
    a_larger ? log_diff_exp(a, b) : log_diff_exp(b, a),
    a_larger ? a_sign : b_sign);
}

}  // namespace math
}  // namespace stan

#endif
