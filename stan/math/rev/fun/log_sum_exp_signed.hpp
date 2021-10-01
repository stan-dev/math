#ifndef STAN_MATH_REV_FUN_LOG_SUM_EXP_SIGNED_HPP
#define STAN_MATH_REV_FUN_LOG_SUM_EXP_SIGNED_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/core/typedefs.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/inv_logit.hpp>
#include <stan/math/prim/fun/log_sum_exp.hpp>
#include <cmath>
#include <vector>

namespace stan {
namespace math {

/**
 * Returns the log sum of exponentials of the input, with the sign of
 *   the exponential provided (i.e., whether to add or subtract)
 *
 * @tparam T1 A type inheriting from EigenBase with scalar type var
 * @tparam T2 A type inheriting from EigenBase with scalar type int
 * @param v input
 * @param signs sign of exponentiated input
 */
template <typename T1, typename T2, require_eigen_st<is_var, T1>* = nullptr,
          require_eigen_st<std::is_integral, T2>* = nullptr,
          require_not_var_matrix_t<T1>* = nullptr>
inline var log_sum_exp_signed(const T1& v, const T2& signs) {
  arena_t<decltype(v)> arena_v = v;
  arena_t<decltype(v.val())> arena_v_val = arena_v.val();
  arena_t<T2> arena_signs = signs;
  var res = log_sum_exp_signed(arena_v.val(), signs);

  reverse_pass_callback([arena_v, arena_v_val, arena_signs, res]() mutable {
    arena_v.adj()
        += res.adj() *
              (arena_v_val.array().val() - res.val()).exp()
                                                    .matrix()
                                                    .cwiseProduct(arena_signs);
  });

  return res;
}

/**
 * Returns the log sum of exponentials of the input, with the sign of
 *   the exponential provided (i.e., whether to add or subtract)
 *
 * @tparam T1 A `var_value` with an input vector or matrix
 * @tparam T2 A type inheriting from EigenBase with scalar type int
 * @param x input
 * @param signs sign of exponentiated input
 */
template <typename T1, typename T2, require_var_matrix_t<T1>* = nullptr,
          require_eigen_st<std::is_integral, T2>* = nullptr>
inline var log_sum_exp_signed(const T1& x, const T2& signs) {
  return make_callback_vari(log_sum_exp_signed(x.val(), signs),
                            [x, signs](const auto& res) mutable {
    x.adj() += res.adj() *
                (x.val().array().val() - res.val()).exp().matrix()
                                                         .cwiseProduct(signs);
  });
}

/**
 * Returns the log sum of exponentials of the input, with the sign of
 *   the exponential provided (i.e., whether to add or subtract)
 *
 * @tparam T1 Type of input vector or matrix
 * @tparam T2 A type inheriting from EigenBase with scalar type int
 * @param x matrix
 * @param signs sign of exponentiated input
 */
template <typename T1, typename T2,
          require_std_vector_st<is_var, T1>* = nullptr,
          require_std_vector_st<std::is_integral, T2>* = nullptr>
inline auto log_sum_exp_signed(const T1& x, const T2& signs) {
  return apply_vector_unary<T1>::reduce(
      x, [&](const auto& v) {
          Eigen::Map<const Eigen::VectorXi> int_vec_map(signs.data(),
                                                        signs.size());
          return log_sum_exp_signed(v, int_vec_map); });
}

}  // namespace math
}  // namespace stan
#endif
