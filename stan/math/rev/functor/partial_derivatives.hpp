#ifndef STAN_MATH_REV_FUNCTOR_PARTIAL_DERIVATIVES_HPP
#define STAN_MATH_REV_FUNCTOR_PARTIAL_DERIVATIVES_HPP

#ifdef STAN_MODEL_FVAR_VAR
#include <stan/math/mix/functor/hessian.hpp>
#endif
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/functor/gradient.hpp>
#include <stan/math/rev/functor/finite_diff_hessian_auto.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Given a function and arguments, return the value and partial derivatives.
 *
 * @tparam T Argument type
 * @tparam F Function type
 * @param f Function
 * @param[in] x Argument vector
 * @param[in] n Index of argument with which to take derivative
 * @param[out] fx Value of function applied to argument
 * @param[out] dfx_dxn Value of partial derivative
 */
template <typename F, typename... TArgs,
          require_all_st_arithmetic<TArgs...>* = nullptr>
inline auto partial_derivatives(const F& func, const TArgs&... args) {
  std::vector<double> serialised_args = serialize<double>(args...);

  auto serial_functor = [&](const auto& v) {
    auto v_deserializer = to_deserializer(v);
    return func(v_deserializer.read(args)...);
  };
  double rtn_value;
  std::vector<double> grad;
  gradient(serial_functor, serialised_args, rtn_value, grad);

  auto grad_deserializer = to_deserializer(grad);
  return std::make_tuple(rtn_value,
                         std::make_tuple(grad_deserializer.read(args)...));
}

template <typename F, typename... TArgs,
          require_any_st_var<TArgs...>* = nullptr>
inline auto partial_derivatives(const F& func, const TArgs&... args) {
  std::vector<double> serialised_args = serialize<double>(value_of(args)...);

  auto serial_functor = [&](const auto& v) {
    auto v_deserializer = to_deserializer(v);
    return func(v_deserializer.read(args)...);
  };

  double rtn_value;
  std::vector<double> grad;
  Eigen::MatrixXd hess;

#ifdef STAN_MODEL_FVAR_VAR
  hessian(serial_functor, serialised_args, rtn_value, grad, hess);
#else
  finite_diff_hessian_auto(serial_functor, serialised_args, rtn_value, grad, hess);
#endif

  auto arena_args = std::make_tuple(stan::math::to_arena(args)...);

  var rtn = rtn_value;
  reverse_pass_callback([rtn, arena_args, grad, hess]() mutable {

  });
  auto grad_deserializer = to_deserializer(grad);
  return std::make_tuple(rtn_value,
                         std::make_tuple(grad_deserializer.read(args)...));
}


//STAN_MODEL_FVAR_VAR

}  // namespace math
}  // namespace stan
#endif
