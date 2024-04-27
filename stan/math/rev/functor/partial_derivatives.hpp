#ifndef STAN_MATH_REV_FUNCTOR_PARTIAL_DERIVATIVES_HPP
#define STAN_MATH_REV_FUNCTOR_PARTIAL_DERIVATIVES_HPP

#ifdef STAN_MODEL_FVAR_VAR
#include <stan/math/mix/functor/hessian.hpp>
#endif
#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/functor/apply_scalar_ternary.hpp>
#include <stan/math/rev/functor/apply_scalar_unary.hpp>
#include <stan/math/rev/functor/gradient.hpp>
#include <stan/math/rev/functor/finite_diff_hessian_auto.hpp>
#include <vector>

namespace stan {
namespace math {
namespace internal {
  template <typename T, require_var_t<T>* = nullptr>
  double get_adj(const T& x) {
    return x.adj();
  }
  template <typename T, require_arithmetic_t<T>* = nullptr>
  double get_adj(const T& x) {
    return 0;
  }

  template <typename T, require_container_t<T>* = nullptr>
  auto get_adj(const T& x) {
  return apply_vector_unary<T>::apply(x,
    [&](const auto& v) {
      return v.unaryExpr([&](const auto& y) { return get_adj(y); });
    });
  }
}

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
  internal::finite_diff_hessian_auto(serial_functor, serialised_args,
                                      rtn_value, grad, hess);
#endif

  auto grad_deserializer = to_deserializer(grad);
  auto arena_args = std::make_tuple(stan::math::to_arena(args)...);
  var rtn = rtn_value;
  auto rtn_grad = std::make_tuple(
    arena_t<promote_scalar_t<var, plain_type_t<TArgs>>>(
      grad_deserializer.read(args))...
  );

  reverse_pass_callback([rtn, arena_args, rtn_grad, hess]() mutable {
    std::vector<double> grad_adjs = math::apply(
      [&](auto&&... grad_args) {
        return serialize<double>(internal::get_adj(grad_args)...);
      }, std::forward<decltype(rtn_grad)>(rtn_grad));
    Eigen::VectorXd hess_grad = hess * as_column_vector_or_scalar(grad_adjs);
    auto hess_deserial = to_deserializer(hess_grad);

    math::for_each([&](auto&& curr_arg, auto&& arg_grad) {
      auto arg_grad_hess
        = hess_deserial.read(plain_type_t<decltype(curr_arg)>(curr_arg));

      if (!is_constant_all<decltype(curr_arg)>::value) {
        apply_scalar_ternary(
          [rtn](auto&& arg, auto&& grad, auto&& arg_hess_cap){
            forward_as<promote_scalar_t<var, decltype(arg)>>(arg).adj()
              += rtn.adj() * grad.val() + arg_hess_cap;
            return 0;
          },
          std::forward<decltype(curr_arg)>(curr_arg),
          std::forward<decltype(arg_grad)>(arg_grad),
          std::forward<decltype(arg_grad_hess)>(arg_grad_hess)
        );
      }
    }, std::forward<decltype(arena_args)>(arena_args),
      std::forward<decltype(rtn_grad)>(rtn_grad));
  });

  using return_t = std::tuple<promote_scalar_t<var, plain_type_t<TArgs>>...>;
  return std::make_tuple(rtn_value, return_t(rtn_grad));
}


//STAN_MODEL_FVAR_VAR

}  // namespace math
}  // namespace stan
#endif
