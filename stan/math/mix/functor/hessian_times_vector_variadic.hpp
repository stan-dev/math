#ifndef STAN_MATH_MIX_FUNCTOR_HESSIAN_TIMES_VECTOR_VARIADIC_HPP
#define STAN_MATH_MIX_FUNCTOR_HESSIAN_TIMES_VECTOR_VARIADIC_HPP

#include <stan/math/mix/functor/internal/autodiff_utils.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/eval.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/apply.hpp>
#include <stan/math/prim/functor/conditional_copy_and_promote.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/to_arena.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/fwd/core.hpp>
#include <ostream>
#include <tuple>

namespace stan {
namespace math {

/**
 * Return the Hessian of a scalar functor times a vector.
 *
 * The functor must be callable as `f(theta, msgs, args...)` with an
 * `fvar<var>` `theta`. This overload is used when none of the inputs contain
 * reverse-mode variables. Forward-mode inputs are not supported.
 *
 * @tparam F type of functor
 * @tparam ThetaVec column vector type
 * @tparam VVec column vector type
 * @tparam Args types of additional arguments
 * @param f functor returning a scalar
 * @param theta argument with respect to which the Hessian is taken
 * @param v vector multiplied by the Hessian
 * @param msgs output stream for functor messages
 * @param args additional arguments passed to the functor
 * @return `H(theta) * v`
 * @throw std::invalid_argument if `theta` and `v` have different sizes
 */
template <typename F, typename ThetaVec, typename VVec, typename... Args,
          require_all_col_vector_t<ThetaVec, VVec>* = nullptr,
          require_all_t<is_vt_not_complex<ThetaVec>, is_vt_not_complex<VVec>>*
          = nullptr,
          require_not_t<
              internal::contains_fvar<ThetaVec, VVec, Args...>>* = nullptr,
          require_not_t<is_any_var_scalar<ThetaVec, VVec, Args...>>* = nullptr>
inline Eigen::VectorXd hessian_times_vector(const F& f, const ThetaVec& theta,
                                            const VVec& v, std::ostream* msgs,
                                            const Args&... args) {
  check_size_match("hessian_times_vector", "theta", theta.size(), "v",
                   v.size());
  const Eigen::Index n = theta.size();
  Eigen::VectorXd hv(n);
  if (n == 0) {
    return hv;
  }

  nested_rev_autodiff nested;
  Eigen::Matrix<fvar<var>, Eigen::Dynamic, 1> theta_fvar(n);
  for (Eigen::Index i = 0; i < n; ++i) {
    theta_fvar.coeffRef(i) = fvar<var>(theta.coeff(i), v.coeff(i));
  }
  fvar<var> fx = f(theta_fvar, msgs, args...);
  grad(fx.d_.vi_);
  for (Eigen::Index i = 0; i < n; ++i) {
    hv.coeffRef(i) = theta_fvar.coeff(i).val_.adj();
  }
  return hv;
}

/**
 * Return a differentiable Hessian-vector product for a scalar functor.
 *
 * For output adjoints `w`, the reverse pass propagates the third-order
 * contraction `grad(v' * H(theta) * w)` to `theta` and `args...`, and
 * `H(theta) * w` to `v`. The functor must be callable with both `fvar<var>`
 * and `fvar<fvar<var>>` inputs. Forward-mode inputs to this function are not
 * supported.
 *
 * @tparam F type of functor
 * @tparam ThetaVec column vector type
 * @tparam VVec column vector type
 * @tparam Args types of additional arguments
 * @param f functor returning a scalar
 * @param theta argument with respect to which the Hessian is taken
 * @param v vector multiplied by the Hessian
 * @param msgs output stream for functor messages
 * @param args additional arguments passed to the functor
 * @return `H(theta) * v`
 * @throw std::invalid_argument if `theta` and `v` have different sizes
 */
template <typename F, typename ThetaVec, typename VVec, typename... Args,
          require_all_col_vector_t<ThetaVec, VVec>* = nullptr,
          require_all_t<is_vt_not_complex<ThetaVec>, is_vt_not_complex<VVec>>*
          = nullptr,
          require_not_t<
              internal::contains_fvar<ThetaVec, VVec, Args...>>* = nullptr,
          require_t<is_any_var_scalar<ThetaVec, VVec, Args...>>* = nullptr>
inline Eigen::Matrix<var, Eigen::Dynamic, 1> hessian_times_vector(
    const F& f, const ThetaVec& theta, const VVec& v, std::ostream* msgs,
    const Args&... args) {
  using Eigen::Index;
  check_size_match("hessian_times_vector", "theta", theta.size(), "v",
                   v.size());
  const Index n = theta.size();
  if (n == 0) {
    return Eigen::Matrix<var, Eigen::Dynamic, 1>(0);
  }

  auto theta_arena = to_arena(theta);
  auto v_arena = to_arena(v);
  auto* args_arena = make_chainable_ptr(std::make_tuple(eval(args)...));
  Eigen::VectorXd hv(n);
  {
    nested_rev_autodiff nested;
    Eigen::Matrix<fvar<var>, Eigen::Dynamic, 1> theta_fvar(n);
    for (Index i = 0; i < n; ++i) {
      theta_fvar.coeffRef(i)
          = fvar<var>(value_of(theta_arena.coeff(i)),
                      value_of(v_arena.coeff(i)));
    }
    fvar<var> fx = math::apply(
        [&](const auto&... inner_args) {
          return f(theta_fvar, msgs, value_of(inner_args)...);
        },
        *args_arena);
    grad(fx.d_.vi_);
    for (Index i = 0; i < n; ++i) {
      hv.coeffRef(i) = theta_fvar.coeff(i).val_.adj();
    }
  }

  arena_t<Eigen::Matrix<var, Eigen::Dynamic, 1>> res = hv;
  reverse_pass_callback(
      [f, msgs, theta_arena, v_arena, args_arena, res]() mutable {
        const Index n = theta_arena.size();
        nested_rev_autodiff nested;
        Eigen::Matrix<fvar<fvar<var>>, Eigen::Dynamic, 1> theta_ffvar(n);
        for (Index i = 0; i < n; ++i) {
          theta_ffvar.coeffRef(i) = fvar<fvar<var>>(
              fvar<var>(theta_arena.coeff(i), v_arena.coeff(i)),
              fvar<var>(res.adj().coeff(i), 0.0));
        }
        auto args_matvar = internal::from_var_value_rec(*args_arena);
        auto args_ffvar
            = internal::shallow_copy_vargs<fvar<fvar<var>>>(args_matvar);
        fvar<fvar<var>> fx = math::apply(
            [&](const auto&... inner_args) {
              return f(theta_ffvar, msgs, inner_args...);
            },
            args_ffvar);
        grad(fx.d_.d_.vi_);
      });
  return res;
}

}  // namespace math
}  // namespace stan
#endif
