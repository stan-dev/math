#ifndef STAN_MATH_MIX_FUNCTOR_LAPLACE_MARGINAL_DENSITY_HPP
#define STAN_MATH_MIX_FUNCTOR_LAPLACE_MARGINAL_DENSITY_HPP
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/mix/functor/laplace_likelihood.hpp>
#include <stan/math/mix/functor/laplace_marginal_density_estimator.hpp>
#include <stan/math/rev/functor/conditional_copy_and_promote.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun.hpp>
#include <stan/math/rev/functor.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/functor/iter_tuple_nested.hpp>
#include <unsupported/Eigen/MatrixFunctions>
#include <cmath>

/**
 * @file
 * Reference for calculations of marginal and its gradients:
 * Margossian et al (2020), https://arxiv.org/abs/2004.12550
 * and Margossian (2023), https://arxiv.org/pdf/2306.14976
 */

namespace stan {
namespace math {

/**
 * For a latent Gaussian model with global parameters phi, latent
 * variables theta, and observations y, this function computes
 * an approximation of the log marginal density, p(y | phi).
 * This is done by marginalizing out theta, using a Laplace
 * approximation. The latter is obtained by finding the mode,
 * using a custom Newton method, and the Hessian of the likelihood.
 *
 * The convergence criterion for the Newton/Wolfe loop is a small change in
 * the optimization objective (not the final Laplace-corrected log marginal
 * density). The user controls the tolerance (i.e. threshold under which the
 * change is deemed small enough) and maximum number of steps.
 *
 * Wrapper for when the hyperparameters are passed as a double.
 *
 * @tparam LLFun Type with a valid `operator(ThetaVec, InnerLLTupleArgs)`
 * where `InnerLLTupleArgs` are the elements of `LLTupleArgs`
 * @tparam LLTupleArgs A tuple whose elements follow the types required for
 * `LLFun`
 * \laplace_common_template_args
 * @param[in] ll_fun A log likelihood functor
 * @param[in] ll_args Tuple containing parameters for `LLFun`
 * \laplace_common_args
 * @param[in] options A set of options for tuning the solver
 * \msg_arg
 * @return the log marginal density, p(y | phi)
 */
template <
    typename LLFun, typename LLTupleArgs, typename CovarFun, typename CovarArgs,
    bool InitTheta,
    require_t<is_all_arithmetic_scalar<CovarArgs, LLTupleArgs>>* = nullptr>
inline auto laplace_marginal_density(LLFun&& ll_fun, LLTupleArgs&& ll_args,
                                     CovarFun&& covariance_function,
                                     CovarArgs&& covar_args,
                                     const laplace_options<InitTheta>& options,
                                     std::ostream* msgs) {
  Eigen::MatrixXd covariance = stan::math::apply(
      [msgs, &covariance_function](auto&&... args) {
        return covariance_function(std::forward<decltype(args)>(args)..., msgs);
      },
      std::forward<CovarArgs>(covar_args));
  return internal::laplace_marginal_density_est(
             std::forward<LLFun>(ll_fun), std::forward<LLTupleArgs>(ll_args),
             std::move(covariance), options, msgs)
      .lmd;
}

namespace internal {

template <bool ZeroInput = false, typename Output, typename Input,
          require_tuple_t<Output>* = nullptr, require_tuple_t<Input>* = nullptr,
          require_t<std::bool_constant<
              std::tuple_size_v<std::decay_t<Output>> == 0>>* = nullptr,
          require_t<std::bool_constant<
              std::tuple_size_v<std::decay_t<Input>> == 0>>* = nullptr>
inline constexpr void copy_compute_s2(Output&& output, Input&& input) {}

/**
 * Copies the adjoints from the input to the output, scaling them by 0.5.
 * @tparam ZeroInput If true, the adjoints of the input will be set to zero
 * @tparam Output A tuple or type where all scalar types are `arithmetic` types
 * @tparam Input A tuple or type where all scalar types are `var` types
 * @param output The output to which the adjoints will be added
 * @param input The input from which the adjoints will be collected
 */
template <bool ZeroInput = false, typename Output, typename Input,
          require_t<is_all_arithmetic_scalar<Output>>* = nullptr,
          require_t<is_any_var_scalar<Input>>* = nullptr>
inline constexpr void copy_compute_s2(Output&& output, Input&& input) {
  if constexpr (is_tuple_v<Output> && is_tuple_v<Input>) {
    static_assert(
        std::tuple_size<std::decay_t<Output>>::value
            == std::tuple_size<std::decay_t<Input>>::value,
        "INTERNAL ERROR:(laplace_marginal_lpdf) copy_compute_s2 called on "
        "tuples of different sizes. This is an internal error, please report "
        "it: "
        "https://github.com/stan-dev/math/issues");
  }
  return iter_tuple_nested(
      [](auto&& output_i, auto&& input_i) {
        using output_i_t = std::decay_t<decltype(output_i)>;
        if constexpr (is_std_vector_v<output_i_t>) {
          using dbl_map_t = Eigen::Map<Eigen::Matrix<double, -1, 1>>;
          using var_map_t = Eigen::Map<Eigen::Matrix<var, -1, 1>>;
          var_map_t input_map(input_i.data(), input_i.size());
          dbl_map_t(output_i.data(), output_i.size()).array()
              += 0.5 * input_map.adj().array();
          if constexpr (ZeroInput) {
            input_map.adj().setZero();
          }
        } else if constexpr (is_eigen_v<output_i_t>) {
          output_i.array() += 0.5 * input_i.adj().array();
          if constexpr (ZeroInput) {
            input_i.adj().setZero();
          }
        } else if constexpr (is_stan_scalar_v<output_i_t>) {
          output_i += (0.5 * input_i.adj());
          if constexpr (ZeroInput) {
            input_i.adj() = 0;
          }
        } else {
          static_assert(
              sizeof(std::decay_t<output_i_t>*) == 0,
              "INTERNAL ERROR:(laplace_marginal_lpdf) copy_compute_s2 was "
              "not able to deduce the actions needed for the given type. "
              "This is an internal error, please report it: "
              "https://github.com/stan-dev/math/issues");
        }
      },
      std::forward<Output>(output), std::forward<Input>(input));
}

}  // namespace internal

/**
 * For a latent Gaussian model with global parameters phi, latent
 * variables theta, and observations y, this function computes
 * an approximation of the log marginal density, p(y | phi).
 * This is done by marginalizing out theta, using a Laplace
 * approximation. The latter is obtained by finding the mode,
 * using a custom Newton method, and the Hessian of the likelihood.
 *
 * The convergence criterion for the Newton/Wolfe loop is a small change in
 * the optimization objective (not the final Laplace-corrected log marginal
 * density). The user controls the tolerance (i.e. threshold under which the
 * change is deemed small enough) and maximum number of steps.
 *
 * Wrapper for when the global parameter is passed as a double.
 *
 * @tparam LLFun Type with a valid `operator(ThetaVec,  InnerLLTupleArgs)`
 * where `InnerLLTupleArgs` are the elements of `LLTupleArgs`
 * @tparam LLTupleArgs A tuple whose elements follow the types required for
 * `LLFun`
 * \laplace_common_template_args
 * @param[in] ll_fun A log likelihood functor
 * @param[in] ll_args Tuple containing parameters for `LLFun`
 * \laplace_common_args
 * @param[in] options A set of options for tuning the solver
 * \msg_arg
 * @return the log marginal density, p(y | phi)
 */
template <typename LLFun, typename LLTupleArgs, typename CovarFun,
          typename CovarArgs, bool InitTheta,
          require_t<is_any_var_scalar<LLTupleArgs, CovarArgs>>* = nullptr>
inline auto laplace_marginal_density(LLFun&& ll_fun, LLTupleArgs&& ll_args,
                                     CovarFun&& covariance_function,
                                     CovarArgs&& covar_args,
                                     const laplace_options<InitTheta>& options,
                                     std::ostream* msgs) {
  auto covar_args_refs = to_ref(std::forward<CovarArgs>(covar_args));
  auto ll_args_refs = to_ref(std::forward<LLTupleArgs>(ll_args));
  // Solver 1, 2, 3
  constexpr bool ll_args_contain_var = is_any_var_scalar<LLTupleArgs>::value;
  auto partial_parm = internal::make_zeroed_arena(ll_args_refs);
  auto covar_args_adj = internal::make_zeroed_arena(covar_args_refs);
  double lmd = 0.0;
  {
    nested_rev_autodiff nested;
    // Make one hard copy here
    auto ll_args_copy = internal::deep_copy_vargs<var>(ll_args_refs);
    auto covar_args_copy = internal::deep_copy_vargs<var>(covar_args_refs);
    auto covariance = stan::math::apply(
        [&covariance_function, &msgs](auto&&... args) {
          if constexpr (is_any_var_scalar_v<decltype(args)...>) {
            return to_var_value(covariance_function(args..., msgs));
          } else {
            return covariance_function(args..., msgs);
          }
        },
        covar_args_copy);
    decltype(auto) covariance_val = value_of(covariance);
    decltype(auto) ll_args_vals = value_of(ll_args_copy);
    auto md_est = internal::laplace_marginal_density_est(
        ll_fun, ll_args_vals, covariance_val, options, msgs);
    auto ll_args_filter = internal::filter_var_scalar_types(ll_args_copy);
    // tuple of references to var types
    // Solver 1, 2
    const bool solver_1_or_2
        = md_est.solver_used == 1 || md_est.solver_used == 2;
    arena_t<Eigen::MatrixXd> R(md_est.theta.size() * solver_1_or_2,
                               md_est.theta.size() * solver_1_or_2);
    // Solver 3
    arena_t<Eigen::MatrixXd> LU_solve_covariance(
        covariance.rows() * (md_est.solver_used == 3),
        covariance.cols() * (md_est.solver_used == 3));
    // Solver 1, 2, 3
    arena_t<Eigen::VectorXd> s2(md_est.theta.size());
    using stan::math::internal::ZeroOut;
    if (md_est.solver_used == 1) {
      if (options.hessian_block_size == 1) {
        arena_t<Eigen::MatrixXd> tmp = md_est.W_r.toDense();
        md_est.L.template triangularView<Eigen::Lower>().solveInPlace(tmp);
        R.noalias() = tmp.transpose() * tmp;
        arena_t<Eigen::MatrixXd> C
            = md_est.L.template triangularView<Eigen::Lower>().solve(
                md_est.W_r * covariance_val);
        if constexpr (ll_args_contain_var) {
          arena_t<Eigen::MatrixXd> A = covariance_val - C.transpose() * C;
          auto s2_tmp = laplace_likelihood::compute_s2(
              ll_fun, md_est.theta, A, options.hessian_block_size, ll_args_copy,
              msgs);
          s2.deep_copy(s2_tmp);
          internal::copy_compute_s2<ZeroOut>(partial_parm, ll_args_filter);
        } else {
          s2.deep_copy(
              (0.5
               * (covariance_val.diagonal() - (C.transpose() * C).diagonal())
                     .cwiseProduct(laplace_likelihood::third_diff(
                         ll_fun, md_est.theta, ll_args_vals, msgs))));
        }

      } else {
        arena_t<Eigen::MatrixXd> tmp = md_est.W_r.toDense();
        md_est.L.template triangularView<Eigen::Lower>().solveInPlace(tmp);
        R.noalias() = tmp.transpose() * tmp;
        arena_t<Eigen::MatrixXd> C
            = md_est.L.template triangularView<Eigen::Lower>().solve(
                md_est.W_r * covariance_val);
        arena_t<Eigen::MatrixXd> A = covariance_val - C.transpose() * C;
        auto s2_tmp = laplace_likelihood::compute_s2(ll_fun, md_est.theta, A,
                                                     options.hessian_block_size,
                                                     ll_args_copy, msgs);
        s2.deep_copy(s2_tmp);
        internal::copy_compute_s2<ZeroOut>(partial_parm, ll_args_filter);
      }
    } else if (md_est.solver_used == 2) {
      R = md_est.W_r
          - md_est.W_r * md_est.K_root
                * md_est.L.transpose()
                      .template triangularView<Eigen::Upper>()
                      .solve(
                          md_est.L.template triangularView<Eigen::Lower>()
                              .solve(md_est.K_root.transpose() * md_est.W_r));

      arena_t<Eigen::MatrixXd> C
          = md_est.L.template triangularView<Eigen::Lower>().solve(
              md_est.K_root.transpose());
      auto s2_tmp = laplace_likelihood::compute_s2(
          ll_fun, md_est.theta, (C.transpose() * C).eval(),
          options.hessian_block_size, ll_args_copy, msgs);
      s2.deep_copy(s2_tmp);
      internal::copy_compute_s2<ZeroOut>(partial_parm, ll_args_filter);
    } else {  // options.solver with LU decomposition
      LU_solve_covariance = md_est.LU.solve(covariance_val);
      auto I_minus_BinvKW
          = Eigen::MatrixXd::Identity(md_est.W_r.rows(), md_est.W_r.cols())
            - LU_solve_covariance * md_est.W_r;
      R = md_est.W_r * I_minus_BinvKW;  // == W - W B^{-1} K W
      arena_t<Eigen::MatrixXd> A
          = covariance_val - covariance_val * md_est.W_r * LU_solve_covariance;
      auto s2_tmp = laplace_likelihood::compute_s2(ll_fun, md_est.theta, A,
                                                   options.hessian_block_size,
                                                   ll_args_copy, msgs);
      s2.deep_copy(s2_tmp);
      internal::copy_compute_s2<ZeroOut>(partial_parm, ll_args_filter);
    }
    if constexpr (is_any_var_scalar_v<scalar_type_t<CovarArgs>>) {
      arena_t<Eigen::MatrixXd> K_adj_arena
          = 0.5 * md_est.a * md_est.a.transpose() - 0.5 * R
            + s2 * md_est.theta_grad.transpose()
            - (R * (covariance.val() * s2)) * md_est.theta_grad.transpose();
      var Z = make_callback_var(
          0.0, [covariance, K_adj_arena](auto&& vi) mutable {
            covariance.adj().array() += vi.adj() * K_adj_arena.array();
          });
      grad(Z.vi_);
      auto covar_args_filter
          = internal::filter_var_scalar_types(covar_args_copy);
      internal::collect_adjoints(covar_args_adj, covar_args_filter);
    }
    if constexpr (ll_args_contain_var) {
      laplace_likelihood::ll_arg_grad(ll_fun, md_est.theta, ll_args_copy, msgs);
      internal::collect_adjoints<ZeroOut>(partial_parm, ll_args_filter);
      arena_t<Eigen::VectorXd> v;
      if (md_est.solver_used == 1 || md_est.solver_used == 2) {
        v = covariance_val * s2 - covariance_val * R * covariance_val * s2;
      } else {
        v = LU_solve_covariance * s2;
      }
      laplace_likelihood::diff_eta_implicit(ll_fun, v, md_est.theta,
                                            ll_args_copy, msgs);
      internal::collect_adjoints<ZeroOut>(partial_parm, ll_args_filter);
    }
    lmd = md_est.lmd;
  }
  var ret(lmd);
  if constexpr (is_any_var_scalar_v<CovarArgs>) {
    auto covar_args_filter = internal::filter_var_scalar_types(covar_args_refs);
    internal::reverse_pass_collect_adjoints(ret, covar_args_filter,
                                            std::move(covar_args_adj));
  }
  if constexpr (ll_args_contain_var) {
    auto ll_args_filter = internal::filter_var_scalar_types(ll_args_refs);
    internal::reverse_pass_collect_adjoints(ret, ll_args_filter,
                                            std::move(partial_parm));
  }
  return ret;
}

}  // namespace math
}  // namespace stan

#endif
