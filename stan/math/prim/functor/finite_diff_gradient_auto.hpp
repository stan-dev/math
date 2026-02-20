#ifndef STAN_MATH_PRIM_FUNCTOR_FINITE_DIFF_GRADIENT_AUTO_HPP
#define STAN_MATH_PRIM_FUNCTOR_FINITE_DIFF_GRADIENT_AUTO_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/finite_diff_stepsize.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <stan/math/prim/functor/apply.hpp>
#include <stan/math/prim/functor/filter_map.hpp>
#include <stan/math/prim/functor/for_each.hpp>
#include <stan/math/prim/functor/make_holder_tuple.hpp>
#include <stan/math/prim/meta/contains_autodiff.hpp>
#include <stan/math/prim/meta/filtered_tuple_indices.hpp>
#include <cmath>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace stan {
namespace math {
namespace internal {

/**
 * Construct a holder tuple containing elements selected by `Is...`.
 *
 * Returned elements preserve reference/value category through
 * `make_holder_tuple`, so this helper can be used with tuples of references.
 *
 * @tparam Tuple tuple-like input type
 * @tparam Is selected index pack
 * @param[in] tuple input tuple
 * @param[in] selected_idxs index sequence selecting tuple entries
 * @return holder tuple containing selected entries
 * @throw None.
 */
template <typename Tuple, std::size_t... Is>
inline auto tuple_subset(Tuple&& tuple,
                         [[maybe_unused]] std::index_sequence<Is...> selected_idxs) {
  return make_holder_tuple(std::get<Is>(std::forward<Tuple>(tuple))...);
}

/**
 * Extract a double-valued coordinate used for finite-difference step size
 * computation from scalar-like autodiff stacks.
 *
 * @tparam T scalar type
 * @param[in] x scalar value
 * @return value used to compute finite-difference step size
 * @throw None.
 */
template <typename T>
inline double finite_diff_stepsize_value(const T& x) {
  if constexpr (is_var<std::decay_t<T>>::value) {
    return x.val();
  } else if constexpr (is_fvar<std::decay_t<T>>::value) {
    return finite_diff_stepsize_value(x.val_);
  } else {
    return value_of_rec(x);
  }
}


/**
 * Recursively traverses `arg_mask`, `arg_work`, and `grad` and invokes `fn` on
 * each perturbable coordinate pair `(arg_work_coord, grad_coord)`.
 *
 * `arg_mask` controls whether a scalar coordinate is perturbable, while
 * `arg_work` is the mutable coordinate that gets perturbed.
 *
 * Supported containers are:
 * - Stan scalars
 * - Eigen dense objects
 * - `std::vector` (possibly nested)
 * - `std::tuple` (including tuple-of-tuples)
 *
 * The `arg_mask` and `arg_work` structures must match exactly. The `grad`
 * structure may either:
 * - match top-level tuple arity with `arg_mask`, or
 * - be compacted to only autodiff-containing top-level tuple entries.
 *
 * @tparam ArgMask traversal mask type
 * @tparam ArgWork mutable argument type
 * @tparam Grad gradient output type
 * @tparam F callable compatible with `(arg_work_coord, grad_coord)`
 * @param[in] arg_mask mask selecting perturbable coordinates
 * @param[in,out] arg_work mutable argument coordinates
 * @param[out] grad gradient output coordinates
 * @param[in] fn callback applied for each perturbable coordinate
 * @return None.
 * @throw None. Fails compilation with `static_assert` for unsupported
 * container categories or mismatched tuple arity.
 */
template <typename ArgMask, typename ArgWork, typename Grad, typename F>
inline void for_each_leaf_pair(ArgMask&& arg_mask, ArgWork&& arg_work,
                               Grad&& grad, F&& fn) {
  using arg_mask_t = std::decay_t<ArgMask>;
  if constexpr (is_stan_scalar_v<arg_mask_t>) {
    if constexpr (contains_autodiff_v<arg_mask_t>) {
      fn(arg_work, grad);
    }
  } else if constexpr (is_eigen_v<arg_mask_t>) {
    for (Eigen::Index i = 0; i < arg_mask.size(); ++i) {
      for_each_leaf_pair(arg_mask.coeffRef(i), arg_work.coeffRef(i),
                         grad.coeffRef(i), fn);
    }
  } else if constexpr (is_std_vector_v<arg_mask_t>) {
    for (std::size_t i = 0; i < arg_mask.size(); ++i) {
      for_each_leaf_pair(arg_mask[i], arg_work[i], grad[i], fn);
    }
  } else if constexpr (is_tuple_v<arg_mask_t>) {
    constexpr auto arg_mask_size = std::tuple_size<std::decay_t<ArgMask>>::value;
    constexpr auto arg_work_size = std::tuple_size<std::decay_t<ArgWork>>::value;
    static_assert(arg_mask_size == arg_work_size,
                  "Tuple size mismatch in finite_diff_gradient_auto traversal "
                  "between arg_mask and arg_work.");
    constexpr auto grad_size = std::tuple_size<std::decay_t<Grad>>::value;
    if constexpr (arg_mask_size == grad_size) {
      stan::math::for_each(
          [&fn](auto&& arg_mask_i, auto&& arg_work_i, auto&& grad_i) {
            for_each_leaf_pair(arg_mask_i, arg_work_i, grad_i, fn);
          },
          std::forward<ArgMask>(arg_mask), std::forward<ArgWork>(arg_work),
          std::forward<Grad>(grad));
    } else {
      using selected_idxs_t
          = filtered_tuple_indices_t<contains_autodiff, std::decay_t<ArgMask>>;
      auto arg_mask_ad
          = tuple_subset(std::forward<ArgMask>(arg_mask), selected_idxs_t{});
      auto arg_work_ad
          = tuple_subset(std::forward<ArgWork>(arg_work), selected_idxs_t{});
      static_assert(
          std::tuple_size<std::decay_t<decltype(arg_mask_ad)>>::value
                  == std::tuple_size<std::decay_t<Grad>>::value
              && std::tuple_size<std::decay_t<decltype(arg_work_ad)>>::value
                     == std::tuple_size<std::decay_t<Grad>>::value,
          "Tuple size mismatch in finite_diff_gradient_auto traversal");
      stan::math::for_each(
          [&fn](auto&& arg_mask_i, auto&& arg_work_i, auto&& grad_i) {
            for_each_leaf_pair(arg_mask_i, arg_work_i, grad_i, fn);
          },
          arg_mask_ad, arg_work_ad, std::forward<Grad>(grad));
    }
  } else {
    static_assert(sizeof(std::decay_t<arg_mask_t>*) == 0,
                  "Unsupported container in finite_diff_gradient_auto.");
  }
}

/**
 * Apply the six-point finite-difference stencil to every perturbable
 * coordinate in one top-level argument and write gradient coordinates.
 *
 * Coordinates in `arg_work` are restored to their original value after each
 * stencil evaluation.
 *
 * @tparam ScalarT scalar type of function value and gradients
 * @tparam F callable type
 * @tparam ArgMask mask argument type
 * @tparam ArgWork mutable working argument type
 * @tparam NthGrad gradient container type for this top-level argument
 * @tparam TupleArgsWork mutable tuple of all working arguments passed to `f`
 * @param[in] f function being finite-differenced
 * @param[in] arg_mask mask argument used to choose perturbable coordinates
 * @param[in,out] arg_work mutable working argument coordinates
 * @param[out] grad gradient output container for this top-level argument
 * @param[in,out] args_work mutable tuple of all working arguments passed to
 * `f`
 * @return None.
 * @throw Any exception thrown by `f`.
 */
template <typename ScalarT, typename F, typename ArgMask, typename ArgWork,
          typename NthGrad, typename TupleArgsWork>
inline void finite_diff_gradient_auto_impl(const F& f, ArgMask&& arg_mask,
                                           ArgWork&& arg_work, NthGrad&& grad,
                                           TupleArgsWork&& args_work) {
  static constexpr int h_scale[6] = {3, 2, 1, -3, -2, -1};
  static constexpr int mults[6] = {1, -9, 45, -1, 9, -45};

  for_each_leaf_pair(
      std::forward<ArgMask>(arg_mask), std::forward<ArgWork>(arg_work),
      std::forward<NthGrad>(grad), [&f, &args_work](auto& arg_coord, auto& grad_coord) {
        const auto orig = arg_coord;
        const double h = finite_diff_stepsize(finite_diff_stepsize_value(orig));
        ScalarT delta_f = 0;
        for (int j = 0; j < 6; ++j) {
          arg_coord = orig + h * h_scale[j];
          delta_f += stan::math::apply(
                         [&f](auto&&... args_i) { return f(args_i...); },
                         args_work)
                     * mults[j];
        }
        arg_coord = orig;
        grad_coord = delta_f / (60 * h);
      });
}

}  // namespace internal

/**
 * Calculate the value and the gradient of the specified function
 * at the specified argument using finite difference.
 *
 * <p>The functor must implement
 *
 * <code>
 * double operator()(const Eigen::Matrix<double, -, 1>&) const;
 * </code>
 *
 * <p>Error of derivative in dimension `i` should be on the should be on
 * order of `epsilon(i)^6`, where `epsilon(i) = sqrt(delta) * abs(x(i))`
 * for input `x` at dimension `i`.
 *
 * The reference for this algorithm is:
 *
 * <br />Robert de Levie. 2009. An improved numerical approximation
 * for the first derivative. Journal of Chemical Sciences 121(5), page
 * 3.
 *
 * <p>The reference for automatically setting the difference is this
 * section of the Wikipedia,
 *
 * <br /><a
 * href="https://en.wikipedia.org/wiki/Numerical_differentiation#Practical_considerations_using_floating-point_arithmetic">Numerical
 * differentiation: practical considerations using floating point
 * arithmetic</a>.
 *
 * <p>Evaluating this function involves 6 calls to the function being
 * differentiated for each dimension in the input, plus one global
 * evaluation.  All evaluations will be for double-precision inputs.
 *
 * @tparam F Type of function
 * @param[in] f function
 * @param[in] x argument to function
 * @param[out] fx function applied to argument
 * @param[out] grad_fx gradient of function at argument
 * @return None.
 * @throw Any exception thrown by `f`.
 */
template <typename F, typename VectorT, typename GradVectorT,
          typename ScalarT = return_type_t<VectorT>,
          require_vector_t<VectorT>* = nullptr>
inline void finite_diff_gradient_auto(const F& f, VectorT&& x, ScalarT& fx,
                                      GradVectorT& grad_fx) {
  using EigT = Eigen::Matrix<ScalarT, -1, 1>;
  static constexpr int h_scale[6] = {3, 2, 1, -3, -2, -1};
  static constexpr int mults[6] = {1, -9, 45, -1, 9, -45};

  fx = f(x);
  grad_fx.resize(x.size());
  Eigen::Map<EigT> grad_map(grad_fx.data(), grad_fx.size());

  grad_map = EigT::NullaryExpr(x.size(), [&f, &x](Eigen::Index i) {
    double h = finite_diff_stepsize(value_of_rec(x[i]));
    ScalarT delta_f = 0;
    for (int j = 0; j < 6; ++j) {
      auto x_temp = EigT::NullaryExpr(x.size(), [&x, &i, &h, &j](Eigen::Index k) {
        return k == i ? x[i] + h * h_scale[j] : x[k];
      });
      delta_f += f(std::move(x_temp)) * mults[j];
    }
    return delta_f / (60 * h);
  });
}

/**
 * Calculate the function value and gradients for heterogeneous tuple arguments
 * while preserving argument scalar types at the call boundary.
 *
 * This overload expects `grads` to mirror the top-level structure of `args`.
 * Non-autodiff leaves are skipped in the finite-difference pass.
 *
 * Pattern used by this overload:
 * 1. Evaluate `f(args...)` once for `fx`.
 * 2. Zip top-level `args` and `grads`.
 * 3. Traverse nested leaves with `for_each_leaf_pair`, perturbing one
 *    coordinate at a time and restoring it after each stencil evaluation.
 *
 * @tparam F callable type
 * @tparam ScalarT scalar type used for function value and gradient storage
 * @tparam TupleArgs tuple of input arguments
 * @tparam TupleGrads tuple of gradient containers
 * @param[in] f function to differentiate
 * @param[out] fx function value at `args`
 * @param[in] args tuple of function arguments
 * @param[out] grads tuple of gradient containers
 * @return None.
 * @throw Any exception thrown by `f`.
 */
template <typename F, typename ScalarT, typename TupleArgs, typename TupleGrads,
          require_tuple_t<TupleArgs>* = nullptr,
          require_tuple_t<TupleGrads>* = nullptr>
inline void finite_diff_gradient_auto(const F& f, ScalarT& fx, TupleArgs&& args,
                                      TupleGrads&& grads) {
  fx = stan::math::apply([&f](auto&&... args_i) { return f(args_i...); }, args);
  stan::math::for_each(
      [&f, &args](auto&& arg_i, auto&& grad_i) {
        if constexpr (!std::is_integral_v<scalar_type_t<decltype(arg_i)>>) {
          internal::finite_diff_gradient_auto_impl<ScalarT>(f, arg_i, arg_i,
                                                            grad_i, args);
        }
      },
      args, grads);
}

/**
 * Tuple overload used for compact gradient tuples with an explicit mask tuple.
 *
 * This overload is used by `stan::math::finite_diff` where:
 * - `args_mask` preserves original autodiff type information
 * - `args_work` is the value-only work tuple passed to function evaluations
 * - `grads` stores only autodiff top-level gradient containers
 *
 * @tparam F callable type
 * @tparam ScalarT scalar type used for function value and gradients
 * @tparam TupleArgsMask tuple used only as perturbation mask
 * @tparam TupleArgsWork mutable tuple used for function evaluation
 * @tparam TupleGrads compact tuple of gradient containers
 * @param[in] f function to differentiate
 * @param[out] fx function value at `args_work`
 * @param[in] args_mask tuple carrying autodiff mask information
 * @param[in,out] args_work mutable value tuple used for finite-difference
 * evaluations
 * @param[out] grads compact tuple of gradient containers
 * @return None.
 * @throw Any exception thrown by `f`.
 */
template <typename F, typename ScalarT, typename TupleArgsMask,
          typename TupleArgsWork, typename TupleGrads,
          require_tuple_t<TupleArgsMask>* = nullptr,
          require_tuple_t<TupleArgsWork>* = nullptr,
          require_tuple_t<TupleGrads>* = nullptr>
inline void finite_diff_gradient_auto(
    const F& f, ScalarT& fx, TupleArgsMask&& args_mask, TupleArgsWork&& args_work,
    TupleGrads&& grads) {
  fx = stan::math::apply([&f](auto&&... args_i) { return f(args_i...); },
                         args_work);
  using selected_idxs_t
      = filtered_tuple_indices_t<contains_autodiff, std::decay_t<TupleArgsMask>>;
  auto args_mask_ad = internal::tuple_subset(args_mask, selected_idxs_t{});
  auto args_work_ad = internal::tuple_subset(args_work, selected_idxs_t{});
  static_assert(std::tuple_size<std::decay_t<decltype(args_mask_ad)>>::value
                        == std::tuple_size<std::decay_t<decltype(args_work_ad)>>::value
                    && std::tuple_size<std::decay_t<decltype(args_mask_ad)>>::value
                           == std::tuple_size<std::decay_t<TupleGrads>>::value,
                "Tuple size mismatch in finite_diff_gradient_auto.");
  stan::math::for_each(
      [&f, &args_work](auto&& arg_mask_i, auto&& arg_work_i, auto&& grad_i) {
        internal::finite_diff_gradient_auto_impl<ScalarT>(
            f, std::forward<decltype(arg_mask_i)>(arg_mask_i),
            std::forward<decltype(arg_work_i)>(arg_work_i),
            std::forward<decltype(grad_i)>(grad_i), args_work);
      },
      args_mask_ad, args_work_ad, std::forward<TupleGrads>(grads));
}

}  // namespace math
}  // namespace stan
#endif
