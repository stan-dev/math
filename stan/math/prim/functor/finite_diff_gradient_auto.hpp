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
 * Type-dependent false helper for `static_assert` in unreachable branches.
 */
template <typename...>
struct dependent_false : std::false_type {};

/**
 * Filter a tuple and return a tuple containing references to the elements
 * whose types satisfy `contains_autodiff`.
 *
 * This mirrors the `filter_var_scalar_types` pattern in rev.
 *
 * @tparam T input tuple or nested tuple-like structure
 * @param t input to filter
 * @return tuple of references to autodiff-containing entries
 */
template <typename T>
inline constexpr decltype(auto) filter_autodiff_types(T&& t) {
  return stan::math::filter_map<contains_autodiff>(
      [](auto&& arg) -> decltype(auto) {
        using arg_t = std::decay_t<decltype(arg)>;
        if constexpr (is_tuple_v<arg_t>) {
          return filter_autodiff_types(std::forward<decltype(arg)>(arg));
        } else {
          return std::forward<decltype(arg)>(arg);
        }
      },
      std::forward<T>(t));
}

template <typename Tuple, std::size_t... Is>
inline auto tuple_subset(Tuple&& tuple, std::index_sequence<Is...>) {
  return make_holder_tuple(std::get<Is>(std::forward<Tuple>(tuple))...);
}

/**
 * Extract a double-valued coordinate used for finite-difference step size
 * computation from scalar-like autodiff stacks.
 *
 * @tparam T scalar type
 * @param x scalar value
 * @return value used to compute finite-difference step size
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
 * @tparam ArgMask traversal mask type
 * @tparam ArgWork mutable argument type
 * @tparam Grad gradient output type
 * @tparam F callable compatible with `(arg_work_coord, grad_coord)`
 * @param arg_mask mask selecting perturbable coordinates
 * @param arg_work mutable argument coordinates
 * @param grad gradient output coordinates
 * @param fn callback applied for each perturbable coordinate
 */
template <typename ArgMask, typename ArgWork, typename Grad, typename F>
inline void for_each_coordinate_pair_mut(ArgMask&& arg_mask, ArgWork&& arg_work,
                                         Grad&& grad, F&& fn) {
  using arg_mask_t = std::decay_t<ArgMask>;
  if constexpr (is_stan_scalar_v<arg_mask_t>) {
    if constexpr (contains_autodiff_v<arg_mask_t>) {
      fn(arg_work, grad);
    }
  } else if constexpr (is_eigen_v<arg_mask_t>) {
    for (Eigen::Index i = 0; i < arg_mask.size(); ++i) {
      for_each_coordinate_pair_mut(arg_mask.coeffRef(i), arg_work.coeffRef(i),
                                   grad.coeffRef(i), fn);
    }
  } else if constexpr (is_std_vector_v<arg_mask_t>) {
    for (std::size_t i = 0; i < arg_mask.size(); ++i) {
      for_each_coordinate_pair_mut(arg_mask[i], arg_work[i], grad[i], fn);
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
            for_each_coordinate_pair_mut(arg_mask_i, arg_work_i, grad_i, fn);
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
            for_each_coordinate_pair_mut(arg_mask_i, arg_work_i, grad_i, fn);
          },
          std::move(arg_mask_ad), std::move(arg_work_ad),
          std::forward<Grad>(grad));
    }
  } else {
    static_assert(dependent_false<arg_mask_t>::value,
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
 * @param f function being finite-differenced
 * @param arg_mask mask argument used to choose perturbable coordinates
 * @param arg_work mutable working argument coordinates
 * @param grad gradient output container for this top-level argument
 * @param args_work mutable tuple of all working arguments passed to `f`
 */
template <typename ScalarT, typename F, typename ArgMask, typename ArgWork,
          typename NthGrad, typename TupleArgsWork>
inline void finite_diff_gradient_auto_impl(const F& f, ArgMask&& arg_mask,
                                           ArgWork&& arg_work, NthGrad&& grad,
                                           TupleArgsWork&& args_work) {
  static constexpr int h_scale[6] = {3, 2, 1, -3, -2, -1};
  static constexpr int mults[6] = {1, -9, 45, -1, 9, -45};

  for_each_coordinate_pair_mut(
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
 * Filter a tuple and return a tuple with references to entries with autodiff
 * scalar type.
 *
 * @tparam T Possibly a tuple, std::vector, Eigen type, or scalar
 * @param[in] t Input to filter
 * @return Filtered input with only autodiff scalar entries
 */
template <typename T>
inline constexpr decltype(auto) filter_ad_scalar_types(T&& t) {
  return internal::filter_autodiff_types(std::forward<T>(t));
}

/**
 * Calculate the function value and gradients for heterogeneous tuple arguments
 * while preserving argument scalar types at the call boundary.
 *
 * This overload expects `grads` to mirror the top-level structure of `args`.
 * Non-autodiff top-level entries are skipped.
 *
 * @tparam F callable type
 * @tparam ScalarT scalar type used for function value and gradient storage
 * @tparam TupleArgs tuple of input arguments
 * @tparam TupleGrads tuple of gradient containers
 * @param[in] f function to differentiate
 * @param[out] fx function value at `args`
 * @param[in] args tuple of function arguments
 * @param[out] grads tuple of gradient containers
 */
template <typename F, typename ScalarT, typename TupleArgs, typename TupleGrads,
          require_tuple_t<TupleArgs>* = nullptr,
          require_tuple_t<TupleGrads>* = nullptr>
inline void finite_diff_gradient_auto(const F& f, ScalarT& fx, TupleArgs&& args,
                                      TupleGrads&& grads) {
  fx = stan::math::apply([&f](auto&&... args_i) { return f(args_i...); }, args);
  stan::math::for_each(
      [&f, &args](auto&& arg_i, auto&& grad_i) {
        if constexpr (contains_autodiff_v<std::decay_t<decltype(arg_i)>>) {
          internal::finite_diff_gradient_auto_impl<ScalarT>(f, arg_i, arg_i,
                                                            grad_i, args);
        }
      },
      args, grads);
}

/**
 * Indexed tuple overload used for compact gradient tuples. `Idxs` select
 * autodiff-containing top-level entries from both `args_mask` and `args_work`,
 * while `GradIdxs` selects matching gradient entries from `grads`.
 *
 * @tparam F callable type
 * @tparam ScalarT scalar type used for function value and gradients
 * @tparam TupleArgsMask tuple used only as perturbation mask
 * @tparam TupleArgsWork mutable tuple used for function evaluation
 * @tparam TupleGrads compact tuple of gradient containers
 */
template <typename F, typename ScalarT, typename TupleArgsMask,
          typename TupleArgsWork, typename TupleGrads, std::size_t... Idxs,
          std::size_t... GradIdxs, require_tuple_t<TupleArgsMask>* = nullptr,
          require_tuple_t<TupleArgsWork>* = nullptr,
          require_tuple_t<TupleGrads>* = nullptr>
inline void finite_diff_gradient_auto(
    const F& f, ScalarT& fx, TupleArgsMask&& args_mask, TupleArgsWork&& args_work,
    TupleGrads&& grads, std::index_sequence<Idxs...> /* idxs */,
    std::index_sequence<GradIdxs...> /* grad_idxs */) {
  static_assert(sizeof...(Idxs) == sizeof...(GradIdxs),
                "Index sequence size mismatch in finite_diff_gradient_auto.");
  fx = stan::math::apply([&f](auto&&... args_i) { return f(args_i...); },
                         args_work);
  using Swallow = int[];
  static_cast<void>(Swallow{
      (static_cast<void>(internal::finite_diff_gradient_auto_impl<ScalarT>(
           f, std::get<Idxs>(args_mask), std::get<Idxs>(args_work),
           std::get<GradIdxs>(grads), args_work)),
       0)...});
}

}  // namespace math
}  // namespace stan
#endif
