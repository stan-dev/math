#ifndef STAN_MATH_PRIM_FUNCTOR_FINITE_DIFF_GRADIENT_AUTO_HPP
#define STAN_MATH_PRIM_FUNCTOR_FINITE_DIFF_GRADIENT_AUTO_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/finite_diff_stepsize.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <stan/math/prim/meta/conditional_copy_and_promote.hpp>
#include <stan/math/prim/functor/apply.hpp>
#include <stan/math/prim/functor/filter_map.hpp>
#include <stan/math/prim/functor/for_each.hpp>
#include <stan/math/prim/functor/make_holder_tuple.hpp>
#include <cmath>
#include <functional>
#include <tuple>
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
 * Trait indicating whether a type contains at least one autodiff scalar.
 *
 * Scalars are checked directly, while container specializations recurse into
 * contained value types.
 */
template <typename T, typename = void>
struct contains_autodiff : bool_constant<is_autodiff_v<std::decay_t<T>>> {};

template <typename... Args>
struct contains_autodiff<std::tuple<Args...>, void>
    : bool_constant<(contains_autodiff<std::decay_t<Args>>::value || ...)> {};

template <typename T, typename... VecArgs>
struct contains_autodiff<std::vector<T, VecArgs...>, void>
    : bool_constant<contains_autodiff<std::decay_t<T>>::value> {};

template <typename T>
inline constexpr bool contains_autodiff_v
    = contains_autodiff<std::decay_t<T>>::value;

template <typename T>
struct contains_autodiff_decayed : contains_autodiff<std::decay_t<T>> {};

template <typename T, typename = void>
struct count_autodiff_args;

template <typename... Args>
struct count_autodiff_args<std::tuple<Args...>, void>
    : std::integral_constant<
          std::size_t,
          (static_cast<std::size_t>(contains_autodiff_v<Args>) + ... + 0)> {};

template <typename TupleArgs>
inline constexpr std::size_t count_autodiff_args_v
    = count_autodiff_args<std::decay_t<TupleArgs>>::value;

template <typename T, template <class...> class Cond, std::size_t pos>
struct indices_of_if_impl {
  static constexpr std::size_t value = (Cond<T>::value ? pos : -1);
};

/**
 * Wrapper that carries a top-level argument by reference and exposes a
 * potentially different mask type for filtering.
 *
 * @tparam MaskT type used for filter decisions
 * @tparam ArgT referenced argument type
 */
template <typename MaskT, typename ArgT>
struct masked_arg_ref {
  using arg_t = std::remove_reference_t<ArgT>;
  std::reference_wrapper<arg_t> arg_;

  explicit masked_arg_ref(arg_t& arg) : arg_(arg) {}
  inline arg_t& get() const { return arg_.get(); }
};

/**
 * Specialization so `filter_map<contains_autodiff>` can filter wrapped
 * top-level arguments based on the explicit mask type.
 *
 * @tparam MaskT mask type used for autodiff detection
 * @tparam ArgT referenced argument type
 */
template <typename MaskT, typename ArgT>
struct contains_autodiff<masked_arg_ref<MaskT, ArgT>, void>
    : contains_autodiff<std::decay_t<MaskT>> {};

/**
 * Trait to detect wrapped top-level arguments.
 */
template <typename T>
struct is_masked_arg_ref : std::false_type {};

template <typename MaskT, typename ArgT>
struct is_masked_arg_ref<masked_arg_ref<MaskT, ArgT>> : std::true_type {};

/**
 * Internal implementation for creating a tuple of masked top-level argument
 * wrappers.
 *
 * @tparam TupleMask tuple type used for autodiff mask decisions
 * @tparam TupleArgs tuple instance to filter
 * @tparam Is compile-time tuple indices
 * @param args tuple instance to filter
 * @return tuple of wrapped top-level arguments carrying explicit mask types
 */
template <typename TupleMask, typename TupleArgs, std::size_t... Is>
inline auto make_masked_arg_ref_tuple_impl(
    TupleArgs&& args, std::index_sequence<Is...>) {
  return make_holder_tuple(
      masked_arg_ref<std::tuple_element_t<Is, std::decay_t<TupleMask>>,
                     decltype(std::get<Is>(std::forward<TupleArgs>(args)))>(
          std::get<Is>(std::forward<TupleArgs>(args)))...);
}

/**
 * Create a tuple of masked top-level argument wrappers.
 *
 * @tparam TupleMask tuple type used for autodiff mask decisions
 * @tparam TupleArgs tuple instance to wrap
 * @param args tuple instance to wrap
 * @return tuple of wrapped top-level arguments carrying explicit mask types
 */
template <typename TupleMask, typename TupleArgs>
inline auto make_masked_arg_ref_tuple(TupleArgs&& args) {
  return make_masked_arg_ref_tuple_impl<TupleMask>(
      std::forward<TupleArgs>(args),
      std::make_index_sequence<std::tuple_size<std::decay_t<TupleMask>>::value>{
      });
}

/**
 * Filter a tuple and return a tuple containing references to the elements
 * whose types satisfy `contains_autodiff`.
 *
 * This mirrors the `filter_var_scalar_types` pattern in rev.
 *
 * @tparam T input tuple or nested tuple-like structure
 * @param t input to filter
 * @return tuple of references to autodiff-containing top-level arguments
 */
template <typename T>
inline constexpr decltype(auto) filter_autodiff_types(T&& t) {
  return stan::math::filter_map<contains_autodiff>(
      [](auto&& arg) -> decltype(auto) {
        using arg_t = std::decay_t<decltype(arg)>;
        if constexpr (is_tuple_v<arg_t>) {
          return filter_autodiff_types(std::forward<decltype(arg)>(arg));
        } else if constexpr (is_masked_arg_ref<arg_t>::value) {
          return std::forward<decltype(arg)>(arg).get();
        } else {
          return std::forward<decltype(arg)>(arg);
        }
      },
      std::forward<T>(t));
}

/**
 * Filter top-level arguments using an explicit tuple mask.
 *
 * @tparam TupleMask tuple type used for autodiff mask decisions
 * @tparam TupleArgs tuple instance to filter
 * @param args tuple instance to filter
 * @return tuple of references to autodiff-containing top-level arguments
 */
template <typename TupleMask, typename TupleArgs>
inline constexpr decltype(auto) filter_autodiff_types_with_mask(
    TupleArgs&& args) {
  return filter_autodiff_types(
      make_masked_arg_ref_tuple<TupleMask>(std::forward<TupleArgs>(args)));
}

/**
 * Builds a zero-initialized gradient container matching `arg` shape for
 * autodiff leaves, while preserving non-autodiff leaves by value.
 *
 * @tparam ScalarT scalar type
 * @tparam Arg argument/container type used as the shape template
 * @param arg shape template argument
 * @return gradient container aligned with `arg`
 */
template <typename Arg>
inline auto zeroed_container(const Arg& arg) {
  using ArgDec = std::decay_t<Arg>;
  using PartialType = partials_type_t<scalar_type_t<ArgDec>>;
  if constexpr (is_tuple_v<ArgDec>) {
    return stan::math::apply(
        [](const auto&... inner_args) {
          return std::make_tuple(zeroed_container(inner_args)...);
        }, arg);
  } else if constexpr (is_std_vector_v<ArgDec>) {
    std::vector<PartialType> ret;
    ret.reserve(arg.size());
    for (const auto& el : arg) {
      ret.push_back(zeroed_container(el));
    }
    return ret;
  } else if constexpr (is_eigen_v<ArgDec>) {
      using eigen_t = promote_scalar_t<PartialType, plain_type_t<ArgDec>>;
      return eigen_t::Zero(arg.rows(), arg.cols()).eval();
  } else if constexpr (is_stan_scalar_v<ArgDec>) {
      return PartialType(0);
  } else {
    static_assert(
        dependent_false<ArgDec>::value,
        "Unsupported gradient container in finite_diff_gradient_auto.");
  }
}

/**
 * Extracts a double-valued coordinate for finite-difference step-size
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
 * Recursively traverses `arg_mask`, `arg_work`, and `grad` with matching
 * structure and invokes `fn` on each autodiff scalar coordinate pair
 * `(arg_work_coord, grad_coord)`.
 *
 * Non-autodiff scalar leaves in `arg_mask` are skipped.
 *
 * @tparam ArgMask mask argument type (uses original scalar categories)
 * @tparam ArgWork mutable working argument type
 * @tparam Grad gradient container type
 * @tparam F callable with signature compatible with `(arg_work_coord, grad_coord)`
 * @param arg_mask mask argument used to decide perturbable coordinates
 * @param arg_work mutable working argument coordinates
 * @param grad gradient output coordinates
 * @param fn callback applied at each perturbable coordinate
 */
template <typename ArgMask, typename Grad, typename F>
inline void for_each_coordinate_pair_mut(ArgMask&& arg_mask, Grad&& grad, F&& fn) {
  using arg_mask_t = std::decay_t<ArgMask>;
  if constexpr (is_stan_scalar_v<arg_mask_t>) {
    if constexpr (is_autodiff_v<arg_mask_t>) {
      fn(arg_mask, grad);
    }
  } else if constexpr (is_eigen_v<arg_mask_t>) {
    for (Eigen::Index i = 0; i < arg_mask.size(); ++i) {
      for_each_coordinate_pair_mut(arg_mask.coeffRef(i),
                                   grad.coeffRef(i), fn);
    }
  } else if constexpr (is_std_vector_v<arg_mask_t>) {
    for (std::size_t i = 0; i < arg_mask.size(); ++i) {
      for_each_coordinate_pair_mut(arg_mask[i], grad[i], fn);
    }
  } else if constexpr (is_tuple_v<arg_mask_t>) {
    static_assert(std::tuple_size<std::decay_t<ArgMask>>::value
                      == std::tuple_size<std::decay_t<Grad>>::value,
                  "Tuple size mismatch in finite_diff_gradient_auto traversal");
    stan::math::for_each(
        [&fn](auto&& arg_mask_i, auto&& grad_i) {
          for_each_coordinate_pair_mut(arg_mask_i, grad_i, fn);
        },
        std::forward<ArgMask>(arg_mask),
        std::forward<Grad>(grad));
  } else {
    static_assert(dependent_false<arg_mask_t>::value,
                  "Unsupported container in finite_diff_gradient_auto.");
  }
}

/**
 * Applies the six-point finite-difference stencil to every perturbable
 * coordinate in one autodiff top-level argument and writes the corresponding
 * gradient coordinate.
 *
 * Coordinates are always restored to their original value after each stencil.
 *
 * @tparam ScalarT scalar type of function value and gradients
 * @tparam F callable type
 * @tparam ArgMask mask argument type (original categories)
 * @tparam ArgWork mutable working argument type
 * @tparam NthGrad gradient container type for this top-level argument
 * @tparam TupleArgsWork mutable tuple of all working arguments
 * @param f function being finite-differenced
 * @param arg_mask mask argument used to choose perturbable coordinates
 * @param arg_work mutable working argument coordinates
 * @param grad gradient output container for this top-level argument
 * @param args_work mutable tuple of all working arguments passed to `f`
 */
template <typename ScalarT, typename F, typename ArgMask,
          typename NthGrad, typename TupleArgsWork>
inline void finite_diff_gradient_auto_impl(const F& f, ArgMask&& arg_mask,
                                           NthGrad&& grad,
                                           TupleArgsWork&& args_work) {
  static constexpr int h_scale[6] = {3, 2, 1, -3, -2, -1};
  static constexpr int mults[6] = {1, -9, 45, -1, 9, -45};

  for_each_coordinate_pair_mut(
      std::forward<ArgMask>(arg_mask), std::forward<NthGrad>(grad),
      [&f, &args_work](auto& arg_coord, auto& grad_coord) {
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
      auto x_temp
          = EigT::NullaryExpr(x.size(), [&x, &i, &h, &j](Eigen::Index k) {
              return k == i ? x[i] + h * h_scale[j] : x[k];
            });
      delta_f += f(std::move(x_temp)) * mults[j];
    }
    return delta_f / (60 * h);
  });
}

/**
 * Filter a tuple and return a tuple with references to the types with an autodiff
 * scalar type.
 * @tparam T Possibly a tuple, std::vector, Eigen type, or scalar
 * @param[in] t Input to filter
 * @return Filtered input with only var scalar types
 */
template <typename T>
inline constexpr decltype(auto) filter_ad_scalar_types(T&& t) {
  return stan::math::filter_map<internal::contains_autodiff>(
      [](auto&& arg) -> decltype(auto) {
        if constexpr (is_tuple_v<std::decay_t<decltype(arg)>>) {
          return filter_ad_scalar_types(std::forward<decltype(arg)>(arg));
        } else {
          return std::forward<decltype(arg)>(arg);
        }
      },
      std::forward<T>(t));
}

/**
 * Calculate the function value and gradients for heterogeneous tuple arguments
 * while preserving non-autodiff argument scalar types at the call boundary.
 *
 * Pattern used by this overload:
 * 1. Build `args_ad` by filtering top-level autodiff-containing entries from
 *    the original `args` tuple.
 * 2. Build one mutable working tuple `args_work` where autodiff-containing
 *    entries are deep-copied and promoted to `ScalarT`, and non-autodiff
 *    entries are forwarded unchanged.
 * 3. Build `args_work_ad` by applying the same top-level filter shape as
 *    `args_ad` (via the `TupleArgs` mask) so both filtered tuples align.
 * 4. Zip `args_ad`, `args_work_ad`, and `grads` and run
 *    `finite_diff_gradient_auto_impl` for each top-level autodiff entry.
 *
 * The split between original mask (`args_ad`) and mutable working values
 * (`args_work_ad`) ensures coordinate traversal only perturbs locations that
 * correspond to autodiff coordinates in the original argument pack, while all
 * function evaluations run against the promoted mutable working tuple.
 *
 * @tparam F callable type
 * @tparam ScalarT scalar type used for function value and gradient storage
 * @tparam TupleArgs tuple of input arguments
 * @tparam TupleGrads tuple of gradient containers aligned to autodiff
 *         top-level arguments
 * @param[in] f function to differentiate
 * @param[out] fx function value at `args`
 * @param[in] args tuple of function arguments
 * @param[out] grads compact tuple of gradient containers
 */
template <typename F, typename ScalarT, typename TupleArgs, typename TupleGrads,
          require_tuple_t<TupleArgs>* = nullptr,
          require_tuple_t<TupleGrads>* = nullptr,
          require_t<bool_constant<
              (internal::count_autodiff_args_v<TupleArgs> > 0)>>* = nullptr>
inline void finite_diff_gradient_auto(const F& f, ScalarT& fx, TupleArgs&& args,
                                      TupleGrads&& grads) {
  fx = stan::math::apply([&f](auto&&... args_i) { return f(args_i...); },
                         args);
  stan::math::for_each([&f, &args, &grads](auto&& arg, auto&& grad) {
    if (!std::is_integral_v<scalar_type_t<decltype(grad)>> && !std::is_integral_v<scalar_type_t<decltype(arg)>>) {
      internal::finite_diff_gradient_auto_impl<ScalarT>(
          f, arg, grad, args);
    }
  }, args, grads);
 }


/**
 * Calculate the function value and gradients for heterogeneous tuple arguments
 * while preserving non-autodiff argument scalar types at the call boundary.
 *
 * Pattern used by this overload:
 * 1. Build `args_ad` by filtering top-level autodiff-containing entries from
 *    the original `args` tuple.
 * 2. Build one mutable working tuple `args_work` where autodiff-containing
 *    entries are deep-copied and promoted to `ScalarT`, and non-autodiff
 *    entries are forwarded unchanged.
 * 3. Build `args_work_ad` by applying the same top-level filter shape as
 *    `args_ad` (via the `TupleArgs` mask) so both filtered tuples align.
 * 4. Zip `args_ad`, `args_work_ad`, and `grads` and run
 *    `finite_diff_gradient_auto_impl` for each top-level autodiff entry.
 *
 * The split between original mask (`args_ad`) and mutable working values
 * (`args_work_ad`) ensures coordinate traversal only perturbs locations that
 * correspond to autodiff coordinates in the original argument pack, while all
 * function evaluations run against the promoted mutable working tuple.
 *
 * @tparam F callable type
 * @tparam ScalarT scalar type used for function value and gradient storage
 * @tparam TupleArgs tuple of input arguments
 * @tparam TupleGrads tuple of gradient containers aligned to autodiff
 *         top-level arguments
 * @param[in] f function to differentiate
 * @param[out] fx function value at `args`
 * @param[in] args tuple of function arguments
 * @param[out] grads compact tuple of gradient containers
 */
template <typename F, typename ScalarT, typename TupleArgs, typename TupleGrads,
          require_tuple_t<TupleArgs>* = nullptr,
          require_tuple_t<TupleGrads>* = nullptr,
          std::size_t... Idxs,
          std::size_t... GradIdxs,
          require_t<bool_constant<
              (internal::count_autodiff_args_v<TupleArgs> > 0)>>* = nullptr>
inline void finite_diff_gradient_auto(const F& f, ScalarT& fx, TupleArgs&& args,
                                      TupleGrads&& grads,
                                      std::index_sequence<Idxs...> /* idxs */,
                                      std::index_sequence<GradIdxs...> /* grad_idxs */) {
  fx = stan::math::apply([&f](auto&&... args_i) { return f(args_i...); },
                         args);
   using Swallow = int[];
  static_cast<void>(Swallow{(
    static_cast<void>(internal::finite_diff_gradient_auto_impl<ScalarT>(
        f, std::get<Idxs>(args),
        std::get<GradIdxs>(grads), args)), 0)...});


 }
}  // namespace math
}  // namespace stan
#endif
