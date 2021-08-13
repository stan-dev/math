#ifndef STAN_MATH_PRIM_FUNCTOR_CLOSURE_ADAPTER_HPP
#define STAN_MATH_PRIM_FUNCTOR_CLOSURE_ADAPTER_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/functor/apply.hpp>
#include <ostream>

namespace stan {
namespace math {
namespace internal {

/**
 * A closure that wraps a C++ lambda and captures values.
 */
template <bool Ref, typename F, typename... Ts>
struct base_closure {
  using return_scalar_t_ = return_type_t<Ts...>;
  /*The base closure with `Ts` as the non-expression partials of `Ts`*/
  using partials_closure_t_
      = base_closure<false, F, decltype(eval(value_of(std::declval<Ts>())))...>;
  using Base_ = base_closure<false, F, Ts...>;
  std::decay_t<F> f_;
  std::tuple<closure_return_type_t<Ts, Ref>...> captures_;
  template <typename FF, require_same_t<FF, F>* = nullptr, typename... Args>
  explicit base_closure(FF&& f, Args&&... args)
      : f_(std::forward<FF>(f)), captures_(std::forward<Args>(args)...) {}

  template <typename... Args>
  auto operator()(std::ostream* msgs, const Args&... args) const {
    return apply([this, msgs, &args...](
                     const auto&... s) { return f_(s..., args..., msgs); },
                 captures_);
  }
};

/**
 * A closure that takes rng argument.
 */
template <bool Ref, typename F, typename... Ts>
struct closure_rng {
  using return_scalar_t_ = double;
  using partials_closure_t_ = closure_rng<false, F, Ts...>;
  using Base_ = closure_rng<false, F, Ts...>;
  std::decay_t<F> f_;
  std::tuple<closure_return_type_t<Ts, Ref>...> captures_;

  template <typename FF, require_same_t<FF, F>* = nullptr, typename... Args>
  explicit closure_rng(FF&& f, Args&&... args)
      : f_(std::forward<FF>(f)), captures_(std::forward<Args>(args)...) {}

  template <typename Rng, typename... Args>
  auto operator()(Rng& rng, std::ostream* msgs, const Args&... args) const {
    return apply(
        [this, &rng, msgs, &args...](const auto&... s) {
          return this->f_(s..., args..., rng, msgs);
        },
        captures_);
  }
};

/**
 * A closure that can be called with `propto` template argument.
 */
template <bool Propto, bool Ref, typename F, typename... Ts>
struct closure_lpdf {
  using return_scalar_t_ = return_type_t<Ts...>;
  using partials_closure_t_ = closure_lpdf<Propto, false, F, Ts...>;
  using Base_ = closure_lpdf<Propto, false, F, Ts...>;
  std::decay_t<F> f_;
  std::tuple<closure_return_type_t<Ts, Ref>...> captures_;

  template <typename FF, require_same_t<FF, F>* = nullptr, typename... Args>
  explicit closure_lpdf(FF&& f, Args&&... args)
      : f_(std::forward<FF>(f)), captures_(std::forward<Args>(args)...) {}

  template <bool propto>
  auto with_propto() {
    return apply(
        [this](const auto&... args) {
          return closure_lpdf < Propto && propto, true, F,
                 Ts... > (this->f_, args...);
        },
        captures_);
  }

  template <bool propto = false, typename... Args>
  auto operator()(std::ostream* msgs, const Args&... args) const {
    return apply(
        [this, msgs, &args...](const auto&... s) {
          return this->f_.template operator()<propto>(s..., args..., msgs);
        },
        captures_);
  }
};

/**
 * A closure that accesses logprob accumulator.
 */
template <bool Propto, bool Ref, typename F, typename... Ts>
struct closure_lp {
  using return_scalar_t_ = return_type_t<Ts...>;
  using partials_closure_t_ = closure_lp<Propto, true, F, Ts...>;
  using Base_ = closure_lp<Propto, true, F, Ts...>;
  std::decay_t<F> f_;
  std::tuple<closure_return_type_t<Ts, Ref>...> captures_;

  template <typename FF, require_same_t<FF, F>* = nullptr, typename... Args>
  explicit closure_lp(FF&& f, Args&&... args)
      : f_(std::forward<FF>(f)), captures_(std::forward<Args>(args)...) {}

  template <bool propto = false, typename T_lp, typename T_lp_accum,
            typename... Args>
  auto operator()(T_lp& lp, T_lp_accum& lp_accum, std::ostream* msgs,
                  const Args&... args) const {
    return apply(
        [this, &lp, &lp_accum, msgs, &args...](const auto&... s) {
          return f_.template operator()<propto>(s..., args..., lp, lp_accum,
                                                msgs);
        },
        captures_);
  }
};

}  // namespace internal

/**
 * Higher-order functor suitable for calling a closure inside variadic ODE
 * solvers.
 */
struct ode_closure_adapter {
  template <typename F, typename T0, typename T1, typename... Args>
  auto operator()(const T0& t, const T1& y, std::ostream* msgs, F&& f,
                  Args&&... args) const {
    return std::forward<F>(f)(msgs, t, y, std::forward<Args>(args)...);
  }
};

struct integrate_ode_closure_adapter {
  template <typename F, typename T0, typename T1, typename... Args>
  auto operator()(const T0& t, const T1& y, std::ostream* msgs, F&& f,
                  Args&&... args) const {
    return to_vector(std::forward<F>(f)(msgs, t, to_array_1d(y),
                                        std::forward<Args>(args)...));
  }
};

/**
 * Create a closure from a C++ lambda and captures.
 */
template <typename F, typename... Args>
auto from_lambda(F&& f, Args&&... args) {
  return internal::base_closure<true, F, Args...>(std::forward<F>(f),
                                                  std::forward<Args>(args)...);
}

/**
 * Create a closure from an rng functor.
 */
template <typename F, typename... Args>
auto rng_from_lambda(F&& f, Args&&... args) {
  return internal::closure_rng<true, F, Args...>(std::forward<F>(f),
                                                 std::forward<Args>(args)...);
}

/**
 * Create a closure from an lpdf functor.
 */
template <bool propto, typename F, typename... Args>
auto lpdf_from_lambda(F&& f, Args&&... args) {
  return internal::closure_lpdf<propto, true, F, Args...>(
      std::forward<F>(f), std::forward<Args>(args)...);
}

/**
 * Create a closure from a functor that needs access to logprob accumulator.
 */
template <bool Propto, typename F, typename... Args>
auto lp_from_lambda(F&& f, Args&&... args) {
  return internal::closure_lp<Propto, true, F, Args...>(
      std::forward<F>(f), std::forward<Args>(args)...);
}

/**
 * Higher-order functor that invokes a closure inside a reduce_sum call.
 */
struct reduce_sum_closure_adapter {
  template <typename F, typename T, typename... Args>
  auto operator()(const std::vector<T>& sub_slice, std::size_t start,
                  std::size_t end, std::ostream* msgs, F&& f,
                  Args&&... args) const {
    return std::forward<F>(f)(msgs, sub_slice, start + error_index::value,
                              end + error_index::value,
                              std::forward<Args>(args)...);
  }
};

}  // namespace math
}  // namespace stan

#endif
