#ifndef STAN_MATH_PRIM_FUNCTOR_CLOSURE_ADAPTER_HPP
#define STAN_MATH_PRIM_FUNCTOR_CLOSURE_ADAPTER_HPP

#include <stan/math/prim/meta/error_index.hpp>
#include <stan/math/prim/meta/return_type.hpp>
#include <stan/math/prim/meta/is_stan_closure.hpp>
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
  using captured_scalar_t__ = return_type_t<Ts...>;
  using ValueOf__
      = base_closure<false, F, decltype(eval(value_of(std::declval<Ts>())))...>;
  using CopyOf__ = base_closure<false, F, Ts...>;
  F f_;
  std::tuple<capture_type_t<Ts, Ref>...> captures_;

  explicit base_closure(const F& f, const Ts&... args)
      : f_(f), captures_(args...) {}

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
  using captured_scalar_t__ = double;
  using ValueOf__ = closure_rng<false, F, Ts...>;
  using CopyOf__ = closure_rng<false, F, Ts...>;
  F f_;
  std::tuple<capture_type_t<Ts, Ref>...> captures_;

  explicit closure_rng(const F& f, const Ts&... args)
      : f_(f), captures_(args...) {}

  template <typename Rng, typename... Args>
  auto operator()(Rng& rng, std::ostream* msgs, const Args&... args) const {
    return apply([this, &rng, msgs, &args...](
                     const auto&... s) { return f_(s..., args..., rng, msgs); },
                 captures_);
  }
};

/**
 * A closure that can be called with `propto` template argument.
 */
template <bool Propto, bool Ref, typename F, typename... Ts>
struct closure_lpdf {
  using captured_scalar_t__ = return_type_t<Ts...>;
  using ValueOf__ = closure_lpdf<Propto, false, F, Ts...>;
  using CopyOf__ = closure_lpdf<Propto, false, F, Ts...>;
  F f_;
  std::tuple<capture_type_t<Ts, Ref>...> captures_;

  explicit closure_lpdf(const F& f, const Ts&... args)
      : f_(f), captures_(args...) {}

  template <bool propto>
  auto with_propto() {
    return apply(
        [this](const auto&... args) {
          return closure_lpdf < Propto && propto, true, F,
                 Ts... > (f_, args...);
        },
        captures_);
  }

  template <bool propto = false, typename... Args>
  auto operator()(std::ostream* msgs, const Args&... args) const {
    return apply(
        [this, msgs, &args...](const auto&... s) {
          return f_.template operator()<propto>(s..., args..., msgs);
        },
        captures_);
  }
};

/**
 * A closure that accesses logprob accumulator.
 */
template <bool Propto, bool Ref, typename F, typename... Ts>
struct closure_lp {
  using captured_scalar_t__ = return_type_t<Ts...>;
  using ValueOf__ = closure_lp<Propto, true, F, Ts...>;
  using CopyOf__ = closure_lp<Propto, true, F, Ts...>;
  F f_;
  std::tuple<capture_type_t<Ts, Ref>...> captures_;

  explicit closure_lp(const F& f, const Ts&... args)
      : f_(f), captures_(args...) {}

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
  auto operator()(const T0& t, const T1& y, std::ostream* msgs, const F& f,
                  Args... args) const {
    return f(msgs, t, y, args...);
  }
};

struct integrate_ode_closure_adapter {
  template <typename F, typename T0, typename T1, typename... Args>
  auto operator()(const T0& t, const T1& y, std::ostream* msgs, const F& f,
                  Args... args) const {
    return to_vector(f(msgs, t, to_array_1d(y), args...));
  }
};

/**
 * Create a closure from a C++ lambda and captures.
 */
template <typename F, typename... Ts>
auto from_lambda(const F& f, const Ts&... a) {
  return internal::base_closure<true, F, Ts...>(f, a...);
}

/**
 * Create a closure from an rng functor.
 */
template <typename F, typename... Ts>
auto rng_from_lambda(const F& f, const Ts&... a) {
  return internal::closure_rng<true, F, Ts...>(f, a...);
}

/**
 * Create a closure from an lpdf functor.
 */
template <bool propto, typename F, typename... Ts>
auto lpdf_from_lambda(const F& f, const Ts&... a) {
  return internal::closure_lpdf<propto, true, F, Ts...>(f, a...);
}

/**
 * Create a closure from a functor that needs access to logprob accumulator.
 */
template <bool Propto, typename F, typename... Ts>
auto lp_from_lambda(const F& f, const Ts&... args) {
  return internal::closure_lp<Propto, true, F, Ts...>(f, args...);
}

/**
 * Higher-order functor that invokes a closure inside a reduce_sum call.
 */
struct reduce_sum_closure_adapter {
  template <typename F, typename T, typename... Args>
  auto operator()(const std::vector<T>& sub_slice, std::size_t start,
                  std::size_t end, std::ostream* msgs, const F& f,
                  Args... args) const {
    return f(msgs, sub_slice, start + error_index::value,
             end + error_index::value, args...);
  }
};

}  // namespace math
}  // namespace stan

#endif
