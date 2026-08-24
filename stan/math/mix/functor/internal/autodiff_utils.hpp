#ifndef STAN_MATH_MIX_FUNCTOR_INTERNAL_AUTODIFF_UTILS_HPP
#define STAN_MATH_MIX_FUNCTOR_INTERNAL_AUTODIFF_UTILS_HPP

#include <stan/math/fwd/meta/is_fvar.hpp>
#include <stan/math/prim/functor/apply.hpp>
#include <stan/math/prim/functor/make_holder_tuple.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/rev/fun/from_var_value.hpp>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace stan {
namespace math {
namespace internal {

template <typename T>
struct contains_fvar_impl
    : is_fvar<scalar_type_t<std::decay_t<T>>> {};

template <typename... Types>
struct contains_fvar_impl<std::tuple<Types...>>
    : std::disjunction<contains_fvar_impl<std::decay_t<Types>>...> {};

template <typename T, typename... VecArgs>
struct contains_fvar_impl<std::vector<T, VecArgs...>>
    : contains_fvar_impl<std::decay_t<T>> {};

template <typename... Types>
struct contains_fvar
    : std::disjunction<contains_fvar_impl<std::decay_t<Types>>...> {};

template <typename Tuple, require_tuple_t<Tuple>* = nullptr>
inline auto from_var_value_rec(Tuple&& arg);

template <typename Vec, require_std_vector_t<Vec>* = nullptr>
inline auto from_var_value_rec(Vec&& arg);

template <typename T, require_var_matrix_t<T>* = nullptr>
inline auto from_var_value_rec(T&& arg);

template <typename T,
          require_all_not_t<is_tuple<std::decay_t<T>>,
                            is_std_vector<std::decay_t<T>>,
                            is_var_matrix<std::decay_t<T>>>* = nullptr>
inline decltype(auto) from_var_value_rec(T&& arg);

/**
 * Convert reverse-mode matrix values to elementwise variables recursively.
 *
 * This lets mixed-mode functors promote both matrix-of-vars and
 * `var_value<Eigen>` arguments to nested forward-mode types while preserving
 * their reverse-mode connections.
 */
template <typename Tuple, require_tuple_t<Tuple>*>
inline auto from_var_value_rec(Tuple&& arg) {
  return math::apply(
      [](auto&&... inner_args) {
        return make_holder_tuple(from_var_value_rec(
            std::forward<decltype(inner_args)>(inner_args))...);
      },
      std::forward<Tuple>(arg));
}

template <typename Vec, require_std_vector_t<Vec>*>
inline auto from_var_value_rec(Vec&& arg) {
  using value_t = typename std::decay_t<Vec>::value_type;
  using promoted_t = std::decay_t<decltype(
      from_var_value_rec(std::declval<const value_t&>()))>;
  std::vector<promoted_t> res;
  res.reserve(arg.size());
  for (auto&& elem : arg) {
    res.emplace_back(from_var_value_rec(
        std::forward<decltype(elem)>(elem)));
  }
  return res;
}

template <typename T, require_var_matrix_t<T>*>
inline auto from_var_value_rec(T&& arg) {
  return from_var_value(std::forward<T>(arg));
}

template <typename T,
          require_all_not_t<is_tuple<std::decay_t<T>>,
                            is_std_vector<std::decay_t<T>>,
                            is_var_matrix<std::decay_t<T>>>*>
inline decltype(auto) from_var_value_rec(T&& arg) {
  return std::forward<T>(arg);
}

}  // namespace internal
}  // namespace math
}  // namespace stan

#endif
