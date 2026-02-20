#ifndef STAN_MATH_PRIM_FUN_ZEROED_CONTAINER_HPP
#define STAN_MATH_PRIM_FUN_ZEROED_CONTAINER_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/functor/filter_map.hpp>
#include <stan/math/prim/meta.hpp>
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
struct zeroed_container_dependent_false : std::false_type {};

/**
 * Build a zero-initialized container matching the shape of `arg`.
 *
 * Scalars produce zero partials scalars, Eigen produces zero matrices of the
 * same dimensions, and std::vector recurses by element.
 */
template <typename Arg>
inline auto zeroed_container_leaf(const Arg& arg) {
  using arg_t = std::decay_t<Arg>;
  using partial_t = partials_type_t<scalar_type_t<arg_t>>;
  if constexpr (is_std_vector_v<arg_t>) {
    using elem_t = std::decay_t<decltype(
        zeroed_container_leaf(std::declval<const value_type_t<arg_t>&>()))>;
    std::vector<elem_t> ret;
    ret.reserve(arg.size());
    for (const auto& element : arg) {
      ret.push_back(zeroed_container_leaf(element));
    }
    return ret;
  } else if constexpr (is_eigen_v<arg_t>) {
    using eig_t = promote_scalar_t<partial_t, plain_type_t<arg_t>>;
    return eig_t::Zero(arg.rows(), arg.cols()).eval();
  } else if constexpr (is_stan_scalar_v<arg_t>) {
    return partial_t(0);
  } else {
    static_assert(
        zeroed_container_dependent_false<arg_t>::value,
        "Unsupported container in zeroed_container.");
  }
}

}  // namespace internal

/**
 * Build a compact tuple of zero-initialized containers for entries in `args`
 * whose types satisfy `Predicate`.
 *
 * Tuple entries that do not satisfy `Predicate` are filtered out. Matching
 * tuple entries recurse so nested tuple output remains compact as well.
 *
 * @tparam Predicate unary type-trait predicate
 * @tparam TupleArgs tuple input type
 * @param args tuple of arguments used for shape and filtering
 * @return tuple of zero-initialized containers for filtered entries
 */
template <template <class...> class Predicate, typename TupleArgs,
          require_tuple_t<TupleArgs>* = nullptr>
inline auto zeroed_container(TupleArgs&& args) {
  return stan::math::filter_map<Predicate>(
      [](auto&& arg) {
        using arg_t = std::decay_t<decltype(arg)>;
        if constexpr (is_tuple_v<arg_t>) {
          return zeroed_container<Predicate>(std::forward<decltype(arg)>(arg));
        } else {
          return internal::zeroed_container_leaf(arg);
        }
      },
      std::forward<TupleArgs>(args));
}

}  // namespace math
}  // namespace stan

#endif
