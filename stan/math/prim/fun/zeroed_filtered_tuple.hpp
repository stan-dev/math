#ifndef STAN_MATH_PRIM_FUN_ZEROED_FILTERED_TUPLE_HPP
#define STAN_MATH_PRIM_FUN_ZEROED_FILTERED_TUPLE_HPP

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
 * Build a zero-initialized container matching the shape of `arg`.
 *
 * Scalars produce zero partials scalars, Eigen produces zero matrices of the
 * same dimensions, and std::vector recurses by element.
 *
 * @tparam Arg input container/scalar type
 * @param[in] arg input value used only for shape
 * @return zero-initialized container with the same structure and dimensions
 * as `arg`
 * @throw `std::bad_alloc` if allocation of nested `std::vector` storage fails.
 * @throw None otherwise. Unsupported categories fail at compile time via
 * `static_assert`.
 */
template <typename Arg>
inline auto zeroed_filtered_tuple_leaf(const Arg& arg) {
  using arg_t = std::decay_t<Arg>;
  using partial_t = partials_type_t<scalar_type_t<arg_t>>;
  if constexpr (is_std_vector_v<arg_t>) {
    using elem_t = std::decay_t<decltype(
        zeroed_filtered_tuple_leaf(std::declval<const value_type_t<arg_t>&>()))>;
    std::vector<elem_t> ret;
    ret.reserve(arg.size());
    for (const auto& element : arg) {
      ret.push_back(zeroed_filtered_tuple_leaf(element));
    }
    return ret;
  } else if constexpr (is_eigen_v<arg_t>) {
    using eig_t = promote_scalar_t<partial_t, plain_type_t<arg_t>>;
    return eig_t::Zero(arg.rows(), arg.cols()).eval();
  } else if constexpr (is_stan_scalar_v<arg_t>) {
    return partial_t(0);
  } else {
    static_assert(
        sizeof(std::decay_t<Arg>*) == 0,
        "Unsupported container in zeroed_filtered_tuple.");
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
 * @param[in] args tuple of arguments used for shape and filtering
 * @return tuple of zero-initialized containers for filtered entries
 * @throw `std::bad_alloc` if nested vector allocation fails while constructing
 * the returned containers.
 * @throw None otherwise. Unsupported categories fail at compile time via
 * `static_assert`.
 */
template <template <class...> class Predicate, typename TupleArgs,
          require_tuple_t<TupleArgs>* = nullptr>
inline auto zeroed_filtered_tuple(TupleArgs&& args) {
  return stan::math::filter_map<Predicate>(
      [](auto&& arg) {
        using arg_t = std::decay_t<decltype(arg)>;
        if constexpr (is_tuple_v<arg_t>) {
          return zeroed_filtered_tuple<Predicate>(std::forward<decltype(arg)>(arg));
        } else {
          return internal::zeroed_filtered_tuple_leaf(arg);
        }
      },
      std::forward<TupleArgs>(args));
}

}  // namespace math
}  // namespace stan

#endif
