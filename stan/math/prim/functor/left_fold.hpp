#ifndef STAN_MATH_PRIM_SCAL_FUNCTOR_LEFT_FOLD_HPP
#define STAN_MATH_PRIM_SCAL_FUNCTOR_LEFT_FOLD_HPP

#include <stan/math/prim/meta.hpp>
#include <functional>
#include <tuple>
#include <utility>

namespace stan {
namespace math {

/**
 * Performs a left fold on a tuple with a given operation
 * @tparam index The index in this recursion for accessing tuple elemenTypes
 * @tparam Op An operation to perform on the tuple elemenTypes
 * @tparam Types types of tuple elements
 * @param op operation to perform on individual tuple elements
 * @param t A tuple to iterate over
 */
template <size_t index, class Op, class... Types>
constexpr auto left_fold(Op op, const std::tuple<Types...>& t) {
  if (index == sizeof...(Types) - 1) {
    return std::get<index>(t);
  } else {
    return op(std::get<index>(t), left_fold<1 + index>(op, t));
  }
}

/**
 * Sum the elements of a tuple
 * @tparam Types types of the tuple's elements
 * @param t tuple whose elements are to be summed
 */
template <typename... Types, require_stan_scalar_t<Types>...>
constexpr auto sum(const std::tuple<Types...>& t) {
  return left_fold<0>(std::plus<>{}, t);
}

}  // namespace math
}  // namespace stan
#endif
