#ifndef STAN_MATH_PRIM_FUN_ALL_HPP
#define STAN_MATH_PRIM_FUN_ALL_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/functor/for_each.hpp>
#include <algorithm>

namespace stan {
namespace math {

/**
 * Return true if all values in the input are true.
 *
 * Overload for a single integral input
 *
 * @tparam T Any type convertible to `bool`
 * @param x integral input
 * @return The input unchanged
 */
template <typename T, require_t<std::is_convertible<T, bool>>* = nullptr>
constexpr inline bool all(T x) {
  return x;
}

/**
 * Return true if all values in the input are true.
 *
 * Overload for Eigen types
 *
 * @tparam ContainerT A type derived from `Eigen::EigenBase` that has an
 *                      `integral` scalar type
 * @param x Eigen object of boolean inputs
 * @return Boolean indicating whether all elements are true
 */
template <typename ContainerT,
          require_eigen_st<std::is_integral, ContainerT>* = nullptr>
inline bool all(const ContainerT& x) {
  return x.all();
}

// Forward-declaration for correct resolution of all(std::vector<std::tuple>)
template <typename... Types>
inline bool all(const std::tuple<Types...>& x);

/**
 * Return true if all values in the input are true.
 *
 * Overload for a std::vector/nested inputs. The Eigen::Map/apply_vector_unary
 * approach cannot be used as std::vector<bool> types do not have a .data()
 * member and are not always stored contiguously.
 *
 * @tparam InnerT Type within std::vector
 * @param x Nested container of boolean inputs
 * @return Boolean indicating whether all elements are true
 */
template <typename InnerT>
inline bool all(const std::vector<InnerT>& x) {
  return std::all_of(x.begin(), x.end(), [](const auto& i) { return all(i); });
}

/**
 * Return true if all values in the input are true.
 *
 * Overload for a tuple input.
 *
 * @tparam Types of items within tuple
 * @param x Tuple of boolean scalar-type elements
 * @return Boolean indicating whether all elements are true
 */
template <typename... Types>
inline bool all(const std::tuple<Types...>& x) {
  bool all_true = true;
  math::for_each(
      [&all_true](const auto& i) {
        all_true = all_true && all(i);
        return;
      },
      x);
  return all_true;
}

}  // namespace math
}  // namespace stan

#endif
