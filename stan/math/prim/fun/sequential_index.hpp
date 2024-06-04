#ifndef STAN_MATH_PRIM_FUN_SEQUENTIAL_INDEX_HPP
#define STAN_MATH_PRIM_FUN_SEQUENTIAL_INDEX_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/meta/is_tuple.hpp>
#include <stan/math/prim/fun/num_elements.hpp>
#include <stan/math/prim/functor/apply_at.hpp>

namespace stan {
namespace math {

/**
 * Utility function for indexing arbitrary types as sequential values, for use
 * as both lvalues and rvalues.
 *
 * Base template for scalars where no indexing is needed.
 *
 * @tparam Type of input scalar
 * @param x Input scalar
 * @return Input scalar unchanged
 */
template <typename T, require_stan_scalar_t<T>* = nullptr>
inline decltype(auto) sequential_index(size_t /* i */, T&& x) {
  return std::forward<T>(x);
}

/**
 * Utility function for indexing arbitrary types as sequential values, for use
 * as both lvalues and rvalues.
 *
 * Template for non-nested std::vectors
 *
 * @tparam Type of non-nested std::vector
 * @param i Index of desired value
 * @param x Input vector
 * @return Value at desired index in container
 */
template <typename T, require_std_vector_vt<is_stan_scalar, T>* = nullptr>
inline decltype(auto) sequential_index(size_t i, T&& x) {
  return x[i];
}

/**
 * Utility function for indexing arbitrary types as sequential values, for use
 * as both lvalues and rvalues.
 *
 * Template for Eigen types
 *
 * @tparam Type of Eigen input
 * @param i Index of desired value
 * @param x Input Eigen object
 * @return Value at desired index in container
 */
template <typename T, require_eigen_t<T>* = nullptr>
inline decltype(auto) sequential_index(size_t i, T&& x) {
  return x.coeffRef(i);
}

/**
 * Utility function for indexing arbitrary types as sequential values, for use
 * as both lvalues and rvalues.
 *
 * Template for nested std::vectors
 *
 * @tparam Type of nested std::vector
 * @param i Index of desired value
 * @param x Input vector
 * @return Value at desired index in container (recursively extracted)
 */
template <typename T, require_std_vector_vt<is_container, T>* = nullptr>
inline decltype(auto) sequential_index(size_t i, T&& x) {
  size_t inner_idx = i;
  size_t elem = 0;
  for (auto&& x_val : x) {
    size_t num_elems = math::num_elements(x_val);
    if (inner_idx <= (num_elems - 1)) {
      break;
    }
    elem++;
    inner_idx -= num_elems;
  }
  return sequential_index(inner_idx, std::forward<decltype(x[elem])>(x[elem]));
}

/**
 * Utility function for indexing arbitrary types as sequential values, for use
 * as both lvalues and rvalues.
 *
 * Template for tuples.
 *
 * @tparam Type of tuple
 * @param i Index of desired value
 * @param x Input tuple
 * @return Value at desired index in tuple (recursively extracted if needed)
 */
template <typename T, math::require_tuple_t<T>* = nullptr>
inline decltype(auto) sequential_index(size_t i, T&& x) {
  size_t inner_idx = i;
  size_t elem = 0;

  auto num_functor = [](auto&& arg) { return math::num_elements(arg); };
  for (size_t j = 0; j < std::tuple_size<std::decay_t<T>>{}; j++) {
    size_t num_elems = math::apply_at(num_functor, j, std::forward<T>(x));
    if (inner_idx <= (num_elems - 1)) {
      break;
    }
    elem++;
    inner_idx -= num_elems;
  }

  auto index_func = [inner_idx](auto&& t_elem) -> decltype(auto) {
    return sequential_index(inner_idx, std::forward<decltype(t_elem)>(t_elem));
  };
  return math::apply_at(index_func, elem, std::forward<T>(x));
}
}  // namespace math
}  // namespace stan
#endif
