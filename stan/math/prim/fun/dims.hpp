#ifndef STAN_MATH_PRIM_FUN_DIMS_HPP
#define STAN_MATH_PRIM_FUN_DIMS_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Pushes dimensions of given argument into given result vector.
 *
 * For a scalar that is a no-op.
 * @tparam type of scalar
 * @param x argument
 * @param result result
 */
template <typename T, require_stan_scalar_t<T>* = nullptr>
inline void dims(const T& x, std::vector<int>& result) {}

/**
 * Pushes dimensions of given argument into given result vector.
 *
 * For an Eigen type those are the numbers of rows and columns.
 * @param x argument
 * @param result result
 */
template <typename T, require_matrix_t<T>* = nullptr>
inline void dims(const T& x, std::vector<int>& result) {
  result.push_back(x.rows());
  result.push_back(x.cols());
}

/**
 * Pushes dimensions of given argument into given result vector.
 *
 * For a `std::vector` that is its size and dimensions of its elements.
 * @tparam type of scalar
 * @tparam Alloc type of allocator
 * @param x argument
 * @param result result
 */
template <typename T, typename Alloc>
inline void dims(const std::vector<T, Alloc>& x, std::vector<int>& result) {
  result.push_back(x.size());
  if (x.size() > 0) {
    dims(x[0], result);
  }
}

/**
 * Determines dimensions of an argument.
 * @param x argument
 * @return vector of sizes in each of arguemnt's dimensions
 */
template <typename T>
inline std::vector<int> dims(const T& x) {
  std::vector<int> result;
  dims(x, result);
  return result;
}

}  // namespace math
}  // namespace stan
#endif
