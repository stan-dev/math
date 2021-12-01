#ifndef STAN_MATH_PRIM_FUN_AS_COLUMN_VECTOR_OR_SCALAR_HPP
#define STAN_MATH_PRIM_FUN_AS_COLUMN_VECTOR_OR_SCALAR_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <vector>

namespace stan {
namespace math {
namespace internal {
template <typename T, typename S, typename Enable>
class empty_broadcast_array;
}
/**
 * no-op that passes the scalar
 *
 * @tparam T Type of scalar element.
 * @param a Specified scalar.
 * @return the scalar.
 */
template <typename T, require_stan_scalar_t<T>* = nullptr>
inline T as_column_vector_or_scalar(const T& a) {
  return a;
}

/**
 * No-op used when working with operands and partials.
 * This is not implimented so it cannot be invoked and only exists so the
 * compiler can resolve it's output type.
 */
template <typename T, typename S>
internal::empty_broadcast_array<T, S, void>& as_column_vector_or_scalar(
    internal::empty_broadcast_array<T, S, void>& a);

/**
 * no-op that returns a column vector.
 *
 * @tparam T Type inheriting from `EigenBase` with dynamic compile time rows
 *  and fixed column of 1.
 * @param a Specified vector.
 * @return Same vector.
 */
template <typename T, require_eigen_col_vector_t<T>* = nullptr>
inline T&& as_column_vector_or_scalar(T&& a) {
  return std::forward<T>(a);
}

/**
 * Converts a row vector to an eigen column vector. For row vectors this returns
 *  a `Transpose<Eigen::Matrix<T, 1, -1>>`.
 *
 * @tparam T Type inheriting from `EigenBase` with dynamic compile time columns
 * and fixed row of 1.
 * @param a Specified vector.
 * @return Transposed vector.
 */
template <typename T, require_eigen_row_vector_t<T>* = nullptr,
          require_not_eigen_col_vector_t<T>* = nullptr>
inline auto as_column_vector_or_scalar(T&& a) {
  return make_holder([](auto& x) { return x.transpose(); }, std::forward<T>(a));
}

/**
 * Converts `std::vector` to a column vector.
 *
 * @tparam T `std::vector` type.
 * @param a Specified vector.
 * @return input converted to a column vector.
 */
template <typename T, require_std_vector_t<T>* = nullptr>
inline auto as_column_vector_or_scalar(T&& a) {
  using plain_vector = Eigen::Matrix<value_type_t<T>, Eigen::Dynamic, 1>;
  using optionally_const_vector
      = std::conditional_t<std::is_const<std::remove_reference_t<T>>::value,
                           const plain_vector, plain_vector>;
  using T_map = Eigen::Map<optionally_const_vector>;
  return make_holder([](auto& x) { return T_map(x.data(), x.size()); },
                     std::forward<T>(a));
}

}  // namespace math
}  // namespace stan

#endif
