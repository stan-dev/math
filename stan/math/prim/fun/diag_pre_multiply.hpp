#ifndef STAN_MATH_PRIM_FUN_DIAG_PRE_MULTIPLY_HPP
#define STAN_MATH_PRIM_FUN_DIAG_PRE_MULTIPLY_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/to_ref.hpp>

namespace stan {
namespace math {

/**
 * Return the product of the diagonal matrix formed from the vector
 * or row_vector and a matrix.
 *
 * @tparam T1 type of the vector/row_vector
 * @tparam T2 type of the matrix
 * @param m1 input vector/row_vector
 * @param m2 input matrix
 *
 * @return product of the diagonal matrix formed from the
 * vector or row_vector and a matrix.
 */
template <typename T1, typename T2, require_eigen_vector_t<T1>* = nullptr,
          require_eigen_t<T2>* = nullptr,
          require_all_not_st_var<T1, T2>* = nullptr>
inline auto diag_pre_multiply(T1&& m1, T2&& m2) {
  auto&& m1_ref = to_ref(std::forward<T1>(m1));
  auto&& m2_ref = to_ref(std::forward<T2>(m2));
  check_size_match("diag_pre_multiply", "m1.size()", m1_ref.size(), "m2.rows()",
                   m2_ref.rows());
  return make_holder(
      [](auto&& m1_, auto&& m2_) {
        return std::forward<decltype(m1_)>(m1_).asDiagonal()
               * std::forward<decltype(m2_)>(m2_);
      },
      std::forward<decltype(m1_ref)>(m1_ref),
      std::forward<decltype(m2_ref)>(m2_ref));
}

}  // namespace math
}  // namespace stan
#endif
