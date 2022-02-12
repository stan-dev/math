#ifndef STAN_MATH_REV_FUN_DIAG_PRE_MULTIPLY_HPP
#define STAN_MATH_REV_FUN_DIAG_PRE_MULTIPLY_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/rev/core.hpp>

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
template <typename T1, typename T2, require_vector_t<T1>* = nullptr,
          require_matrix_t<T2>* = nullptr,
          require_any_st_var<T1, T2>* = nullptr>
auto diag_pre_multiply(const T1& m1, const T2& m2) {
  check_size_match("diag_pre_multiply", "m1.size()", m1.size(), "m2.rows()",
                   m2.rows());
  check_size_match("diag_pre_multiply", "m1.size()", m1.size(), "m2.rows()",
                   m2.rows());
  auto val_fun = [&](auto&& x, auto&& y) { return x.asDiagonal() * y; };
  auto grad_fun_m1 = [&](auto&& val, auto&& adj, auto&& x, auto&& y) {
    return y.cwiseProduct(adj).rowwise().sum();
  };
  auto grad_fun_m2 = [&](auto&& val, auto&& adj, auto&& x, auto&& y) {
    return x.asDiagonal() * adj;
  };
  return user_gradients(std::forward_as_tuple(m1, m2),
                        std::forward<decltype(val_fun)>(val_fun),
                        std::forward_as_tuple(grad_fun_m1, grad_fun_m2));
}

}  // namespace math
}  // namespace stan

#endif
