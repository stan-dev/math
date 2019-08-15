#ifndef STAN_MATH_PRIM_MAT_FUN_QUAD_FORM_HPP
#define STAN_MATH_PRIM_MAT_FUN_QUAD_FORM_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/mat/err/check_multiplicable.hpp>
#include <stan/math/prim/mat/err/check_square.hpp>

namespace stan {
namespace math {
/**
 * Compute B^T A B
 **/
template <typename T1, typename T2, typename = enable_if_all_eigen<T1, T2>, std::enable_if_t<T2::ColsAtCompileTime != 1>* = nullptr>
inline auto quad_form(const T1& A, const T2& B) {
  check_square("quad_form", "A", A);
  check_multiplicable("quad_form", "A", A, "B", B);
  return B.transpose() * A * B;
}

template <typename T1, typename T2, typename = enable_if_all_eigen<T1, T2>, std::enable_if_t<T2::ColsAtCompileTime == 1>* = nullptr>
inline auto quad_form(const T1& A,
                   const T2& B) {
  check_square("quad_form", "A", A);
  check_multiplicable("quad_form", "A", A, "B", B);
  return B.dot(A * B);
}

template <typename T1, typename T2, typename = enable_if_all_eigen<T1, T2>, std::enable_if_t<T1::ColsAtCompileTime == 1>* = nullptr>
inline auto quad_form(const T1& A,
                   const T2& B) {
  check_square("quad_form", "A", A);
  check_multiplicable("quad_form", "A", A, "B", B);
  return A.dot(B * A);
}

}  // namespace math
}  // namespace stan

#endif
