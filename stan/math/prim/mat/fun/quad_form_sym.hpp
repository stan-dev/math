#ifndef STAN_MATH_PRIM_MAT_FUN_QUAD_FORM_SYM_HPP
#define STAN_MATH_PRIM_MAT_FUN_QUAD_FORM_SYM_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/mat/err/check_multiplicable.hpp>
#include <stan/math/prim/mat/err/check_square.hpp>
#include <stan/math/prim/mat/err/check_symmetric.hpp>

namespace stan {
namespace math {

template <typename T1, typename T2, typename = enable_if_all_eigen<T1, T2>,
          std::enable_if_t<T2::ColsAtCompileTime != 1>* = nullptr>
inline auto quad_form_sym(const T1& A, const T2& B) {
  check_square("quad_form_sym", "A", A);
  check_multiplicable("quad_form_sym", "A", A, "B", B);
  check_symmetric("quad_form_sym", "A", A);
  using scalar_val = return_type_t<typename T1::Scalar, typename T2::Scalar>;
  auto ret(B.transpose() * A * B);
  return scalar_val(0.5) * (ret + ret.transpose());
}

template <typename T1, typename T2, typename = enable_if_all_eigen<T1, T2>,
          std::enable_if_t<T2::ColsAtCompileTime == 1>* = nullptr>
inline auto quad_form_sym(const T1& A, const T2& B) {
  check_square("quad_form_sym", "A", A);
  check_multiplicable("quad_form_sym", "A", A, "B", B);
  check_symmetric("quad_form_sym", "A", A);
  return B.dot(A * B);
}

template <typename T1, typename T2, typename = enable_if_all_eigen<T1, T2>,
          std::enable_if_t<T1::ColsAtCompileTime == 1>* = nullptr>
inline auto quad_form_sym(const T1& A, const T2& B) {
  check_square("quad_form_sym", "A", A);
  check_multiplicable("quad_form_sym", "A", A, "B", B);
  check_symmetric("quad_form_sym", "A", A);
  return A.dot(B * A);
}

}  // namespace math
}  // namespace stan
#endif
