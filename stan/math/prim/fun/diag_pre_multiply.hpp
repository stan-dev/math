#ifndef STAN_MATH_PRIM_FUN_DIAG_PRE_MULTIPLY_HPP
#define STAN_MATH_PRIM_FUN_DIAG_PRE_MULTIPLY_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

template <typename T1, typename T2, require_eigen_vector_t<T1>* = nullptr,
          require_eigen_t<T2>* = nullptr>
auto diag_pre_multiply(const T1& m1, const T2& m2) {
  check_size_match("diag_pre_multiply", "m1.size()", m1.size(), "m2.rows()",
                   m2.rows());
  return (m1.asDiagonal() * m2).eval();
}

}  // namespace math
}  // namespace stan
#endif
