#ifndef STAN_MATH_PRIM_FUN_DIAG_POST_MULTIPLY_HPP
#define STAN_MATH_PRIM_FUN_DIAG_POST_MULTIPLY_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

template <typename T1, typename T2, require_eigen_t<T1>* = nullptr,
          require_eigen_vector_t<T2>* = nullptr>
auto diag_post_multiply(const T1& m1, const T2& m2) {
  check_size_match("diag_post_multiply", "m2.size()", m2.size(), "m1.cols()",
                   m1.cols());
  return (m1 * m2.asDiagonal()).eval();
}

}  // namespace math
}  // namespace stan
#endif
