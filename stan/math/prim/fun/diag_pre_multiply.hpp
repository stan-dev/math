#ifndef STAN_MATH_PRIM_FUN_DIAG_PRE_MULTIPLY_HPP
#define STAN_MATH_PRIM_FUN_DIAG_PRE_MULTIPLY_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

template <typename EigMat1, typename EigMat2,
 typename = require_all_eigen_t<EigMat1, EigMat2>>
auto diag_pre_multiply(EigMat1&& m1, EigMat2&& m2) {
  check_vector("diag_pre_multiply", "m1", m1);
  check_size_match("diag_pre_multiply", "m1.size()", m1.size(), "m2.rows()",
                   m2.rows());
  return m1.asDiagonal() * m2;
}

}  // namespace math
}  // namespace stan
#endif
