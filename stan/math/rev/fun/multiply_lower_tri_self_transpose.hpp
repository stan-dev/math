#ifndef STAN_MATH_REV_FUN_MULTIPLY_LOWER_TRI_SELF_TRANSPOSE_HPP
#define STAN_MATH_REV_FUN_MULTIPLY_LOWER_TRI_SELF_TRANSPOSE_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/dot_product.hpp>
#include <stan/math/rev/fun/dot_self.hpp>
#include <stan/math/rev/fun/typedefs.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

inline matrix_v multiply_lower_tri_self_transpose(const matrix_v& L) {
  // check_square("multiply_lower_tri_self_transpose",
  // L, "L", (double*)0);
  int K = L.rows();
  int J = L.cols();
  matrix_v LLt(K, K);
  if (K == 0) {
    return LLt;
  }
  // if (K == 1) {
  //   LLt(0, 0) = L(0, 0) * L(0, 0);
  //   return LLt;
  // }
  int Knz;
  if (K >= J) {
    Knz = (K - J) * J + (J * (J + 1)) / 2;
  } else {  // if (K < J)
    Knz = (K * (K + 1)) / 2;
  }
  vari** vs = reinterpret_cast<vari**>(
      ChainableStack::instance_->memalloc_.alloc(Knz * sizeof(vari*)));
  int pos = 0;
  for (int m = 0; m < K; ++m) {
    for (int n = 0; n < ((J < (m + 1)) ? J : (m + 1)); ++n) {
      vs[pos++] = L(m, n).vi_;
    }
  }
  for (int m = 0, mpos = 0; m < K; ++m, mpos += (J < m) ? J : m) {
    LLt.coeffRef(m, m) = var(
        new internal::dot_self_vari(vs + mpos, (J < (m + 1)) ? J : (m + 1)));
    for (int n = 0, npos = 0; n < m; ++n, npos += (J < n) ? J : n) {
      LLt.coeffRef(m, n) = LLt.coeffRef(n, m)
          = dot_product(L.row(m).head((J < (n + 1)) ? J : (n + 1)),
                        L.row(n).head((J < (n + 1)) ? J : (n + 1)));
    }
  }
  return LLt;
}

}  // namespace math
}  // namespace stan
#endif
