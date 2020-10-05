#ifndef STAN_MATH_REV_FUN_MDIVIDE_LEFT_LDLT_HPP
#define STAN_MATH_REV_FUN_MDIVIDE_LEFT_LDLT_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/LDLT_alloc.hpp>
#include <stan/math/rev/fun/LDLT_factor.hpp>
#include <stan/math/rev/fun/typedefs.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <memory>

namespace stan {
namespace math {

/**
 * Returns the solution of the system Ax=b given an LDLT_factor of A
 *
 * @tparam R1 number of rows in the LDLT_factor, can be Eigen::Dynamic
 * @tparam C1 number of columns in the LDLT_factor, can be Eigen::Dynamic
 *
 * @param A LDLT_factor
 * @param B Right hand side matrix or vector.
 * @return x = A^-1 B, solution of the linear system.
 * @throws std::domain_error if rows of B don't match the size of A.
 */
template <int R, int C, typename T1, typename T2,
          require_eigen_t<T2>* = nullptr, require_any_st_var<T1, T2>* = nullptr>
inline auto mdivide_left_ldlt(const LDLT_factor<T1, R, C>& A, const T2& B) {
  using ret_type = Eigen::Matrix<var, R, T2::ColsAtCompileTime>;
  
  check_multiplicable("mdivide_left_ldlt", "A", A, "B", B);

  const auto& B_ref = to_ref(B);

  if (A.cols() == 0) {
    return ret_type(0, B.cols());
  }

  arena_t<promote_scalar_t<var, T2>> arena_B;

  if (!is_constant<T2>::value) {
    arena_B = B_ref;
  }

  arena_t<ret_type> res = A.solve(value_of(B_ref));

  reverse_pass_callback([A, arena_B, res]() mutable {
    promote_scalar_t<double, T2> adjB = A.solve(res.adj());

    if (!is_constant<T1>::value)
      forward_as<LDLT_factor<var, R, C>>(A).alloc_->arena_A_.adj()
          -= adjB * res.val().transpose().eval();

    if (!is_constant<T2>::value)
      arena_B.adj() += adjB;
  });

  return ret_type(res);
}

}  // namespace math
}  // namespace stan
#endif
