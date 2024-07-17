#ifndef STAN_MATH_REV_FUN_MDIVIDE_LEFT_SPD_HPP
#define STAN_MATH_REV_FUN_MDIVIDE_LEFT_SPD_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/core/typedefs.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <vector>

namespace stan {
namespace math {
/**
 * Returns the solution of the system Ax=B where A is symmetric positive
 * definite.
 *
 *
 * @tparam T1 type of the first matrix
 * @tparam T2 type of the right-hand side matrix or vector
 *
 * @param A Matrix.
 * @param B Right hand side matrix or vector.
 * @return x = A^-1 B, solution of the linear system.
 * @throws std::domain_error if A is not square or B does not have
 * as many rows as A has columns.
 */
template <typename T1, typename T2, require_all_matrix_t<T1, T2> * = nullptr,
  require_any_st_var<T1, T2>* = nullptr>
inline auto mdivide_left_spd(T1 &&A, T2 &&B) {
  using ret_val_type = plain_type_t<decltype(value_of(A) * value_of(B))>;
  using ret_type = return_var_matrix_t<ret_val_type, T1, T2>;
  if (A.size() == 0) {
    return arena_t<ret_type>(ret_val_type(0, B.cols()));
  }
  check_multiplicable("mdivide_left_spd", "A", A, "B", B);
  if constexpr (is_autodiffable_v<T1, T2>) {
    arena_t<T1> arena_A = std::forward<T1>(A);
    check_symmetric("mdivide_left_spd", "A", arena_A.val());
    check_not_nan("mdivide_left_spd", "A", arena_A.val());
    auto A_llt = arena_A.val().llt();
    check_pos_definite("mdivide_left_spd", "A", A_llt);
    arena_t<Eigen::MatrixXd> arena_A_llt = A_llt.matrixL();
    arena_t<T2> arena_B = std::forward<T2>(B);
    arena_t<ret_type> res = A_llt.solve(arena_B.val());
    reverse_pass_callback([arena_A, arena_B, arena_A_llt, res]() mutable {
      arena_t<decltype(res.adj().eval())> adjB = res.adj().eval();
      arena_A_llt.template triangularView<Eigen::Lower>().solveInPlace(adjB);
      arena_A_llt.template triangularView<Eigen::Lower>()
          .transpose()
          .solveInPlace(adjB);

      arena_A.adj() -= adjB * res.val_op().transpose();
      arena_B.adj() += adjB;
    });
    return res;
  } else if constexpr (is_autodiffable_v<T1>) {
    arena_t<T1> arena_A = std::forward<T1>(A);
    check_symmetric("mdivide_left_spd", "A", arena_A.val());
    check_not_nan("mdivide_left_spd", "A", arena_A.val());
    auto A_llt = arena_A.val().llt();
    check_pos_definite("mdivide_left_spd", "A", A_llt);
    arena_t<Eigen::MatrixXd> arena_A_llt = A_llt.matrixL();
    arena_t<ret_type> res = A_llt.solve(B);
    reverse_pass_callback([arena_A, arena_A_llt, res]() mutable {
      arena_t<decltype(res.adj().eval())> adjB = res.adj().eval();
      arena_A_llt.template triangularView<Eigen::Lower>().solveInPlace(adjB);
      arena_A_llt.template triangularView<Eigen::Lower>()
          .transpose()
          .solveInPlace(adjB);

      arena_A.adj() -= adjB * res.val().transpose().eval();
    });
    return res;
  } else {
    auto&& A_ref = to_ref(A);
    check_symmetric("mdivide_left_spd", "A", A_ref);
    check_not_nan("mdivide_left_spd", "A", A_ref);
    auto A_llt = A_ref.llt();
    check_pos_definite("mdivide_left_spd", "A", A_llt);
    arena_t<Eigen::MatrixXd> arena_A_llt = A_llt.matrixL();
    arena_t<T2> arena_B = std::forward<T2>(B);
    arena_t<ret_type> res = A_llt.solve(arena_B.val());
    reverse_pass_callback([arena_B, arena_A_llt, res]() mutable {
      arena_t<decltype(res.adj().eval())> adjB = res.adj().eval();
      arena_A_llt.template triangularView<Eigen::Lower>().solveInPlace(adjB);
      arena_A_llt.template triangularView<Eigen::Lower>()
          .transpose()
          .solveInPlace(adjB);
      arena_B.adj() += adjB;
    });
    return res;
  }
}

}  // namespace math
}  // namespace stan
#endif
