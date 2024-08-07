#ifndef STAN_MATH_REV_FUN_MDIVIDE_LEFT_TRI_HPP
#define STAN_MATH_REV_FUN_MDIVIDE_LEFT_TRI_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/core/typedefs.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/typedefs.hpp>

namespace stan {
namespace math {

/**
 * Returns the solution of the system Ax=B when A is triangular.
 *
 *
 * @tparam TriView Specifies whether A is upper (Eigen::Upper)
 * or lower triangular (Eigen::Lower).
 * @tparam T1 type of the triangular matrix
 * @tparam T2 type of the right-hand side matrix or vector
 *
 * @param A Triangular matrix.
 * @param B Right hand side matrix or vector.
 * @return x = A^-1 B, solution of the linear system.
 * @throws std::domain_error if A is not square or B does not have
 * as many rows as A has columns.
 */
template <Eigen::UpLoType TriView, typename T1, typename T2,
          require_all_matrix_t<T1, T2>* = nullptr,
          require_any_st_var<T1, T2>* = nullptr>
inline auto mdivide_left_tri(T1&& A, T2&& B) {
  using ret_val_type = plain_type_t<decltype(value_of(A) * value_of(B))>;
  using ret_type = return_var_matrix_t<ret_val_type, T1, T2>;
  if (A.size() == 0) {
    return arena_t<ret_type>(ret_val_type(0, B.cols()));
  }

  check_square("mdivide_left_tri", "A", A);
  check_multiplicable("mdivide_left_tri", "A", A, "B", B);

  arena_t<T1> arena_A = std::forward<T1>(A);
  if constexpr (is_autodiffable_v<T1, T2>) {
    arena_t<T2> arena_B = std::forward<T2>(B);
    auto arena_A_val = to_arena(arena_A.val());
    arena_t<ret_type> res
        = arena_A_val.template triangularView<TriView>().solve(arena_B.val());
    reverse_pass_callback([arena_A, arena_B, arena_A_val, res]() mutable {
      arena_t<promote_scalar_t<double, T2>> adjB
          = arena_A_val.template triangularView<TriView>().transpose().solve(
              res.adj());
      arena_B.adj() += adjB;
      arena_A.adj() -= (adjB * res.val().transpose().eval())
                           .template triangularView<TriView>();
    });
    return res;
  } else if constexpr (is_autodiffable_v<T1>) {
    auto arena_A_val = to_arena(std::forward<T1>(A).val());
    arena_t<ret_type> res
        = arena_A_val.template triangularView<TriView>().solve(B);
    reverse_pass_callback([arena_A, arena_A_val, res]() mutable {
      arena_A.adj()
          -= (arena_A_val.template triangularView<TriView>().transpose().solve(
                  res.adj())
              * res.val().transpose().eval())
                 .template triangularView<TriView>();
    });
    return res;
  } else {
    arena_t<T2> arena_B = std::forward<T2>(B);
    arena_t<ret_type> res
        = arena_A.template triangularView<TriView>().solve(arena_B.val());
    reverse_pass_callback([arena_A, arena_B, res]() mutable {
      arena_B.adj()
          += arena_A.template triangularView<TriView>().transpose().solve(
              res.adj());
    });
    return res;
  }
}

}  // namespace math
}  // namespace stan
#endif
