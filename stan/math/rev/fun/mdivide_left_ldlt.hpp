#ifndef STAN_MATH_REV_FUN_MDIVIDE_LEFT_LDLT_HPP
#define STAN_MATH_REV_FUN_MDIVIDE_LEFT_LDLT_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/LDLT_factor.hpp>
#include <stan/math/rev/core/typedefs.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <memory>

namespace stan {
namespace math {

/**
 * Returns the solution of the system Ax=b given an LDLT_factor of A
 *
 * @tparam T type of B
 * @param A LDLT_factor
 * @param B Right hand side matrix or vector.
 * @return x = A^-1 B, solution of the linear system.
 * @throws std::domain_error if rows of B don't match the size of A.
 */
template <typename T1, typename T2, require_all_matrix_t<T1, T2>* = nullptr,
          require_any_st_var<T1, T2>* = nullptr>
inline auto mdivide_left_ldlt(LDLT_factor<T1>& A, T2&& B) {
  using ret_val_type = Eigen::Matrix<double, Eigen::Dynamic,
                                     std::decay_t<T2>::ColsAtCompileTime>;
  using ret_type = promote_var_matrix_t<ret_val_type, T1, T2>;

  check_multiplicable("mdivide_left_ldlt", "A", A.matrix().val(), "B", B);

  if (A.matrix().size() == 0) {
    return arena_t<ret_type>(ret_val_type(0, B.cols()));
  }

  if constexpr (!is_constant_v<T1> && !is_constant_v<T2>) {
    arena_t<T2> arena_B = std::forward<T2>(B);
    arena_t<T1> arena_A = A.matrix();
    arena_t<ret_type> res = A.ldlt().solve(arena_B.val());
    const auto* ldlt_ptr = make_chainable_ptr(A.ldlt());

    reverse_pass_callback([arena_A, arena_B, ldlt_ptr, res]() mutable {
      promote_scalar_t<double, T2> adjB = ldlt_ptr->solve(res.adj());

      arena_A.adj() -= adjB * res.val_op().transpose();
      arena_B.adj() += adjB;
    });

    return res;
  } else if constexpr (!is_constant_v<T1>) {
    arena_t<T1> arena_A = A.matrix();
    arena_t<ret_type> res = A.ldlt().solve(std::forward<T2>(B));
    const auto* ldlt_ptr = make_chainable_ptr(A.ldlt());

    reverse_pass_callback([arena_A, ldlt_ptr, res]() mutable {
      arena_A.adj() -= ldlt_ptr->solve(res.adj()) * res.val_op().transpose();
    });

    return res;
  } else {
    arena_t<T2> arena_B = std::forward<T2>(B);
    arena_t<ret_type> res = A.ldlt().solve(arena_B.val());
    const auto* ldlt_ptr = make_chainable_ptr(A.ldlt());

    reverse_pass_callback([arena_B, ldlt_ptr, res]() mutable {
      arena_B.adj() += ldlt_ptr->solve(res.adj());
    });

    return res;
  }
}

}  // namespace math
}  // namespace stan
#endif
