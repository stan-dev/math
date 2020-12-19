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
 * @param A LDLT_factor2
 * @param B Right hand side matrix or vector.
 * @return x = A^-1 B, solution of the linear system.
 * @throws std::domain_error if rows of B don't match the size of A.
 */
template <typename T1, bool alloc_in_arena, typename T2,
	  require_all_matrix_t<T1, T2>* = nullptr,
          require_any_st_var<T1, T2>* = nullptr>
inline auto mdivide_left_ldlt(const LDLT_factor2<T1, alloc_in_arena>& A, const T2& B) {
  using ret_val_type = Eigen::Matrix<double, Eigen::Dynamic, T2::ColsAtCompileTime>;
  using ret_type = promote_var_matrix_t<ret_val_type, T1,  T2>;

  check_multiplicable("mdivide_left_ldlt", "A", A.matrix().val(), "B", B);

  if (A.matrix().size() == 0) {
    return ret_type(ret_val_type(0, B.cols()));
  }

  if (!is_constant<T1>::value && !is_constant<T2>::value) {
    arena_t<promote_scalar_t<var, T2>> arena_B = B;
    arena_t<ret_type> res = A.ldlt().solve(arena_B.val());

    reverse_pass_callback([A, arena_B, res]() mutable {
      promote_scalar_t<double, T2> adjB = A.ldlt().solve(res.adj());

      forward_as<promote_scalar_t<var, T1>>(A.matrix()).adj() -= adjB * res.val().transpose().eval();
      arena_B.adj() += adjB;
    });

    return ret_type(res);
  } else if (!is_constant<T1>::value) {
    arena_t<ret_type> res = A.ldlt().solve(value_of(B));

    reverse_pass_callback([A, res]() mutable {
      promote_scalar_t<double, T2> adjB = A.ldlt().solve(res.adj());

      forward_as<promote_scalar_t<var, T1>>(A.matrix()).adj() -= adjB * res.val().transpose().eval();
    });

    return ret_type(res);
  } else {
    arena_t<promote_scalar_t<var, T2>> arena_B = B;
    arena_t<ret_type> res = A.ldlt().solve(arena_B.val());

    reverse_pass_callback([A, arena_B, res]() mutable {
      promote_scalar_t<double, T2> adjB = A.ldlt().solve(res.adj());

      arena_B.adj() += adjB;
    });

    return ret_type(res);
  }
}

/**
 * Returns the solution of the system Ax=b given an LDLT_factor of A
 *
 * @tparam T type of B
 * @tparam R rowtype of A
 * @tparam C coltype of A
 * @param A LDLT_factor
 * @param B Right hand side matrix or vector.
 * @return x = A^-1 B, solution of the linear system.
 * @throws std::domain_error if rows of B don't match the size of A.
 */
/*template <typename T, int R, int C,
          require_var_matrix_t<T>* = nullptr>
inline auto mdivide_left_ldlt(const LDLT_factor<double, R, C>& A, const T& B) {
  using ret_val_type = Eigen::Matrix<double, R, T::ColsAtCompileTime>;
  using ret_type = var_value<ret_val_type>;

  check_multiplicable("mdivide_left_ldlt", "A", A, "B", B.val());

  if (A.rows() == 0) {
    return ret_type(ret_val_type(0, B.cols()));
  }

  arena_t<ret_type> res = A.solve(B.val());

  reverse_pass_callback([A, B, res]() mutable {
    B.adj() += A.solve(res.adj());
  });

  return ret_type(res);
  }*/

}  // namespace math
}  // namespace stan
#endif
