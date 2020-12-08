#ifndef STAN_MATH_REV_FUN_MDIVIDE_LEFT_HPP
#define STAN_MATH_REV_FUN_MDIVIDE_LEFT_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/core/typedefs.hpp>
#include <stan/math/rev/core/chainable_object.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Return the solution `X` of `AX = B`.
 *
 * A must be a square matrix, but B can be a matrix or a vector
 *
 * @tparam T1 type of first matrix
 * @tparam T2 type of second matrix
 *
 * @param[in] A square matrix
 * @param[in] B right hand side
 * @return solution of AX = B
 */
template <typename T1, typename T2, require_all_matrix_t<T1, T2>* = nullptr,
          require_any_st_var<T1, T2>* = nullptr>
inline auto mdivide_left(const T1& A, const T2& B) {
  using ret_val_type = plain_type_t<decltype(value_of(A) * value_of(B))>;
  using ret_type = promote_var_matrix_t<ret_val_type, T1, T2>;

  check_square("mdivide_left", "A", A);
  check_multiplicable("mdivide_left", "A", A, "B", B);

  if (A.size() == 0) {
    return ret_type(ret_val_type(0, B.cols()));
  }

  if (!is_constant<T1>::value && !is_constant<T2>::value) {
    arena_t<promote_scalar_t<var, T1>> arena_A = A;
    arena_t<promote_scalar_t<var, T2>> arena_B = B;

    auto hqr_A_ptr = make_chainable_ptr(arena_A.val().householderQr());
    arena_t<ret_type> res = hqr_A_ptr->solve(arena_B.val());
    reverse_pass_callback([arena_A, arena_B, hqr_A_ptr, res]() mutable {
      promote_scalar_t<double, T2> adjB
          = hqr_A_ptr->householderQ()
            * hqr_A_ptr->matrixQR()
                  .template triangularView<Eigen::Upper>()
                  .transpose()
                  .solve(res.adj());
      arena_A.adj() -= adjB * res.val_op().transpose();
      arena_B.adj() += adjB;
    });

    return ret_type(res);
  } else if (!is_constant<T2>::value) {
    arena_t<promote_scalar_t<var, T2>> arena_B = B;

    auto hqr_A_ptr = make_chainable_ptr(value_of(A).householderQr());
    arena_t<ret_type> res = hqr_A_ptr->solve(arena_B.val());
    reverse_pass_callback([arena_B, hqr_A_ptr, res]() mutable {
      arena_B.adj() += hqr_A_ptr->householderQ()
                       * hqr_A_ptr->matrixQR()
                             .template triangularView<Eigen::Upper>()
                             .transpose()
                             .solve(res.adj());
    });
    return ret_type(res);
  } else {
    arena_t<promote_scalar_t<var, T1>> arena_A = A;

    auto hqr_A_ptr = make_chainable_ptr(arena_A.val().householderQr());
    arena_t<ret_type> res = hqr_A_ptr->solve(value_of(B));
    reverse_pass_callback([arena_A, hqr_A_ptr, res]() mutable {
      arena_A.adj() -= hqr_A_ptr->householderQ()
                       * hqr_A_ptr->matrixQR()
                             .template triangularView<Eigen::Upper>()
                             .transpose()
                             .solve(res.adj())
                       * res.val_op().transpose();
    });
    return ret_type(res);
  }
}

}  // namespace math
}  // namespace stan
#endif
