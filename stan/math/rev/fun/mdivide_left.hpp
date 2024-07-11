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
inline auto mdivide_left(T1&& A, T2&& B) {
  using ret_val_type = plain_type_t<decltype(value_of(A) * value_of(B))>;
  using ret_type = promote_var_matrix_t<ret_val_type, T1, T2>;

  check_square("mdivide_left", "A", A);
  check_multiplicable("mdivide_left", "A", A, "B", B);

  if (A.size() == 0) {
    return arena_t<ret_type>(ret_val_type(0, B.cols()));
  }

  if constexpr (!is_constant_v<T1> && !is_constant_v<T2>) {
    arena_t<T1> arena_A = std::forward<T1>(A);
    arena_t<T2> arena_B = std::forward<T2>(B);

    auto hqr_A_ptr = make_chainable_ptr(arena_A.val().householderQr());
    arena_t<ret_type> res = hqr_A_ptr->solve(arena_B.val());
    reverse_pass_callback([arena_A, arena_B, hqr_A_ptr, res]() mutable {
      using T2_t = std::decay_t<T2>;
      arena_t<Eigen::Matrix<double, T2_t::RowsAtCompileTime, T2_t::ColsAtCompileTime>> adjB
          = hqr_A_ptr->householderQ()
            * hqr_A_ptr->matrixQR()
                  .template triangularView<Eigen::Upper>()
                  .transpose()
                  .solve(res.adj());
      arena_A.adj() -= adjB * res.val_op().transpose();
      arena_B.adj() += adjB;
    });

    return res;
  } else if constexpr (!is_constant_v<T2>) {
    arena_t<T2> arena_B = std::forward<T2>(B);

    auto hqr_A_ptr = make_chainable_ptr(value_of(A).householderQr());
    arena_t<ret_type> res = hqr_A_ptr->solve(arena_B.val());
    reverse_pass_callback([arena_B, hqr_A_ptr, res]() mutable {
      arena_B.adj() += hqr_A_ptr->householderQ()
                       * hqr_A_ptr->matrixQR()
                             .template triangularView<Eigen::Upper>()
                             .transpose()
                             .solve(res.adj());
    });
    return res;
  } else {
    arena_t<T1> arena_A = std::forward<T1>(A);

    auto hqr_A_ptr = make_chainable_ptr(arena_A.val().householderQr());
    arena_t<ret_type> res = hqr_A_ptr->solve(B);
    reverse_pass_callback([arena_A, hqr_A_ptr, res]() mutable {
      arena_A.adj() -= hqr_A_ptr->householderQ()
                       * hqr_A_ptr->matrixQR()
                             .template triangularView<Eigen::Upper>()
                             .transpose()
                             .solve(res.adj())
                       * res.val_op().transpose();
    });
    return res;
  }
}

}  // namespace math
}  // namespace stan
#endif
