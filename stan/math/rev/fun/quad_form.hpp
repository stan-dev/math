#ifndef STAN_MATH_REV_FUN_QUAD_FORM_HPP
#define STAN_MATH_REV_FUN_QUAD_FORM_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/to_var_value.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/quad_form.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <type_traits>

namespace stan {
namespace math {

/**
 * Return the quadratic form \f$ B^T A B \f$.
 *
 * Symmetry of the resulting matrix is not guaranteed due to numerical
 * precision.
 *
 * @tparam Mat1 type of the first (square) matrix
 * @tparam Mat2 type of the second matrix
 *
 * @param A square matrix
 * @param B second matrix
 * @param symmetric indicates whether the output should be made symmetric
 * @return The quadratic form
 * @throws std::invalid_argument if A is not square, or if A cannot be
 * multiplied by B
 */
template <typename Mat1, typename Mat2,
          require_any_rev_matrix_t<Mat1, Mat2>* = nullptr>
inline auto quad_form(const Mat1& A, const Mat2& B, bool symmetric = false) {
  check_square("quad_form", "A", A);
  check_multiplicable("quad_form", "A", A, "B", B);

  using val_return_t = plain_type_t<decltype(value_of(B).transpose().eval()
                                  * value_of(A) * value_of(B).eval())>;
  using return_t = promote_var_matrix_t<val_return_t, Mat1, Mat2>;

  if (!is_constant<Mat1>::value && !is_constant<Mat2>::value) {
    arena_t<promote_scalar_t<var, Mat1>> arena_A(A);
    arena_t<promote_scalar_t<var, Mat2>> arena_B(B);

    check_not_nan("multiply", "A", arena_A.val());
    check_not_nan("multiply", "B", arena_B.val());

    arena_t<val_return_t> res_val(arena_B.val_op().transpose() * arena_A.val_op() * arena_B.val_op());
    if (symmetric) {
      res_val += (res_val.transpose()).eval();
    }

    arena_t<return_t> res(res_val);
    reverse_pass_callback([arena_A, arena_B, res]() mutable {
      auto C_adj_B_t = (res.adj_op() * arena_B.val_op().transpose()).eval();
      if (is_var_matrix<Mat1>::value) {
        arena_A.adj().noalias() += arena_B.val_op() * C_adj_B_t;
      } else {
        arena_A.adj() += arena_B.val_op() * C_adj_B_t;
      }
      if (is_var_matrix<Mat2>::value) {
        arena_B.adj().noalias() += arena_A.val_op() * C_adj_B_t.transpose()
               + arena_A.val_op().transpose() * arena_B.val_op() * res.adj_op();
      } else {
        arena_B.adj() += arena_A.val_op() * C_adj_B_t.transpose()
               + arena_A.val_op().transpose() * arena_B.val_op() * res.adj_op();
      }
    });
    return return_t(res);
  } else if (!is_constant<Mat2>::value) {
    arena_t<promote_scalar_t<double, Mat1>> arena_A(value_of(A));
    arena_t<promote_scalar_t<var, Mat2>> arena_B(B);

    check_not_nan("multiply", "A", arena_A);
    check_not_nan("multiply", "B", arena_B.val());
    arena_t<val_return_t> res_val(arena_B.val_op().transpose() * arena_A * arena_B.val_op());
    if (symmetric) {
      res_val += (res_val.transpose()).eval();
    }

    arena_t<return_t> res(res_val);
    reverse_pass_callback([arena_A, arena_B, res]() mutable {
      if (is_var_matrix<Mat2>::value) {
        arena_B.adj().noalias() += arena_A * (res.adj_op() * arena_B.val_op().transpose()).transpose()
               + arena_A.transpose() * arena_B.val_op() * res.adj_op();
      } else {
        arena_B.adj() += arena_A * (res.adj_op() * arena_B.val_op().transpose()).transpose()
               + arena_A.transpose() * arena_B.val_op() * res.adj_op();
      }
    });
    return return_t(res);
  } else if (!is_constant<Mat1>::value) {
    arena_t<promote_scalar_t<var, Mat1>> arena_A(A);
    arena_t<promote_scalar_t<double, Mat2>> arena_B(value_of(B));

    check_not_nan("multiply", "A", arena_A.val());
    check_not_nan("multiply", "B", arena_B);

    arena_t<val_return_t> res_val(arena_B.transpose() * arena_A.val_op() * arena_B);
    if (symmetric) {
      res_val += (res_val.transpose()).eval();
    }

    arena_t<return_t> res(res_val);
    reverse_pass_callback([arena_A, arena_B, res]() mutable {
      if (is_var_matrix<Mat1>::value) {
        arena_A.adj().noalias() += arena_B * (res.adj_op() * arena_B.transpose());
      } else {
        arena_A.adj() += arena_B * (res.adj_op() * arena_B.transpose());
      }
    });
    return return_t(res);
  }

}


}  // namespace math
}  // namespace stan
#endif
