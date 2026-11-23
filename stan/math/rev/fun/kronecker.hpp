#ifndef STAN_MATH_REV_FUN_KRONECKER_HPP
#define STAN_MATH_REV_FUN_KRONECKER_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/kronecker.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * Return the Kronecker product of two matrices, at least one of which
 * holds `var`s.
 *
 * Column \f$c = q j + l\f$ of the adjoint \f$G\f$ of the result, read as a
 * \f$p \times m\f$ matrix \f$G_c\f$, gives
 * \f$\bar{A}_{\cdot j} \mathrel{+}= G_c^\top B_{\cdot l}\f$ and
 * \f$\bar{B}_{\cdot l} \mathrel{+}= G_c A_{\cdot j}\f$.
 *
 * @tparam T1 type of the first matrix
 * @tparam T2 type of the second matrix
 * @param A first matrix (m x n)
 * @param B second matrix (p x q)
 * @return Kronecker product of `A` and `B` (mp x nq)
 */
template <typename T1, typename T2, require_all_matrix_t<T1, T2>* = nullptr,
          require_return_type_t<is_var, T1, T2>* = nullptr>
inline auto kronecker(T1&& A, T2&& B) {
  using ret_type = return_var_matrix_t<Eigen::MatrixXd, T1, T2>;
  arena_t<T1> arena_A(std::forward<T1>(A));
  arena_t<T2> arena_B(std::forward<T2>(B));
  arena_t<ret_type> res = [&]() {
    if constexpr (is_var_matrix<ret_type>::value) {
      // hands the heap value to the chainable stack; var_value would copy it
      return arena_t<Eigen::MatrixXd>(
          kronecker(value_of(arena_A), value_of(arena_B)));
    } else {
      return kronecker(value_of(arena_A), value_of(arena_B));
    }
  }();
  reverse_pass_callback([arena_A, arena_B, res]() mutable {
    const Eigen::Index m = arena_A.rows();
    const Eigen::Index p = arena_B.rows();
    const Eigen::Index q = arena_B.cols();
    decltype(auto) res_adj = to_ref(res.adj_op());
    decltype(auto) A_val = to_ref(value_of(arena_A));
    decltype(auto) B_val = to_ref(value_of(arena_B));
    Eigen::MatrixXd A_adj = Eigen::MatrixXd::Zero(m, arena_A.cols());
    Eigen::MatrixXd B_adj = Eigen::MatrixXd::Zero(p, q);
    for (Eigen::Index c = 0; c < res_adj.cols(); ++c) {
      const auto res_adj_c = res_adj.col(c).reshaped(p, m);
      if constexpr (is_autodiff_v<T1>) {
        A_adj.col(c / q).noalias() += res_adj_c.transpose() * B_val.col(c % q);
      }
      if constexpr (is_autodiff_v<T2>) {
        B_adj.col(c % q).noalias() += res_adj_c * A_val.col(c / q);
      }
    }
    if constexpr (is_autodiff_v<T1>) {
      arena_A.adj() += A_adj;
    }
    if constexpr (is_autodiff_v<T2>) {
      arena_B.adj() += B_adj;
    }
  });
  return ret_type(res);
}

}  // namespace math
}  // namespace stan

#endif
