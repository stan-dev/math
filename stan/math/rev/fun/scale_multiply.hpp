#ifndef STAN_MATH_REV_FUN_SCALE_MULTIPLY_HPP
#define STAN_MATH_REV_FUN_SCALE_MULTIPLY_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/to_arena.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Return the product of a scalar and two matrices, `s * A * B`, where at
 * least one argument has `var` scalars. The whole product is one node on
 * the reverse-mode stack.
 *
 * @tparam Scal type of the scalar
 * @tparam Mat1 type of the first matrix or expression
 * @tparam Mat2 type of the second matrix or expression
 * @param s scalar
 * @param A first matrix
 * @param B second matrix
 * @return `s * A * B`
 * @throw std::invalid_argument if the number of columns of `A` does not
 * match the number of rows of `B`
 */
template <typename Scal, typename Mat1, typename Mat2,
          require_stan_scalar_t<Scal>* = nullptr,
          require_all_eigen_t<Mat1, Mat2>* = nullptr,
          require_any_st_var<Scal, Mat1, Mat2>* = nullptr,
          require_all_not_complex_t<Scal, value_type_t<Mat1>,
                                    value_type_t<Mat2>>* = nullptr>
inline auto scale_multiply(Scal&& s, Mat1&& A, Mat2&& B) {
  check_multiplicable("scale_multiply", "A", A, "B", B);
  arena_t<Mat1> arena_A(std::forward<Mat1>(A));
  arena_t<Mat2> arena_B(std::forward<Mat2>(B));
  using return_t = return_var_matrix_t<decltype(arena_A * arena_B), Mat1, Mat2>;
  arena_t<return_t> res = value_of(s) * (value_of(arena_A) * value_of(arena_B));
  reverse_pass_callback([s, arena_A, arena_B, res]() mutable {
    // Eigen's matrix product cannot read the adjoint view directly
    const auto res_adj = res.adj().eval();
    if constexpr (is_autodiff_v<Scal>) {
      s.adj() += (value_of(arena_A).transpose() * res_adj)
                     .cwiseProduct(value_of(arena_B))
                     .sum();
    }
    if constexpr (is_autodiff_v<Mat1>) {
      arena_A.adj() += value_of(s) * res_adj * value_of(arena_B).transpose();
    }
    if constexpr (is_autodiff_v<Mat2>) {
      arena_B.adj() += value_of(s) * value_of(arena_A).transpose() * res_adj;
    }
  });
  return res;
}

}  // namespace math
}  // namespace stan

#endif
