#ifndef STAN_MATH_PRIM_FUN_SCALE_MULTIPLY_HPP
#define STAN_MATH_PRIM_FUN_SCALE_MULTIPLY_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Return the product of a scalar and two matrices, `s * A * B`.
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
          require_all_st_arithmetic<Scal, Mat1, Mat2>* = nullptr>
inline auto scale_multiply(Scal&& s, Mat1&& A, Mat2&& B) {
  check_multiplicable("scale_multiply", "A", A, "B", B);
  return make_holder(
      [s](auto&& A_, auto&& B_) {
        return s
               * (std::forward<decltype(A_)>(A_)
                  * std::forward<decltype(B_)>(B_));
      },
      std::forward<Mat1>(A), std::forward<Mat2>(B));
}

}  // namespace math
}  // namespace stan

#endif
