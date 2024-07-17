#ifndef STAN_MATH_REV_FUN_ATAN2_HPP
#define STAN_MATH_REV_FUN_ATAN2_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the principal value of the arc tangent, in radians, of
 * the first variable divided by the second (cmath).
 *
 * The partial derivatives are defined by
 *
 * \f$ \frac{\partial}{\partial x} \arctan \frac{x}{y} = \frac{y}{x^2 + y^2}\f$,
 * and
 *
 * \f$ \frac{\partial}{\partial y} \arctan \frac{x}{y} = \frac{-x}{x^2 +
 * y^2}\f$.
 *
 * @param a Numerator variable.
 * @param b Denominator variable.
 * @return The arc tangent of the fraction, in radians.
 */
inline var atan2(const var& a, const var& b) {
  return make_callback_var(
      std::atan2(a.val(), b.val()), [a, b](const auto& vi) mutable {
        double a_sq_plus_b_sq = (a.val() * a.val()) + (b.val() * b.val());
        a.adj() += vi.adj_ * b.val() / a_sq_plus_b_sq;
        b.adj() += -vi.adj_ * a.val() / a_sq_plus_b_sq;
      });
}

/**
 * Return the principal value of the arc tangent, in radians, of
 * the first variable divided by the second scalar (cmath).
 *
 * The derivative with respect to the variable is
 *
 * \f$ \frac{d}{d x} \arctan \frac{x}{c} = \frac{c}{x^2 + c^2}\f$.
 *
 * @param a Numerator variable.
 * @param b Denominator scalar.
 * @return The arc tangent of the fraction, in radians.
 */
inline var atan2(const var& a, double b) {
  return make_callback_var(
      std::atan2(a.val(), b), [a, b](const auto& vi) mutable {
        double a_sq_plus_b_sq = (a.val() * a.val()) + (b * b);
        a.adj() += vi.adj_ * b / a_sq_plus_b_sq;
      });
}

/**
 * Return the principal value of the arc tangent, in radians, of
 * the first scalar divided by the second variable (cmath).
 *
 * The derivative with respect to the variable is
 *
 * \f$ \frac{\partial}{\partial y} \arctan \frac{c}{y} = \frac{-c}{c^2 +
 y^2}\f$.
 *
 *
   \f[
   \mbox{atan2}(x, y) =
   \begin{cases}
     \arctan\left(\frac{x}{y}\right) & \mbox{if } -\infty\leq x \leq \infty,
 -\infty\leq y \leq \infty \\[6pt] \textrm{NaN} & \mbox{if } x = \textrm{NaN or
 } y = \textrm{NaN} \end{cases} \f]

   \f[
   \frac{\partial\, \mbox{atan2}(x, y)}{\partial x} =
   \begin{cases}
     \frac{y}{x^2+y^2} & \mbox{if } -\infty\leq x\leq \infty, -\infty\leq y \leq
 \infty \\[6pt] \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{atan2}(x, y)}{\partial y} =
   \begin{cases}
     -\frac{x}{x^2+y^2} & \mbox{if } -\infty\leq x\leq \infty, -\infty\leq y
 \leq \infty \\[6pt] \textrm{NaN} & \mbox{if } x = \textrm{NaN or } y =
 \textrm{NaN} \end{cases} \f]
 *
 * @param a Numerator scalar.
 * @param b Denominator variable.
 * @return The arc tangent of the fraction, in radians.
 */
inline var atan2(double a, const var& b) {
  return make_callback_var(
      std::atan2(a, b.val()), [a, b](const auto& vi) mutable {
        double a_sq_plus_b_sq = (a * a) + (b.val() * b.val());
        b.adj() += -vi.adj_ * a / a_sq_plus_b_sq;
      });
}

template <typename Mat1, typename Mat2,
          require_any_var_matrix_t<Mat1, Mat2>* = nullptr,
          require_all_matrix_t<Mat1, Mat2>* = nullptr>
inline auto atan2(Mat1&& a, Mat2&& b) {
  arena_t<Mat1> arena_a = std::forward<Mat1>(a);
  arena_t<Mat2> arena_b = std::forward<Mat2>(b);
    auto a_sq_plus_b_sq
        = to_arena(value_of(arena_a).array().square()
                   + value_of(arena_b).array().square());
    return make_callback_var(
        atan2(arena_a.val(), arena_b.val()),
        [arena_a, arena_b, a_sq_plus_b_sq](auto& vi) mutable {
          if constexpr (is_autodiffable_v<Mat1>) {
            arena_a.adj().array()
                += vi.adj().array() * value_of(arena_b).array() / a_sq_plus_b_sq;
          }
          if constexpr (is_autodiffable_v<Mat2>) {
            arena_b.adj().array()
                += -vi.adj().array() * value_of(arena_a).array() / a_sq_plus_b_sq;
          }
        });
}

template <typename Scalar, typename VarMat,
          require_var_matrix_t<VarMat>* = nullptr,
          require_stan_scalar_t<Scalar>* = nullptr>
inline auto atan2(Scalar a, VarMat&& b) {
  arena_t<VarMat> arena_b = std::forward<VarMat>(b);
  auto a_sq_plus_b_sq = to_arena(
      square(value_of(a)) + (value_of(arena_b).array().square()));
    return make_callback_var(
        atan2(value_of(a), value_of(arena_b)),
        [a, arena_b, a_sq_plus_b_sq](auto& vi) mutable {
          if constexpr (is_autodiffable_v<Scalar>) {
            a.adj() += (vi.adj().array() * value_of(arena_b).array() / a_sq_plus_b_sq)
                          .sum();
          }
          if constexpr (is_autodiffable_v<VarMat>) {
            arena_b.adj().array() += -vi.adj().array() * value_of(a) / a_sq_plus_b_sq;
          }
        });
}

template <typename VarMat, typename Scalar,
          require_var_matrix_t<VarMat>* = nullptr,
          require_stan_scalar_t<Scalar>* = nullptr>
inline auto atan2(VarMat&& a, Scalar b) {
  arena_t<VarMat> arena_a = std::forward<VarMat>(a);
  auto a_sq_plus_b_sq = to_arena(value_of(arena_a).array().square() + square(value_of(b)));
    return make_callback_var(
        atan2(value_of(arena_a), value_of(b)),
        [arena_a, b, a_sq_plus_b_sq](auto& vi) mutable {
        if constexpr (is_autodiffable_v<VarMat>) {
          arena_a.adj().array() += vi.adj().array() * value_of(b) / a_sq_plus_b_sq;
        }
        if constexpr (is_autodiffable_v<Scalar>) {
          b.adj()
              += -(vi.adj().array() * value_of(arena_a).array() / a_sq_plus_b_sq)
                      .sum();
        }
        });
}

}  // namespace math
}  // namespace stan
#endif
