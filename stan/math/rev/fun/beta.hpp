#ifndef STAN_MATH_REV_FUN_BETA_HPP
#define STAN_MATH_REV_FUN_BETA_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/beta.hpp>
#include <stan/math/prim/fun/digamma.hpp>

namespace stan {
namespace math {

/**
 * Returns the beta function and gradients for two var inputs.
 *
   \f[
     \mathrm{beta}(a,b) = \left(B\left(a,b\right)\right)
   \f]

   \f[
    \frac{\partial }{\partial a} = \left(\psi^{\left(0\right)}\left(a\right)
                                      - \psi^{\left(0\right)}
                                      \left(a + b\right)\right)
                                    * \mathrm{beta}(a,b)
   \f]

   \f[
    \frac{\partial }{\partial b} = \left(\psi^{\left(0\right)}\left(b\right)
                                      - \psi^{\left(0\right)}
                                      \left(a + b\right)\right)
                                    * \mathrm{beta}(a,b)
   \f]
 *
 * @param a var Argument
 * @param b var Argument
 * @return Result of beta function
 */
inline var beta(const var& a, const var& b) {
  double digamma_ab = digamma(a.val() + b.val());
  double digamma_a = digamma(a.val()) - digamma_ab;
  double digamma_b = digamma(b.val()) - digamma_ab;
  return make_callback_var(beta(a.val(), b.val()),
                           [a, b, digamma_a, digamma_b](auto& vi) mutable {
                             const double adj_val = vi.adj() * vi.val();
                             a.adj() += adj_val * digamma_a;
                             b.adj() += adj_val * digamma_b;
                           });
}

/**
 * Returns the beta function and gradient for first var input.
 *
   \f[
     \mathrm{beta}(a,b) = \left(B\left(a,b\right)\right)
   \f]

   \f[
    \frac{\partial }{\partial a} = \left(\psi^{\left(0\right)}\left(a\right)
                                      - \psi^{\left(0\right)}
                                      \left(a + b\right)\right)
                                    * \mathrm{beta}(a,b)
   \f]
 *
 * @param a var Argument
 * @param b double Argument
 * @return Result of beta function
 */
inline var beta(const var& a, double b) {
  auto digamma_ab = digamma(a.val()) - digamma(a.val() + b);
  return make_callback_var(beta(a.val(), b),
                           [a, b, digamma_ab](auto& vi) mutable {
                             a.adj() += vi.adj() * digamma_ab * vi.val();
                           });
}

/**
 * Returns the beta function and gradient for second var input.
 *
   \f[
     \mathrm{beta}(a,b) = \left(B\left(a,b\right)\right)
   \f]

   \f[
    \frac{\partial }{\partial b} = \left(\psi^{\left(0\right)}\left(b\right)
                                      - \psi^{\left(0\right)}
                                      \left(a + b\right)\right)
                                    * \mathrm{beta}(a,b)
   \f]
 *
 * @param a double Argument
 * @param b var Argument
 * @return Result of beta function
 */
inline var beta(double a, const var& b) {
  auto beta_val = beta(a, b.val());
  auto digamma_ab = (digamma(b.val()) - digamma(a + b.val())) * beta_val;
  return make_callback_var(beta_val, [a, b, digamma_ab](auto& vi) mutable {
    b.adj() += vi.adj() * digamma_ab;
  });
}

template <typename Mat1, typename Mat2,
          require_any_var_matrix_t<Mat1, Mat2>* = nullptr,
          require_all_matrix_t<Mat1, Mat2>* = nullptr>
inline auto beta(const Mat1& a, const Mat2& b) {
  if (!is_constant<Mat1>::value && !is_constant<Mat2>::value) {
    arena_t<promote_scalar_t<var, Mat1>> arena_a = a;
    arena_t<promote_scalar_t<var, Mat2>> arena_b = b;
    auto beta_val = beta(arena_a.val(), arena_b.val());
    auto digamma_ab
        = to_arena(digamma(arena_a.val().array() + arena_b.val().array()));
    return make_callback_var(
        beta(arena_a.val(), arena_b.val()),
        [arena_a, arena_b, digamma_ab](auto& vi) mutable {
          const auto adj_val = (vi.adj().array() * vi.val().array()).eval();
          arena_a.adj().array()
              += adj_val * (digamma(arena_a.val().array()) - digamma_ab);
          arena_b.adj().array()
              += adj_val * (digamma(arena_b.val().array()) - digamma_ab);
        });
  } else if (!is_constant<Mat1>::value) {
    arena_t<promote_scalar_t<var, Mat1>> arena_a = a;
    arena_t<promote_scalar_t<double, Mat2>> arena_b = value_of(b);
    auto digamma_ab
        = to_arena(digamma(arena_a.val()).array()
                   - digamma(arena_a.val().array() + arena_b.array()));
    return make_callback_var(beta(arena_a.val(), arena_b),
                             [arena_a, arena_b, digamma_ab](auto& vi) mutable {
                               arena_a.adj().array() += vi.adj().array()
                                                        * digamma_ab
                                                        * vi.val().array();
                             });
  } else if (!is_constant<Mat2>::value) {
    arena_t<promote_scalar_t<double, Mat1>> arena_a = value_of(a);
    arena_t<promote_scalar_t<var, Mat2>> arena_b = b;
    auto beta_val = beta(arena_a, arena_b.val());
    auto digamma_ab
        = to_arena((digamma(arena_b.val()).array()
                    - digamma(arena_a.array() + arena_b.val().array()))
                   * beta_val.array());
    return make_callback_var(
        beta_val, [arena_a, arena_b, digamma_ab](auto& vi) mutable {
          arena_b.adj().array() += vi.adj().array() * digamma_ab.array();
        });
  }
}

template <typename Scalar, typename VarMat,
          require_var_matrix_t<VarMat>* = nullptr,
          require_stan_scalar_t<Scalar>* = nullptr>
inline auto beta(const Scalar& a, const VarMat& b) {
  if (!is_constant<Scalar>::value && !is_constant<VarMat>::value) {
    var arena_a = a;
    arena_t<promote_scalar_t<var, VarMat>> arena_b = b;
    auto beta_val = beta(arena_a.val(), arena_b.val());
    auto digamma_ab = to_arena(digamma(arena_a.val() + arena_b.val().array()));
    return make_callback_var(
        beta(arena_a.val(), arena_b.val()),
        [arena_a, arena_b, digamma_ab](auto& vi) mutable {
          const auto adj_val = (vi.adj().array() * vi.val().array()).eval();
          arena_a.adj()
              += (adj_val * (digamma(arena_a.val()) - digamma_ab)).sum();
          arena_b.adj().array()
              += adj_val * (digamma(arena_b.val().array()) - digamma_ab);
        });
  } else if (!is_constant<Scalar>::value) {
    var arena_a = a;
    arena_t<promote_scalar_t<double, VarMat>> arena_b = value_of(b);
    auto digamma_ab = to_arena(digamma(arena_a.val())
                               - digamma(arena_a.val() + arena_b.array()));
    return make_callback_var(
        beta(arena_a.val(), arena_b),
        [arena_a, arena_b, digamma_ab](auto& vi) mutable {
          arena_a.adj()
              += (vi.adj().array() * digamma_ab * vi.val().array()).sum();
        });
  } else if (!is_constant<VarMat>::value) {
    double arena_a = value_of(a);
    arena_t<promote_scalar_t<var, VarMat>> arena_b = b;
    auto beta_val = beta(arena_a, arena_b.val());
    auto digamma_ab = to_arena((digamma(arena_b.val()).array()
                                - digamma(arena_a + arena_b.val().array()))
                               * beta_val.array());
    return make_callback_var(
        beta_val, [arena_a, arena_b, digamma_ab](auto& vi) mutable {
          arena_b.adj().array() += vi.adj().array() * digamma_ab.array();
        });
  }
}

template <typename VarMat, typename Scalar,
          require_var_matrix_t<VarMat>* = nullptr,
          require_stan_scalar_t<Scalar>* = nullptr>
inline auto beta(const VarMat& a, const Scalar& b) {
  if (!is_constant<VarMat>::value && !is_constant<Scalar>::value) {
    arena_t<promote_scalar_t<var, VarMat>> arena_a = a;
    var arena_b = b;
    auto beta_val = beta(arena_a.val(), arena_b.val());
    auto digamma_ab = to_arena(digamma(arena_a.val().array() + arena_b.val()));
    return make_callback_var(
        beta(arena_a.val(), arena_b.val()),
        [arena_a, arena_b, digamma_ab](auto& vi) mutable {
          const auto adj_val = (vi.adj().array() * vi.val().array()).eval();
          arena_a.adj().array()
              += adj_val * (digamma(arena_a.val().array()) - digamma_ab);
          arena_b.adj()
              += (adj_val * (digamma(arena_b.val()) - digamma_ab)).sum();
        });
  } else if (!is_constant<VarMat>::value) {
    arena_t<promote_scalar_t<var, VarMat>> arena_a = a;
    double arena_b = value_of(b);
    auto digamma_ab = to_arena(digamma(arena_a.val()).array()
                               - digamma(arena_a.val().array() + arena_b));
    return make_callback_var(beta(arena_a.val(), arena_b),
                             [arena_a, arena_b, digamma_ab](auto& vi) mutable {
                               arena_a.adj().array() += vi.adj().array()
                                                        * digamma_ab
                                                        * vi.val().array();
                             });
  } else if (!is_constant<Scalar>::value) {
    arena_t<promote_scalar_t<double, VarMat>> arena_a = value_of(a);
    var arena_b = b;
    auto beta_val = beta(arena_a, arena_b.val());
    auto digamma_ab = to_arena(
        (digamma(arena_b.val()) - digamma(arena_a.array() + arena_b.val()))
        * beta_val.array());
    return make_callback_var(
        beta_val, [arena_a, arena_b, digamma_ab](auto& vi) mutable {
          arena_b.adj() += (vi.adj().array() * digamma_ab.array()).sum();
        });
  }
}

}  // namespace math
}  // namespace stan
#endif
