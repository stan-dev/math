#ifndef STAN_MATH_FWD_FUN_EXP_HPP
#define STAN_MATH_FWD_FUN_EXP_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <cmath>

namespace stan {
namespace math {
/**
 * Return the exponentiation of the specified variable (cmath).
 *
 * The derivative is defined by
 *
   \f[
   \frac{\partial\, \mbox{exp}(x)}{\partial x} =
   \begin{cases}
     e^x & \mbox{if } -\infty\leq x\leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @tparam T Inner type of the fvar
 * @param x fvar<T> variable
 * @return Exponentiated variable.
 */
template <typename T>
inline fvar<T> exp(const fvar<T>& x) {
  using std::exp;
  T exp_val = exp(x.val_);
  return fvar<T>(exp_val, x.d_ * exp_val);
}

/**
 * Return the exponential of each variable in a container.
 *
 * @tparam Container Type of container
 * @param x Container of fvar
 * @return Exponential of each variable in container.
 */
template <typename Container,
          require_container_st<is_container, is_fvar, Container>...>
inline auto exp(const Container& x) {
  return apply_vector_unary<Container>::apply(
      x, [](const auto& v) {
        using T_plain = plain_type_t<decltype(v)>;
        const Eigen::Ref<const T_plain>& v_ref = v;

        T_plain result(v_ref.rows(), v_ref.cols());
        result.val() = exp(v_ref.val());
        result.d() = v_ref.d().cwiseProduct(result.val());

        return result;
});
}

}  // namespace math
}  // namespace stan
#endif
