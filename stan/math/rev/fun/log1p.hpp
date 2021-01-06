#ifndef STAN_MATH_REV_FUN_LOG1P_HPP
#define STAN_MATH_REV_FUN_LOG1P_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/log1p.hpp>

namespace stan {
namespace math {

namespace internal {
class log1p_vari : public op_v_vari {
 public:
  explicit log1p_vari(vari* avi) : op_v_vari(log1p(avi->val_), avi) {}
  void chain() { avi_->adj_ += adj_ / (1 + avi_->val_); }
};
}  // namespace internal

/**
 * The log (1 + x) function for variables (C99).
 *
 * The derivative is given by
 *
 * \f$\frac{d}{dx} \log (1 + x) = \frac{1}{1 + x}\f$.
 *
 * @param a The variable.
 * @return The log of 1 plus the variable.
 */
inline var log1p(const var& a) { return var(new internal::log1p_vari(a.vi_)); }

/**
 * Return the elementwise log of (1 + x)
 *
 * @tparam T type of input
 * @param x input
 * @return Elementwise log(1 + x)
 */
template <typename T, require_var_matrix_t<T>* = nullptr>
inline auto log1p(const T& x) {
  check_greater_or_equal("log1p", "x", x.val(), -1.0);
  return make_callback_var(
      x.val().array().log1p().matrix(), [x](const auto& vi) {
        x.adj().array() += vi.adj().array() / (1 + x.val().array());
      });
}

}  // namespace math
}  // namespace stan
#endif
