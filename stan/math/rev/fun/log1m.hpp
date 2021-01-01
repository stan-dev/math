#ifndef STAN_MATH_REV_FUN_LOG1M_HPP
#define STAN_MATH_REV_FUN_LOG1M_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/log1m.hpp>

namespace stan {
namespace math {

namespace internal {
class log1m_vari : public op_v_vari {
 public:
  explicit log1m_vari(vari* avi) : op_v_vari(log1m(avi->val_), avi) {}
  void chain() { avi_->adj_ += adj_ / (avi_->val_ - 1); }
};
}  // namespace internal

/**
 * The log (1 - x) function for variables.
 *
 * The derivative is given by
 *
 * \f$\frac{d}{dx} \log (1 - x) = -\frac{1}{1 - x}\f$.
 *
 * @param a The variable.
 * @return The variable representing log of 1 minus the variable.
 */
inline var log1m(const var& a) { return var(new internal::log1m_vari(a.vi_)); }

/**
 * Return the elementwise log of 1 - x
 *
 * @tparam T type of x
 * @param x argument
 * @return elementwise log of 1 - x
 */
template <typename T, require_var_matrix_t<T>* = nullptr>
inline auto log1m(const T& x) {
  return make_callback_var(
      stan::math::log1m(x.val()), [x](const auto& vi) mutable {
        x.adj() += (vi.adj().array() / (x.val().array() - 1.0)).matrix();
      });
}

}  // namespace math
}  // namespace stan
#endif
