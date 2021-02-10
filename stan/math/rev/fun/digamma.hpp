#ifndef STAN_MATH_REV_FUN_DIGAMMA_HPP
#define STAN_MATH_REV_FUN_DIGAMMA_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/trigamma.hpp>
#include <stan/math/prim/fun/digamma.hpp>

namespace stan {
namespace math {

namespace internal {
class digamma_vari : public op_v_vari {
 public:
  explicit digamma_vari(vari* avi) : op_v_vari(digamma(avi->val_), avi) {}
  void chain() { avi_->adj_ += adj_ * trigamma(avi_->val_); }
};
}  // namespace internal

inline var digamma(const var& a) {
  return make_callback_var(digamma(a.val()), [a](auto& vi) {
    a.adj() += vi.adj() * trigamma(a.val());
  });
}

template <typename T, require_var_matrix_t<T>* = nullptr>
inline auto digamma(const T& a) {
  return make_callback_var(
      a.val()
          .array()
          .unaryExpr([](auto& x) { return digamma(x); })
          .matrix()
          .eval(),
      [a](auto& vi) mutable {
        a.adj().array()
            += vi.adj().array()
               * a.val().array().unaryExpr([](auto& x) { return trigamma(x); });
      });
}

}  // namespace math
}  // namespace stan
#endif
