#ifndef STAN_MATH_REV_FUNCTOR_REVERSE_PASS_CALLBACK_HPP
#define STAN_MATH_REV_FUNCTOR_REVERSE_PASS_CALLBACK_HPP

#include <stan/math/rev/core/vari.hpp>

namespace stan {
namespace math {
namespace internal {

template <typename F>
struct reverse_pass_callback_vari : public vari_base {
  F rev_functor_;

  explicit reverse_pass_callback_vari(F&& rev_functor)
      : rev_functor_(std::forward<F>(rev_functor)) {
    ChainableStack::instance_->var_stack_.push_back(this);
  }

  inline void chain() final { rev_functor_(); }
  inline void set_zero_adjoint() final {}
  inline void init_dependent() {}
};

}  // namespace internal

/**
 * Puts a callback on the autodiff stack to be called in reverse pass.
 *
 * The intended use case is for the callable to ba a lambda function that
 * captures any arguments it needs to work with. All captured values must be
 * trivially destructible or they will leak memory! `to_AD_stack()` function can
 * be used to ensure that.
 *
 * @tparam F type of callable
 * @param functor funtor or other callable to call in the reverse pass
 */
template <typename F>
inline void reverse_pass_callback(F&& functor) {
  new internal::reverse_pass_callback_vari<F>(std::forward<F>(functor));
}

}  // namespace math
}  // namespace stan

#endif
