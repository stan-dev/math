#ifndef STAN_MATH_REV_CORE_PRECOMP_V_VARI_HPP
#define STAN_MATH_REV_CORE_PRECOMP_V_VARI_HPP

#include <stan/math/rev/core/vari.hpp>
#include <stan/math/rev/core/v_vari.hpp>

namespace stan {
namespace math {

// use for single precomputed partials
template <typename VariVal, typename Vari>
class precomp_v_vari : public op_vari<VariVal, Vari*> {
 protected:
  using Scalar = typename Vari::Scalar;
  using op_vari<VariVal, Vari*>::avi;
  Scalar da_;

 public:
  precomp_v_vari(Scalar val, Vari* avi, Scalar da)
      : op_vari<VariVal, Vari*>(val, avi), da_(da) {}
  void chain() { avi()->adj_ += this->adj_ * da_; }
};

}  // namespace math
}  // namespace stan
#endif
