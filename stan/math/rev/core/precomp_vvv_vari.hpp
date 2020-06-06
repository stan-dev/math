#ifndef STAN_MATH_REV_CORE_PRECOMP_VVV_VARI_HPP
#define STAN_MATH_REV_CORE_PRECOMP_VVV_VARI_HPP

#include <stan/math/rev/core/vari.hpp>
#include <stan/math/rev/core/op_vari.hpp>

namespace stan {
namespace math {

// use for single precomputed partials
template <typename VariVal, typename Vari1, typename Vari2, typename Vari3>
class precomp_vvv_vari final : public op_vari<VariVal, Vari1*, Vari2*, Vari3*> {
  using op_vari<VariVal, Vari1*, Vari2*, Vari3*>::avi;
  using op_vari<VariVal, Vari1*, Vari2*, Vari3*>::bvi;
  using op_vari<VariVal, Vari1*, Vari2*, Vari3*>::cvi;
  using Scalar1 = typename Vari1::Scalar;
  using Scalar2 = typename Vari2::Scalar;
  using Scalar3 = typename Vari3::Scalar;

 protected:
  Scalar1 da_;
  Scalar2 db_;
  Scalar3 dc_;

 public:
  precomp_vvv_vari(VariVal val, Vari1* avi, Vari2* bvi, Vari3* cvi, Scalar1 da,
                   Scalar2 db, Scalar3 dc)
      : op_vari<VariVal, Vari1*, Vari2*, Vari3*>(val, avi, bvi, cvi),
        da_(da),
        db_(db),
        dc_(dc) {}
  void chain() {
    avi()->adj_ += this->adj_ * da_;
    bvi()->adj_ += this->adj_ * db_;
    cvi()->adj_ += this->adj_ * dc_;
  }
};

}  // namespace math
}  // namespace stan
#endif
