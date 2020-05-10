#ifndef STAN_MATH_REV_CORE_PRECOMP_VV_VARI_HPP
#define STAN_MATH_REV_CORE_PRECOMP_VV_VARI_HPP

#include <stan/math/rev/core/vari.hpp>
#include <stan/math/rev/core/op_vari.hpp>

namespace stan {
namespace math {

// use for single precomputed partials
template <typename VariVal, typename Vari1, typename Vari2>
class precomp_vv_vari : public op_vari<VariVal, Vari1*, Vari2*> {
  using op_vari<VariVal, Vari1*, Vari2*>::avi;
  using op_vari<VariVal, Vari1*, Vari2*>::bvi;
  using Scalar1 = typename Vari1::Scalar;
  using Scalar2 = typename Vari2::Scalar;

 protected:
  Scalar1 da_;
  Scalar2 db_;

 public:
  precomp_vv_vari(VariVal val, Vari1* avi, Vari2* bvi, Scalar1 da, Scalar2 db)
      : op_vari<VariVal, Vari1*, Vari2*>(val, avi, bvi), da_(da), db_(db) {}
  void chain() {
    avi()->adj_ += this->adj_ * da_;
    bvi()->adj_ += this->adj_ * db_;
  }
};

}  // namespace math
}  // namespace stan
#endif
