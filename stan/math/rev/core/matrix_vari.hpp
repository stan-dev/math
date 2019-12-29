#ifndef STAN_MATH_REV_CORE_MATRIX_VARI_HPP
#define STAN_MATH_REV_CORE_MATRIX_VARI_HPP

#include <stan/math/rev/mat/fun/Eigen_NumTraits.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/core/vari.hpp>
#include <stan/math/rev/mat/fun/typedefs.hpp>

namespace stan {
namespace math {

class op_matrix_vari : public vari {
 protected:
  const size_t size_;
  vari** vis_;

 public:
  template <typename Derived, require_t<is_var<typename Derived::Scalar>>...>
  op_matrix_vari(double f, const Eigen::MatrixBase<Derived>& vs)
      : vari(f), size_(vs.size()) {
    vis_ = ChainableStack::instance_->memalloc_.alloc_array<vari*>(size_);
    Eigen::Map<Eigen::Matrix<vari*, -1, -1>>(vis_, vs.rows(), vs.cols())
        = vs.vi();
  }
  vari* operator[](size_t n) const { return vis_[n]; }
  size_t size() { return size_; }
};

}  // namespace math
}  // namespace stan
#endif
