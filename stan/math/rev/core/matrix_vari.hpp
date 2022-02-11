#ifndef STAN_MATH_REV_CORE_MATRIX_VARI_HPP
#define STAN_MATH_REV_CORE_MATRIX_VARI_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/core/Eigen_NumTraits.hpp>
#include <stan/math/rev/core/typedefs.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/core/vari.hpp>
#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

class op_matrix_vari : public vari {
 protected:
  const size_t size_;
  vari** vis_;

 public:
  template <typename T, require_eigen_vt<is_var, T>* = nullptr>
  op_matrix_vari(double f, const T& vs) : vari(f), size_(vs.size()) {
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
