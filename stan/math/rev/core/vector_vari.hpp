#ifndef STAN_MATH_REV_CORE_VECTOR_VARI_HPP
#define STAN_MATH_REV_CORE_VECTOR_VARI_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/core/vari.hpp>
#include <vector>

namespace stan {
namespace math {

class op_vector_vari : public vari {
 protected:
  const size_t size_;
  vari** vis_;

 public:
  template <typename Arith, typename VecVar,
            require_arithmetic_t<Arith>* = nullptr,
            require_vector_like_vt<is_var, VecVar>* = nullptr>
  op_vector_vari(Arith f, VecVar&& vs) : vari(f), size_(vs.size()) {
    vis_ = reinterpret_cast<vari**>(operator new(sizeof(vari*) * vs.size()));
    for (size_t i = 0; i < vs.size(); ++i) {
      vis_[i] = vs[i].vi_;
    }
  }
  vari* operator[](size_t n) const { return vis_[n]; }
  size_t size() { return size_; }
};

}  // namespace math
}  // namespace stan
#endif
