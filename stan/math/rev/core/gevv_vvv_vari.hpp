#ifndef STAN_MATH_REV_CORE_GEVV_VVV_VARI_HPP
#define STAN_MATH_REV_CORE_GEVV_VVV_VARI_HPP

#include <stan/math/rev/core/vari.hpp>
#include <stan/math/rev/core/var.hpp>
#include <stan/math/rev/core/chainablestack.hpp>

namespace stan {
namespace math {

class gevv_vvv_vari : public vari {
 protected:
  vari* alpha_;
  vari** v1_;
  vari** v2_;
  double dotval_;
  size_t length_;
  inline static double eval_gevv(const var* alpha, const var* v1, int stride1,
                                 const var* v2, int stride2, size_t length,
                                 double* dotprod) {
    double result = 0;
    for (size_t i = 0; i < length; i++) {
      result += v1[i * stride1].vi_->val_ * v2[i * stride2].vi_->val_;
    }
    *dotprod = result;
    return alpha->vi_->val_ * result;
  }

 public:
  gevv_vvv_vari(const var* alpha, const var* v1, int stride1, const var* v2,
                int stride2, size_t length)
      : vari(eval_gevv(alpha, v1, stride1, v2, stride2, length, &dotval_)),
        length_(length) {
    alpha_ = alpha->vi_;
    v1_ = ChainableStack::instance_->memalloc_.alloc_array<vari*>(2 * length_);
    v2_ = v1_ + length_;
    for (size_t i = 0; i < length_; i++) {
      v1_[i] = v1[i * stride1].vi_;
    }
    for (size_t i = 0; i < length_; i++) {
      v2_[i] = v2[i * stride2].vi_;
    }
  }
  virtual ~gevv_vvv_vari() {}
  void chain() {
    const double adj_alpha = adj_ * alpha_->val_;
    for (size_t i = 0; i < length_; i++) {
      v1_[i]->adj_ += adj_alpha * v2_[i]->val_;
      v2_[i]->adj_ += adj_alpha * v1_[i]->val_;
    }
    alpha_->adj_ += adj_ * dotval_;
  }
};

class gevv_vvd_vari : public vari {
 protected:
  vari* alpha_;
  vari** v1_;
  double dotval_;
  size_t length_;
  inline static double eval_gevd(const var* alpha, const var* v1, int stride1,
                                 const double* v2, int stride2, size_t length,
                                 double* dotprod) {
    double result = 0;
    for (size_t i = 0; i < length; i++) {
      result += v1[i * stride1].vi_->val_ * v2[i * stride2];
    }
    *dotprod = result;
    return alpha->vi_->val_ * result;
  }

 public:
  gevv_vvd_vari(const var* alpha, const var* v1, int stride1, const double* v2,
                int stride2, size_t length)
      : vari(eval_gevd(alpha, v1, stride1, v2, stride2, length, &dotval_)),
        length_(length) {
    alpha_ = alpha->vi_;
    v1_ = ChainableStack::instance_->memalloc_.alloc_array<vari*>(length_);
    for (size_t i = 0; i < length_; i++) {
      v1_[i] = v1[i * stride1].vi_;
    }
  }
  virtual ~gevv_vvd_vari() {}
  void chain() {
    const double adj_alpha = adj_ * alpha_->val_;
    for (size_t i = 0; i < length_; i++) {
      v1_[i]->adj_ += adj_alpha;
    }
    alpha_->adj_ += adj_ * dotval_;
  }
};

class gevv_vdv_vari : public vari {
 protected:
  vari* alpha_;
  vari** v2_;
  double dotval_;
  size_t length_;
  inline static double eval_gedv(const var* alpha, const double* v1, int stride1,
                                 const var* v2, int stride2, size_t length,
                                 double* dotprod) {
    double result = 0;
    for (size_t i = 0; i < length; i++) {
      result += v1[i * stride1] * v2[i * stride2].vi_->val_;
    }
    *dotprod = result;
    return alpha->vi_->val_ * result;
  }

 public:
  gevv_vdv_vari(const var* alpha, const double* v1, int stride1, const var* v2,
                int stride2, size_t length)
      : vari(eval_gedv(alpha, v1, stride1, v2, stride2, length, &dotval_)),
        length_(length) {
    alpha_ = alpha->vi_;
    v2_ = ChainableStack::instance_->memalloc_.alloc_array<vari*>(length_);
    for (size_t i = 0; i < length_; i++) {
      v2_[i] = v2[i * stride2].vi_;
    }
  }
  virtual ~gevv_vdv_vari() {}
  void chain() {
    const double adj_alpha = adj_ * alpha_->val_;
    for (size_t i = 0; i < length_; i++) {
      v2_[i]->adj_ += adj_alpha;
    }
    alpha_->adj_ += adj_ * dotval_;
  }
};

}  // namespace math
}  // namespace stan
#endif
