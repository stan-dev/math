#ifndef STAN_MATH_REV_FUN_INC_BETA_HPP
#define STAN_MATH_REV_FUN_INC_BETA_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/scal/fun/inc_beta_dda.hpp>
#include <stan/math/prim/scal/fun/inc_beta_ddb.hpp>
#include <stan/math/prim/scal/fun/inc_beta_ddz.hpp>
#include <stan/math/prim/scal/fun/digamma.hpp>

namespace stan {
namespace math {

namespace internal {

class inc_beta_vvv_vari : public op_vvv_vari {
 public:
  inc_beta_vvv_vari(vari* avi, vari* bvi, vari* cvi)
      : op_vvv_vari(inc_beta(avi->val_, bvi->val_, cvi->val_), avi, bvi, cvi) {}
  void chain() {
    const double& ad = avi_->val_;
    const double& bd = bvi_->val_;
    const double& cd = cvi_->val_;
    const double digamma_abd = digamma(ad + bd);
    avi_->adj_ += adj_ * inc_beta_dda(ad, bd, cd, digamma(ad), digamma_abd);
    bvi_->adj_ += adj_ * inc_beta_ddb(ad, bd, cd, digamma(bd), digamma_abd);
    cvi_->adj_ += adj_ * inc_beta_ddz(ad, bd, cd);
  }
};

class inc_beta_dvv_vari : public op_dvv_vari {
 public:
  inc_beta_dvv_vari(double ad, vari* bvi, vari* cvi)
      : op_dvv_vari(inc_beta(ad, bvi->val_, cvi->val_), ad, bvi, cvi) {}
  void chain() {
    const double& bd = bvi_->val_;
    const double& cd = cvi_->val_;
    bvi_->adj_
        += adj_ * inc_beta_ddb(ad_, bd, cd, digamma(bd), digamma(ad_ + bd));
    cvi_->adj_ += adj_ * inc_beta_ddz(ad_, bd, cd);
  }
};

class inc_beta_vdv_vari : public op_vdv_vari {
 public:
  inc_beta_vdv_vari(vari* avi, double bd, vari* cvi)
      : op_vdv_vari(inc_beta(avi->val_, bd, cvi->val_), avi, bd, cvi) {}
  void chain() {
    const double& ad = avi_->val_;
    const double& cd = cvi_->val_;
    avi_->adj_
        += adj_ * inc_beta_dda(ad, bd_, cd, digamma(ad), digamma(ad + bd_));
    cvi_->adj_ += adj_ * inc_beta_ddz(ad, bd_, cd);
  }
};

class inc_beta_vvd_vari : public op_vvd_vari {
 public:
  inc_beta_vvd_vari(vari* avi, vari* bvi, const double cd)
      : op_vvd_vari(inc_beta(avi->val_, bvi->val_, cd), avi, bvi, cd) {}
  void chain() {
    const double& ad = avi_->val_;
    const double& bd = bvi_->val_;
    const double digamma_abd = digamma(ad + bd);
    avi_->adj_ += adj_ * inc_beta_dda(ad, bd, cd_, digamma(ad), digamma_abd);
    bvi_->adj_ += adj_ * inc_beta_ddb(ad, bd, cd_, digamma(bd), digamma_abd);
  }
};

class inc_beta_ddv_vari : public op_ddv_vari {
 public:
  inc_beta_ddv_vari(const double ad, const double bd, vari* cvi)
      : op_ddv_vari(inc_beta(ad, bd, cvi->val_), ad, bd, cvi) {}
  void chain() {
    const double& cd = cvi_->val_;
    cvi_->adj_ += adj_ * inc_beta_ddz(ad_, bd_, cd);
  }
};

class inc_beta_vdd_vari : public op_vdd_vari {
 public:
  inc_beta_vdd_vari(vari* avi, double bd, double cd)
      : op_vdd_vari(inc_beta(avi->val_, bd, cd), avi, bd, cd) {}
  void chain() {
    const double& ad = avi_->val_;
    avi_->adj_
        += adj_ * inc_beta_dda(ad, bd_, cd_, digamma(ad), digamma(ad + bd_));
  }
};

class inc_beta_dvd_vari : public op_dvd_vari {
 public:
  inc_beta_dvd_vari(double ad, vari* bvi, double cd)
      : op_dvd_vari(inc_beta(ad, bvi->val_, cd), ad, bvi, cd) {}
  void chain() {
    const double& bd = bvi_->val_;
    bvi_->adj_
        += adj_ * inc_beta_ddb(ad_, bd, cd_, digamma(bd), digamma(ad_ + bd));
  }
};

}  // namespace internal

inline var inc_beta(const var& a, const var& b, const var& c) {
  return var(new internal::inc_beta_vvv_vari(a.vi_, b.vi_, c.vi_));
}

inline var inc_beta(const double a, const var& b, const var& c) {
  return var(new internal::inc_beta_dvv_vari(a, b.vi_, c.vi_));
}

inline var inc_beta(const var& a, const double b, const var& c) {
  return var(new internal::inc_beta_vdv_vari(a.vi_, b, c.vi_));
}

inline var inc_beta(const var& a, const var& b, const double c) {
  return var(new internal::inc_beta_vvd_vari(a.vi_, b.vi_, c));
}

inline var inc_beta(const double a, const double b, const var& c) {
  return var(new internal::inc_beta_ddv_vari(a, b, c.vi_));
}

inline var inc_beta(const var& a, const double b, const double c) {
  return var(new internal::inc_beta_vdd_vari(a.vi_, b, c));
}

inline var inc_beta(const double a, const var& b, const double c) {
  return var(new internal::inc_beta_dvd_vari(a, b.vi_, c));
}

}  // namespace math
}  // namespace stan
#endif
