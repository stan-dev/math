#ifndef STAN_MATH_REV_CORE_PROPOGATE_STATIC_MATRIX_HPP
#define STAN_MATH_REV_CORE_PROPOGATE_STATIC_MATRIX_HPP

#include <stan/math/rev/core/vari.hpp>
#include <stan/math/prim/meta.hpp>
#include <utility>
#include <vector>

namespace stan {
namespace math {

template <int R, int C>
class from_static_vari : public vari_base {
  size_t N_;
  vari_value<Eigen::Matrix<double, R, C>>* input_vi_;
  vari** output_vis_;

 public:
  from_static_vari(size_t N, vari_value<Eigen::Matrix<double, R, C>>* input_vi,
                   vari** output_vis)
      : N_(N), input_vi_(input_vi), output_vis_(output_vis) {}
  void chain() {
    /*input_vi_->adj_ +=
      Eigen::Map<Eigen::Matrix<vari*, R, C>>(output_vis_,
					     input_vi_->adj_.rows(),
					     input_vi_->adj_.cols()).adj_();*/
    for (size_t n = 0; n < N_; ++n) {
      input_vi_->adj_(n) += output_vis_[n]->adj_;
    }
  }
};

template <int R, int C>
class to_static_vari : public vari_base {
  Eigen::Index N_;
  vari** input_vis_;
  vari_value<Eigen::Matrix<double, R, C>>* output_vi_;

 public:
  to_static_vari(size_t N, vari** input_vis,
                 vari_value<Eigen::Matrix<double, R, C>>* output_vi)
      : N_(N), input_vis_(input_vis), output_vi_(output_vi) {}
  void chain() {
    /*Eigen::Map<Eigen::Matrix<vari*, R, C>>(input_vis_,
					   output_vi_->adj_.rows(),
					   output_vi_->adj_.cols()).adj_() +=
					   output_vi_->adj_;*/
    for (size_t n = 0; n < N_; ++n) {
      input_vis_[n]->adj_ += output_vi_->adj_(n);
    }
  }
};

}  // namespace math
}  // namespace stan
#endif
