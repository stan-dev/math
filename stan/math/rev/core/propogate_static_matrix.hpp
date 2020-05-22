#ifndef STAN_MATH_REV_CORE_PROPOGATE_STATIC_MATRIX_HPP
#define STAN_MATH_REV_CORE_PROPOGATE_STATIC_MATRIX_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/rev/meta.hpp>

#include <utility>
#include <vector>

namespace stan {
  namespace math {
    template <typename T>
    class dynamic_to_static_vari : public vari_base {
      Eigen::Index N_;
      vari_value<double>* output_;
      vari_value<T>* input_;
    public:
      dynamic_to_static_vari(vari_value<T>* input, vari_value<double>* output, Eigen::Index N) :
        input_(input), N_(N),
          output_(output) {
      }
      void chain() {
        for(size_t n = 0; n < N_; ++n) {
          input_->adj_(n) += output_->adj_;
          ++output_;
        }
        output_ -= N_;
      }
    };

    template <typename T>
    class static_to_dynamic_vari : public vari_base {
      Eigen::Index N_;
      vari_value<T>* output_;
      vari_value<double>* input_;
    public:
      static_to_dynamic_vari(vari_value<double>* input, vari_value<T>* output, Eigen::Index N) :
        input_(input), N_(N),
          output_(output) {
      }
      void chain() {
        for(size_t n = 0; n < N_; ++n) {
          input_->adj_ += output_->adj_(n);
          --input_;
        }
        input_ += N_;
      }
    };
  }
}
#endif
