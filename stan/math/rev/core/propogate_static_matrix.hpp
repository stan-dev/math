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
      vari_value<double>** output_;
      vari_value<T>* input_;
    public:
      template <typename TT>
      dynamic_to_static_vari(vari_value<T>* input, TT&& output, Eigen::Index N) :
        input_(input), N_(N),
          output_(ChainableStack::instance_->memalloc_.alloc_array<vari*>(N_)) {
        using matrix_vi = Eigen::Matrix<vari*, Eigen::Dynamic, Eigen::Dynamic>;
        Eigen::Map<matrix_vi>(output_, output.rows(), output.cols()) = output.vi();
      }
      void chain() {
        for(size_t n = 0; n < N_; ++n) {
          input_->adj_(n) += output_[n]->adj_;
        }
      }
    };

    template <typename T>
    class static_to_dynamic_vari : public vari_base {
      Eigen::Index N_;
      vari_value<T>* output_;
      vari_value<double>** input_;
    public:
      template <typename TT>
      static_to_dynamic_vari(TT&& input, vari_value<T>* output, Eigen::Index N) :
      N_(N),
      input_(ChainableStack::instance_->memalloc_.alloc_array<vari*>(N_)),
      output_(output) {
        using matrix_vi = Eigen::Matrix<vari*, Eigen::Dynamic, Eigen::Dynamic>;
        Eigen::Map<matrix_vi>(input_, input.rows(), input.cols()) = input.vi();
      }
      void chain() {
        for(size_t n = 0; n < N_; ++n) {
          input_[n]->adj_ += output_->adj_(n);
        }
      }
    };
  }
}
#endif
