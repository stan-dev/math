#ifndef STAN_MATH_PRIM_MAT_META_SCALAR_SEQ_VIEW_HPP
#define STAN_MATH_PRIM_MAT_META_SCALAR_SEQ_VIEW_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/scal/meta/scalar_seq_view.hpp>

namespace stan {
  /**
   * This specialization allows scalar_seq_view to wrap matrices; required
   * because matrices are indexed by () instead of [].
   *
   * @tparam T scalar type inside the Matrix.
   */
  template <typename T>
  class scalar_seq_view<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>, T> {
    typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> matrix_t;

  public:
    explicit scalar_seq_view(const matrix_t& m) : m_(m) {}

    const T& operator[] (int i) const {
      return m_(i);
    }

    int size() const {
      return m_.size();
    }

  private:
    const matrix_t& m_;
  };
}

#endif
