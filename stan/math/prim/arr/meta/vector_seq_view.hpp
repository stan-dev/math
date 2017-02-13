#ifndef STAN_MATH_PRIM_ARR_META_VECTOR_SEQ_VIEW_HPP
#define STAN_MATH_PRIM_ARR_META_VECTOR_SEQ_VIEW_HPP

#include <stan/math/prim/scal/meta/vector_seq_view.hpp>
#include <vector>

namespace stan {
  /**
   * This specialization will allow indexing into std::vectors of underlying
   * items of type T. See base template class doc for more details.
   *
   * @tparam T type stored in the std::vector
   */
  template <typename T>
  class vector_seq_view<std::vector<T> > {
  public:
    explicit vector_seq_view(const std::vector<T>& c) : c_(c) {}

    const T& operator[](int i) const {
      return c_[i];
    }

    int size() const {
      return c_.size();
    }

  private:
    const std::vector<T>& c_;
  };
}

#endif
