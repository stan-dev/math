#ifndef STAN_MATH_PRIM_SCAL_META_SCALAR_SEQ_VIEW_WRITABLE_HPP
#define STAN_MATH_PRIM_SCAL_META_SCALAR_SEQ_VIEW_WRITABLE_HPP

#include <stan/math/prim/scal/meta/scalar_type.hpp>

namespace stan {
  /**
   * scalar_seq_view_writable provides a uniform sequence-like wrapper around
   * either a scalar or a sequence of scalars. It performs exactly the same as
   * scalar_seq_view but without the const variables
   *
   * @tparam C the container type; will be the scalar type if wrapping a scalar
   * @tparam T the scalar type
   */
  template <typename C, typename T = typename scalar_type<C>::type>
  class scalar_seq_view_writable {
  public:
    explicit scalar_seq_view_writable(C& c) : c_(c) {}

    /**
     * Segfaults if out of bounds.
     * @param i index
     * @return the element at the specified position in the container
     */
    T& operator[](int i) {
      return c_[i];
    }

    int size() {
      return c_.size();
    }

  private:
    C& c_;
  };

  /**
   * This specialization handles wrapping a scalar as if it were a sequence.
   *
   * @tparam T the scalar type
   */
  template <typename T>
  class scalar_seq_view_writable<T, T> {
  public:
    explicit scalar_seq_view_writable(T& t) : t_(t) {}

    T& operator[](int /* i */) {
      return t_;
    }

    int size() {
      return 1;
    }

  private:
    T& t_;
  };
}
#endif
