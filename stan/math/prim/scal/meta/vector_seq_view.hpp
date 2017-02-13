#ifndef STAN_MATH_PRIM_SCAL_META_VECTOR_SEQ_VIEW_HPP
#define STAN_MATH_PRIM_SCAL_META_VECTOR_SEQ_VIEW_HPP

namespace stan {
  /**
   * vector_seq_view provides a sequence-like wrapper around either a single
   * item or a std::vector of items. It is similar to scalar_seq_view but is
   * intended to wrap Eigen vectors instead of scalars. In the base form, it
   * assumes that a single item has been passed in and exposes a sequence
   * interface to that single item. It has just one specialization, for
   * std::vector<T>. In that specialization, operator[] indexes into the vector
   * instead of always returning the item.
   *
   * @tparam T the type to be wrapped
   */
  template <typename T>
  class vector_seq_view {
  public:
    explicit vector_seq_view(const T& t) : t_(t) {}

    const T& operator[](int /* i */) const {
      return t_;
    }

    int size() const {
      return 1;
    }

  private:
    const T& t_;

  };
}

#endif
