#ifndef STAN_MATH_PRIM_MAT_META_VECTOR_SEQ_VIEW_HPP
#define STAN_MATH_PRIM_MAT_META_VECTOR_SEQ_VIEW_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <vector>

namespace stan {
  using std::vector;
  using Eigen::Matrix;
  using Eigen::Dynamic;

  /**
   * Provides a low-cost wrapper for situations where you either need an Eigen
   * Vector or RowVector or a std::vector of them and you want to be agnostic
   * between those two options. This is similar to scalar_seq_view but instead
   * of being a sequence-like view over a scalar or seq of scalars, it's a
   * sequence-like view over a vector or seq of vectors. Notably this version
   * only allows std::vectors as the sequence type, since we would have
   * difficulty figuring out which contained type was the container otherwise.
   *
   * @tparam T the underlying scalar type inside the underlying Vector
   */
  template <typename T>
  class vector_seq_view { };

  template <typename T>
  class vector_seq_view<Matrix<T, Dynamic, 1> > {
    public:
      explicit vector_seq_view(const Matrix<T, Dynamic, 1>& m): m_(m) {}
      int size() const {
        return 1;
      }
      Matrix<T, Dynamic, 1> operator[](int /* i */) const {
        return m_;
      }
    private:
      const Matrix<T, Dynamic, 1>& m_;
  };

  template <typename T>
  class vector_seq_view<Matrix<T, 1, Dynamic> > {
    public:
      explicit vector_seq_view(const Matrix<T, 1, Dynamic>& m): m_(m) {}
      int size() const {
        return 1;
      }
      Matrix<T, 1, Dynamic> operator[](int /* i */) const {
        return m_;
      }
    private:
      const Matrix<T, 1, Dynamic>& m_;
  };

  template <typename T>
  class vector_seq_view<vector<Matrix<T, Dynamic, 1> > > {
    public:
      explicit vector_seq_view(const vector<Matrix<T, Dynamic, 1> >& v):
        v_(v) {}
      int size() const {
        return v_.size();
      }
      Matrix<T, Dynamic, 1> operator[](int i) const {
        return v_[i];
      }
    private:
      const vector<Matrix<T, Dynamic, 1> >& v_;
  };

  template <typename T>
  class vector_seq_view<vector<Matrix<T, 1, Dynamic> > > {
    public:
      explicit vector_seq_view(const vector<Matrix<T, 1, Dynamic> >& v):
        v_(v) {}
      int size() const {
        return v_.size();
      }
      Matrix<T, 1, Dynamic> operator[](int i) const {
        return v_[i];
      }
    private:
      const vector<Matrix<T, 1, Dynamic> >& v_;
  };
}  // namespace stan

#endif
