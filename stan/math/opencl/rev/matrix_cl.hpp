#ifndef STAN_MATH_OPENCL_REV_MATRIX_CL_HPP
#define STAN_MATH_OPENCL_REV_MATRIX_CL_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <stan/math/prim/fun/vec_concat.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/rev/matrix_cl_vari.hpp>
#include <CL/cl2.hpp>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

namespace stan {
namespace math {

template <typename T>
class matrix_cl<T, require_var_t<T>> {
 private:
  /**
   * cl::Buffer provides functionality for working with the OpenCL buffer.
   * An OpenCL buffer allocates the memory in the device that
   * is provided by the context.
   */
  const int rows_{0};
  const int cols_{0};
  matrix_cl_view view_{matrix_cl_view::Entire};

 public:
  matrix_cl<vari>* vi_;
  using Scalar = T;
  using type = T;
  inline int rows() const { return rows_; }

  inline int cols() const { return cols_; }

  inline int size() const { return rows_ * cols_; }

  inline auto& val() const { return vi_->val_; }
  inline auto& adj() const { return vi_->adj_; }
  inline auto& vi() {
    return vi_;
  }
  matrix_cl() {}

  // Forward declare the methods that work in place on the matrix
  template <matrix_cl_view matrix_view = matrix_cl_view::Entire>
  inline void zeros();
  template <matrix_cl_view matrix_view = matrix_cl_view::Entire>
  inline void zeros_strict_tri();
  template <TriangularMapCL triangular_map = TriangularMapCL::LowerToUpper>
  inline void triangular_transpose();

  inline void sub_block(const matrix_cl<T, require_var_t<T>>& A, size_t A_i,
                        size_t A_j, size_t this_i, size_t this_j, size_t nrows,
                        size_t ncols);

  inline const matrix_cl_view& view() const { return view_; }

  inline void view(const matrix_cl_view& view) {
    view_ = view;
    vi_->view(view);
  }

  template <typename Mat, require_eigen_st<is_var, Mat>...>
  matrix_cl(Mat&& A,
                     matrix_cl_view partial_view = matrix_cl_view::Entire)
      : rows_(A.rows()),
        cols_(A.cols()),
        vi_(new matrix_cl<vari>(A.eval(), partial_view)),
        view_(partial_view) {}

  matrix_cl(const matrix_vi& A,
                     matrix_cl_view partial_view = matrix_cl_view::Entire)
      : rows_(A.rows()),
        cols_(A.cols()),
        vi_(new matrix_cl<vari>(A.eval(), partial_view)),
        view_(partial_view) {}

  matrix_cl(vari** A, const int& R, const int& C,
                     matrix_cl_view partial_view = matrix_cl_view::Entire)
      : rows_(R),
        cols_(C),
        view_(partial_view),
        vi_(new matrix_cl<vari>(A, R, C, partial_view)) {  }

  matrix_cl(const std::vector<var>& A, const int& R, const int& C,
                     matrix_cl_view partial_view = matrix_cl_view::Entire)
      : rows_(R),
        cols_(C),
        view_(partial_view),
        vi_(new matrix_cl<vari>(A, R, C, partial_view)) { }

  matrix_cl(const int& rows, const int& cols,
                     matrix_cl_view partial_view = matrix_cl_view::Entire)
      : rows_(rows),
        cols_(cols),
        vi_(new matrix_cl<vari>(rows, cols, partial_view)),
        view_(partial_view) {}

  matrix_cl(matrix_cl<vari>* A) : vi_(A), view_(A->view()), rows_(A->rows()), cols_(A->cols()){}
  static inline void* operator new(size_t nbytes) {
    return ChainableStack::instance_->memalloc_.alloc(nbytes);
  }
    // We need to write a stack allocator for this but this is just to clean up after tests
    ~matrix_cl() {
      delete vi_;
    }

};


inline matrix_cl<var> operator+(const matrix_cl<var>& a, const matrix_cl<var>& b) {
  auto* value = new internal::add_vari<matrix_cl<vari>, decltype(a.vi_), decltype(b.vi_)>(a.vi_, b.vi_);
  return matrix_cl<var>(value);
}

}  // namespace math
}  // namespace stan

#endif
#endif
