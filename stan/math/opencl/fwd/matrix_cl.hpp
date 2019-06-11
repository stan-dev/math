#ifndef STAN_MATH_OPENCL_FWD_MATRIX_CL_HPP
#define STAN_MATH_OPENCL_FWD_MATRIX_CL_HPP
#ifdef STAN_OPENCL
#include <stan/math/opencl/opencl_context.hpp>
#include <stan/math/opencl/constants.hpp>
#include <stan/math/opencl/err/check_opencl.hpp>
#include <stan/math/opencl/buffer_types.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/scal/err/check_size_match.hpp>
#include <stan/math/prim/scal/err/domain_error.hpp>
#include <stan/math/prim/arr/fun/vec_concat.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/opencl/rev/matrix_cl.hpp>
#include <CL/cl.hpp>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

namespace stan {
namespace math {

template <typename T>
class fvar;

template <>
template <typename T>
class matrix_cl<fvar<T>>;

template <typename T>
class matrix_cl<fvar<T>> {
 private:
  /**
   * cl::Buffer provides functionality for working with the OpenCL buffer.
   * An OpenCL buffer allocates the memory in the device that
   * is provided by the context.
   */
  const int rows_;
  const int cols_;
  mutable matrix_cl<T> val_;
  mutable matrix_cl<T> d_;

 public:
  int rows() const { return rows_; }

  int cols() const { return cols_; }

  int size() const { return rows_ * cols_; }

  matrix_cl<T>& val() const {return val_;}
  matrix_cl<T>& d() const {return d_;}
  explicit matrix_cl() : rows_(0), cols_(0) {}
  template <int R, int C>
  explicit matrix_cl(const Eigen::Matrix<fvar<T>, R, C>& A)
      : val_(A.val().eval()), d_(A.d().eval()), rows_(A.rows()), cols_(A.cols()) {}

  explicit matrix_cl(const int& rows, const int& cols) :
  rows_(rows), cols_(cols) {}

  matrix_cl<fvar<T>> operator=(const matrix_cl<fvar<T>>& A) {
    val_ = A.val();
    d_ = A.d();
    return *this;
  }
};

}
}

#endif
#endif
