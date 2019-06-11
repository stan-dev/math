#ifndef STAN_MATH_OPENCL_REV_MATRIX_CL_HPP
#define STAN_MATH_OPENCL_REV_MATRIX_CL_HPP
#ifdef STAN_OPENCL
#include <stan/math/mix/mat.hpp>
#include <stan/math/prim/mat/fun/typedefs.hpp>
#include <stan/math/rev/mat/fun/typedefs.hpp>
#include <stan/math/fwd/mat/fun/typedefs.hpp>
#include <stan/math/mix/mat.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/rev/mat/fun/typedefs.hpp>
#include <stan/math/opencl/opencl_context.hpp>
#include <stan/math/opencl/constants.hpp>
#include <stan/math/opencl/err/check_opencl.hpp>
#include <stan/math/opencl/buffer_types.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/scal/err/check_size_match.hpp>
#include <stan/math/prim/scal/err/domain_error.hpp>
#include <stan/math/prim/arr/fun/vec_concat.hpp>
#include <stan/math/rev/core.hpp>
#include <CL/cl.hpp>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

namespace stan {
namespace math {

class var;
template <>
class matrix_cl<var> {
 private:
  /**
   * cl::Buffer provides functionality for working with the OpenCL buffer.
   * An OpenCL buffer allocates the memory in the device that
   * is provided by the context.
   */
  const int rows_;
  const int cols_;
  mutable matrix_cl<double> val_;
  mutable matrix_cl<double> adj_;
  mutable TriangularViewCL triangular_view_;
 public:
  int rows() const { return rows_; }

  int cols() const { return cols_; }

  int size() const { return rows_ * cols_; }

  matrix_cl<double>& val() const {return val_;}
  matrix_cl<double>& adj() const {return adj_;}
  explicit matrix_cl() : rows_(0), cols_(0) {}

  template <TriangularViewCL triangular_view = TriangularViewCL::Entire>
  void zeros();
  template <TriangularMapCL triangular_map = TriangularMapCL::LowerToUpper>
  void triangular_transpose();
  template <TriangularViewCL triangular_view = TriangularViewCL::Entire>
  void sub_block(const matrix_cl& A, size_t A_i, size_t A_j, size_t this_i,
                 size_t this_j, size_t nrows, size_t ncols);

  template <int R, int C>
  explicit matrix_cl(const Eigen::Matrix<var, R, C>& A)
      : rows_(A.rows()), cols_(A.cols()),
      val_(A.val().eval()),
      adj_(A.adj().eval()) {}

  explicit matrix_cl(const matrix_vi& A)
      : rows_(A.rows()), cols_(A.cols()),
      val_(A.val().eval()),
      adj_(A.adj().eval()) {}

  explicit matrix_cl(vari** A, const int& R, const int& C,
    TriangularViewCL triangular_view = TriangularViewCL::Entire) :
    rows_(R), cols_(C), triangular_view_(triangular_view),
    val_(R, C), adj_(R, C) {
    cl::Context& ctx = opencl_context.context();
    cl::CommandQueue& queue = opencl_context.queue();
    cl::Event transfer_event1;
    cl::Event transfer_event2;
    if (triangular_view == TriangularViewCL::Entire) {
      const matrix_vi foo = Eigen::Map<const matrix_vi>(A, R, C);
      val_.buffer() = cl::Buffer(ctx, CL_MEM_READ_WRITE, sizeof(double) * R * C);
      adj_.buffer() = cl::Buffer(ctx, CL_MEM_READ_WRITE, sizeof(double) * R * C);
      queue.enqueueWriteBuffer(val_.buffer(), CL_FALSE, 0,
                               sizeof(double) * R * C,
                               foo.val().eval().data(), NULL,
                               &transfer_event1);
      queue.enqueueWriteBuffer(adj_.buffer(), CL_FALSE, 0,
                              sizeof(double) * R * C,
                              foo.adj().eval().data(), NULL,
                              &transfer_event2);
    } else {
      const int packed_size = R * (R + 1) / 2;
      val_.buffer() = cl::Buffer(ctx, CL_MEM_READ_WRITE, sizeof(double) * packed_size);
      adj_.buffer() = cl::Buffer(ctx, CL_MEM_READ_WRITE, sizeof(double) * packed_size);
      const matrix_vi foo = Eigen::Map<const matrix_vi>(A, 1, packed_size);
      queue.enqueueWriteBuffer(val_.buffer(), CL_FALSE, 0,
                               sizeof(double) * packed_size,
                               foo.val().eval().data(), NULL,
                               &transfer_event1);
      queue.enqueueWriteBuffer(adj_.buffer(), CL_FALSE, 0,
                              sizeof(double) * packed_size,
                              foo.adj().eval().data(), NULL,
                              &transfer_event2);
    }
    val_.add_write_event(transfer_event1);
    adj_.add_write_event(transfer_event2);
  }

  explicit matrix_cl(const int& rows, const int& cols) :
  rows_(rows), cols_(cols), val_(rows, cols), adj_(rows, cols) {}

  matrix_cl<var> operator=(const matrix_cl<var>& A) {
    val_ = A.val();
    adj_ = A.adj();
    return *this;
  }
};

}
}

#endif
#endif
