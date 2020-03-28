#ifndef STAN_MATH_OPENCL_REV_MATRIX_CL_VARI_HPP
#define STAN_MATH_OPENCL_REV_MATRIX_CL_VARI_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <stan/math/prim/fun/vec_concat.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <CL/cl2.hpp>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

namespace stan {
namespace math {

  template <typename T>
  class matrix_cl<T, require_var_t<T>> ;

template <typename T>
class matrix_cl<T, require_t<std::is_same<vari, std::decay_t<T>>>>  {
 private:
   friend class matrix_cl<var>;
  /**
   * cl::Buffer provides functionality for working with the OpenCL buffer.
   * An OpenCL buffer allocates the memory in the device that
   * is provided by the context.
   */
  const int rows_{0};
  const int cols_{0};
  matrix_cl_view view_{matrix_cl_view::Entire};

 public:
  mutable matrix_cl<double> val_;
  mutable matrix_cl<double> adj_;
  using Scalar = T;
  using type = T;
  inline int rows() const { return rows_; }

  inline int cols() const { return cols_; }

  inline int size() const { return rows_ * cols_; }

  inline matrix_cl<double>& val() const { return val_; }
  inline matrix_cl<double>& adj() const { return adj_; }
  matrix_cl() : rows_(0), cols_(0), val_(0, 0), adj_(0, 0) {}

  // Forward declare the methods that work in place on the matrix
  template <matrix_cl_view matrix_view = matrix_cl_view::Entire>
  inline void zeros();
  template <matrix_cl_view matrix_view = matrix_cl_view::Entire>
  inline void zeros_strict_tri();
  template <TriangularMapCL triangular_map = TriangularMapCL::LowerToUpper>
  inline void triangular_transpose();

  inline void sub_block(const matrix_cl<T, internal::require_vari_t<T>>& A, size_t A_i,
                        size_t A_j, size_t this_i, size_t this_j, size_t nrows,
                        size_t ncols);

  inline const matrix_cl_view& view() const { return view_; }

  inline void view(const matrix_cl_view& view) {
    view_ = view;
    val_.view(view);
    adj_.view(view);
  }

  template <typename Mat, require_eigen_st<is_var, Mat>...>
  matrix_cl(Mat&& A,
                     matrix_cl_view partial_view = matrix_cl_view::Entire)
      : rows_(A.rows()),
        cols_(A.cols()),
        val_(A.val().eval(), partial_view),
        adj_(A.adj().eval(), partial_view),
        view_(partial_view) {}

  matrix_cl(const matrix_vi& A,
                     matrix_cl_view partial_view = matrix_cl_view::Entire)
      : rows_(A.rows()),
        cols_(A.cols()),
        val_(A.val().eval(), partial_view),
        adj_(A.adj().eval(), partial_view),
        view_(partial_view) {}

  matrix_cl(vari** A, const int& R, const int& C,
                     matrix_cl_view partial_view = matrix_cl_view::Entire)
      : rows_(R),
        cols_(C),
        view_(partial_view),
        val_(R, C, partial_view),
        adj_(R, C, partial_view) {
    cl::Context& ctx = opencl_context.context();
    cl::CommandQueue& queue = opencl_context.queue();
    cl::Event transfer_event1;
    cl::Event transfer_event2;
    if (view_ == matrix_cl_view::Entire) {
      const matrix_vi foo = Eigen::Map<const matrix_vi>(A, R, C);
      val_.buffer()
          = cl::Buffer(ctx, CL_MEM_READ_WRITE, sizeof(double) * R * C);
      adj_.buffer()
          = cl::Buffer(ctx, CL_MEM_READ_WRITE, sizeof(double) * R * C);
      queue.enqueueWriteBuffer(val_.buffer(), CL_FALSE, 0,
                               sizeof(double) * R * C, foo.val().eval().data(),
                               NULL, &transfer_event1);
      queue.enqueueWriteBuffer(adj_.buffer(), CL_FALSE, 0,
                               sizeof(double) * R * C, foo.adj().eval().data(),
                               NULL, &transfer_event2);
    } else {
      const int packed_size = R * (R + 1) / 2;
      val_.buffer()
          = cl::Buffer(ctx, CL_MEM_READ_WRITE, sizeof(double) * packed_size);
      adj_.buffer()
          = cl::Buffer(ctx, CL_MEM_READ_WRITE, sizeof(double) * packed_size);
      const matrix_vi foo = Eigen::Map<const matrix_vi>(A, 1, packed_size);
      queue.enqueueWriteBuffer(val_.buffer(), CL_FALSE, 0,
                               sizeof(double) * packed_size,
                               foo.val().eval().data(), NULL, &transfer_event1);
      queue.enqueueWriteBuffer(adj_.buffer(), CL_FALSE, 0,
                               sizeof(double) * packed_size,
                               foo.adj().eval().data(), NULL, &transfer_event2);
    }
    val_.add_write_event(transfer_event1);
    adj_.add_write_event(transfer_event2);
    ChainableStack::instance_->var_stack_.push_back(this);
  }
  template <typename TT, require_t<std::is_assignable<matrix_cl<double>, std::decay_t<TT>>>* = nullptr>
  matrix_cl(TT&& A) : val_(A.eval()), adj_(A.rows(), A.cols()) {
    adj_.template zeros<matrix_cl_view::Entire>();
    ChainableStack::instance_->var_stack_.push_back(this);
    }

  explicit matrix_cl(const std::vector<var>& A, const int& R, const int& C,
                     matrix_cl_view partial_view = matrix_cl_view::Entire)
      : rows_(R),
        cols_(C),
        view_(partial_view),
        val_(R, C, partial_view),
        adj_(R, C, partial_view) {
    cl::Context& ctx = opencl_context.context();
    cl::CommandQueue& queue = opencl_context.queue();
    cl::Event transfer_event1;
    cl::Event transfer_event2;
    if (view_ == matrix_cl_view::Entire) {
      const matrix_v foo = Eigen::Map<const matrix_v>(A.data(), R, C);
      val_.buffer()
          = cl::Buffer(ctx, CL_MEM_READ_WRITE, sizeof(double) * R * C);
      adj_.buffer()
          = cl::Buffer(ctx, CL_MEM_READ_WRITE, sizeof(double) * R * C);
      queue.enqueueWriteBuffer(val_.buffer(), CL_FALSE, 0,
                               sizeof(double) * R * C, foo.val().eval().data(),
                               NULL, &transfer_event1);
      queue.enqueueWriteBuffer(adj_.buffer(), CL_FALSE, 0,
                               sizeof(double) * R * C, foo.adj().eval().data(),
                               NULL, &transfer_event2);
    } else {
      const int packed_size = R * (R + 1) / 2;
      val_.buffer()
          = cl::Buffer(ctx, CL_MEM_READ_WRITE, sizeof(double) * packed_size);
      adj_.buffer()
          = cl::Buffer(ctx, CL_MEM_READ_WRITE, sizeof(double) * packed_size);
      const matrix_v foo = Eigen::Map<const matrix_v>(A.data(), 1, packed_size);
      queue.enqueueWriteBuffer(val_.buffer(), CL_FALSE, 0,
                               sizeof(double) * packed_size,
                               foo.val().eval().data(), NULL, &transfer_event1);
      queue.enqueueWriteBuffer(adj_.buffer(), CL_FALSE, 0,
                               sizeof(double) * packed_size,
                               foo.adj().eval().data(), NULL, &transfer_event2);
    }
    val_.add_write_event(transfer_event1);
    adj_.add_write_event(transfer_event2);
    ChainableStack::instance_->var_stack_.push_back(this);
  }

  explicit matrix_cl(const int& rows, const int& cols,
                     matrix_cl_view partial_view = matrix_cl_view::Entire)
      : rows_(rows),
        cols_(cols),
        val_(rows, cols),
        adj_(rows, cols),
        view_(partial_view) {
          ChainableStack::instance_->var_stack_.push_back(this);
        }
        static inline void* operator new(size_t nbytes) {
          return ChainableStack::instance_->memalloc_.alloc(nbytes);
        }
        static inline void operator delete(void* /* ignore arg */) { /* no op */
        }
};
}
}
#endif
#endif
