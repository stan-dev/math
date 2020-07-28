#ifndef STAN_MATH_OPENCL_REV_COPY_HPP
#define STAN_MATH_OPENCL_REV_COPY_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/vec_concat.hpp>
#include <stan/math/opencl/rev/opencl.hpp>

#include <CL/cl2.hpp>
#include <iostream>
#include <vector>
#include <algorithm>
#include <type_traits>

namespace stan {
namespace math {

namespace internal {
template <typename T_arg_adj, require_eigen_t<T_arg_adj>* = nullptr>
class op_copy_to_cl_vari final
    : public vari_value<matrix_cl<value_type_t<T_arg_adj>>> {
  T_arg_adj arg_adj_;

 public:
  template <typename T_arg_val, require_eigen_t<T_arg_val>* = nullptr,
            require_vt_same<T_arg_val, T_arg_adj>* = nullptr>
  explicit op_copy_to_cl_vari(const T_arg_val& val, T_arg_adj adj)
      : vari_value<matrix_cl<value_type_t<T_arg_adj>>>(to_matrix_cl(val)),
        arg_adj_(adj) {}

  virtual void chain() {
    arg_adj_ += from_matrix_cl<std::decay_t<T_arg_adj>::RowsAtCompileTime,
                               std::decay_t<T_arg_adj>::ColsAtCompileTime>(
        this->adj_);
  }
};
}  // namespace internal

/** \ingroup opencl
 * Copies the source var containing Eigen matrices to destination var that has
 * data stored on the OpenCL device.
 *
 * @tparam T type of the Eigen matrix
 * @param a source Eigen matrix
 * @return var with a copy of the data on the OpenCL device
 */
template <typename T>
inline var_value<matrix_cl<value_type_t<T>>> to_matrix_cl(
    const var_value<T>& a) {
  return new internal::op_copy_to_cl_vari<decltype(a.vi_->adj_)>(a.val(),
                                                                 a.vi_->adj_);
}

namespace internal {
template <typename T, int Rows, int Cols,
          require_all_kernel_expressions_t<T>* = nullptr>
class op_copy_from_cl_vari final
    : public vari_value<Eigen::Matrix<value_type_t<T>, Rows, Cols>> {
  vari_value<T>& a_;

 public:
  explicit op_copy_from_cl_vari(vari_value<T>& a)
      : vari_value<Eigen::Matrix<value_type_t<T>, Rows, Cols>>(
            from_matrix_cl<Rows, Cols>(a.val_)),
        a_(a) {}

  virtual void chain() { a_.adj_ = a_.adj_ + to_matrix_cl(this->adj_); }
};
}  // namespace internal

/** \ingroup opencl
 * Copies the source var that has data stored on the OpenCL device to
 * destination var containing Eigen matrices.
 *
 * @tparam Rows number of compile time rows of the destination matrix
 * @tparam Rows number of compile time columns of the destination matrix
 * @tparam T type of the matrix or expression on the OpenCL device
 * @param a source matrix_cl or expression
 * @return var with a copy of the data on the host
 */
template <int Rows = Eigen::Dynamic, int Cols = Eigen::Dynamic, typename T,
          require_all_kernel_expressions_t<T>* = nullptr>
inline var_value<Eigen::Matrix<value_type_t<T>, Rows, Cols>> from_matrix_cl(
    const var_value<T>& a) {
  return new internal::op_copy_from_cl_vari<T, Rows, Cols>(*a.vi_);
}

/** \ingroup opencl
 * Copies the source Eigen matrix of vars to
 * the destination matrix that is stored
 * on the OpenCL device.
 *
 * @tparam R Compile time rows of the Eigen matrix
 * @tparam C Compile time columns of the Eigen matrix
 * @param src source Eigen matrix
 * @return matrix_cl with a copy of the data in the source matrix
 */
template <typename T, require_eigen_vt<is_var, T>* = nullptr>
inline var_value<matrix_cl<value_type_t<value_type_t<T>>>> to_matrix_cl(
    const T& src) {
  // the matrix can go out of scope before chain() is called. So we store a map
  // to the data
  var* src_array
      = ChainableStack::instance_->memalloc_.alloc_array<var>(src.size());
  Eigen::Map<plain_type_t<T>> src_stacked(src_array, src.rows(), src.cols());
  src_stacked = src;
  return new internal::op_copy_to_cl_vari<decltype(src_stacked.adj())>(
      src_stacked.val(), src_stacked.adj());
}

/** \ingroup opencl
 * Copies the adjoint of the source matrix of vars that is stored
 * on the OpenCL device to the destination Eigen
 * matrix.
 *
 * @param src source matrix of vars on the OpenCL device
 * @return Eigen matrix with a copy of the adjoints in the source matrix
 */
template <int R = Eigen::Dynamic, int C = Eigen::Dynamic, typename T,
          require_var_t<T>* = nullptr>
inline Eigen::Matrix<double, R, C> from_matrix_cl(const matrix_cl<T>& src) try {
  Eigen::Matrix<double, R, C> dst(src.rows(), src.cols());
  if (src.size() == 0) {
    return dst;
  }
  dst = from_matrix_cl(src.adj());
  return dst;
} catch (const cl::Error& e) {
  check_opencl_error("copy var (OpenCL)->Eigen", e);
  Eigen::Matrix<double, R, C> dst(src.rows(), src.cols());
  return dst;
}

/** \ingroup opencl
 * Packs the flat triangular matrix on the OpenCL device and
 * copies the adjoint values to the std::vector.
 *
 * @param src the flat triangular source matrix on the OpenCL device
 * @return the packed std::vector of adjoint values.
 * @throw <code>std::invalid_argument</code> if the matrix is not triangular
 */
template <typename T, typename = require_var_t<T>>
inline std::vector<double> packed_copy(const matrix_cl<T>& src) try {
  check_triangular("var packed_copy", "src", src);
  const int packed_size = src.rows() * (src.rows() + 1) / 2;
  std::vector<double> dst(packed_size);
  if (dst.size() == 0) {
    return dst;
  }
  dst = packed_copy(src.adj());
  return dst;
} catch (const cl::Error& e) {
  check_opencl_error("packed_copy var (OpenCL->std::vector)", e);
  std::vector<double> dst(0);
  return dst;
}

/** \ingroup opencl
 * Copies the packed triangular matrix from
 * the source std::vector to an OpenCL buffer and
 * unpacks it to a flat matrix on the OpenCL device.
 *
 * @tparam matrix_view the triangularity of the source matrix
 * @param src the packed source std::vector
 * @param rows the number of rows in the flat matrix
 * @return the destination flat matrix on the OpenCL device
 * @throw <code>std::invalid_argument</code> if the
 * size of the vector does not match the expected size
 * for the packed triangular matrix
 */
template <matrix_cl_view matrix_view>
inline matrix_cl<var> packed_copy(vari** src, int rows) try {
  matrix_cl<var> dst(rows, rows, matrix_view);
  if (dst.size() == 0) {
    return dst;
  }
  const Eigen::Matrix<vari*, -1, 1> foo
      = Eigen::Map<const Eigen::Matrix<vari*, -1, 1>>(src,
                                                      rows * (rows + 1) / 2, 1);
  dst.val() = packed_copy<matrix_view>(foo.val().eval(), rows);
  dst.adj() = packed_copy<matrix_view>(foo.adj().eval(), rows);
  return dst;
} catch (const cl::Error& e) {
  check_opencl_error("packed_copy var (std::vector->OpenCL)", e);
  matrix_cl<var> dst(rows, rows, matrix_view);
  return dst;
}

/** \ingroup opencl
 * Copies the source matrix to the
 * destination matrix. Both matrices
 * are stored on the OpenCL device.
 *
 * @tparam T An arithmetic type to pass the value from the OpenCL matrix to.
 * @param src the source matrix
 * @return matrix_cl with copies of values in the source matrix
 * @throw <code>std::invalid_argument</code> if the
 * matrices do not have matching dimensions
 */
template <typename T, require_var_t<T>* = nullptr>
inline matrix_cl<T> copy_cl(const matrix_cl<T>& src) {
  matrix_cl<T> dst(src.rows(), src.cols(), src.view());
  if (src.size() == 0) {
    return dst;
  }
  try {
    dst.val() = copy_cl(src.val());
    dst.adj() = copy_cl(src.adj());
  } catch (const cl::Error& e) {
    check_opencl_error("copy_cl var (OpenCL)->(OpenCL)", e);
  }
  return dst;
}

}  // namespace math
}  // namespace stan
#endif
#endif
