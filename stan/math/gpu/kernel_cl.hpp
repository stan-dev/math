#ifndef STAN_MATH_GPU_CL_HPP
#define STAN_MATH_GPU_CL_HPP
#ifdef STAN_OPENCL
#include <stan/math/gpu/opencl_context.hpp>
#include <stan/math/gpu/constants.hpp>
#include <CL/cl.hpp>
#include <string>
#include <algorithm>
#include <map>
#include <vector>

namespace stan {
namespace math {
class kernel_cl;  // Declare here so kernel_cl_base can see the class
namespace internal {
class kernel_cl_base {
  friend class stan::math::kernel_cl;

 private:
  std::string copy_matrix =
#include <stan/math/gpu/kernels/copy_matrix.cl>  // NOLINT
      ;                                          // NOLINT
  std::string transpose_matrix =
#include <stan/math/gpu/kernels/transpose_matrix.cl>  // NOLINT
      ;                                               // NOLINT
  std::string zeros_matrix =
#include <stan/math/gpu/kernels/zeros_matrix.cl>  // NOLINT
      ;                                           // NOLINT
  std::string identity_matrix =
#include <stan/math/gpu/kernels/identity_matrix.cl>  // NOLINT
      ;                                              // NOLINT
  std::string copy_triangular_matrix =
#include <stan/math/gpu/kernels/copy_triangular_matrix.cl>  // NOLINT
      ;                                                     // NOLINT
  std::string copy_triangular_transposed_matrix =
#include <stan/math/gpu/kernels/triangular_transpose.cl>  // NOLINT
      ;                                                   // NOLINT
  std::string copy_submatrix =
#include <stan/math/gpu/kernels/sub_block.cl>  // NOLINT
      ;                                        // NOLINT
  std::string check_nan =
#include <stan/math/gpu/kernels/check_nan.cl>  // NOLINT
      ;                                        // NOLINT
  std::string check_diagonal_zeros =
#include <stan/math/gpu/kernels/check_diagonal_zeros.cl>  // NOLINT
      ;                                                   // NOLINT
  std::string check_symmetric =
#include <stan/math/gpu/kernels/check_symmetric.cl>  // NOLINT
      ;                                              // NOLINT
  std::string subtract_symmetric =
#include <stan/math/gpu/kernels/subtract_matrix.cl>  // NOLINT
      ;                                              // NOLINT
  std::string add_symmetric =
#include <stan/math/gpu/kernels/add_matrix.cl>  // NOLINT
      ;                                         // NOLINT
  kernel_cl_base() {
    std::string kernel_opts = "";
    std::string kernel_source = helpers;
    for (auto kernel_info : kernel_table) {
      kernel_source += kernel_info.second;
      for (auto comp_opts : base_opts) {
        kernel_opts += std::string(" -D") + comp_opts.first + "="
                       + std::to_string(comp_opts.second);
      }
    }
    try {
      cl::Program::Sources source(
          1,
          std::make_pair(kernel_source.c_str(), strlen(kernel_source.c_str())));
      cl::Program program_ = cl::Program(opencl_context.context(), source);
      program_.build({opencl_context.device()}, kernel_opts.c_str());

      cl_int err = CL_SUCCESS;
      // Iterate over all the kernels and create the cl::Kernel objects
      for (auto kernel_info : kernel_table) {
        kernels[kernel_info.first]
            = cl::Kernel(program_, kernel_info.first, &err);
      }
    } catch (const cl::Error& e) {
      check_opencl_error("Kernel Compilation", e);
    }
  }

 protected:
  std::string helpers =                      // Helper macros for the kernels.
#include <stan/math/gpu/kernels/helpers.cl>  // NOLINT
      ;                                      // NOLINT
  // Holds Default parameter values for each Kernel.
  typedef std::map<const char*, int> map_base_opts;
  map_base_opts base_opts
      = {{"LOWER", static_cast<int>(TriangularViewGPU::Lower)},
         {"UPPER", static_cast<int>(TriangularViewGPU::Upper)},
         {"ENTIRE", static_cast<int>(TriangularViewGPU::Entire)},
         {"UPPER_TO_LOWER", static_cast<int>(TriangularMapGPU::UpperToLower)},
         {"LOWER_TO_UPPER", static_cast<int>(TriangularMapGPU::LowerToUpper)}};

  typedef std::map<const char*, std::string> map_table;
  /**
   * Map of a kernel name (first) and it's meta information (second).
   *
   * Kernel  | Description
   * ------- | -------------
   * add | C = A + B
   * copy  | Copy matrix A on the GPU to matrix B.
   * copy_triangular | Copy the lower or upper triangular matrix of A to B.
   * copy_triangular_transposed | Copy the transpose lower/upper triangular.
   * identity | Create a NxM identity matrix.
   * is_nan | Check if a matrix on the GPU contains nan values.
   * is_zero_on_diagonal | Check if a matrix has zeros on the diagonal.
   * is_symmetric | Check if a matrix is symmetric.
   * sub_block | Copy a subset of matrix A into B.
   * subtract | C = A - B
   * transpose  | Take the transpose of a matrix.
   * zeros | Make a lower, upper, or full matrix of zeros.
   */
  const map_table kernel_table
      = {{"add", add_symmetric},
         {"copy", copy_matrix},
         {"copy_triangular", copy_triangular_matrix},
         {"copy_triangular_transposed", copy_triangular_transposed_matrix},
         {"identity", identity_matrix},
         {"is_nan", check_nan},
         {"is_symmetric", check_symmetric},
         {"is_zero_on_diagonal", check_diagonal_zeros},
         {"sub_block", copy_submatrix},
         {"subtract", subtract_symmetric},
         {"transpose", transpose_matrix},
         {"zeros", zeros_matrix}};
  typedef std::map<const char*, cl::Kernel> map_kern;
  map_kern kernels;  // The compiled kernels

  static kernel_cl_base& getInstance() {
    static kernel_cl_base instance_;
    return instance_;
  }
  kernel_cl_base(kernel_cl_base const&) = delete;
  void operator=(kernel_cl_base const&) = delete;
};
}  // namespace internal

/**
 * The adapter class that can access the <code>kernel_cl_base</code> class
 *
 *  This class operates as an API to the <code>kernel_cl_base</code>
 * class. Kernals can be called
 */
class kernel_cl {
 private:
  /**
   * return all compiled kernels.
   */
  inline internal::kernel_cl_base::map_kern kernels() {
    return internal::kernel_cl_base::getInstance().kernels;
  }

  /**
   * Terminating function for recursively setting arguments
   *  in an OpenCL kernel.
   *
   * @param k An OpenCL kernel.
   * @param i The <code>i</code>th argument to the kernel.
   * @note This function definition serves to end the recursive call for
   * <code>set_args()</code>
   */
  inline void recursive_args(cl::Kernel& k, int i) {}

  /**
   * Used in <code>set_args()</code> to add arguments to an OpenCL
   * kernel.
   *
   * @param kernel An OpenCL kernel.
   * @param i the position of the argument to the OpenCL kernel.
   * @param first_arg The first argument to go into the OpenCL kernel.
   * @param extra_args The remaining arguments to go into the OpenCL kernel.
   * @tparam T the type of <code>first_arg</code>.
   * @tparam Args The types of <code>extra_args</code>.
   */
  template <typename T, typename... Args>
  inline void recursive_args(cl::Kernel& kernel, int i, const T& first_arg,
                             const Args&... extra_args) {
    kernel.setArg(i, first_arg);
    this->recursive_args(kernel, i + 1, extra_args...);
  }

  /**
   * Adds arguments to an OpenCL kernel.
   *
   * @param kernel An OpenCL kernel.
   * @param args The arguments to mote to the OpenCL kernel.
   * @tparam Args The types of <code>extra_args</code>.
   */
  template <typename... Args>
  inline void set_args(cl::Kernel& kernel, const Args&... args) {
    this->recursive_args(kernel, 0, args...);
  }

  /**
   * return a kernel with no parameters set
   * @param kernel_name The name of a kernel specified in the
   *  <code>kernel_table</code>. See the docs of <code>kernel_table</code> for
   *  all available kernels.
   */
  inline auto get_kernel(const char* kernel_name) {
    return internal::kernel_cl_base::getInstance().kernels[(kernel_name)];
  }

 public:
  /**
   * API for accessing and returning OpenCL kernels.
   */
  kernel_cl() = default;

  /**
   * Returns kernel with parameters filled for addition of two matrices.
   *
   * @param[out] C Buffer of the output <code>matrix_gpu</code>.
   * @param[in] A Buffer of a <code>matrix_gpu</code>.
   * @param[in] B Buffer of a <code>matrix_gpu</code>.
   * @param rows Number of rows for matrix A.
   * @param cols Number of rows for matrix B.
   *
   * @note This kernel uses the helper macros available in helpers.cl.
   */
  inline auto add(cl::Buffer A, cl::Buffer B, cl::Buffer C, int rows,
                  int cols) {
    auto kern = this->get_kernel("add");
    this->set_args(kern, A, B, C, rows, cols);
    return kern;
  }

  /**
   * Returns kernel with parameters for copy one matrix on the GPU to another.
   * @param[in] A Buffor of the <code>matrix_gpu</code> to copy.
   * @param[out] B Buffer of the <code>matrix_gpu</code> to copy A to.
   * @param rows The number of rows in A.
   * @param cols The number of cols in A.
   *
   * @note Kernel used in math/gpu/matrix_gpu.hpp.
   *  This kernel uses the helper macros available in helpers.cl.
   */
  inline auto copy(cl::Buffer A, cl::Buffer B, int rows, int cols) {
    auto kern = this->get_kernel("copy");
    this->set_args(kern, A, B, rows, cols);
    return kern;
  }

  /**
   * Returns kernel for copying a lower/upper triangular of a matrix to it's
   * upper/lower.
   *
   * @param[in,out] A Buffer of a <code>matrix_gpu</code>.
   * @param rows The number of rows in A.
   * @param cols The number of cols in A.
   * @tparam triangular_map Specifies if the copy is
   * lower-to-upper or upper-to-lower triangular. The value
   * must be of type TriangularMap.
   *
   * @note Used in mat/gpu/copy_triangular_transposed.hpp.
   *  This kernel uses the helper macros available in helpers.cl.
   */
  template <TriangularMapGPU triangular_map = TriangularMapGPU::LowerToUpper>
  inline auto copy_triangular_transposed(cl::Buffer A, int rows, int cols) {
    auto kern = this->get_kernel("copy_triangular_transposed");
    this->set_args(kern, A, rows, cols, static_cast<int>(triangular_map));
    return kern;
  }
  /**
   * Returns the kernel for copying the lower or upper
   * triangular of the source matrix to
   * the destination matrix.
   * Both matrices are stored on the GPU.
   *
   * @param[out] A Buffer of the output <code>matrix_gpu</code> to copy
   * triangular to.
   * @param[in] B Buffer of the <code>matrix_gpu</code> to copy the triangular
   * from.
   * @param rows The number of rows of B.
   * @param cols The number of cols of B.
   * @tparam triangular_map int to describe
   * which part of the matrix to copy:
   * TriangularViewGPU::Lower - copies the lower triangular
   * TriangularViewGPU::Upper - copes the upper triangular
   *
   * @note Used in math/gpu/copy_triangular_opencl.hpp.
   *  This kernel uses the helper macros available in helpers.cl.
   */
  template <TriangularViewGPU triangular_view = TriangularViewGPU::Entire>
  inline auto copy_triangular(cl::Buffer A, cl::Buffer B, int rows, int cols) {
    auto kern = this->get_kernel("copy_triangular");
    this->set_args(kern, A, B, rows, cols, static_cast<int>(triangular_view));
    return kern;
  }

  /**
   * Returns the kernel to make an identity matrix on the GPU
   *
   * @param[in,out] A Buffer of the identity matrix output.
   * @param rows The number of rows for A.
   * @param cols The number of cols for A.
   *
   * @note Used in math/gpu/identity_opencl.hpp.
   *  This kernel uses the helper macros available in helpers.cl.
   */
  inline auto identity(cl::Buffer A, int rows, int cols) {
    auto kern = this->get_kernel("identity");
    this->set_args(kern, A, rows, cols);
    return kern;
  }

  /**
   * Returns the kernel to check if the <code>matrix_gpu</code> has NaN values.
   *
   * @param[in] A Buffer of the <code>matrix_gpu</code> to check.
   * @param rows The number of rows in matrix A.
   * @param cols The number of columns in matrix A.
   * @param[out] flag the flag to be written to if any diagonal is zero.
   *
   * @note Kernel for stan/math/gpu/err/check_nan.hpp.
   *  This kernel uses the helper macros available in helpers.cl.
   */
  inline auto is_nan(cl::Buffer A, int rows, int cols, cl::Buffer flag) {
    auto kern = this->get_kernel("is_nan");
    this->set_args(kern, A, rows, cols, flag);
    return kern;
  }

  /**
   * Returns kernel to check if the <code>matrix_gpu</code> is symmetric
   *
   * @param[in] A Buffer of the <code>matrix_gpu</code> to check.
   * @param[in] rows The number of rows in matrix A.
   * @param[in] cols The number of columns in matrix A.
   * @param[out] flag the flag to be written to if any diagonal is zero.
   * @param[in] tol The numeric tolerance when comparing doubles for symmetry.
   *
   * @note Kernel for stan/math/gpu/err/check_symmetric.hpp.
   *  This kernel uses the helper macros available in helpers.cl.
   */
  inline auto is_symmetric(cl::Buffer A, int rows, int cols, cl::Buffer flag,
                           double tol) {
    auto kern = this->get_kernel("is_symmetric");
    this->set_args(kern, A, rows, cols, flag, tol);
    return kern;
  }

  /**
   * Returns kernel to check if the <code>matrix_gpu</code> has zeros on the
   * diagonal
   *
   * @param[in] A Buffer of the Matrix to check.
   * @param rows The number of rows for A.
   * @param cols The number of cols of A.
   * @param[out] flag the flag to be written to if any diagonal is zero.
   *
   * @note Kernel for stan/math/gpu/err/check_diagonal_zeros.hpp.
   *  This kernel uses the helper macros available in helpers.cl.
   */
  inline auto is_zero_on_diagonal(cl::Buffer A, int rows, int cols,
                                  cl::Buffer flag) {
    auto kern = this->get_kernel("is_zero_on_diagonal");
    this->set_args(kern, A, rows, cols, flag);
    return kern;
  }

  /**
   * Returns kernel to copy a submatrix of the source matrix to
   * the destination matrix. The submatrix to copy
   * starts at (0, 0)
   * and is of size size_rows x size_cols.
   * The submatrix is copied to the
   * destination matrix starting at
   * (dst_offset_rows, dst_offset_cols)
   *
   * @param[in] src Buffer of the source matrix.
   * @param[out] dst Buffer of the destination submatrix.
   * @param src_i the offset row in A
   * @param src_j the offset column in A
   * @param this_i the offset row for the matrix to be subset into
   * @param this_j the offset col for the matrix to be subset into
   * @param nrows the number of rows in the submatrix
   * @param ncols the number of columns in the submatrix
   * @param src_rows The number of rows for matrix A
   * @param src_cols The number of cols for matrix A
   * @param dst_rows The number of rows for matrix B
   * @param dst_cols The number of rows for matrix B
   *
   * @note used in math/gpu/matrix_gpu.hpp.
   *  This kernel uses the helper macros available in helpers.cl.
   *
   */
  template <typename T>
  inline auto sub_block(const cl::Buffer& src, const cl::Buffer& dst, T src_i,
                        T src_j, T this_i, T this_j, T nrows, T ncols,
                        T src_rows, T src_cols, T dst_rows, T dst_cols) {
    auto kern = this->get_kernel("sub_block");
    this->set_args(kern, src, dst, src_i, src_j, this_i, this_j, nrows, ncols,
                   src_rows, src_cols, dst_rows, dst_cols);
    return kern;
  }

  /**
   * Returns kernel for matrix subtraction on the GPU Subtracts the second
   * matrix from the first matrix and stores the result in the third matrix
   * (C=A-B).
   *
   * @param[out] C Buffer of the output matrix.
   * @param[in] B Buffer of RHS input matrix.
   * @param[in] A Buffer LHS input matrix.
   * @param rows The number of rows for matrix A.
   * @param cols The number of columns for matrix A.
   *
   * @note Used in math/gpu/subtract_opencl.hpp
   *  This kernel uses the helper macros available in helpers.cl.
   */
  inline auto subtract(cl::Buffer A, cl::Buffer B, cl::Buffer C, int rows,
                       int cols) {
    auto kern = this->get_kernel("subtract");
    this->set_args(kern, A, B, C, rows, cols);
    return kern;
  }

  /**
   * Returns kernel to take the transpose of a matrix on the GPU.
   *
   * @param[out] B Buffer of the output matrix to hold transpose of A.
   * @param[in] A Buffer of the input matrix to transpose into B.
   * @param rows The number of rows for A.
   * @param cols The number of columns for A.
   *
   * @note This kernel uses the helper macros available in helpers.cl.
   */
  inline auto transpose(cl::Buffer A, cl::Buffer B, int rows, int cols) {
    auto kern = this->get_kernel("transpose");
    this->set_args(kern, A, B, rows, cols);
    return kern;
  }

  /**
   * Returns kernel to stores zeros in a matrix on the GPU.
   * Supports writing zeroes to the lower and upper triangular or
   * the whole matrix.
   *
   * @param[out] A Buffer of a matrix
   * @param rows Number of rows for matrix A
   * @param cols Number of columns for matrix A
   * @tparam triangular_view optional parameter that describes where to assign
   * zeros: LOWER - lower triangular UPPER - upper triangular if the part
   * parameter is not specified, zeros are assigned to the whole matrix.
   *
   * @note  This kernel uses the helper macros available in helpers.cl.
   */
  template <TriangularViewGPU triangular_view = TriangularViewGPU::Entire>
  inline auto zeros(cl::Buffer A, int rows, int cols) {
    auto kern = this->get_kernel("zeros");
    this->set_args(kern, A, rows, cols, static_cast<int>(triangular_view));
    return kern;
  }
};
kernel_cl kernel_cl;
}  // namespace math
}  // namespace stan

#endif
#endif
