#ifndef STAN_MATH_GPU_CL_HPP
#define STAN_MATH_GPU_CL_HPP
#ifdef STAN_OPENCL
#include <stan/math/gpu/opencl_context.hpp>
#include <stan/math/gpu/constants.hpp>
#include <CL/cl.hpp>
#include <string>
#include <map>
#include <vector>

namespace stan {
namespace math {
class kernel_cl; // Declare here so kernel_cl_base can see the class
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
  std::string helpers = // Helper macros for the kernels.
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
	 * copy  | Copy matrix A on the GPU to matrix B.
	 * transpose  | Take the transpose of a matrix.
	 * zeros | Make a lower, upper, or full matrix of zeros.
	 * identity | Create a NxM identity matrix.
	 * copy_triangular | Copy the lower or upper triangular matrix of A to B.
	 * copy_triangular_transposed | Copy the transpose lower/upper triangular.
	 * copy_submatrix | Copy a subset of matrix A into B.
	 * add | C = A + B
	 * subtract | C = A - B
	 * is_nan | Check if a matrix on the GPU contains nan values.
	 * is_zero_on_diagonal | Check if a matrix has zeros on the diagonal.
	 * is_symmetric | Check if a matrix is symmetric.
	 */
  const map_table kernel_table
      = {{"copy", copy_matrix},
         {"transpose", transpose_matrix},
         {"zeros", zeros_matrix},
         {"identity", identity_matrix},
         {"copy_triangular", copy_triangular_matrix},
         {"copy_triangular_transposed", copy_triangular_transposed_matrix},
         {"copy_submatrix", copy_submatrix},
         {"add", add_symmetric},
         {"subtract", subtract_symmetric},
         {"is_nan", check_nan},
         {"is_zero_on_diagonal", check_diagonal_zeros},
         {"is_symmetric", check_symmetric}};
  typedef std::map<const char*, cl::Kernel> map_kern;
  map_kern kernels;  // The compiled kernels

  static kernel_cl_base& getInstance() {
    static kernel_cl_base instance_;
    return instance_;
  }
  kernel_cl_base(kernel_cl_base const&) = delete;
  void operator=(kernel_cl_base const&) = delete;
};
}

/**
 * The adapter class that can access the <code>kernel_cl_base</code> class
 *
 *  This class operates as an API to the <code>kernel_cl_base</code>
 * class. Kernals can be called
 */
class kernel_cl {
 public:
  cl::Kernel compiled_;
  /**
   * Returns a kernel
   * @param kernel_name The name of a kernel specified in the
	 *  <code>kernel_table</code>. See the docs of <code>kernel_table</code> for
	 *  all available kernels.
   * @throw std::domain_error
   */
  explicit kernel_cl(const char* kernel_name) {
    if (internal::kernel_cl_base::getInstance().kernels.count((kernel_name)) == 0) {
      domain_error("compiling kernels", kernel_name, " kernel does not exist",
                   "");
    } else {
      compiled_ = internal::kernel_cl_base::getInstance().kernels[(kernel_name)];
    }
  }

  /**
   * return all compiled kernels.
   */
  inline internal::kernel_cl_base::map_kern kernels() {
    return internal::kernel_cl_base::getInstance().kernels;
  }

	/**
   * Adds arguments to an OpenCL kernel.
   *
   * @param args The arguments to mote to the OpenCL kernel.
   * @tparam Args The types of <code>extra_args</code>.
   * @note Comes from:
   * simpleopencl.blogspot.com/2013/04/calling-kernels-with-large-number-of.html
   */
	template <typename... Args>
	inline void set_args(const Args&... args) {
		this->recursive_args(this->compiled_, 0, args...);
	}

private:
  /**
   * Terminating function for recursively setting arguments in an OpenCL kernel.
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
   * @note Comes from:
   * simpleopencl.blogspot.com/2013/04/calling-kernels-with-large-number-of.html
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
   * @note Comes from:
   * simpleopencl.blogspot.com/2013/04/calling-kernels-with-large-number-of.html
   */
  template <typename... Args>
  inline void set_args(cl::Kernel& kernel, const Args&... args) {
    this->recursive_args(kernel, 0, args...);
  }
};
}  // namespace math
}  // namespace stan

#endif
#endif
