#ifndef STAN_MATH_OPENCL_REV_DESERIALIZER_HPP
#define STAN_MATH_OPENCL_REV_DESERIALIZER_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>

namespace stan::math {


// Example function to compute the padded offset.
// alignment_in_doubles is the alignment requirement in units of doubles.
inline size_t compute_padded_offset(size_t current_offset, size_t alignment_in_doubles) {
  size_t remainder = current_offset % alignment_in_doubles;
  if (remainder == 0)
    return current_offset;
  else
    return current_offset + (alignment_in_doubles - remainder);
}

/**
 * Create a `var_value<matrix_cl<double>>` with the correct padding and values
 * for the subbuffers.
 *
 * @tparam Model The model class which has methods `get_param_names` and `get_dims`.
 * @param model A pointer to a Stan model.
 * @param values The Eigen vector whose values will be placed in the correct positions
 * for the values of the `var_value<matrix_cl<double>>`
 * @return a var_value<matrix_cl<double>> whose buffer is allocated with proper padding.
 */
inline auto make_deserializer_data(const std::vector<std::vector<size_t>>& dims, const Eigen::VectorXd& values) {

  /*
   * `CL_DEVICE_MEM_BASE_ADDR_ALIGN` is the alignment requirement (in bits) for
   * sub-buffer offsets. The minimum value
   * is the size (in bits) of the largest OpenCL built-in data type supported
   * by the device (long16 in FULL profile, long16 or int16 in EMBEDDED profile)
   * for devices that are not of type CL_DEVICE_TYPE_CUSTOM.
   * The standard range for this is from 128 to 4096 bits (16 to 512 bytes).
   */
  const size_t alignment_in_doubles = (opencl_context.device()[0].getInfo<CL_DEVICE_MEM_BASE_ADDR_ALIGN>() / 8)/sizeof(double);

  // STEP 3. Compute padded offsets for each parameter.
  size_t total_doubles = 0;
  std::vector<size_t> offsets;  // will hold the starting offset (in doubles) for each parameter
  for (const auto& dim : dims) {
    // Compute the number of elements for this parameter.
    size_t param_size = 1;
    for (size_t d : dim) {
      param_size *= d;
    }
    // Ensure the starting offset is aligned.
    total_doubles = compute_padded_offset(total_doubles, alignment_in_doubles);
    offsets.push_back(total_doubles);
    total_doubles += param_size;
  }
  // Now total_doubles is the size of our padded global buffer.

  // STEP 4. Allocate a padded host vector and copy in the values.
  Eigen::VectorXd padded_data(total_doubles);
  padded_data.setZero();  // initialize to zero so padding regions are defined
  size_t input_pos = 0;
  for (size_t i = 0; i < dims.size(); ++i) {
    size_t param_size = 1;
    for (size_t d : dims[i])
      param_size *= d;
    // Ensure we have enough input values.
    assert(input_pos + param_size <= static_cast<size_t>(values.size()));
    // Copy the parameter's data from the input vector into the correct subregion.
    padded_data.segment(offsets[i], param_size) = values.segment(input_pos, param_size);
    input_pos += param_size;
  }
  // Optionally, you might also want to check that input_pos == values.size()

  // STEP 5. Create a matrix_cl from the padded_data.
  // We assume that the matrix_cl constructor can accept an Eigen vector
  // and will allocate a buffer on the GPU (with SVM or via OpenCL buffers).
  stan::math::matrix_cl<double> gpu_buffer(padded_data);

  // STEP 6. Wrap the matrix_cl in a var_value and return it.
  return stan::math::var_value<stan::math::matrix_cl<double>>(std::move(gpu_buffer));
}

template <typename Mat>
class deserializer_cl {
  std::reference_wrapper<Mat> mat_;
  std::size_t alignment_in_doubles_{64};
  std::size_t current_offset_{0};
  deserializer_cl(Mat& mat, std::size_t alignment_in_doubles)
      : mat_(mat), alignment_in_doubles_(alignment_in_doubles), current_offset_(0) {}

  explicit deserializer_cl(Mat& mat)
      : mat_(mat),
      alignment_in_doubles_((opencl_context.device()[0].getInfo<CL_DEVICE_MEM_BASE_ADDR_ALIGN>() / 8)/sizeof(double)),
      current_offset_(0) {}
  inline void align_offset() {
    const size_t remainder = current_offset_ % alignment_in_doubles_;
    if (remainder == 0)
      return;
    else
      current_offset_ = current_offset_ + (alignment_in_doubles_ - remainder);
  }
}

}
#endif
#endif
