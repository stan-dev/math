#ifndef STAN_MATH_OPENCL_BUFFER_TYPES_HPP
#define STAN_MATH_OPENCL_BUFFER_TYPES_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/constants.hpp>
#include <CL/cl.hpp>

namespace stan {
namespace math {
namespace opencl_kernels {

// A read buffer will only add events to the read_write stack.
struct read_buffer {
  typedef cl::Buffer buffer;
};

// Write buffers will add events to the write event stack.
struct write_buffer {
  typedef cl::Buffer buffer;
};

namespace internal {

/**
 * meta template struct for changing read/write buffer argument types to
 * cl::Buffer types.
 * @tparam T A template typename that for cases of non-read/write buffer cases
 * will return a typedef holding only it's original type. For read and write
 * buffers this will return the a cl::Buffer.
 */
template <typename T = cl::Buffer>
struct to_buffer {
  typedef T type;
};

template <>
struct to_buffer<read_buffer> {
  typedef cl::Buffer type;
};

template <>
struct to_buffer<write_buffer> {
  typedef cl::Buffer type;
};

/**
 * meta template struct for changing cl::Buffer argument types to matrix_cl.
 * @tparam T A template typename that for cases of cl::Buffer will return a
 * typedef with a matrix_cl type. Otherwise will return a typedef with the
 * input's type.
 */
template <typename T = cl::Buffer>
struct to_matrix {
  typedef T type;
};

template <>
struct to_matrix<cl::Buffer> {
  typedef matrix_cl type;
};

/**
 * meta template struct for changing read and write buffer argument types to
 * matrix_cl. typedef with a matrix_cl type. Otherwise will return a typedef
 * with the input's type.
 */
template <>
struct to_matrix<read_buffer> {
  typedef matrix_cl type;
};

template <>
struct to_matrix<write_buffer> {
  typedef matrix_cl type;
};

}  // namespace internal
}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan

#endif
#endif
