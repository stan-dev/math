#ifndef STAN_MATH_OPENCL_BUFFER_TYPES_HPP
#define STAN_MATH_OPENCL_BUFFER_TYPES_HPP
#ifdef STAN_OPENCL

namespace stan {
namespace math {
namespace opencl_kernels {
struct read_buffer {
  typedef cl::Buffer buffer;
  static const eventTypeCL event_type = eventTypeCL::read;
};

struct write_buffer {
  typedef cl::Buffer buffer;
  static const eventTypeCL event_type = eventTypeCL::write;
};

struct read_write_buffer {
  typedef cl::Buffer buffer;
  static const eventTypeCL event_type = eventTypeCL::read_write;
};

namespace internal {

/**
 * meta template struct for changing cl::Buffer argument types to matrix_cl.
 * @tparam T A template typename that for cases of cl::Buffer will return a
 * typedef with a matrix_cl type. Otherwise will return a typedef with the
 * input's type.
 */
template <typename T = cl::Buffer>
struct to_buffer {
  typedef T type;
};

/**
 * meta template struct for changing cl::Buffer argument types to matrix_cl.
 * typedef with a matrix_cl type. Otherwise will return a typedef with the
 * input's type.
 */
template <>
struct to_buffer<read_buffer> {
  typedef cl::Buffer type;
};

template <>
struct to_buffer<write_buffer> {
  typedef cl::Buffer type;
};

template <>
struct to_buffer<read_write_buffer> {
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
 * meta template struct for changing cl::Buffer argument types to matrix_cl.
 * typedef with a matrix_cl type. Otherwise will return a typedef with the
 * input's type.
 */
template <>
struct to_matrix<read_buffer> {
  typedef matrix_cl type;
};

template <>
struct to_matrix<write_buffer> {
  typedef matrix_cl type;
};

template <>
struct to_matrix<read_write_buffer> {
  typedef matrix_cl type;
};

}  // namespace internal
}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan

#endif
#endif
