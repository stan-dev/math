#ifndef STAN_MATH_OPENCL_SCALAR_CL_HPP
#define STAN_MATH_OPENCL_SCALAR_CL_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/opencl/copy.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <CL/opencl.hpp>
#include <tbb/concurrent_vector.h>
#include <utility>

namespace stan {
namespace math {
namespace opencl {

/** \ingroup opencl
 * A double that lives on the OpenCL device.
 *
 * The value is stored in a 1x1 `matrix_cl<double>` so it can be read and
 * written by kernels without being transferred to the host. There are no
 * implicit conversions to or from host values: use the explicit
 * constructor to upload a value and `opencl::to_host()` to read it back.
 *
 * Copies duplicate the value on the device. Moves transfer ownership of the
 * device buffer.
 */
template <>
class ScalarCl<double> {
 private:
  matrix_cl<double> buf_;

 public:
  using value_type = double;

  /**
   * Constructs a device scalar holding zero. No host data is transferred.
   */
  ScalarCl() : buf_(1, 1) { buf_.setZero(); }

  /**
   * Uploads a host value to the device.
   *
   * The value is passed to `matrix_cl` as an rvalue so the write blocks until
   * it completes; an lvalue would enqueue a non-blocking write reading from
   * host memory that may no longer exist when the write executes.
   * @param value value to upload
   */
  explicit ScalarCl(double value)
      : buf_(std::move(value), matrix_cl_view::Entire) {}

  ScalarCl(const ScalarCl& other) = default;
  ScalarCl(ScalarCl&& other) = default;
  ScalarCl& operator=(const ScalarCl& other) = default;
  ScalarCl& operator=(ScalarCl&& other) = default;

  /**
   * @return the backing 1x1 `matrix_cl`
   */
  inline const matrix_cl<double>& matrix() const noexcept { return buf_; }
  /**
   * @return the backing 1x1 `matrix_cl`
   */
  inline matrix_cl<double>& matrix() noexcept { return buf_; }

  /**
   * @return the OpenCL buffer holding the value
   */
  inline const cl::Buffer& buffer() const noexcept { return buf_.buffer(); }
  /**
   * @return the OpenCL buffer holding the value
   */
  inline cl::Buffer& buffer() noexcept { return buf_.buffer(); }

  /**
   * @return events of all operations writing to the value
   */
  inline const tbb::concurrent_vector<cl::Event>& write_events() const {
    return buf_.write_events();
  }
  /**
   * @return events of all operations reading the value
   */
  inline const tbb::concurrent_vector<cl::Event>& read_events() const {
    return buf_.read_events();
  }
  /**
   * @return events of all operations reading or writing the value
   */
  inline tbb::concurrent_vector<cl::Event> read_write_events() const {
    return buf_.read_write_events();
  }
  /**
   * Adds an event of an operation reading the value.
   * @param new_event event to add
   */
  inline void add_read_event(cl::Event new_event) const {
    buf_.add_read_event(std::move(new_event));
  }
  /**
   * Adds an event of an operation writing the value.
   * @param new_event event to add
   */
  inline void add_write_event(cl::Event new_event) const {
    buf_.add_write_event(std::move(new_event));
  }
  /**
   * Adds an event of an operation reading and writing the value.
   * @param new_event event to add
   */
  inline void add_read_write_event(cl::Event new_event) const {
    buf_.add_read_write_event(std::move(new_event));
  }
};

/** \ingroup opencl
 * Copies a device scalar to the host. Blocks until all writes to the scalar
 * have finished and the value has been read.
 * @param x device scalar
 * @return host copy of the value
 */
inline double to_host(const ScalarCl<double>& x) {
  return from_matrix_cl<double>(x.matrix());
}

}  // namespace opencl
}  // namespace math
}  // namespace stan

#endif
#endif
