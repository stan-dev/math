#ifndef STAN_MATH_OPENCL_CACHE_COPY_HPP
#define STAN_MATH_OPENCL_CACHE_COPY_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/opencl_context.hpp>
#include <stan/math/opencl/kernel_cl.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <CL/cl.hpp>

namespace stan {
namespace math {
namespace internal {
template <int R, int C>
inline void cache_copy(cl::Buffer dst, const Eigen::Matrix<double, R, C>& src) {
  cl::Context& ctx = opencl_context.context();
  cl::CommandQueue queue = opencl_context.queue();
#ifdef STAN_OPENCL_CACHE
  if (src.opencl_buffer_() != NULL) {
    queue.enqueueCopyBuffer(src.opencl_buffer_, dst, 0, 0,
                            sizeof(double) * src.size());
  } else {
    try {
      src.opencl_buffer_
          = cl::Buffer(ctx, CL_MEM_READ_WRITE, sizeof(double) * src.size());
      /**
       * Writes the contents of src to the OpenCL buffer
       * starting at the offset 0
       * CL_TRUE denotes that the call is blocking
       * We do not want to execute any further kernels
       * on the device until we are sure that the data is transferred)
       */
      cl::Event copy_event;
      auto source_size = sizeof(double) * src.size();
      queue.enqueueWriteBuffer(dst, CL_TRUE, 0, source_size, src.data());
      queue.enqueueCopyBuffer(dst, src.opencl_buffer_, 0, 0, source_size, NULL,
                              &copy_event);
      copy_event.wait();
    } catch (const cl::Error& e) {
      check_opencl_error("copy Eigen->GPU", e);
    }
  }
#else
  queue.enqueueWriteBuffer(dst, CL_TRUE, 0, sizeof(double) * src.size(),
                           src.data());
#endif
}
}  // namespace internal
}  // namespace math
}  // namespace stan
#endif
#endif
