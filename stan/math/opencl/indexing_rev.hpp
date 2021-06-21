#ifndef STAN_MATH_OPENCL_INDEXING_REV_HPP
#define STAN_MATH_OPENCL_INDEXING_REV_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/kernels/indexing_rev.hpp>
#include <stan/math/opencl/kernels/add.hpp>

namespace stan {
namespace math {

/**
 * Performs reverse pass for indexing operation on the OpenCL device. Depending
 * on the size of indexed matrix and the amount of local memory available on the
 * device selects the best kernel to use for the operation.
 *
 * @param[in,out] adj adjoint of the argument to indexing
 * @param idx indices
 * @param res adjoint of the result of the indexing operation
 */
void indexing_rev(matrix_cl<double>& adj, const matrix_cl<int>& idx,
                  const matrix_cl<double>& res) {
  int local_mem_size
      = opencl_context.device()[0].getInfo<CL_DEVICE_LOCAL_MEM_SIZE>();
  int preferred_work_groups
      = opencl_context.device()[0].getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>();
  int local_size = 64;
  int n_threads = preferred_work_groups * 16 * local_size;

  try {
    if (local_mem_size > sizeof(double) * adj.size() * local_size * 2) {
      stan::math::opencl_kernels::indexing_rev_local_independent(
          cl::NDRange(n_threads), cl::NDRange(local_size), adj, idx, res,
          cl::Local(sizeof(double) * adj.size() * local_size), res.size(),
          adj.size());
    } else if (local_mem_size > sizeof(double) * adj.size()) {
      stan::math::opencl_kernels::indexing_rev_local_atomic(
          cl::NDRange(n_threads), cl::NDRange(local_size), adj, idx, res,
          cl::Local(sizeof(double) * adj.size()), res.size(), adj.size());
    } else {
      stan::math::opencl_kernels::indexing_rev_global_atomic(
          cl::NDRange(n_threads), cl::NDRange(local_size), adj, idx, res,
          res.size());
    }
  } catch (cl::Error& e) {
    check_opencl_error("indexing reverse pass", e);
  }
}

}  // namespace math
}  // namespace stan

#endif
#endif
