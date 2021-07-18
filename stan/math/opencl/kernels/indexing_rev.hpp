#ifndef STAN_MATH_OPENCL_KERNELS_INDEXING_REV_HPP
#define STAN_MATH_OPENCL_KERNELS_INDEXING_REV_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_cl.hpp>
#include <stan/math/opencl/buffer_types.hpp>
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/opencl/kernels/device_functions/atomic_add_double.hpp>
#include <string>

namespace stan {
namespace math {
namespace opencl_kernels {

// \cond
static const std::string indexing_rev_global_atomic_kernel_code = STRINGIFY(
    // \endcond
    /** \ingroup opencl_kernels
     *
     * Increments adjoint of the indexing operation argument given the indices
     * and adjoints of the indexing result.
     *
     * This kernel uses global atomics and is the fastest for large sizes of the
     * indexed matrix.
     *
     * @param[in,out] adj adjoint to increment
     * @param index
     * @param res adjoint of the result of indexing
     * @param batch_size Number of matrices in the batch.
     * @note Code is a <code>const char*</code> held in
     * <code>add_batch_kernel_code.</code>
     */
    __kernel void indexing_rev(__global double* adj, const __global int* index,
                               const __global double* res, int size) {
      const int gid = get_global_id(0);
      const int gsize = get_global_size(0);
      for (int i = gid; i < size; i += gsize) {
        atomic_add_double(adj + index[i], res[i]);
      }
    }
    // \cond
);
// \endcond

/** \ingroup opencl_kernels
 * See the docs for \link kernels/add.hpp add_batch() \endlink
 */
const kernel_cl<in_out_buffer, in_buffer, in_buffer, int>
    indexing_rev_global_atomic("indexing_rev",
                               {atomic_add_double_device_function,
                                indexing_rev_global_atomic_kernel_code});

// \cond
static const std::string indexing_rev_local_atomic_kernel_code = STRINGIFY(
    // \endcond
    /** \ingroup opencl_kernels
     *
     * Increments adjoint of the indexing operation argument given the indices
     * and adjoints of the indexing result.
     *
     * This kernel uses local atomics and is the fastest for medium sizes of the
     * indexed matrix (and will not work for large sizes).
     *
     * @param[in,out] adj adjoint to increment
     * @param index
     * @param res adjoint of the result of indexing
     * @param batch_size Number of matrices in the batch.
     * @note Code is a <code>const char*</code> held in
     * <code>add_batch_kernel_code.</code>
     */
    __kernel void indexing_rev(__global double* adj, const __global int* index,
                               const __global double* res,
                               __local double* adj_loc, int index_size,
                               int adj_size) {
      const int gid = get_global_id(0);
      const int lid = get_local_id(0);
      const int gsize = get_global_size(0);
      const int lsize = get_local_size(0);
      for (int i = lid; i < adj_size; i += lsize) {
        adj_loc[i] = 0;
      }
      barrier(CLK_LOCAL_MEM_FENCE);
      for (int i = gid; i < index_size; i += gsize) {
        local_atomic_add_double(adj_loc + index[i], res[i]);
      }
      barrier(CLK_LOCAL_MEM_FENCE);
      for (int i = lid; i < adj_size; i += lsize) {
        atomic_add_double(adj + i, adj_loc[i]);
      }
    }
    // \cond
);
// \endcond

/** \ingroup opencl_kernels
 * See the docs for \link kernels/add.hpp add_batch() \endlink
 */
const kernel_cl<in_out_buffer, in_buffer, in_buffer, cl::LocalSpaceArg, int,
                int>
    indexing_rev_local_atomic("indexing_rev",
                              {atomic_add_double_device_function,
                               indexing_rev_local_atomic_kernel_code});

// \cond
static const std::string indexing_rev_local_independent_kernel_code = STRINGIFY(
    // \endcond
    /** \ingroup opencl_kernels
     *
     * Increments adjoint of the indexing operation argument given the indices
     * and adjoints of the indexing result.
     *
     * This kernel makes each thread build its own copy of the adjoints before
     * combining them. It is the fastest (and only works for) small size of the
     * indexed matrix.
     *
     * @param[in,out] adj adjoint to increment
     * @param index
     * @param res adjoint of the result of indexing
     * @param batch_size Number of matrices in the batch.
     * @note Code is a <code>const char*</code> held in
     * <code>add_batch_kernel_code.</code>
     */
    __kernel void indexing_rev(__global double* adj, const __global int* index,
                               const __global double* res,
                               __local double* adj_loc, int index_size,
                               int adj_size) {
      const int gid = get_global_id(0);
      const int lid = get_local_id(0);
      const int gsize = get_global_size(0);
      const int lsize = get_local_size(0);
      for (int i = lid; i < adj_size * lsize; i += lsize) {
        adj_loc[i] = 0;
      }
      barrier(CLK_LOCAL_MEM_FENCE);
      for (int i = gid; i < index_size; i += gsize) {
        adj_loc[index[i] + lid * adj_size] += res[i];
      }
      barrier(CLK_LOCAL_MEM_FENCE);
      for (int i = lid; i < adj_size; i += lsize) {
        double p = adj_loc[i + adj_size];
        for (int j = 2; j < lsize; j++) {
          p += adj_loc[i + j * adj_size];
        }
        adj_loc[i] += p;
      }
      barrier(CLK_LOCAL_MEM_FENCE);
      for (int i = lid; i < adj_size; i += lsize) {
        atomic_add_double(adj + i, adj_loc[i]);
      }
    }
    // \cond
);
// \endcond

/** \ingroup opencl_kernels
 * See the docs for \link kernels/add.hpp add_batch() \endlink
 */
const kernel_cl<in_out_buffer, in_buffer, in_buffer, cl::LocalSpaceArg, int,
                int>
    indexing_rev_local_independent(
        "indexing_rev", {atomic_add_double_device_function,
                         indexing_rev_local_independent_kernel_code});
}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan
#endif
#endif
