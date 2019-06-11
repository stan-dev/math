#ifndef STAN_MATH_OPENCL_SUB_BLOCK_HPP
#define STAN_MATH_OPENCL_SUB_BLOCK_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/opencl_context.hpp>
#include <stan/math/opencl/constants.hpp>
#include <stan/math/opencl/kernels/sub_block.hpp>
#include <stan/math/prim/scal/err/domain_error.hpp>
#include <stan/math/opencl/rev/matrix_cl.hpp>
#include <CL/cl.hpp>
#include <vector>

namespace stan {
namespace math {

  template <TriangularViewCL triangular_view>
  inline void matrix_cl<var>::sub_block(const matrix_cl<var>& A, size_t A_i, size_t A_j,
                                   size_t this_i, size_t this_j, size_t nrows,
                                   size_t ncols) try {
    this->val().sub_block<triangular_view>(A.val(), A_i, A_j, this_i, this_j, nrows, ncols);
    this->adj().sub_block<triangular_view>(A.adj(), A_i, A_j, this_i, this_j, nrows, ncols);
  } catch (const cl::Error& e) {
   check_opencl_error("copy_submatrix", e);
  }

  }  // namespace math
  }  // namespace stan

  #endif
  #endif
