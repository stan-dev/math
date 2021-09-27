#ifndef STAN_MATH_OPENCL_REV_CHOLESKY_DECOMPOSE_HPP
#define STAN_MATH_OPENCL_REV_CHOLESKY_DECOMPOSE_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/prim/cholesky_decompose.hpp>
#include <stan/math/opencl/prim/symmetrize_from_lower_tri.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * Returns the lower-triangular Cholesky factor (i.e., matrix
 * square root) of the specified square, symmetric reverse mode matrix on the
 * OpenCL device. The return value \f$L\f$ will be a lower-triangular matrix
 * such that the original matrix \f$A\f$ is given by <p>\f$A = L \times L^T\f$.
 * @param A Input square matrix
 * @return Square root of matrix.
 * @throw std::domain_error if m is not a symmetric matrix or
 *   if m is not positive definite (if m has more than 0 elements)
 */
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
inline var_value<matrix_cl<double>> cholesky_decompose(const var_value<T>& A) {
  check_cl("cholesky_decompose (OpenCL)", "A", A.val(), "not NaN")
      = !isnan(A.val());

  return make_callback_var(
      cholesky_decompose(A.val()),
      [A](vari_value<matrix_cl<double>>& L_A) mutable {
        int M_ = A.rows();
        int block_size
            = M_ / opencl_context.tuning_opts().cholesky_rev_block_partition;
        block_size = std::max(block_size, 8);
        block_size = std::min(
            block_size,
            opencl_context.tuning_opts().cholesky_rev_min_block_size);
        matrix_cl<double> A_adj = L_A.adj();
        for (int k = M_; k > 0; k -= block_size) {
          const int j = std::max(0, k - block_size);
          const int k_j_ind = k - j;
          const int m_k_ind = M_ - k;

          auto&& R_val = block_zero_based(L_A.val(), j, 0, k_j_ind, j);
          auto&& R_adj = block_zero_based(A_adj, j, 0, k_j_ind, j);
          matrix_cl<double> D_val
              = block_zero_based(L_A.val(), j, j, k_j_ind, k_j_ind);
          matrix_cl<double> D_adj
              = block_zero_based(A_adj, j, j, k_j_ind, k_j_ind);
          auto&& B_val = block_zero_based(L_A.val(), k, 0, m_k_ind, j);
          auto&& B_adj = block_zero_based(A_adj, k, 0, m_k_ind, j);
          auto&& C_val = block_zero_based(L_A.val(), k, j, m_k_ind, k_j_ind);
          auto&& C_adj = block_zero_based(A_adj, k, j, m_k_ind, k_j_ind);

          C_adj = C_adj * tri_inverse(D_val);
          B_adj = B_adj - C_adj * R_val;
          D_adj = D_adj - transpose(C_adj) * C_val;

          D_adj = symmetrize_from_lower_tri(transpose(D_val) * D_adj);
          D_val = transpose(tri_inverse(D_val));
          D_adj = symmetrize_from_lower_tri(D_val * transpose(D_val * D_adj));

          R_adj = R_adj - transpose(C_adj) * B_val - D_adj * R_val;
          diagonal(D_adj) = diagonal(D_adj) * 0.5;

          block_zero_based(A_adj, j, j, k_j_ind, k_j_ind) = D_adj;
        }
        A_adj.view(matrix_cl_view::Lower);
        A.adj() += A_adj;
      });
}

}  // namespace math
}  // namespace stan

#endif
#endif
