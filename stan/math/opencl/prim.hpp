#ifndef STAN_MATH_OPENCL_PRIM_HPP
#define STAN_MATH_OPENCL_PRIM_HPP
#ifdef STAN_OPENCL

/**
 * \defgroup opencl OpenCL
 * Stan's OpenCL backend allows for computation to be executed in parallel
 *  on a GPU or in multithreaded CPUs. It is meant to easily conform with Eigen
 * such that you can create and read from a `matrix_cl` by doing
 *
 *```cpp
 * Eigen::MatrixXd A_eig = Eigen::MatrixXd::Random(10, 10);
 * matrix_cl<double> A(A_eig);
 * matrix_cl<double> B = to_matrix_cl(A_eig);
 * matrix_cl<double> C = cholesky_decompose(A * B);
 * // Read back to eigen matrix.
 * Eigen::MatrixXd C_eig = from_matrix_cl(C);
 *
 * // Also for vectors and raw pointers of pointers
 * std::vector<var> A_vec(10, 0);
 * matrix_cl<var> B_var(A_vec, 10, 1);
 *
 * vari** A_vari= // fill
 * matrix_cl<var> B_vari(A_vari, 10, 1);
 *
 *```
 *
 * Execution is performed in async and Kernel operations are compounded and
 * compiled Just In Time. This allows for a low amount of overhead when passing
 * data to and from the device and executing computations.
 *
 * For more details see the paper on Arvix.
 * https://arxiv.org/pdf/1907.01063.pdf
 */

/**
 * \ingroup opencl
 * \defgroup error_checks_opencl Error Checks
 */

/**
 * \ingroup opencl
 * \defgroup kernel_executor_opencl Kernel Executor
 * The kernel executor allows OpenCL kernels to be executed in async. GPUs
 * have the capability to perform reads, writes, and computation at the same
 * time. In order to maximize the throughput to the device the Kernel
 * Executor assigns matrices read and write events. Write events are blocking
 * in that no further operations can be completed until the write event
 * is finished. However, read events can happen in async together.
 * Take the following for example
 *
 *```cpp
 * matrix_cl<double> A = from_eigen_cl(A_eig);
 * matrix_cl<double> B = from_eigen_cl(B_eig);
 * matrix_cl<double> C = A * B;
 * matrix_cl<double> D = A + B;
 * matrix_cl<double> E = C + D;
 *```
 * In the above, When `A` and `B` are created from the Eigen matrices `A_eig`
 *and `B_eig`. they are both assigned write events to their write event stack.
 * `C` and `D` depend on `A` and `B` while `E` depends on
 * `C` and `D`. When executing `C`'s operation, `A` and `B` are assigned
 * events to their read event stack while `C` is assigned an event to it's write
 *event stack. Once `A` and `B` have finished their write event the kernel to
 *compute `C` can begin. The excution to create `D` also waits for the write
 *events of `A` and `B`, but does not have to wait for the execution of `C` to
 *finish. Executing `E` requires waiting for for the write events of both `C`
 *and `D`.
 *
 */

/**
 * \ingroup opencl
 * \defgroup opencl_kernels Custom OpenCL kernels
 */

/**
 * \ingroup opencl
 * \defgroup prim_fun_opencl OpenCL overloads of stan/math/prim functions
 */
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/opencl_context.hpp>
#include <stan/math/opencl/matrix_cl.hpp>

#include <stan/math/opencl/scalar_type.hpp>
#include <stan/math/opencl/copy.hpp>
#include <stan/math/opencl/cholesky_decompose.hpp>
#include <stan/math/opencl/is_constant.hpp>
#include <stan/math/opencl/tri_inverse.hpp>
#include <stan/math/opencl/multiply_transpose.hpp>
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/opencl/pinned_matrix.hpp>
#include <stan/math/opencl/plain_type.hpp>
#include <stan/math/opencl/ref_type_for_opencl.hpp>
#include <stan/math/opencl/ref_type.hpp>
#include <stan/math/opencl/to_ref_for_opencl.hpp>
#include <stan/math/opencl/value_type.hpp>
#include <stan/math/opencl/zeros_strict_tri.hpp>
#include <stan/math/opencl/qr_decomposition.hpp>

#include <stan/math/opencl/prim_constraint.hpp>

#include <stan/math/opencl/prim/add_diag.hpp>
#include <stan/math/opencl/prim/append_array.hpp>
#include <stan/math/opencl/prim/bernoulli_cdf.hpp>
#include <stan/math/opencl/prim/bernoulli_lccdf.hpp>
#include <stan/math/opencl/prim/bernoulli_lcdf.hpp>
#include <stan/math/opencl/prim/bernoulli_lpmf.hpp>
#include <stan/math/opencl/prim/bernoulli_logit_lpmf.hpp>
#include <stan/math/opencl/prim/bernoulli_logit_glm_lpmf.hpp>
#include <stan/math/opencl/prim/beta_binomial_lpmf.hpp>
#include <stan/math/opencl/prim/beta_lpdf.hpp>
#include <stan/math/opencl/prim/beta_proportion_lpdf.hpp>
#include <stan/math/opencl/prim/binomial_logit_lpmf.hpp>
#include <stan/math/opencl/prim/binomial_lpmf.hpp>
#include <stan/math/opencl/prim/binomial_logit_glm_lpmf.hpp>
#include <stan/math/opencl/prim/block.hpp>
#include <stan/math/opencl/prim/categorical_logit_glm_lpmf.hpp>
#include <stan/math/opencl/prim/cauchy_cdf.hpp>
#include <stan/math/opencl/prim/cauchy_lccdf.hpp>
#include <stan/math/opencl/prim/cauchy_lcdf.hpp>
#include <stan/math/opencl/prim/cauchy_lpdf.hpp>
#include <stan/math/opencl/prim/chi_square_lpdf.hpp>
#include <stan/math/opencl/prim/cholesky_decompose.hpp>
#include <stan/math/opencl/prim/col.hpp>
#include <stan/math/opencl/prim/cols.hpp>
#include <stan/math/opencl/prim/columns_dot_product.hpp>
#include <stan/math/opencl/prim/columns_dot_self.hpp>
#include <stan/math/opencl/prim/crossprod.hpp>
#include <stan/math/opencl/prim/cumulative_sum.hpp>
#include <stan/math/opencl/prim/diag_matrix.hpp>
#include <stan/math/opencl/prim/diag_pre_multiply.hpp>
#include <stan/math/opencl/prim/diag_post_multiply.hpp>
#include <stan/math/opencl/prim/dims.hpp>
#include <stan/math/opencl/prim/dirichlet_lpdf.hpp>
#include <stan/math/opencl/prim/distance.hpp>
#include <stan/math/opencl/prim/divide.hpp>
#include <stan/math/opencl/prim/divide_columns.hpp>
#include <stan/math/opencl/prim/dot_product.hpp>
#include <stan/math/opencl/prim/dot_self.hpp>
#include <stan/math/opencl/prim/double_exponential_cdf.hpp>
#include <stan/math/opencl/prim/double_exponential_lccdf.hpp>
#include <stan/math/opencl/prim/double_exponential_lcdf.hpp>
#include <stan/math/opencl/prim/double_exponential_lpdf.hpp>
#include <stan/math/opencl/prim/eigenvalues_sym.hpp>
#include <stan/math/opencl/prim/eigenvectors_sym.hpp>
#include <stan/math/opencl/prim/exp_mod_normal_cdf.hpp>
#include <stan/math/opencl/prim/exp_mod_normal_lccdf.hpp>
#include <stan/math/opencl/prim/exp_mod_normal_lcdf.hpp>
#include <stan/math/opencl/prim/exp_mod_normal_lpdf.hpp>
#include <stan/math/opencl/prim/exponential_cdf.hpp>
#include <stan/math/opencl/prim/exponential_lccdf.hpp>
#include <stan/math/opencl/prim/exponential_lcdf.hpp>
#include <stan/math/opencl/prim/exponential_lpdf.hpp>
#include <stan/math/opencl/prim/frechet_cdf.hpp>
#include <stan/math/opencl/prim/frechet_lccdf.hpp>
#include <stan/math/opencl/prim/frechet_lcdf.hpp>
#include <stan/math/opencl/prim/frechet_lpdf.hpp>
#include <stan/math/opencl/prim/gamma_lpdf.hpp>
#include <stan/math/opencl/prim/gp_dot_prod_cov.hpp>
#include <stan/math/opencl/prim/gp_exponential_cov.hpp>
#include <stan/math/opencl/prim/gp_exp_quad_cov.hpp>
#include <stan/math/opencl/prim/gp_matern32_cov.hpp>
#include <stan/math/opencl/prim/gp_matern52_cov.hpp>
#include <stan/math/opencl/prim/gumbel_cdf.hpp>
#include <stan/math/opencl/prim/gumbel_lccdf.hpp>
#include <stan/math/opencl/prim/gumbel_lcdf.hpp>
#include <stan/math/opencl/prim/gumbel_lpdf.hpp>
#include <stan/math/opencl/prim/head.hpp>
#include <stan/math/opencl/prim/identity_matrix.hpp>
#include <stan/math/opencl/prim/inv.hpp>
#include <stan/math/opencl/prim/inv_chi_square_lpdf.hpp>
#include <stan/math/opencl/prim/inv_cloglog.hpp>
#include <stan/math/opencl/prim/inv_gamma_lpdf.hpp>
#include <stan/math/opencl/prim/inv_sqrt.hpp>
#include <stan/math/opencl/prim/log_mix.hpp>
#include <stan/math/opencl/prim/log_softmax.hpp>
#include <stan/math/opencl/prim/logistic_cdf.hpp>
#include <stan/math/opencl/prim/logistic_lccdf.hpp>
#include <stan/math/opencl/prim/logistic_lcdf.hpp>
#include <stan/math/opencl/prim/logistic_lpdf.hpp>
#include <stan/math/opencl/prim/log_sum_exp.hpp>
#include <stan/math/opencl/prim/lognormal_cdf.hpp>
#include <stan/math/opencl/prim/lognormal_lccdf.hpp>
#include <stan/math/opencl/prim/lognormal_lcdf.hpp>
#include <stan/math/opencl/prim/lognormal_lpdf.hpp>
#include <stan/math/opencl/prim/matrix_power.hpp>
#include <stan/math/opencl/prim/mdivide_left_tri_low.hpp>
#include <stan/math/opencl/prim/mdivide_right_tri_low.hpp>
#include <stan/math/opencl/prim/mean.hpp>
#include <stan/math/opencl/prim/multi_normal_cholesky_lpdf.hpp>
#include <stan/math/opencl/prim/multiply_lower_tri_self_transpose.hpp>
#include <stan/math/opencl/prim/neg_binomial_lpmf.hpp>
#include <stan/math/opencl/prim/neg_binomial_2_lpmf.hpp>
#include <stan/math/opencl/prim/neg_binomial_2_log_lpmf.hpp>
#include <stan/math/opencl/prim/neg_binomial_2_log_glm_lpmf.hpp>
#include <stan/math/opencl/prim/normal_id_glm_lpdf.hpp>
#include <stan/math/opencl/prim/normal_cdf.hpp>
#include <stan/math/opencl/prim/normal_lccdf.hpp>
#include <stan/math/opencl/prim/normal_lcdf.hpp>
#include <stan/math/opencl/prim/normal_lpdf.hpp>
#include <stan/math/opencl/prim/num_elements.hpp>
#include <stan/math/opencl/prim/ordered_logistic_glm_lpmf.hpp>
#include <stan/math/opencl/prim/ordered_logistic_lpmf.hpp>
#include <stan/math/opencl/prim/pareto_cdf.hpp>
#include <stan/math/opencl/prim/pareto_lccdf.hpp>
#include <stan/math/opencl/prim/pareto_lcdf.hpp>
#include <stan/math/opencl/prim/pareto_lpdf.hpp>
#include <stan/math/opencl/prim/pareto_type_2_cdf.hpp>
#include <stan/math/opencl/prim/pareto_type_2_lccdf.hpp>
#include <stan/math/opencl/prim/pareto_type_2_lcdf.hpp>
#include <stan/math/opencl/prim/pareto_type_2_lpdf.hpp>
#include <stan/math/opencl/prim/poisson_log_glm_lpmf.hpp>
#include <stan/math/opencl/prim/poisson_log_lpmf.hpp>
#include <stan/math/opencl/prim/poisson_lpmf.hpp>
#include <stan/math/opencl/prim/prod.hpp>
#include <stan/math/opencl/prim/qr_Q.hpp>
#include <stan/math/opencl/prim/qr_R.hpp>
#include <stan/math/opencl/prim/qr_thin_Q.hpp>
#include <stan/math/opencl/prim/qr_thin_R.hpp>
#include <stan/math/opencl/prim/rank.hpp>
#include <stan/math/opencl/prim/rayleigh_cdf.hpp>
#include <stan/math/opencl/prim/rayleigh_lccdf.hpp>
#include <stan/math/opencl/prim/rayleigh_lcdf.hpp>
#include <stan/math/opencl/prim/rayleigh_lpdf.hpp>
#include <stan/math/opencl/prim/rep_array.hpp>
#include <stan/math/opencl/prim/rep_matrix.hpp>
#include <stan/math/opencl/prim/rep_row_vector.hpp>
#include <stan/math/opencl/prim/rep_vector.hpp>
#include <stan/math/opencl/prim/reverse.hpp>
#include <stan/math/opencl/prim/row.hpp>
#include <stan/math/opencl/prim/rows.hpp>
#include <stan/math/opencl/prim/rows_dot_product.hpp>
#include <stan/math/opencl/prim/rows_dot_self.hpp>
#include <stan/math/opencl/prim/scaled_inv_chi_square_lpdf.hpp>
#include <stan/math/opencl/prim/sd.hpp>
#include <stan/math/opencl/prim/segment.hpp>
#include <stan/math/opencl/prim/sign.hpp>
#include <stan/math/opencl/prim/size.hpp>
#include <stan/math/opencl/prim/softmax.hpp>
#include <stan/math/opencl/prim/sort_asc.hpp>
#include <stan/math/opencl/prim/sort_desc.hpp>
#include <stan/math/opencl/prim/squared_distance.hpp>
#include <stan/math/opencl/prim/sub_col.hpp>
#include <stan/math/opencl/prim/sub_row.hpp>
#include <stan/math/opencl/prim/std_normal_cdf.hpp>
#include <stan/math/opencl/prim/std_normal_lccdf.hpp>
#include <stan/math/opencl/prim/std_normal_lcdf.hpp>
#include <stan/math/opencl/prim/std_normal_lpdf.hpp>
#include <stan/math/opencl/prim/student_t_lpdf.hpp>
#include <stan/math/opencl/prim/skew_double_exponential_cdf.hpp>
#include <stan/math/opencl/prim/skew_double_exponential_lcdf.hpp>
#include <stan/math/opencl/prim/skew_double_exponential_lccdf.hpp>
#include <stan/math/opencl/prim/skew_double_exponential_lpdf.hpp>
#include <stan/math/opencl/prim/skew_normal_lpdf.hpp>
#include <stan/math/opencl/prim/sum.hpp>
#include <stan/math/opencl/prim/symmetrize_from_lower_tri.hpp>
#include <stan/math/opencl/prim/symmetrize_from_upper_tri.hpp>
#include <stan/math/opencl/prim/tail.hpp>
#include <stan/math/opencl/prim/tcrossprod.hpp>
#include <stan/math/opencl/prim/to_array_1d.hpp>
#include <stan/math/opencl/prim/to_array_2d.hpp>
#include <stan/math/opencl/prim/to_matrix.hpp>
#include <stan/math/opencl/prim/to_row_vector.hpp>
#include <stan/math/opencl/prim/to_vector.hpp>
#include <stan/math/opencl/prim/trace.hpp>
#include <stan/math/opencl/prim/uniform_cdf.hpp>
#include <stan/math/opencl/prim/uniform_lccdf.hpp>
#include <stan/math/opencl/prim/uniform_lcdf.hpp>
#include <stan/math/opencl/prim/uniform_lpdf.hpp>
#include <stan/math/opencl/prim/variance.hpp>
#include <stan/math/opencl/prim/weibull_cdf.hpp>
#include <stan/math/opencl/prim/weibull_lccdf.hpp>
#include <stan/math/opencl/prim/weibull_lcdf.hpp>
#include <stan/math/opencl/prim/weibull_lpdf.hpp>

#include <stan/math/opencl/err.hpp>

#endif
#endif
