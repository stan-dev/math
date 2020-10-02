#ifndef STAN_MATH_OPENCL_OPENCL
#define STAN_MATH_OPENCL_OPENCL
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

#include <stan/math/opencl/copy.hpp>
#include <stan/math/opencl/copy_triangular.hpp>
#include <stan/math/opencl/cholesky_decompose.hpp>
#include <stan/math/opencl/diagonal_multiply.hpp>
#include <stan/math/opencl/identity.hpp>
#include <stan/math/opencl/is_constant.hpp>
#include <stan/math/opencl/tri_inverse.hpp>
#include <stan/math/opencl/multiply_transpose.hpp>
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/opencl/plain_type.hpp>
#include <stan/math/opencl/scalar_type.hpp>
#include <stan/math/opencl/sub_block.hpp>
#include <stan/math/opencl/triangular_transpose.hpp>
#include <stan/math/opencl/value_type.hpp>
#include <stan/math/opencl/zeros.hpp>

#include <stan/math/opencl/prim/add.hpp>
#include <stan/math/opencl/prim/bernoulli_logit_glm_lpmf.hpp>
#include <stan/math/opencl/prim/categorical_logit_glm_lpmf.hpp>
#include <stan/math/opencl/prim/cholesky_decompose.hpp>
#include <stan/math/opencl/prim/col.hpp>
#include <stan/math/opencl/prim/cols.hpp>
#include <stan/math/opencl/prim/crossprod.hpp>
#include <stan/math/opencl/prim/dims.hpp>
#include <stan/math/opencl/prim/divide.hpp>
#include <stan/math/opencl/prim/divide_columns.hpp>
#include <stan/math/opencl/prim/gp_exp_quad_cov.hpp>
#include <stan/math/opencl/prim/inv.hpp>
#include <stan/math/opencl/prim/inv_cloglog.hpp>
#include <stan/math/opencl/prim/inv_sqrt.hpp>
#include <stan/math/opencl/prim/mdivide_left_tri_low.hpp>
#include <stan/math/opencl/prim/mdivide_right_tri_low.hpp>
#include <stan/math/opencl/prim/neg_binomial_2_log_glm_lpmf.hpp>
#include <stan/math/opencl/prim/normal_id_glm_lpdf.hpp>
#include <stan/math/opencl/prim/ordered_logistic_glm_lpmf.hpp>
#include <stan/math/opencl/prim/poisson_log_glm_lpmf.hpp>
#include <stan/math/opencl/prim/rep_matrix.hpp>
#include <stan/math/opencl/prim/rep_row_vector.hpp>
#include <stan/math/opencl/prim/rep_vector.hpp>
#include <stan/math/opencl/prim/row.hpp>
#include <stan/math/opencl/prim/rows.hpp>
#include <stan/math/opencl/prim/size.hpp>
#include <stan/math/opencl/prim/sum.hpp>
#include <stan/math/opencl/prim/tcrossprod.hpp>

#include <stan/math/opencl/err.hpp>

#endif
#endif
