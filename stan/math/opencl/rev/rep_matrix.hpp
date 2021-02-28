#ifndef STAN_MATH_OPENCL_REV_REP_MATRIX_HPP
#define STAN_MATH_OPENCL_REV_REP_MATRIX_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/opencl/kernels/rep_matrix.hpp>
#include <stan/math/opencl/prim/rep_matrix.hpp>
#include <stan/math/opencl/prim/sum.hpp>

namespace stan {
namespace math {

/** \ingroup opencl
 * Creates a matrix_cl by replicating the given value of arithmetic type.
 *
 * @tparam T type of the result matrix
 * @param A the input value
 * @param n number of rows in the result matrix
 * @param m number of columns in the result matrix
 *
 * @return matrix_cl with replicated value from the input
 *
 * @throw <code>domain_error</code> if the
 * requested dimensions are negative
 *
 */
template <typename T_ret, require_var_vt<is_matrix_cl, T_ret>* = nullptr>
inline var_value<matrix_cl<double>> rep_matrix(const var& A, int n, int m) {
  return make_callback_var(rep_matrix<matrix_cl<double>>(A.val(), n, m),
                           [A](vari_value<matrix_cl<double>>& res) mutable {
                             A.adj() += sum(res.adj());
                           });
}

/** \ingroup opencl
 * Creates a matrix_cl by replicating the input
 * vector or row_vector.  The elements of the
 * vector or row_vector must be of arithmetic type.
 *
 * @tparam T type of elements in the input matrix
 * @param A the input matrix_cl (vector or row_vector)
 * @param m number of rows (if x is a row_vector) or columns
 *  (if x is a vector) in the results matrix
 *
 * @return result matrix with replicated rows or columns
 *
 * @throw <code>domain_error</code> if the
 * requested dimensions are negative
 *
 */
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
inline var_value<matrix_cl<double>> rep_matrix(const var_value<T>& A, int m) {
  return make_callback_var(
      rep_matrix(A.val(), m), [A](vari_value<matrix_cl<double>>& res) mutable {
        if (A.adj().size() != 0) {
          matrix_cl<double> A_adj = std::move(A.adj());
          try {
            opencl_kernels::rep_matrix_rev(
                cl::NDRange(A_adj.rows(), A_adj.cols()), A_adj, res.adj(),
                res.adj().rows(), res.adj().cols(), res.adj().view());
          } catch (const cl::Error& e) {
            check_opencl_error("rep_matrix(rev OpenCL)", e);
          }
          A.adj() = std::move(A_adj);
        }
      });
}

}  // namespace math
}  // namespace stan

#endif
#endif
