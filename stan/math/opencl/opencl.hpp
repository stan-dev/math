#ifndef STAN_MATH_OPENCL_OPENCL
#define STAN_MATH_OPENCL_OPENCL
#ifdef STAN_OPENCL

#include <stan/math/opencl/opencl_context.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/copy.hpp>
#include <stan/math/opencl/copy_triangular.hpp>
#include <stan/math/opencl/cholesky_decompose.hpp>
#include <stan/math/opencl/diagonal_multiply.hpp>
#include <stan/math/opencl/identity.hpp>
#include <stan/math/opencl/is_matrix_cl.hpp>
#include <stan/math/opencl/tri_inverse.hpp>
#include <stan/math/opencl/multiply_transpose.hpp>
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/opencl/prim/rep_matrix.hpp>
#include <stan/math/opencl/prim/rep_row_vector.hpp>
#include <stan/math/opencl/prim/rep_vector.hpp>
#include <stan/math/opencl/scalar_type.hpp>
#include <stan/math/opencl/sub_block.hpp>
#include <stan/math/opencl/triangular_transpose.hpp>
#include <stan/math/opencl/value_type.hpp>
#include <stan/math/opencl/zeros.hpp>

#include <stan/math/opencl/prim/add.hpp>
#include <stan/math/opencl/prim/bernoulli_logit_glm_lpmf.hpp>
#include <stan/math/opencl/prim/categorical_logit_glm_lpmf.hpp>
#include <stan/math/opencl/prim/cholesky_decompose.hpp>
#include <stan/math/opencl/prim/divide_columns.hpp>
#include <stan/math/opencl/prim/gp_exp_quad_cov.hpp>
#include <stan/math/opencl/prim/mdivide_left_tri_low.hpp>
#include <stan/math/opencl/prim/mdivide_right_tri_low.hpp>
#include <stan/math/opencl/prim/multiply.hpp>
#include <stan/math/opencl/prim/neg_binomial_2_log_glm_lpmf.hpp>
#include <stan/math/opencl/prim/normal_id_glm_lpdf.hpp>
#include <stan/math/opencl/prim/ordered_logistic_glm_lpmf.hpp>
#include <stan/math/opencl/prim/poisson_log_glm_lpmf.hpp>
#include <stan/math/opencl/prim/subtract.hpp>
#include <stan/math/opencl/prim/transpose.hpp>

#include <stan/math/opencl/err/check_diagonal_zeros.hpp>
#include <stan/math/opencl/err/check_invalid_matrix_view.hpp>
#include <stan/math/opencl/err/check_matching_dims.hpp>
#include <stan/math/opencl/err/check_nan.hpp>
#include <stan/math/opencl/err/check_opencl.hpp>
#include <stan/math/opencl/err/check_mat_size_one.hpp>
#include <stan/math/opencl/err/check_square.hpp>
#include <stan/math/opencl/err/check_symmetric.hpp>
#include <stan/math/opencl/err/check_vector.hpp>

#endif
#endif
