#ifndef STAN_MATH_OPENCL_REV_HPP
#define STAN_MATH_OPENCL_REV_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/prim.hpp>
#include <stan/math/opencl/rev/as_column_vector_or_scalar.hpp>
#include <stan/math/opencl/rev/acos.hpp>
#include <stan/math/opencl/rev/acosh.hpp>
#include <stan/math/opencl/rev/add.hpp>
#include <stan/math/opencl/rev/add_diag.hpp>
#include <stan/math/opencl/rev/adjoint_results.hpp>
#include <stan/math/opencl/rev/arena_matrix_cl.hpp>
#include <stan/math/opencl/rev/arena_type.hpp>
#include <stan/math/opencl/rev/asin.hpp>
#include <stan/math/opencl/rev/asinh.hpp>
#include <stan/math/opencl/rev/atan.hpp>
#include <stan/math/opencl/rev/atanh.hpp>
#include <stan/math/opencl/rev/beta.hpp>
#include <stan/math/opencl/rev/block.hpp>
#include <stan/math/opencl/rev/cbrt.hpp>
#include <stan/math/opencl/rev/ceil.hpp>
#include <stan/math/opencl/rev/cholesky_decompose.hpp>
#include <stan/math/opencl/rev/columns_dot_product.hpp>
#include <stan/math/opencl/rev/columns_dot_self.hpp>
#include <stan/math/opencl/rev/copy.hpp>
#include <stan/math/opencl/rev/cos.hpp>
#include <stan/math/opencl/rev/cosh.hpp>
#include <stan/math/opencl/rev/crossprod.hpp>
#include <stan/math/opencl/rev/diag_matrix.hpp>
#include <stan/math/opencl/rev/diag_pre_multiply.hpp>
#include <stan/math/opencl/rev/diag_post_multiply.hpp>
#include <stan/math/opencl/rev/diagonal.hpp>
#include <stan/math/opencl/rev/digamma.hpp>
#include <stan/math/opencl/rev/divide.hpp>
#include <stan/math/opencl/rev/dot_product.hpp>
#include <stan/math/opencl/rev/dot_self.hpp>
#include <stan/math/opencl/rev/elt_divide.hpp>
#include <stan/math/opencl/rev/elt_multiply.hpp>
#include <stan/math/opencl/rev/erf.hpp>
#include <stan/math/opencl/rev/erfc.hpp>
#include <stan/math/opencl/rev/exp.hpp>
#include <stan/math/opencl/rev/exp2.hpp>
#include <stan/math/opencl/rev/expm1.hpp>
#include <stan/math/opencl/rev/fabs.hpp>
#include <stan/math/opencl/rev/floor.hpp>
#include <stan/math/opencl/rev/hypot.hpp>
#include <stan/math/opencl/rev/inv.hpp>
#include <stan/math/opencl/rev/inv_cloglog.hpp>
#include <stan/math/opencl/rev/inv_logit.hpp>
#include <stan/math/opencl/rev/inv_Phi.hpp>
#include <stan/math/opencl/rev/inv_sqrt.hpp>
#include <stan/math/opencl/rev/inv_square.hpp>
#include <stan/math/opencl/rev/lbeta.hpp>
#include <stan/math/opencl/rev/ldexp.hpp>
#include <stan/math/opencl/rev/lgamma.hpp>
#include <stan/math/opencl/rev/log.hpp>
#include <stan/math/opencl/rev/log10.hpp>
#include <stan/math/opencl/rev/log1m.hpp>
#include <stan/math/opencl/rev/log1p.hpp>
#include <stan/math/opencl/rev/log2.hpp>
#include <stan/math/opencl/rev/logit.hpp>
#include <stan/math/opencl/rev/log1m_exp.hpp>
#include <stan/math/opencl/rev/log1p_exp.hpp>
#include <stan/math/opencl/rev/log_inv_logit.hpp>
#include <stan/math/opencl/rev/log1m_inv_logit.hpp>
#include <stan/math/opencl/rev/log_inv_logit_diff.hpp>
#include <stan/math/opencl/rev/log_diff_exp.hpp>
#include <stan/math/opencl/rev/mdivide_left_tri_low.hpp>
#include <stan/math/opencl/rev/mdivide_right_tri_low.hpp>
#include <stan/math/opencl/rev/multiply.hpp>
#include <stan/math/opencl/rev/multiply_log.hpp>
#include <stan/math/opencl/rev/operands_and_partials.hpp>
#include <stan/math/opencl/rev/operator_unary_minus.hpp>
#include <stan/math/opencl/rev/operator_unary_plus.hpp>
#include <stan/math/opencl/rev/Phi.hpp>
#include <stan/math/opencl/rev/Phi_approx.hpp>
#include <stan/math/opencl/rev/pow.hpp>
#include <stan/math/opencl/rev/round.hpp>
#include <stan/math/opencl/rev/rows_dot_product.hpp>
#include <stan/math/opencl/rev/rows_dot_self.hpp>
#include <stan/math/opencl/rev/sin.hpp>
#include <stan/math/opencl/rev/sinh.hpp>
#include <stan/math/opencl/rev/sqrt.hpp>
#include <stan/math/opencl/rev/square.hpp>
#include <stan/math/opencl/rev/squared_distance.hpp>
#include <stan/math/opencl/rev/subtract.hpp>
#include <stan/math/opencl/rev/sum.hpp>
#include <stan/math/opencl/rev/tan.hpp>
#include <stan/math/opencl/rev/tanh.hpp>
#include <stan/math/opencl/rev/tcrossprod.hpp>
#include <stan/math/opencl/rev/tgamma.hpp>
#include <stan/math/opencl/rev/transpose.hpp>
#include <stan/math/opencl/rev/to_arena.hpp>
#include <stan/math/opencl/rev/trunc.hpp>
#include <stan/math/rev.hpp>

#endif
#endif
