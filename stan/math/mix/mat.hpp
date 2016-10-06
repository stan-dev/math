#ifndef STAN_MATH_MIX_MAT_HPP
#define STAN_MATH_MIX_MAT_HPP

#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/scal/meta/is_fvar.hpp>
#include <stan/math/fwd/scal/meta/partials_type.hpp>

#include <stan/math/rev/core.hpp>
#include <stan/math/rev/scal/meta/is_var.hpp>
#include <stan/math/rev/scal/meta/partials_type.hpp>

#include <stan/math/prim/mat.hpp>
#include <stan/math/fwd/mat.hpp>
#include <stan/math/rev/mat.hpp>

// FIXME(carpenter): include these in stan.hpp, too, for reverse?
// these are here because they have to come after the general defs
#include <stan/math/prim/mat/fun/acosh.hpp>
#include <stan/math/prim/mat/fun/asinh.hpp>
#include <stan/math/prim/mat/fun/atanh.hpp>
#include <stan/math/prim/mat/fun/exp2.hpp>
#include <stan/math/prim/mat/fun/log1m_inv_logit.hpp>

#include <stan/math/mix/mat/fun/typedefs.hpp>

#include <stan/math/mix/mat/functor/derivative.hpp>
#include <stan/math/mix/mat/functor/finite_diff_grad_hessian.hpp>
#include <stan/math/mix/mat/functor/grad_hessian.hpp>
#include <stan/math/mix/mat/functor/grad_tr_mat_times_hessian.hpp>
#include <stan/math/mix/mat/functor/gradient_dot_vector.hpp>
#include <stan/math/mix/mat/functor/hessian.hpp>
#include <stan/math/mix/mat/functor/hessian_times_vector.hpp>
#include <stan/math/mix/mat/functor/partial_derivative.hpp>

#endif
