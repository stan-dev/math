#ifndef STAN_MATH_MIX_HPP
#define STAN_MATH_MIX_HPP

#include <stan/math/manual_forward_decl.hpp>
/**
 *If this is defined then the script to generate the forward decls may have deleted this file
 */
#ifndef STAN_MATH_FORWARD_DECL_HPP
#include <stan/math/forward_decl.hpp>
#endif
#include <stan/math/mix/meta.hpp>
#include <stan/math/mix/fun.hpp>
#include <stan/math/mix/functor.hpp>

#include <stan/math/rev/core.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/fun.hpp>
#include <stan/math/rev/functor.hpp>
#include <stan/math/rev/prob.hpp>
#include <stan/math/rev/constraint.hpp>

#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/fun.hpp>
#include <stan/math/fwd/functor.hpp>
#include <stan/math/fwd/prob.hpp>
#include <stan/math/fwd/constraint.hpp>

#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math/opencl/rev_constraint.hpp>
#endif

#include <stan/math/prim.hpp>

#endif
