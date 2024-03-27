#ifndef STAN_MATH_REV_HPP
#define STAN_MATH_REV_HPP

#include <stan/math/prim/fun/Eigen.hpp>

#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <stan/math/opencl/rev_constraint.hpp>
#endif

#include <stan/math/rev/constraint.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/fun.hpp>
#include <stan/math/rev/functor.hpp>
#include <stan/math/rev/prob.hpp>

#include <stan/math/prim.hpp>

#endif
