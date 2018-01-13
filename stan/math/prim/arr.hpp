#ifndef STAN_MATH_PRIM_ARR_HPP
#define STAN_MATH_PRIM_ARR_HPP

#include <stan/math/prim/arr/meta/get.hpp>
#include <stan/math/prim/arr/meta/index_type.hpp>
#include <stan/math/prim/arr/meta/is_constant_struct.hpp>
#include <stan/math/prim/arr/meta/is_vector.hpp>
#include <stan/math/prim/arr/meta/length.hpp>
#include <stan/math/prim/arr/meta/scalar_type.hpp>
#include <stan/math/prim/arr/meta/value_type.hpp>
#include <stan/math/prim/arr/meta/VectorBuilderHelper.hpp>

#include <stan/math/prim/arr/err/check_matching_sizes.hpp>
#include <stan/math/prim/arr/err/check_nonzero_size.hpp>
#include <stan/math/prim/arr/err/check_ordered.hpp>

#include <stan/math/prim/arr/fun/array_builder.hpp>
#include <stan/math/prim/arr/fun/common_type.hpp>
#include <stan/math/prim/arr/fun/dot.hpp>
#include <stan/math/prim/arr/fun/dot_self.hpp>
#include <stan/math/prim/arr/fun/fill.hpp>
#include <stan/math/prim/arr/fun/inverse_softmax.hpp>
#include <stan/math/prim/arr/fun/log_sum_exp.hpp>
#include <stan/math/prim/arr/fun/promote_elements.hpp>
#include <stan/math/prim/arr/fun/promote_scalar.hpp>
#include <stan/math/prim/arr/fun/promote_scalar_type.hpp>
#include <stan/math/prim/arr/fun/rep_array.hpp>
#include <stan/math/prim/arr/fun/scaled_add.hpp>
#include <stan/math/prim/arr/fun/sort_asc.hpp>
#include <stan/math/prim/arr/fun/sort_desc.hpp>
#include <stan/math/prim/arr/fun/sub.hpp>
#include <stan/math/prim/arr/fun/sum.hpp>
#include <stan/math/prim/arr/fun/value_of.hpp>
#include <stan/math/prim/arr/fun/value_of_rec.hpp>

#include <stan/math/prim/arr/functor/coupled_ode_observer.hpp>
#include <stan/math/prim/arr/functor/coupled_ode_system.hpp>
#include <stan/math/prim/arr/functor/integrate_ode_rk45.hpp>

#include <stan/math/prim/scal.hpp>

#endif
