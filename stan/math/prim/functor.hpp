#ifndef STAN_MATH_PRIM_FUNCTOR_HPP
#define STAN_MATH_PRIM_FUNCTOR_HPP

#include <stan/math/prim/functor/apply.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>
#include <stan/math/prim/functor/apply_scalar_binary.hpp>
#include <stan/math/prim/functor/apply_vector_unary.hpp>
#include <stan/math/prim/functor/coupled_ode_system.hpp>
#include <stan/math/prim/functor/finite_diff_gradient.hpp>
#include <stan/math/prim/functor/finite_diff_gradient_auto.hpp>
#include <stan/math/prim/functor/for_each.hpp>
#include <stan/math/prim/functor/hcubature.hpp>
#include <stan/math/prim/functor/integrate_1d.hpp>
#include <stan/math/prim/functor/integrate_1d_adapter.hpp>
#include <stan/math/prim/functor/integrate_ode_rk45.hpp>
#include <stan/math/prim/functor/integrate_ode_std_vector_interface_adapter.hpp>
#include <stan/math/prim/functor/ode_ckrk.hpp>
#include <stan/math/prim/functor/ode_rk45.hpp>
#include <stan/math/prim/functor/ode_store_sensitivities.hpp>
#include <stan/math/prim/functor/map_rect.hpp>
#include <stan/math/prim/functor/map_rect_combine.hpp>
#include <stan/math/prim/functor/map_rect_concurrent.hpp>
#include <stan/math/prim/functor/map_rect_reduce.hpp>
#include <stan/math/prim/functor/mpi_cluster.hpp>
#include <stan/math/prim/functor/mpi_command.hpp>
#include <stan/math/prim/functor/mpi_distributed_apply.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <stan/math/prim/functor/reduce_sum.hpp>
#include <stan/math/prim/functor/reduce_sum_static.hpp>

#endif
