#ifndef STAN_MATH_PRIM_FUNCTOR_REDUCE_SUM_STATIC_HPP
#define STAN_MATH_PRIM_FUNCTOR_REDUCE_SUM_STATIC_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/functor/reduce_sum.hpp>
#include <tbb/task_arena.h>
#include <tbb/parallel_reduce.h>
#include <tbb/blocked_range.h>

#include <tuple>
#include <vector>

namespace stan {
namespace math {

/**
 * Call an instance of the function `ReduceFunction` on every element
 *   of an input sequence and sum these terms.
 *
 * This defers to reduce_sum_impl for the appropriate implementation
 *
 * ReduceFunction must define an operator() with the same signature as:
 *   T f(Vec&& vmapped_subset, int start, int end, std::ostream* msgs, Args&&...
 * args)
 *
 * `ReduceFunction` must be default constructible without any arguments
 *
 * grainsize must be greater than or equal to 1
 *
 * If STAN_THREADS is not defined, do all the work with one ReduceFunction call.
 *
 * @tparam ReduceFunction Type of reducer function
 * @tparam ReturnType An arithmetic type
 * @tparam Vec Type of sliced argument
 * @tparam Args Types of shared arguments
 * @param vmapped Sliced arguments used only in some sum terms
 * @param grainsize Suggested grainsize for tbb
 * @param[in, out] msgs The print stream for warning messages
 * @param args Shared arguments used in every sum term
 * @return Sum of terms
 */
template <typename ReduceFunction, typename Vec,
          typename = require_vector_like_t<Vec>,
	  require_stan_closure_t<ReduceFunction>* = nullptr,
	  typename... Args>
auto reduce_sum_static(Vec&& vmapped, int grainsize, std::ostream* msgs,
                       const ReduceFunction& f, Args&&... args) {
  using return_type = return_type_t<ReduceFunction, Vec, Args...>;

  check_positive("reduce_sum", "grainsize", grainsize);

#ifdef STAN_THREADS
  return internal::reduce_sum_impl<ReduceFunction, void, return_type, Vec,
                                   Args...>()(std::forward<Vec>(vmapped), false,
                                              grainsize, msgs, f,
                                              std::forward<Args>(args)...);
#else
  if (vmapped.empty()) {
    return return_type(0);
  }

  return f(std::forward<Vec>(vmapped), 0, vmapped.size() - 1,
	   msgs, std::forward<Args>(args)...);
#endif
}

template <typename ReduceFunction, typename Vec,
          typename = require_vector_like_t<Vec>,
	  require_not_stan_closure_t<ReduceFunction>* = nullptr,
	  typename... Args>
auto reduce_sum_static(Vec&& vmapped, int grainsize, std::ostream* msgs,
                       Args&&... args) {
  ReduceFunction f;
  internal::reduce_sum_closure_adapter<ReduceFunction> cl(f);
  return reduce_sum_static(vmapped, grainsize, msgs, cl, args...);
}

}  // namespace math
}  // namespace stan

#endif
