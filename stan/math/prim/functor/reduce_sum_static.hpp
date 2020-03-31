#ifndef STAN_MATH_PRIM_FUNCTOR_REDUCE_SUM_STATIC_HPP
#define STAN_MATH_PRIM_FUNCTOR_REDUCE_SUM_STATIC_HPP

#include <stan/math/prim/meta.hpp>

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
 * An instance, f, of `ReduceFunction` should have the signature:
 *   T f(int start, int end, Vec&& vmapped_subset, std::ostream* msgs, Args&&...
 * args)
 *
 * `ReduceFunction` must be default constructible without any arguments
 *
 * grainsize must be greater than or equal to 1
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
          typename = require_vector_like_t<Vec>, typename... Args>
auto reduce_sum_static(Vec&& vmapped, int grainsize, std::ostream* msgs,
                       Args&&... args) {
  using return_type = return_type_t<Vec, Args...>;

  check_positive("reduce_sum", "grainsize", grainsize);

  return internal::reduce_sum_impl<ReduceFunction, void, return_type, Vec,
                                   Args...>()(std::forward<Vec>(vmapped), false,
                                              grainsize, msgs,
                                              std::forward<Args>(args)...);
}

}  // namespace math
}  // namespace stan

#endif
