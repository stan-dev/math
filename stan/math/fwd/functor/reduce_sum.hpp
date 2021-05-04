#ifndef STAN_MATH_FWD_FUNCTOR_REDUCE_SUM_HPP
#define STAN_MATH_FWD_FUNCTOR_REDUCE_SUM_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>

#include <algorithm>
#include <tuple>
#include <vector>

namespace stan {
namespace math {

namespace internal {
template <typename ReduceFunction, typename Enable, typename ReturnType,
          typename Vec, typename... Args>
struct reduce_sum_impl;

/**
 * reduce_sum_impl implementation for any autodiff type.
 *
 * @tparam ReduceFunction Type of reducer function
 * @tparam ReturnType An arithmetic type
 * @tparam Vec Type of sliced argument
 * @tparam Args Types of shared arguments
 */
template <typename ReduceFunction, typename Enable, typename ReturnType,
          typename Vec, typename... Args>
struct reduce_sum_impl {
  /**
   * Call an instance of the function `ReduceFunction` on every element
   *   of an input sequence and sum these terms.
   *
   * This specialization is not parallelized and works for any autodiff types.
   *
   * An instance, f, of `ReduceFunction` should have the signature:
   *   T f(Vec&& vmapped_subset, int start, int end, std::ostream* msgs,
   * Args&&... args)
   *
   * `ReduceFunction` must be default constructible without any arguments
   *
   * Each call to `ReduceFunction` is responsible for computing the
   *   start through end terms (inclusive) of the overall sum. All args are
   * passed from this function through to the `ReduceFunction` instances.
   *   However, only the start through end  (inclusive) elements of the vmapped
   * argument are passed to the `ReduceFunction` instances (as the
   * `vmapped_subset` argument).
   *
   * If auto partitioning is true, do the calculation with one
   *  ReduceFunction call. If false, break work into pieces strictly smaller
   *  than grainsize.
   *
   * grainsize must be greater than or equal to 1
   *
   * @param vmapped Vector containing one element per term of sum
   * @param auto_partitioning Work partitioning style (ignored)
   * @param grainsize Suggested grainsize for tbb
   * @param[in, out] msgs The print stream for warning messages
   * @param args Shared arguments used in every sum term
   * @return Summation of all terms
   */
  inline return_type_t<Vec, Args...> operator()(Vec&& vmapped,
                                                bool auto_partitioning,
                                                int grainsize,
                                                std::ostream* msgs,
                                                Args&&... args) const {
    if (vmapped.empty()) {
      return 0.0;
    }

    if (auto_partitioning) {
      return ReduceFunction()(std::forward<Vec>(vmapped), 0, vmapped.size() - 1,
                              msgs, std::forward<Args>(args)...);
    } else {
      return_type_t<Vec, Args...> sum = 0.0;
      for (size_t i = 0; i < (vmapped.size() + grainsize - 1) / grainsize;
           ++i) {
        size_t start = i * grainsize;
        size_t end = std::min((i + 1) * grainsize, vmapped.size()) - 1;

        std::decay_t<Vec> sub_slice;
        sub_slice.reserve(end - start + 1);
        for (size_t i = start; i <= end; ++i) {
          sub_slice.emplace_back(vmapped[i]);
        }

        sum += ReduceFunction()(std::forward<Vec>(sub_slice), start, end, msgs,
                                std::forward<Args>(args)...);
      }
      return sum;
    }
  }
};

}  // namespace internal

}  // namespace math
}  // namespace stan

#endif
