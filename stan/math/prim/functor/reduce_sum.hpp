#ifndef STAN_MATH_PRIM_FUNCTOR_REDUCE_SUM_HPP
#define STAN_MATH_PRIM_FUNCTOR_REDUCE_SUM_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>

#include <tbb/task_arena.h>
#include <tbb/parallel_reduce.h>
#include <tbb/blocked_range.h>

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
 * Specialization of reduce_sum_impl for arithmetic types
 *
 * @tparam ReduceFunction Type of reducer function
 * @tparam ReturnType An arithmetic type
 * @tparam Vec Type of sliced argument
 * @tparam Args Types of shared arguments
 */
template <typename ReduceFunction, typename ReturnType, typename Vec,
          typename... Args>
struct reduce_sum_impl<ReduceFunction, require_arithmetic_t<ReturnType>,
                       ReturnType, Vec, Args...> {
  /**
   * This struct is used by the TBB to accumulate partial
   *  sums over consecutive ranges of the input. To distribute the workload,
   *  the TBB can split larger partial sums into smaller ones in which
   *  case the splitting copy constructor is used. It is designed to
   *  meet the Imperative form requirements of `tbb::parallel_reduce`.
   *
   * @note see link [here](https://tinyurl.com/vp7xw2t) for requirements.
   */
  struct recursive_reducer {
    Vec vmapped_;
    std::stringstream msgs_;
    std::tuple<Args...> args_tuple_;
    return_type_t<Vec, Args...> sum_{0.0};

    recursive_reducer(Vec&& vmapped, std::ostream* msgs, Args&&... args)
        : vmapped_(std::forward<Vec>(vmapped)),
          args_tuple_(std::forward<Args>(args)...) {}

    /**
     * This is the copy operator as required for tbb::parallel_reduce
     *   Imperative form. This requires sum_ be reset to zero since
     *   the newly created reducer is used to accumulate an independent
     *   partial sum.
     */
    recursive_reducer(recursive_reducer& other, tbb::split)
        : vmapped_(other.vmapped_), args_tuple_(other.args_tuple_) {}

    /**
     * Compute the value and of `ReduceFunction` over the range defined by r
     *   and accumulate those in member variable sum_. This function may
     *   be called multiple times per object instantiation (so the sum_
     *   must be accumulated, not just assigned).
     *
     * @param r Range over which to compute `ReduceFunction`
     */
    inline void operator()(const tbb::blocked_range<size_t>& r) {
      if (r.empty()) {
        return;
      }

      std::decay_t<Vec> sub_slice;
      sub_slice.reserve(r.size());
      for (size_t i = r.begin(); i < r.end(); ++i) {
        sub_slice.emplace_back(vmapped_[i]);
      }

      sum_ += apply(
          [&](auto&&... args) {
            return ReduceFunction()(sub_slice, r.begin(), r.end() - 1, &msgs_,
                                    args...);
          },
          args_tuple_);
    }

    /**
     * Join reducers. Accumuluate the value (sum_) of the other reducer.
     *
     * @param rhs Another partial sum
     */
    inline void join(const recursive_reducer& rhs) {
      sum_ += rhs.sum_;
      msgs_ << rhs.msgs_.str();
    }
  };

  /**
   * Call an instance of the function `ReduceFunction` on every element
   *   of an input sequence and sum these terms.
   *
   * This specialization is parallelized using tbb and works only for
   *   arithmetic types.
   *
   * ReduceFunction must define an operator() with the same signature as:
   *   double f(Vec&& vmapped_subset, int start, int end, std::ostream* msgs,
   * Args&&... args)
   *
   * `ReduceFunction` must be default constructible without any arguments
   *
   * Each call to `ReduceFunction` is responsible for computing the
   *   start through end (inclusive) terms of the overall sum. All args are
   * passed from this function through to the `ReduceFunction` instances.
   *   However, only the start through end (inclusive) elements of the vmapped
   * argument are passed to the `ReduceFunction` instances (as the
   * `vmapped_subset` argument).
   *
   * This function distributes computation of the desired sum
   *   over multiple threads by coordinating calls to `ReduceFunction`
   * instances.
   *
   * If auto partitioning is true, break work into pieces automatically,
   *  taking grainsize as a recommended work size. The partitioning is
   *  not deterministic nor is the order guaranteed in which partial
   *  sums are accumulated. Due to floating point imprecisions this will likely
   *  lead to slight differences in the accumulated results between
   *  multiple runs. If false, break work deterministically into pieces smaller
   *  than or equal to grainsize and accumulate all the partial sums
   *  in the same order. This still may not achieve bitwise reproducibility.
   *
   * grainsize must be greater than or equal to 1
   *
   * @param vmapped Vector containing one element per term of sum
   * @param auto_partitioning Work partitioning style
   * @param grainsize Suggested grainsize for tbb
   * @param[in, out] msgs The print stream for warning messages
   * @param args Shared arguments used in every sum term
   * @return Summation of all terms
   */
  inline ReturnType operator()(Vec&& vmapped, bool auto_partitioning,
                               int grainsize, std::ostream* msgs,
                               Args&&... args) const {
    const std::size_t num_terms = vmapped.size();
    if (vmapped.empty()) {
      return 0.0;
    }
    recursive_reducer worker(std::forward<Vec>(vmapped), msgs,
                             std::forward<Args>(args)...);

    if (auto_partitioning) {
      tbb::parallel_reduce(
          tbb::blocked_range<std::size_t>(0, num_terms, grainsize), worker);
    } else {
      tbb::simple_partitioner partitioner;
      tbb::parallel_deterministic_reduce(
          tbb::blocked_range<std::size_t>(0, num_terms, grainsize), worker,
          partitioner);
    }
    if (msgs) {
      *msgs << worker.msgs_.str();
    }

    return worker.sum_;
  }
};

}  // namespace internal

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
 * @tparam ReduceFunction Type of reducer function
 * @tparam ReturnType An arithmetic type
 * @tparam Vec Type of sliced argument
 * @tparam Args Types of shared arguments
 * @param vmapped Vector containing one element per term of sum
 * @param grainsize Suggested grainsize for tbb
 * @param[in, out] msgs The print stream for warning messages
 * @param args Shared arguments used in every sum term
 * @return Sum of terms
 */
template <typename ReduceFunction, typename Vec,
          typename = require_vector_like_t<Vec>, typename... Args>
inline auto reduce_sum(Vec&& vmapped, int grainsize, std::ostream* msgs,
                       Args&&... args) {
  using return_type = return_type_t<Vec, Args...>;

  check_positive("reduce_sum", "grainsize", grainsize);

#ifdef STAN_THREADS
  return internal::reduce_sum_impl<ReduceFunction, void, return_type, Vec,
                                   ref_type_t<Args&&>...>()(
      std::forward<Vec>(vmapped), true, grainsize, msgs,
      std::forward<Args>(args)...);
#else
  if (vmapped.empty()) {
    return return_type(0.0);
  }

  return ReduceFunction()(std::forward<Vec>(vmapped), 0, vmapped.size() - 1,
                          msgs, std::forward<Args>(args)...);
#endif
}

}  // namespace math
}  // namespace stan

#endif
