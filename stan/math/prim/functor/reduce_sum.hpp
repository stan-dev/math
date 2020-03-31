#ifndef STAN_MATH_PRIM_FUNCTOR_REDUCE_SUM_HPP
#define STAN_MATH_PRIM_FUNCTOR_REDUCE_SUM_HPP

#include <stan/math/prim/meta.hpp>

#include <tbb/task_arena.h>
#include <tbb/parallel_reduce.h>
#include <tbb/blocked_range.h>

#include <algorithm>
#include <tuple>
#include <vector>

namespace stan {
namespace math {

namespace internal {

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
   *   T f(int start, int end, Vec&& vmapped_subset, std::ostream* msgs,
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
   * @param vmapped Sliced arguments used only in some sum terms
   * @param auto_partitioning Work partitioning style (ignored)
   * @param grainsize Suggested grainsize for tbb
   * @param[in, out] msgs The print stream for warning messages
   * @param args Shared arguments used in every sum term
   * @return Summation of all terms
   */
  return_type_t<Vec, Args...> operator()(Vec&& vmapped, bool auto_partitioning,
                                         int grainsize, std::ostream* msgs,
                                         Args&&... args) const {
    const std::size_t num_jobs = vmapped.size();

    if (num_jobs == 0) {
      return 0.0;
    }

    if (auto_partitioning) {
      return ReduceFunction()(0, vmapped.size() - 1, std::forward<Vec>(vmapped),
                              msgs, std::forward<Args>(args)...);
    } else {
      return_type_t<Vec, Args...> sum = 0.0;
      for (size_t i = 0; i < (vmapped.size() + grainsize - 1) / grainsize;
           ++i) {
        size_t start = i * grainsize;
        size_t end = std::min((i + 1) * grainsize, vmapped.size()) - 1;

        std::decay_t<Vec> sub_slice;
        sub_slice.reserve(end - start + 1);
        for (int i = start; i <= end; ++i) {
          sub_slice.emplace_back(vmapped[i]);
        }

        sum += ReduceFunction()(start, end, std::forward<Vec>(sub_slice), msgs,
                                std::forward<Args>(args)...);
      }
      return sum;
    }
  }
};

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
   * Internal object meeting the Imperative form requirements of
   * `tbb::parallel_reduce`
   *
   * @note see link [here](https://tinyurl.com/vp7xw2t) for requirements.
   */
  struct recursive_reducer {
    Vec vmapped_;
    std::ostream* msgs_;
    std::tuple<Args...> args_tuple_;
    return_type_t<Vec, Args...> sum_{0.0};

    recursive_reducer(Vec&& vmapped, std::ostream* msgs, Args&&... args)
        : vmapped_(std::forward<Vec>(vmapped)),
          msgs_(msgs),
          args_tuple_(std::forward<Args>(args)...) {}

    /*
     * This is the copy operator as required for tbb::parallel_reduce
     *   Imperative form. This requires sum_ be reset to zero.
     */
    recursive_reducer(recursive_reducer& other, tbb::split)
        : vmapped_(other.vmapped_),
          msgs_(other.msgs_),
          args_tuple_(other.args_tuple_) {}

    /**
     * Compute the value and of `ReduceFunction` over the range defined by r
     *   and accumulate those in member variable sum_. This function may
     *   be called multiple times per object instantiation (so the sum_
     *   must be accumulated, not just assigned).
     *
     * @param r Range over which to compute `ReduceFunction`
     */
    void operator()(const tbb::blocked_range<size_t>& r) {
      if (r.empty()) {
        return;
      }

      std::decay_t<Vec> sub_slice;
      sub_slice.reserve(r.size());
      for (int i = r.begin(); i < r.end(); ++i) {
        sub_slice.emplace_back(vmapped_[i]);
      }

      sum_ += apply(
          [&](auto&&... args) {
            return ReduceFunction()(r.begin(), r.end() - 1, sub_slice, msgs_,
                                    args...);
          },
          this->args_tuple_);
    }

    /**
     * Join reducers. Accumuluate the value (sum_) of the other reducer.
     *
     * @param rhs Another partial sum
     */
    void join(const recursive_reducer& child) { this->sum_ += child.sum_; }
  };

  /**
   * Call an instance of the function `ReduceFunction` on every element
   *   of an input sequence and sum these terms.
   *
   * This specialization is parallelized using tbb and works only for
   *   arithmetic types.
   *
   * An instance, f, of `ReduceFunction` should have the signature:
   *   double f(int start, int end, Vec&& vmapped_subset, std::ostream* msgs,
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
   *  taking grainsize as a recommended work size (this process
   *  is not deterministic). If false, break work deterministically
   *  into pieces smaller than or equal to grainsize. The execution
   *  order is non-deterministic.
   *
   * grainsize must be greater than or equal to 1
   *
   * @param vmapped Sliced arguments used only in some sum terms
   * @param auto_partitioning Work partitioning style
   * @param grainsize Suggested grainsize for tbb
   * @param[in, out] msgs The print stream for warning messages
   * @param args Shared arguments used in every sum term
   * @return Summation of all terms
   */
  ReturnType operator()(Vec&& vmapped, bool auto_partitioning, int grainsize,
                        std::ostream* msgs, Args&&... args) const {
    const std::size_t num_jobs = vmapped.size();
    if (num_jobs == 0) {
      return 0.0;
    }
    recursive_reducer worker(std::forward<Vec>(vmapped), msgs,
                             std::forward<Args>(args)...);

    if (auto_partitioning) {
      tbb::parallel_reduce(
          tbb::blocked_range<std::size_t>(0, num_jobs, grainsize), worker);
    } else {
      tbb::simple_partitioner partitioner;
      tbb::parallel_deterministic_reduce(
          tbb::blocked_range<std::size_t>(0, num_jobs, grainsize), worker,
          partitioner);
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
auto reduce_sum(Vec&& vmapped, int grainsize, std::ostream* msgs,
                Args&&... args) {
  using return_type = return_type_t<Vec, Args...>;

  check_positive("reduce_sum", "grainsize", grainsize);

  return internal::reduce_sum_impl<ReduceFunction, void, return_type, Vec,
                                   Args...>()(std::forward<Vec>(vmapped), true,
                                              grainsize, msgs,
                                              std::forward<Args>(args)...);
}

}  // namespace math
}  // namespace stan

#endif
