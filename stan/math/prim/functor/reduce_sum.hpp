#ifndef STAN_MATH_PRIM_SCAL_FUNCTOR_REDUCE_SUM_HPP
#define STAN_MATH_PRIM_SCAL_FUNCTOR_REDUCE_SUM_HPP

#include <stan/math/prim/meta.hpp>

#include <tbb/task_arena.h>
#include <tbb/parallel_reduce.h>
#include <tbb/blocked_range.h>

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
   *   start through end - 1 terms of the overall sum. All args are passed
   *   from this function through to the `ReduceFunction` instances.
   *   However, only elements start through end - 1 of the vmapped argument are
   *   passed to the `ReduceFunction` instances (as the `vmapped_subset`
   * argument).
   *
   * @param vmapped Sliced arguments used only in some sum terms
   * @param grainsize Suggested grainsize for tbb
   * @param[in, out] msgs The print stream for warning messages
   * @param args Shared arguments used in every sum term
   * @return Summation of all terms
   */
  return_type_t<Vec, Args...> operator()(Vec&& vmapped, std::size_t grainsize,
                                         std::ostream* msgs,
                                         Args&&... args) const {
    const std::size_t num_jobs = vmapped.size();

    if (num_jobs == 0) {
      return 0.0;
    }

    return ReduceFunction()(0, vmapped.size() - 1, std::forward<Vec>(vmapped),
                            msgs, std::forward<Args>(args)...);
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
   *   start through end - 1 terms of the overall sum. All args are passed
   *   from this function through to the `ReduceFunction` instances.
   *   However, only elements start through end - 1 of the vmapped argument are
   *   passed to the `ReduceFunction` instances (as the `vmapped_subset`
   * argument).
   *
   * This function distributes computation of the desired sum
   *   over multiple threads by coordinating calls to `ReduceFunction`
   * instances.
   *
   * @param vmapped Sliced arguments used only in some sum terms
   * @param grainsize Suggested grainsize for tbb
   * @param[in, out] msgs The print stream for warning messages
   * @param args Shared arguments used in every sum term
   * @return Summation of all terms
   */
  ReturnType operator()(Vec&& vmapped, std::size_t grainsize,
                        std::ostream* msgs, Args&&... args) const {
    const std::size_t num_jobs = vmapped.size();
    if (num_jobs == 0) {
      return 0.0;
    }
    recursive_reducer worker(std::forward<Vec>(vmapped), msgs,
                             std::forward<Args>(args)...);

#ifdef STAN_DETERMINISTIC
    tbb::static_partitioner partitioner;
    tbb::parallel_deterministic_reduce(
        tbb::blocked_range<std::size_t>(0, num_jobs, grainsize), worker,
        partitioner);
#else
    tbb::parallel_reduce(
        tbb::blocked_range<std::size_t>(0, num_jobs, grainsize), worker);
#endif

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
 * @tparam ReduceFunction Type of reducer function
 * @tparam ReturnType An arithmetic type
 * @tparam Vec Type of sliced argument
 * @tparam Args Types of shared arguments
 * @return Sum of terms
 */
template <typename ReduceFunction, typename Vec,
          typename = require_vector_like_t<Vec>, typename... Args>
auto reduce_sum(Vec&& vmapped, std::size_t grainsize, std::ostream* msgs,
                Args&&... args) {
  using return_type = return_type_t<Vec, Args...>;
  return internal::reduce_sum_impl<ReduceFunction, void, return_type, Vec,
                                   Args...>()(
      std::forward<Vec>(vmapped), grainsize, msgs, std::forward<Args>(args)...);
}

}  // namespace math
}  // namespace stan

#endif
