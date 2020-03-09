#ifndef STAN_MATH_PRIM_SCAL_FUNCTOR_REDUCE_SUM_HPP
#define STAN_MATH_PRIM_SCAL_FUNCTOR_REDUCE_SUM_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun.hpp>

#include <boost/iterator/counting_iterator.hpp>

#include <tbb/task_arena.h>
#include <tbb/parallel_reduce.h>
#include <tbb/blocked_range.h>

#include <iostream>
#include <iterator>
#include <tuple>
#include <vector>

namespace stan {
namespace math {

namespace internal {

/**
 * Base definition of implimentation called by reduce `reduce_sum`.
 * @tparam ReduceFunction An type with a valid `operator()`
 * @tparam ReturnType The return type of the call to
 * `ReduceFunction::operator()`
 * @tparam Vec type of the first input argument with an `operator[]`
 * @tparam Args Parameter pack holding the types of the object send to
 * `ReduceFunction::operator()`
 */
template <typename ReduceFunction, typename Enable, typename ReturnType,
          typename Vec, typename... Args>
struct reduce_sum_impl {
  template <typename OpVec, typename... OpArgs,
  typename = require_same_t<Vec, OpVec>,
  typename = require_all_t<std::is_same<Args, OpArgs>...>>
  return_type_t<OpVec, OpArgs...> operator()(OpVec&& vmapped, std::size_t grainsize, std::ostream* msgs,
    OpArgs&&... args) const {
    const std::size_t num_jobs = vmapped.size();

    if (num_jobs == 0) {
      return 0.0;      
    }

    return ReduceFunction()(0, vmapped.size() - 1, std::forward<OpVec>(vmapped),
      msgs, std::forward<OpArgs>(args)...);
  }
};

/**
 * Arithmetic definition of implimentation called by reduce `reduce_sum`.
 * @tparam ReduceFunction An type with a valid `operator()`
 * @tparam ReturnType The return type of the call to
 * `ReduceFunction::operator()`
 * @tparam Vec type of the first input argument with an `operator[]`
 * @tparam Args Parameter pack holding the types of the object send to
 * `ReduceFunction::operator()`
 */
template <typename ReduceFunction, typename ReturnType, typename Vec,
          typename... Args>
struct reduce_sum_impl<ReduceFunction, require_arithmetic_t<ReturnType>,
                       ReturnType, Vec, Args...> {
  /**
   * Internal object used in `tbb::parallel_reduce`
   * @note see link [here](https://tinyurl.com/vp7xw2t) for requirements.
   */
  struct recursive_reducer {
    Vec vmapped_;
    std::ostream* msgs_;
    std::tuple<Args...> args_tuple_;
    return_type_t<Vec, Args...> sum_{0.0};

    template <typename OpVec, typename... OpArgs,
    typename = require_same_t<Vec, OpVec>,
    typename = require_all_t<std::is_same<Args, OpArgs>...>>
    recursive_reducer(OpVec&& vmapped, std::ostream* msgs, OpArgs&&... args)
        : vmapped_(std::forward<OpVec>(vmapped)), msgs_(msgs),
          args_tuple_(std::forward<OpArgs>(args)...) {}

    recursive_reducer(recursive_reducer& other, tbb::split)
        : vmapped_(other.vmapped_),
          msgs_(other.msgs_),
          args_tuple_(other.args_tuple_) {}

    void operator()(const tbb::blocked_range<size_t>& r) {
      if (r.empty()) {
        return;
      }

      auto start = this->vmapped_.begin();
      std::advance(start, r.begin());
      auto end = this->vmapped_.begin();
      std::advance(end, r.end());

      const std::decay_t<Vec> sub_slice(start, end);

      sum_ += apply(
          [&](auto&&... args) {
            return ReduceFunction()(r.begin(), r.end() - 1, sub_slice, msgs_,
                                    args...);
          },
          this->args_tuple_);
    }
    /**
     * Join the sum of the reducers
     * @param child A sub-reducer that is aggregated to make the final reduced
     * value.
     * @note side effect of updating `recursive_reducer->sum_`
     */
    void join(const recursive_reducer& child) { this->sum_ += child.sum_; }
  };

  template <typename OpVec, typename... OpArgs,
  typename = require_same_t<Vec, OpVec>,
  typename = require_all_t<std::is_same<Args, OpArgs>...>>
  return_type_t<OpVec, OpArgs...> operator()(OpVec&& vmapped,
                                         std::size_t grainsize,
                                         std::ostream* msgs,
                                         OpArgs&&... args) const {
    const std::size_t num_jobs = vmapped.size();
    if (num_jobs == 0) {
      return 0.0;
    }
    recursive_reducer worker(std::forward<OpVec>(vmapped), msgs,
     std::forward<OpArgs>(args)...);

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

/*
 * Performs a map and reduce sum over a function.
 * @tparam ReduceFunction An type with a valid `operator()`
 * @tparam Vec type of the first input argument with an `operator[]`
 * @tparam Args Parameter pack holding the types of the object send to
 * @note The ReduceFunction is only passed in as type to prohibit
 * that any internal state of the functor is causing trouble. Thus,
 * the functor must be default constructible without any arguments.
 */
template <typename ReduceFunction, typename Vec, typename = require_vector_like_t<Vec>,
          typename... Args>
auto reduce_sum(Vec&& vmapped, std::size_t grainsize, std::ostream* msgs,
                Args&&... args) {
  using return_type = return_type_t<Vec, Args...>;
  return internal::reduce_sum_impl<ReduceFunction, void, return_type, Vec,
                                   Args...>()(std::forward<Vec>(vmapped),
                                              grainsize, msgs,
                                              std::forward<Args>(args)...);
}

}  // namespace math
}  // namespace stan

#endif
