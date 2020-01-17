#ifndef STAN_MATH_PRIM_SCAL_FUNCTOR_REDUCE_SUM_HPP
#define STAN_MATH_PRIM_SCAL_FUNCTOR_REDUCE_SUM_HPP

#include <stan/math/prim/meta/return_type.hpp>
#include <stan/math/prim/fun/typedefs.hpp>

#include <boost/iterator/counting_iterator.hpp>

#include <tbb/task_arena.h>
#include <tbb/parallel_reduce.h>
#include <tbb/blocked_range.h>

#include <iostream>
#include <iterator>
#include <vector>

namespace stan {
namespace math {

namespace internal {

template <typename ReduceFunction, typename Enable, typename ReturnType,
          typename Vec, typename... Args>
struct reduce_sum_impl {
  return_type_t<Vec, Args...> operator()(const Vec& vmapped,
                                         std::size_t grainsize,
                                         const Args&... args) const {
    const std::size_t num_jobs = vmapped.size();

    if (num_jobs == 0)
      return 0.0;

    return ReduceFunction()(0, vmapped.size() - 1, vmapped, args...);
  }
};

template <typename ReduceFunction, typename ReturnType, typename Vec,
          typename... Args>
struct reduce_sum_impl<ReduceFunction, require_arithmetic_t<ReturnType>,
                       ReturnType, Vec, Args...> {
  struct recursive_reducer {
    using vmapped_t = Vec;
    std::tuple<const Args&...> args_tuple_;
    size_t tuple_size_ = sizeof...(Args);
    const vmapped_t& vmapped_;
    return_type_t<Vec, Args...> sum_;

    recursive_reducer(const vmapped_t& vmapped, const Args&... args)
        : vmapped_(vmapped), args_tuple_(args...), sum_(0.0) {}

    recursive_reducer(recursive_reducer& other, tbb::split)
        : vmapped_(other.vmapped_), args_tuple_(other.args_tuple_), sum_(0.0) {}

    void operator()(const tbb::blocked_range<size_t>& r) {
      if (r.empty()) {
        return;
      }

      auto start = vmapped_.begin();
      std::advance(start, r.begin());
      auto end = vmapped_.begin();
      std::advance(end, r.end());

      const vmapped_t sub_slice(start, end);

      sum_ += apply(
          [&](auto&&... args) {
            return ReduceFunction()(r.begin(), r.end() - 1, sub_slice, args...);
          },
          args_tuple_);
    }

    void join(const recursive_reducer& child) { sum_ += child.sum_; }
  };

  return_type_t<Vec, Args...> operator()(const Vec& vmapped,
                                         std::size_t grainsize,
                                         const Args&... args) const {
    const std::size_t num_jobs = vmapped.size();

    if (num_jobs == 0)
      return 0.0;

    recursive_reducer worker(vmapped, args...);

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
 * Note that the ReduceFunction is only passed in as type to prohibit
 * that any internal state of the functor is causing trouble. Thus,
 * the functor must be default constructible without any arguments.
 */
template <typename ReduceFunction, typename Vec, require_vector_like_t<Vec>...,
          typename... Args>
auto reduce_sum(const Vec& vmapped, std::size_t grainsize,
                          const Args&... args) {
  using return_type = return_type_t<Vec, Args...>;
  return internal::reduce_sum_impl<ReduceFunction, void, return_type, Vec,
                                   Args...>()(vmapped, grainsize, args...);
}

}  // namespace math
}  // namespace stan

#endif
