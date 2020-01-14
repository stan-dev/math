#ifndef STAN_MATH_PRIM_SCAL_FUNCTOR_REDUCE_SUM_HPP
#define STAN_MATH_PRIM_SCAL_FUNCTOR_REDUCE_SUM_HPP

#include <stan/math/prim/meta/return_type.hpp>
#include <stan/math/prim/mat/fun/typedefs.hpp>

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


// base definition => compile error
template <typename ReduceFunction, typename Enable = void, typename M, typename T, typename... Args>
struct reduce_sum_impl {};

// todo, double check if I need enable if here
template <typename ReduceFunction, typename M, typename T, typename... Args>
struct reduce_sum_impl<ReduceFunction, M, T, double, Args...> {
  using vmapped_t = std::vector<M>;
  using arg1_t = std::vector<Arg1>;
  using arg2_t = std::vector<Arg2>;
  using arg3_t = std::vector<Arg3>;
  using arg4_t = std::vector<Arg4>;
  std::tuple<std::vector<Args>...> arg_;
  const vmapped_t& vmapped_;
  const arg1_t& arg1_;
  const arg2_t& arg2_;
  const arg3_t& arg3_;
  const arg4_t& arg4_;
  T terms_sum_;
    reduce_sum_impl(const vmapped_t& vmapped, const T& init, Args&&... args)
        : vmapped_(vmapped),
          arg_(std::make_tuple(std::forward<Args>(args)...))
          terms_sum_(value_of(init)) {}

    reduce_sum_impl(reduce_sum_impl& other, tbb::split)
        : vmapped_(other.vmapped_),
        arg_(other.arg_),
        terms_sum_(0.0) {}

    void operator()(const tbb::blocked_range<size_t>& r) {
      if (r.empty()) {
        return;
      }

      auto start = vmapped_.begin();
      std::advance(start, r.begin());
      auto end = vmapped_.begin();
      std::advance(end, r.end());

      const vmapped_t sub_slice(start, end);

      terms_sum_ += ReduceFunction()(r.begin(), r.end() - 1, sub_slice, arg1_,
                                     arg2_, arg3_, arg4_);
    }

    void join(const reduce_sum_impl& child) {
      terms_sum_ += child.terms_sum_;
    }

  // Todo: Better name for this than operator(), collect()?
  template <typename... OtherArgs>
  T operator()(const vmapped_t& vmapped, T init, std::size_t grainsize,
               Args&&... args) const {
    const std::size_t num_jobs = vmapped.size();
    reduce_sum_impl<ReduceFunction, M, T, OtherArgs...> worker(vmapped, init, args...);
    tbb::parallel_reduce(
        tbb::blocked_range<std::size_t>(0, num_jobs, grainsize), worker);
    return std::move(worker.terms_sum_);
  }
};

}  // namespace internal

/*
 * Note that the ReduceFunction is only passed in as type to prohibit
 * that any internal state of the functor is causing trouble. Thus,
 * the functor must be default constructible without any arguments.
 */
template <typename ReduceFunction, typename M, typename T, typename... Args>
constexpr T reduce_sum(const std::vector<M>& vmapped, T init,
                       std::size_t grainsize, Args&&... args) {
  typedef T return_base_t;
  // void here but need to figure out enable_if stuff
  // We do this somewhere in the opencl code
  return internal::reduce_sum_impl<ReduceFunction, M, T, void, Args...>()(
    vmapped, init, grainsize, std::forward<Args>(args)...);
}


}  // namespace math
}  // namespace stan

#endif
