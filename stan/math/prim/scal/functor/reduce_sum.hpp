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
template <class ReduceFunction, class M, class T, class Arg1, class Arg2,
          class Arg3, class Arg4, class T_return_type>
struct reduce_sum_impl {};

template <class ReduceFunction, class M, class T, class Arg1, class Arg2,
          class Arg3, class Arg4>
struct reduce_sum_impl<ReduceFunction, M, T, Arg1, Arg2, Arg3, Arg4, double> {
  using vmapped_t = std::vector<M>;
  using arg1_t = std::vector<Arg1>;
  using arg2_t = std::vector<Arg2>;
  using arg3_t = std::vector<Arg3>;
  using arg4_t = std::vector<Arg4>;

  struct recursive_reducer {
    const vmapped_t& vmapped_;
    const arg1_t& arg1_;
    const arg2_t& arg2_;
    const arg3_t& arg3_;
    const arg4_t& arg4_;
    T terms_sum_;

    recursive_reducer(const vmapped_t& vmapped, const T& init,
                      const arg1_t& arg1, const arg2_t& arg2,
                      const arg3_t& arg3, const arg4_t& arg4, )
        : vmapped_(vmapped),
          arg1_(arg1),
          arg2_(arg2),
          arg3_(arg3),
          arg4_(arg4),
          terms_sum_(value_of(init)) {}

    recursive_reducer(recursive_reducer& other, tbb::split)
        : vmapped_(other.vmapped_),
          arg1_(other.arg1_),
          arg2_(other.arg2_),
          arg3_(other.arg3_),
          arg4_(other.arg4_),
          terms_sum_(0.0) {}

    void operator()(const tbb::blocked_range<size_t>& r) {
      if (r.empty())
        return;

      auto start = vmapped_.begin();
      std::advance(start, r.begin());
      auto end = vmapped_.begin();
      std::advance(end, r.end());

      const vmapped_t sub_slice(start, end);

      terms_sum_ += ReduceFunction()(r.begin(), r.end() - 1, sub_slice, arg1_,
                                     arg2_, arg3_, arg4_);
    }

    void join(const recursive_reducer& child) {
      terms_sum_ += child.terms_sum_;
    }
  };

  T operator()(const vmapped_t& vmapped, T init, std::size_t grainsize,
               const arg1_t& arg1, const arg2_t& arg2, const arg3_t& arg3,
               const arg4_t& arg4) const {
    const std::size_t num_jobs = vmapped.size();
    recursive_reducer worker(vmapped, init, arg1, arg2, arg3, arg4);
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
 *
 */
template <class ReduceFunction, class M, class T, class Arg1, class Arg2,
          class Arg3, class Arg4>
constexpr T reduce_sum(const std::vector<M>& vmapped, T init,
                       std::size_t grainsize, const std::vector<Arg1>& arg1,
                       const std::vector<Arg2>& arg2,
                       const std::vector<Arg3>& arg3,
                       const std::vector<Arg4>& arg4) {
  typedef T return_base_t;
  return internal::reduce_sum_impl<ReduceFunction, M, T, Arg1, Arg2, Arg3, Arg4,
                                   return_base_t>()(vmapped, init, grainsize,
                                                    arg1, arg2, arg3, arg4);
}

}  // namespace math
}  // namespace stan

#endif
