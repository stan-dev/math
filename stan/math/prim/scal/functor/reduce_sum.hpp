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
          class T_return_type>
struct reduce_sum_impl {};

template <class ReduceFunction, class M, class T, class Arg1, class Arg2>
struct reduce_sum_impl<ReduceFunction, M, T, Arg1, Arg2, double> {
  using vmapped_t = std::vector<M>;
  using arg1_t = std::vector<Arg1>;
  using arg2_t = std::vector<Arg2>;

  struct recursive_reducer {
    const vmapped_t& vmapped_;
    const arg1_t& arg1_;
    const arg2_t& arg2_;
    T terms_sum_;

    recursive_reducer(const vmapped_t& vmapped, const T& init,
                      const arg1_t& arg1, const arg2_t& arg2)
        : vmapped_(vmapped),
          arg1_(arg1),
          arg2_(arg2),
          terms_sum_(value_of(init)) {}

    recursive_reducer(recursive_reducer& other, tbb::split)
        : vmapped_(other.vmapped_),
          arg1_(other.arg1_),
          arg2_(other.arg2_),
          terms_sum_(0.0) {}

    void operator()(const tbb::blocked_range<size_t>& r) {
      if (r.empty())
        return;

      auto start = vmapped_.begin();
      std::advance(start, r.begin());
      auto end = vmapped_.begin();
      std::advance(end, r.end());

      vmapped_t sub_slice(start, end);

      terms_sum_
          += ReduceFunction()(r.begin(), r.end() - 1, sub_slice, arg1_, arg2_);
    }

    void join(const recursive_reducer& child) {
      terms_sum_ += child.terms_sum_;
    }
  };

  T operator()(const vmapped_t& vmapped, T init, std::size_t grainsize,
               const arg1_t& arg1, const arg2_t& arg2) const {
    const std::size_t num_jobs = vmapped.size();
    recursive_reducer worker(vmapped, init, arg1, arg2);
    tbb::parallel_reduce(
        tbb::blocked_range<std::size_t>(0, num_jobs, grainsize), worker);
    return std::move(worker.terms_sum_);
  }
};
}  // namespace internal

/*
 * The ReduceFunction is expected to have an operator with the
 * signature
 *
 * inline T operator()(std::size_t start, std::size_t end,
 *                     std::vector<int>& sub_slice, std::vector<T>& lambda,
 *                     const std::vector<int>& gidx) const
 *
 * So the input data is sliced into smaller pieces and given into the
 * reducer. Which exact sub-slice is given to the function is defined
 * by start and end index. These can be used to map data to groups,
 * for example. It would be good to extend this signature in two ways:
 *
 * 1. The large std::vector which is sliced into smaller pieces should
 *    be allowed to be any data (so int, double, vectors, matrices,
 *    ...). In a second step this should also be extended to vars
 *    which will require more considerations, but is beneficial. One
 *    could leave this generalization to a future version.
 *
 * 2. The arguments past sub_slice should be variable. So the
 *    arguments past sub_slice would need to be represented by the
 *    "..." thing in C++. The complication for this is managing all
 *    the calls.
 *
 * So the reduce_sum signature should look like
 *
 * template <class ReduceFunction, class InputIt, class T, ...>
 * constexpr T reduce_sum(InputIt first, InputIt last, T init,
 *                                std::size_t grainsize, ...)
 *
 * which corresponds to a ReduceFunction which looks like
 *
 * inline T operator()(std::size_t start, std::size_t end,
 *                     std::vector<InputIt*>& sub_slice, ...)
 *
 * Note that the ReduceFunction is only passed in as type to prohibit
 * that any internal state of the functor is causing trouble. Thus,
 * the functor must be default constructible without any arguments.
 *
 */
template <class ReduceFunction, class M, class T, class Arg1, class Arg2>
constexpr T reduce_sum(const std::vector<M>& vmapped, T init,
                                std::size_t grainsize,
                                const std::vector<Arg1>& arg1,
                                const std::vector<Arg2>& arg2) {
  typedef T return_base_t;
  return internal::reduce_sum_impl<ReduceFunction, M, T, Arg1, Arg2,
                                            return_base_t>()(
      vmapped, init, grainsize, arg1, arg2);
}

}  // namespace math
}  // namespace stan

#endif
