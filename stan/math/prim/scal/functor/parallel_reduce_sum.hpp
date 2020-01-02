#ifndef STAN_MATH_PRIM_SCAL_FUNCTOR_PARALLEL_REDUCE_SUM_HPP
#define STAN_MATH_PRIM_SCAL_FUNCTOR_PARALLEL_REDUCE_SUM_HPP

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
template <class ReduceFunction, class InputIt, class T, class Arg1,
          class T_return_type>
struct parallel_reduce_sum_impl {};

template <class ReduceFunction, class InputIt, class T, class Arg1>
struct parallel_reduce_sum_impl<ReduceFunction, InputIt, T, Arg1, double> {
  struct recursive_reducer {
    InputIt first_;
    const std::vector<Arg1>& varg1_;
    const std::vector<int>& idata1_;
    T sum_;
    typedef typename std::iterator_traits<InputIt>::value_type elem_t;

    recursive_reducer(InputIt first, const T& init,
                      const std::vector<Arg1>& varg1,
                      const std::vector<int>& idata1)
        : first_(first), varg1_(varg1), idata1_(idata1), sum_(init) {}

    recursive_reducer(recursive_reducer& other, tbb::split)
        : first_(other.first_),
          varg1_(other.varg1_),
          idata1_(other.idata1_),
          sum_(T(0.0)) {}

    void operator()(const tbb::blocked_range<size_t>& r) {
      if (r.empty())
        return;

      auto start = first_;
      std::advance(start, r.begin());
      auto end = first_;
      std::advance(end, r.end());

      std::vector<elem_t> sub_slice(start, end);

      sum_ += ReduceFunction()(r.begin(), r.end() - 1, sub_slice, varg1_,
                               idata1_);
    }

    void join(const recursive_reducer& child) { sum_ += child.sum_; }
  };

  T operator()(InputIt first, InputIt last, T init, std::size_t grainsize,
               const std::vector<Arg1>& varg1,
               const std::vector<int>& idata1) const {
    const std::size_t num_jobs = std::distance(first, last);
    recursive_reducer worker(first, init, varg1, idata1);
    tbb::parallel_reduce(
        tbb::blocked_range<std::size_t>(0, num_jobs, grainsize), worker);
    return std::move(worker.sum_);
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
 * So the parallel_reduce_sum signature should look like
 *
 * template <class ReduceFunction, class InputIt, class T, ...>
 * constexpr T parallel_reduce_sum(InputIt first, InputIt last, T init,
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
template <class ReduceFunction, class InputIt, class T, class Arg1>
constexpr T parallel_reduce_sum(InputIt first, InputIt last, T init,
                                std::size_t grainsize,
                                const std::vector<Arg1>& varg1,
                                const std::vector<int>& idata1) {
  typedef T return_base_t;
  return internal::parallel_reduce_sum_impl<ReduceFunction, InputIt, T, Arg1,
                                            return_base_t>()(
      first, last, init, grainsize, varg1, idata1);
}

}  // namespace math
}  // namespace stan

#endif
