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
    std::vector<Arg1>& varg1_;
    const std::vector<int>& idata1_;
    T sum_;
    typedef typename std::iterator_traits<InputIt>::value_type elem_t;

    recursive_reducer(InputIt first, const T& init, std::vector<Arg1>& varg1,
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
               std::vector<Arg1>& varg1, const std::vector<int>& idata1) const {
    const std::size_t num_jobs = std::distance(first, last);
    recursive_reducer worker(first, init, varg1, idata1);
    tbb::parallel_reduce(
        tbb::blocked_range<std::size_t>(0, num_jobs, grainsize), worker);
    return std::move(worker.sum_);
  }
};
}  // namespace internal

template <class ReduceFunction, class InputIt, class T, class Arg1>
constexpr T parallel_reduce_sum(InputIt first, InputIt last, T init,
                                std::size_t grainsize, std::vector<Arg1>& varg1,
                                const std::vector<int>& idata1) {
  typedef T return_base_t;
  return internal::parallel_reduce_sum_impl<ReduceFunction, InputIt, T, Arg1,
                                            return_base_t>()(
      first, last, init, grainsize, varg1, idata1);
}

}  // namespace math
}  // namespace stan

#endif
