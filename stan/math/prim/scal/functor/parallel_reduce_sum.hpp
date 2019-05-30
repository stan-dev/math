#ifndef STAN_MATH_PRIM_SCAL_FUNCTOR_PARALLEL_REDUCE_SUM_HPP
#define STAN_MATH_PRIM_SCAL_FUNCTOR_PARALLEL_REDUCE_SUM_HPP

#include <stan/math/prim/scal/meta/return_type.hpp>
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
template <class InputIt, class T, class BinaryFunction, class T_return_type>
struct parallel_reduce_sum_impl {};

template <class InputIt, class T, class BinaryFunction>
struct parallel_reduce_sum_impl<InputIt, T, BinaryFunction, double> {
  struct recursive_reducer {
    InputIt first_;
    const T& init_;
    const BinaryFunction& f_;
    T sum_;

    recursive_reducer(InputIt first, const T& init, const BinaryFunction& f)
        : first_(first), init_(init), f_(f), sum_(init) {}

    recursive_reducer(recursive_reducer& other, tbb::split)
        : first_(other.first_), init_(other.init_), f_(other.f_), sum_(init_) {}

    void operator()(const tbb::blocked_range<size_t>& r) {
      if (r.empty())
        return;
      auto start = first_;
      std::advance(start, r.begin());
      auto end = first_;
      std::advance(end, r.end() - 1);
      sum_ += f_(*start, *end);
    }

    void join(const recursive_reducer& child) { sum_ += child.sum_; }
  };

  T operator()(InputIt first, InputIt last, T init, BinaryFunction f,
               std::size_t grainsize) const {
    const std::size_t num_jobs = std::distance(first, last);

    recursive_reducer worker(first, init, f);

    tbb::parallel_reduce(
        tbb::blocked_range<std::size_t>(0, num_jobs, grainsize), worker);

    return std::move(worker.sum_);
  }
};

}  // namespace internal

template <class InputIt, class T, class BinaryFunction>
constexpr T parallel_reduce_sum(InputIt first, InputIt last, T init,
                                BinaryFunction f, std::size_t grainsize = 1) {
  typedef T return_base_t;
  return internal::parallel_reduce_sum_impl<InputIt, T, BinaryFunction,
                                            return_base_t>()(first, last, init,
                                                             f, grainsize);
}

}  // namespace math
}  // namespace stan

#endif
