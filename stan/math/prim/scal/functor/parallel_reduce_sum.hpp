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
/*
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
*/

template <class InputIt, class T, class BinaryFunction>
constexpr T parallel_reduce_sum(InputIt first, InputIt last, T init,
                                BinaryFunction f, std::size_t grainsize = 0) {
  // avoid TBB overhead in case of only 1 core
  const std::size_t num_threads = internal::get_num_threads();
  if (num_threads == 1)
    return f(*first, *(last - 1));

  typedef boost::counting_iterator<std::size_t> count_iter;

  const std::size_t num_terms = std::distance(first, last);

  const std::size_t blocksize
      = grainsize == 0 ? num_terms / num_threads : grainsize;

  const std::size_t num_blocks = num_terms / blocksize;
  // std::size_t extra_terms = num_terms % grainsize;

  // if we allow non-determinism in outputs, then we could code the
  // parallel_for call here directly and take advantage of blocked
  // reduce sweeps (so that we get larger reduces in a single go)

  // TODO: Should we use grainsize of the parallel for loop?

  std::vector<T> partial_sums = parallel_map(
      count_iter(0), count_iter(num_blocks), [&](int block) -> T {
        int start = block * blocksize;
        int end = block == num_blocks - 1 ? num_terms : (block + 1) * blocksize;

        auto start_iter = first;
        std::advance(start_iter, start);
        auto end_iter = first;
        std::advance(end_iter, end - 1);

        return f(*start_iter, *end_iter);
      });

  return init + sum(partial_sums);
}

}  // namespace math
}  // namespace stan

#endif
