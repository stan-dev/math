#ifndef STAN_MATH_REV_SCAL_FUNCTOR_PARALLEL_REDUCE_SUM_HPP
#define STAN_MATH_REV_SCAL_FUNCTOR_PARALLEL_REDUCE_SUM_HPP

#include <stan/math/prim/meta/return_type.hpp>
#include <stan/math/rev/mat/fun/typedefs.hpp>

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

template <class ReduceFunction, class InputIt, class T, class Arg1>
struct parallel_reduce_sum_impl<ReduceFunction, InputIt, T, Arg1, var> {
  using ops_partials_t = operands_and_partials<Arg1>;

  struct recursive_reducer {
    InputIt first_;
    const T& init_;
    Arg1& arg1_;
    std::vector<double> terms_adjs_;
    std::vector<double> terms_vals_;
    typedef typename std::iterator_traits<InputIt>::value_type elem_t;

    recursive_reducer(InputIt first, const T& init, Arg1& arg1)
        : first_(first),
          init_(init),
          arg1_(arg1),
          terms_adjs_(1, 0.0),
          terms_vals_(1, init_.val()) {}

    recursive_reducer(recursive_reducer& other, tbb::split)
        : first_(other.first_),
          init_(other.init_),
          arg1_(other.arg1_),
          terms_adjs_(1, 0.0),
          terms_vals_(1, init_.val()) {}
    // terms_adjs_(other.terms_adjs_), terms_vals_(other.terms_vals_) {}

    void operator()(const tbb::blocked_range<size_t>& r) {
      if (r.empty())
        return;

      auto start = first_;
      std::advance(start, r.begin());
      auto end = first_;
      std::advance(end, r.end());

      try {
        start_nested();

        // elem_t should better be a non-var!!! Otherwise we interfere
        // with the outer AD tree as we mess with the adjoints (and
        // right now we discard that it is a var if it would be)
        std::vector<elem_t> sub_slice(start, end);

        // create a copy of arg1 which is not tied to the out AD tree
        Arg1 arg1(var(new vari(arg1_.val(), false)));

        T sub_sum_v = ReduceFunction()(r.begin(), r.end() - 1, sub_slice, arg1);

        sub_sum_v.grad();

        terms_vals_.push_back(sub_sum_v.val());
        terms_adjs_.push_back(arg1.adj());

      } catch (const std::exception& e) {
        recover_memory_nested();
        throw;
      }
      recover_memory_nested();
    }
    void join(const recursive_reducer& child) {
      terms_adjs_.insert(terms_adjs_.end(), child.terms_adjs_.begin(),
                         child.terms_adjs_.end());
      terms_vals_.insert(terms_vals_.end(), child.terms_vals_.begin(),
                         child.terms_vals_.end());
    }
  };

  T operator()(InputIt first, InputIt last, T init, std::size_t grainsize,
               Arg1& arg1) const {
    const std::size_t num_jobs = std::distance(first, last);
    recursive_reducer worker(first, init, arg1);
    tbb::parallel_reduce(
        tbb::blocked_range<std::size_t>(0, num_jobs, grainsize), worker);
    std::vector<T> terms(worker.terms_adjs_.size());
    ops_partials_t ops_sum(arg1);
    ops_sum.edge1_.partials_[0] = sum(worker.terms_adjs_);
    return ops_sum.build(sum(worker.terms_vals_));
  }
};
}  // namespace internal

}  // namespace math
}  // namespace stan

#endif
