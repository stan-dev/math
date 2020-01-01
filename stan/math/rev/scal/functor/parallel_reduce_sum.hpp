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
  using ops_partials_t
      = operands_and_partials<std::vector<Arg1>, const std::vector<int>>;

  struct recursive_reducer {
    InputIt first_;
    std::vector<Arg1>& varg1_;
    const std::vector<int>& idata1_;
    ops_partials_t terms_partials_;
    double terms_sum_;
    typedef typename std::iterator_traits<InputIt>::value_type elem_t;

    recursive_reducer(InputIt first, const T& init, std::vector<Arg1>& varg1,
                      const std::vector<int>& idata1)
        : first_(first),
          varg1_(varg1),
          idata1_(idata1),
          terms_partials_(varg1_, idata1_),
          terms_sum_(init.val()) {}

    recursive_reducer(recursive_reducer& other, tbb::split)
        : first_(other.first_),
          varg1_(other.varg1_),
          idata1_(other.idata1_),
          terms_partials_(varg1_, idata1_),
          terms_sum_(0.0) {}

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

        // create a deep copy of arg1 which is not tied to the outer
        // AD tree
        // todo: these copies can be put onto a thread-local storage
        // object such that we only create these once per thread which
        // should significantly boost performance as many copies are
        // avoided... this comes at the cost that the adjoints have to
        // be zeroed "manually"
        std::vector<Arg1> varg1;
        varg1.reserve(varg1_.size());

        for (Arg1& elem : varg1_)
          varg1.emplace_back(var(new vari(elem.val(), false)));

        T sub_sum_v = ReduceFunction()(r.begin(), r.end() - 1, sub_slice, varg1,
                                       idata1_);

        sub_sum_v.grad();

        terms_sum_ += sub_sum_v.val();

        for (std::size_t i = 0; i != varg1.size(); ++i)
          terms_partials_.edge1_.partials_[i] += varg1[i].adj();

      } catch (const std::exception& e) {
        recover_memory_nested();
        throw;
      }
      recover_memory_nested();
    }

    void join(const recursive_reducer& child) {
      terms_sum_ += child.terms_sum_;

      for (std::size_t i = 0; i != varg1_.size(); ++i)
        terms_partials_.edge1_.partials_[i]
            += child.terms_partials_.edge1_.partials_[i];
    }
  };

  T operator()(InputIt first, InputIt last, T init, std::size_t grainsize,
               std::vector<Arg1>& varg1, const std::vector<int>& idata) const {
    const std::size_t num_jobs = std::distance(first, last);
    recursive_reducer worker(first, init, varg1, idata);
    tbb::parallel_reduce(
        tbb::blocked_range<std::size_t>(0, num_jobs, grainsize), worker);
    return worker.terms_partials_.build(worker.terms_sum_);
  }
};
}  // namespace internal

}  // namespace math
}  // namespace stan

#endif
