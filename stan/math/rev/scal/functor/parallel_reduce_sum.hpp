#ifndef STAN_MATH_REV_SCAL_FUNCTOR_PARALLEL_REDUCE_SUM_HPP
#define STAN_MATH_REV_SCAL_FUNCTOR_PARALLEL_REDUCE_SUM_HPP

#include <stan/math/prim/scal/functor/parallel_reduce_sum.hpp>
#include <stan/math/prim/mat/fun/typedefs.hpp>

#include <stan/math/rev/core/scoped_chainablestack.hpp>

#include <boost/iterator/counting_iterator.hpp>

#include <tbb/parallel_reduce.h>
#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>

#include <iostream>
#include <iterator>
#include <map>
#include <tuple>
#include <vector>
#include <atomic>

namespace stan {
namespace math {
namespace internal {

template <class InputIt, class T, class BinaryFunction>
struct parallel_reduce_sum_impl<InputIt, T, BinaryFunction, var> {
  typedef ChainableStack::AutodiffStackStorage chainablestack_t;
  typedef tbb::enumerable_thread_specific<ScopedChainableStack>
      tls_scoped_stack_t;

  struct recursive_reducer {
    InputIt first_;
    const BinaryFunction& f_;
    tls_scoped_stack_t& tls_scoped_stack_;
    std::vector<var> sum_terms_;
    /*
    const std::size_t worker_id_;

    static std::size_t get_worker_count() {
      static std::atomic<std::size_t> worker_count{0};
      return worker_count.fetch_add(1);
    }
    */

    recursive_reducer(InputIt first, const T& init, const BinaryFunction& f,
                      tls_scoped_stack_t& tls_scoped_stack)
        : first_(first),
          f_(f),
          tls_scoped_stack_(tls_scoped_stack),
          sum_terms_(1, init)  //,
                               // worker_id_(get_worker_count())
    {}

    recursive_reducer(recursive_reducer& other, tbb::split)
        : first_(other.first_),
          f_(other.f_),
          tls_scoped_stack_(other.tls_scoped_stack_) {
      // std::cout << "Splitting off some work with worker id " << worker_id_ <<
      // "..." << std::endl;
    }

    void operator()(const tbb::blocked_range<size_t>& r) {
      // std::cout << "Summing " << r.begin() << " - " << r.end() << " with
      // worker " << worker_id_ << std::endl;
      tls_scoped_stack_.local().execute([&] {
        auto start = first_;
        std::advance(start, r.begin());
        auto end = first_;
        std::advance(end, r.end());
        sum_terms_.emplace_back(f_(*start, *end));
      });
    }

    void join(recursive_reducer& child) {
      // std::cout << "Joining a child " << child.worker_id_ << " into worker "
      // << worker_id_ << std::endl;
      tls_scoped_stack_.local().execute([&] {
        if (child.sum_terms_.size() == 1) {
          sum_terms_.emplace_back(child.sum_terms_[0]);
        } else {
          sum_terms_.emplace_back(sum(child.sum_terms_));
        }
        child.sum_terms_.clear();
      });
    }
  };

  T operator()(InputIt first, InputIt last, T init, BinaryFunction f) const {
    const std::size_t num_jobs = std::distance(first, last);

    // std::cout << "Running NEW var parallel_reduce_sum implementation ..."
    //          << std::endl;

    // All AD terms are written to thread-local AD tapes which are all
    // stored as part of the parent nochain stacks.
    chainablestack_t& parent_stack = ChainableStack::instance();

    tls_scoped_stack_t child_stacks(
        [&parent_stack]() { return ScopedChainableStack(parent_stack); });

    recursive_reducer worker(first, init, f, child_stacks);

    // static tbb::affinity_partitioner partitioner;
    tbb::auto_partitioner partitioner;
    // tbb::static_partitioner partitioner;
    // it seems that best performance is attained with the simple
    // partititioner and a reasonable grainsize
    // tbb::simple_partitioner partitioner;

    // TODO: make grainsize a parameter??!!!
    tbb::parallel_reduce(tbb::blocked_range<std::size_t>(0, num_jobs, 1),
                         worker, partitioner);

    child_stacks.combine_each(
        [&parent_stack](ScopedChainableStack& child_scoped_stack) {
          child_scoped_stack.append_to_parent();
        });

    return sum(worker.sum_terms_);
  }
};

}  // namespace internal
}  // namespace math
}  // namespace stan

#endif
