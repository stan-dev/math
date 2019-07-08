#ifndef STAN_MATH_REV_SCAL_FUNCTOR_PARALLEL_REDUCE_SUM_HPP
#define STAN_MATH_REV_SCAL_FUNCTOR_PARALLEL_REDUCE_SUM_HPP

#include <stan/math/prim/scal/functor/parallel_reduce_sum.hpp>
#include <stan/math/prim/scal/functor/parallel_for_each.hpp>
#include <stan/math/rev/scal/functor/parallel_for_each.hpp>
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

/*
template <class InputIt, class T, class BinaryFunction>
struct parallel_reduce_sum_impl<InputIt, T, BinaryFunction, var> {
  typedef ChainableStack::AutodiffStackStorage chainablestack_t;

  typedef tbb::enumerable_thread_specific<ScopedChainableStack>
      tls_scoped_stack_t;

  T operator()(InputIt first, InputIt last, T init, BinaryFunction f,
               std::size_t grainsize) const {

    // avoid TBB overhead in case of only 1 core
    if (get_num_threads() == 1)
      return init + f(*first, *(last-1));

    const std::size_t num_jobs = std::distance(first, last);

    // All AD terms are written to thread-local AD tapes which are all
    // stored as part of the parent nochain stacks.
    chainablestack_t& parent_stack = *ChainableStack::instance_;

    tls_scoped_stack_t child_stacks(
        [&parent_stack]() { return ScopedChainableStack(parent_stack); });

    tbb::enumerable_thread_specific<var> partial_sums(var(0.0));

    // static tbb::affinity_partitioner partitioner;
    tbb::auto_partitioner partitioner;
    //tbb::static_partitioner partitioner;
    // it seems that best performance is attained with the simple
    // partititioner and a reasonable grainsize
    // tbb::simple_partitioner partitioner;

    // TODO: make grainsize a parameter??!!!
    tbb::parallel_for(
        tbb::blocked_range<std::size_t>(0, num_jobs, grainsize),
        [&](const tbb::blocked_range<size_t>& r) {
          child_stacks.local().execute([&](){
                                         auto start_iter = first;
                                         std::advance(start_iter, r.begin());
                                         auto end_iter = first;
                                         std::advance(end_iter, r.end()-1);
                                         partial_sums.local() += f(*start_iter,
*end_iter);
                                       });
        }, partitioner);

    child_stacks.combine_each(
        [&parent_stack](ScopedChainableStack& child_scoped_stack) {
          child_scoped_stack.append_to_parent();
        });

    var total = init;
    partial_sums.combine_each([&total](var& term) {
                                total += term;
                              });

    return total;
  }

};

*/

/*
template <class InputIt, class T, class BinaryFunction>
struct parallel_reduce_sum_impl<InputIt, T, BinaryFunction, var> {
  T operator()(InputIt first, InputIt last, T init, BinaryFunction f,
               std::size_t grainsize) const {

    // avoid TBB overhead in case of only 1 core
    if (get_num_threads() == 1)
      return init + f(*first, *(last-1));

    typedef boost::counting_iterator<std::size_t> count_iter;

    const std::size_t num_terms = std::distance(first, last);

    const std::size_t blocksize = grainsize == 0
                                  ? num_terms / get_num_threads()
                                  : grainsize ;

    const std::size_t num_blocks = num_terms / blocksize;
    //std::size_t extra_terms = num_terms % grainsize;

    // if we allow non-determinism in outputs, then we could code the
    // parallel_for call here directly and take advantage of blocked
    // reduce sweeps (so that we get larger reduces in a single go)

    std::vector<T> partial_sums = parallel_map(count_iter(0),
                                               count_iter(num_blocks),
                                               [&](int block) -> T {
                                                 int start = block * blocksize;
                                                 int end = block == num_blocks-1
                                                           ? num_terms
                                                           : (block+1) *
blocksize;

                                                 auto start_iter = first;
                                                 std::advance(start_iter,
start); auto end_iter = first; std::advance(end_iter, end-1);

                                                 return f(*start_iter,
*end_iter);
                                               });

    return init + sum(partial_sums);
  }
};
*/

/*
template <class InputIt, class T, class BinaryFunction>
struct parallel_reduce_sum_impl<InputIt, T, BinaryFunction, var> {
  typedef ChainableStack::AutodiffStackStorage chainablestack_t;
  typedef tbb::enumerable_thread_specific<ScopedChainableStack>
      tls_scoped_stack_t;

  struct recursive_reducer {
    const InputIt first_;
    const BinaryFunction& f_;
    tls_scoped_stack_t& tls_scoped_stack_;
    std::vector<var> sum_terms_;

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
      if (r.empty())
        return;
      // std::cout << "Summing " << r.begin() << " - " << r.end() << " with
      // worker " << worker_id_ << std::endl;
      tls_scoped_stack_.local().execute([&] {
        auto start = first_;
        std::advance(start, r.begin());
        auto end = first_;
        std::advance(end, r.end() - 1);
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

  T operator()(InputIt first, InputIt last, T init, BinaryFunction f,
               std::size_t grainsize) const {
    const std::size_t num_jobs = std::distance(first, last);

    // std::cout << "Running NEW var parallel_reduce_sum implementation ..."
    //          << std::endl;

    // All AD terms are written to thread-local AD tapes which are all
    // stored as part of the parent nochain stacks.
    chainablestack_t& parent_stack = *ChainableStack::instance_;

    tls_scoped_stack_t child_stacks(
        [&parent_stack]() { return ScopedChainableStack(parent_stack); });

    recursive_reducer worker(first, init, f, child_stacks);

    // static tbb::affinity_partitioner partitioner;
    // tbb::auto_partitioner partitioner;
    tbb::static_partitioner partitioner;
    // it seems that best performance is attained with the simple
    // partititioner and a reasonable grainsize
    // tbb::simple_partitioner partitioner;

    // TODO: make grainsize a parameter??!!!
    tbb::parallel_deterministic_reduce(
        tbb::blocked_range<std::size_t>(0, num_jobs, grainsize), worker,
        partitioner);

    child_stacks.combine_each(
        [&parent_stack](ScopedChainableStack& child_scoped_stack) {
          child_scoped_stack.append_to_parent();
        });

    return sum(worker.sum_terms_);
  }
};
*/

/*
template <class InputIt, class T, class BinaryFunction>
struct parallel_reduce_sum_impl<InputIt, T, BinaryFunction, var> {
  typedef ChainableStack::AutodiffStackStorage chainablestack_t;
  typedef tbb::concurrent_vector<ScopedChainableStack>
      scoped_stack_queue_t;

  struct recursive_reducer {
    const InputIt first_;
    BinaryFunction f_;
    chainablestack_t& parent_stack_;
    scoped_stack_queue_t& scoped_stack_queue_;
    ScopedChainableStack& scoped_stack_;
    std::vector<var> sum_terms_;

    recursive_reducer(InputIt first, const T& init, BinaryFunction f,
                      chainablestack_t& parent_stack,
                      scoped_stack_queue_t& scoped_stack_queue)
        : first_(first),
          f_(f),
          parent_stack_(parent_stack),
          scoped_stack_queue_(scoped_stack_queue),
          scoped_stack_(*scoped_stack_queue_.emplace_back(parent_stack_)),
          sum_terms_(1, init)  //,
                               // worker_id_(get_worker_count())
    {}

    recursive_reducer(recursive_reducer& other, tbb::split)
        : first_(other.first_),
          f_(other.f_),
          parent_stack_(other.parent_stack_),
          scoped_stack_queue_(other.scoped_stack_queue_),
          scoped_stack_(*scoped_stack_queue_.emplace_back(parent_stack_)) {
      // std::cout << "Splitting off some work with worker id " << worker_id_ <<
      // "..." << std::endl;
    }

    void operator()(const tbb::blocked_range<size_t>& r) {
      if (r.empty())
        return;
      // std::cout << "Summing " << r.begin() << " - " << r.end() << " with
      // worker " << worker_id_ << std::endl;
      scoped_stack_.execute([&] {
        auto start = first_;
        std::advance(start, r.begin());
        auto end = first_;
        std::advance(end, r.end() - 1);
        sum_terms_.emplace_back(f_(*start, *end));
      });
    }

    void join(recursive_reducer& child) {
      // std::cout << "Joining a child " << child.worker_id_ << " into worker "
      // << worker_id_ << std::endl;
      scoped_stack_.execute([&] {
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
    chainablestack_t& parent_stack = *ChainableStack::instance_;

    scoped_stack_queue_t child_stacks;

    recursive_reducer worker(first, init, f, parent_stack, child_stacks);

    //static tbb::affinity_partitioner partitioner;
    tbb::auto_partitioner partitioner;
    //tbb::static_partitioner partitioner;
    // it seems that best performance is attained with the simple
    // partititioner and a reasonable grainsize
    //tbb::simple_partitioner partitioner;

    // TODO: make grainsize a parameter??!!!
    tbb::parallel_reduce(tbb::blocked_range<std::size_t>(0, num_jobs),
                         worker, partitioner);

    for(std::size_t i=0; i != child_stacks.size(); ++i) {
      child_stacks[i].append_to_parent();
    }

    return sum(worker.sum_terms_);
  }
};
*/

}  // namespace internal
}  // namespace math
}  // namespace stan

#endif
