#ifndef STAN_MATH_REV_SCAL_FUNCTOR_PARALLEL_REDUCE_SUM_HPP
#define STAN_MATH_REV_SCAL_FUNCTOR_PARALLEL_REDUCE_SUM_HPP

#include <stan/math/prim/scal/functor/parallel_reduce_sum.hpp>
#include <stan/math/prim/mat/fun/typedefs.hpp>

#include <stan/math/rev/core/nest_chainablestack.hpp>

#include <boost/iterator/counting_iterator.hpp>

#include <tbb/task_arena.h>
#include <tbb/concurrent_vector.h>
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
  typedef tbb::concurrent_vector<var> partial_sum_buffer_t;
  typedef ChainableStack::AutodiffStackStorage chainablestack_t;
  typedef ChainableStack::AutodiffStackQueue chainablequeue_t;
  typedef std::shared_ptr<chainablestack_t> stack_ptr_t;
  typedef tbb::enumerable_thread_specific<stack_ptr_t> tls_stack_ptr_t;
  typedef tbb::enumerable_thread_specific<ScopedChainableStack>
      tls_scoped_stack_t;

  struct recursive_reducer {
    InputIt first_;
    const BinaryFunction& f_;
    // tls_stack_ptr_t& tls_stack_ptr_;
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
                      // tls_stack_ptr_t& tls_stack_ptr
                      tls_scoped_stack_t& tls_scoped_stack)
        : first_(first),
          f_(f),
          // tls_stack_ptr_(tls_stack_ptr),
          tls_scoped_stack_(tls_scoped_stack),
          sum_terms_(1, init)  //,
                               // worker_id_(get_worker_count())
    {}

    recursive_reducer(recursive_reducer& other, tbb::split)
        : first_(other.first_),
          f_(other.f_),
          // tls_stack_ptr_(other.tls_stack_ptr_)
          tls_scoped_stack_(other.tls_scoped_stack_) {
      // std::cout << "Splitting off some work with worker id " << worker_id_ <<
      // "..." << std::endl;
    }

    void operator()(const tbb::blocked_range<size_t>& r) {
      // std::cout << "Summing " << r.begin() << " - " << r.end() << " with
      // worker " << worker_id_ << std::endl;
      // move current AD stack out
      // of the way in child thread
      chainablequeue_t& local_queue = ChainableStack::queue();

      // TODO?? Add a start_managed_nested(external_ad_tape) ???
      /*
      try {
        start_nested();
        const std::size_t nested_stack_instance = local_queue.current_instance_;
        stack_ptr_t nested_stack
            = local_queue.instance_stack_[nested_stack_instance];
        local_queue.instance_stack_[nested_stack_instance]
            = tls_stack_ptr_.local();
        ChainableStack::instance_
            = local_queue.instance_stack_[nested_stack_instance].get();

        auto start = first_;
        std::advance(start, r.begin());
        auto end = first_;
        std::advance(end, r.end());
        sum_terms_.emplace_back(f_(*start, *end));

        local_queue.instance_stack_[nested_stack_instance] = nested_stack;
        ChainableStack::instance_ = nested_stack.get();
        recover_memory_nested();
      } catch (const std::exception& e) {
        local_queue.instance_stack_[local_queue.current_instance_].reset(
            new chainablestack_t(ChainableStack::queue().stack_id_));
        recover_memory_nested();
        throw;
      }
      */

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

      /*
      chainablequeue_t& local_queue = ChainableStack::queue();

      try {
        start_nested();
        const std::size_t nested_stack_instance = local_queue.current_instance_;
        stack_ptr_t nested_stack
            = local_queue.instance_stack_[nested_stack_instance];
        local_queue.instance_stack_[nested_stack_instance]
            = tls_stack_ptr_.local();
        ChainableStack::instance_
            = local_queue.instance_stack_[nested_stack_instance].get();

        if (child.sum_terms_.size() == 1) {
          sum_terms_.emplace_back(child.sum_terms_[0]);
        } else {
          sum_terms_.emplace_back(sum(child.sum_terms_));
        }
        child.sum_terms_.clear();

        local_queue.instance_stack_[nested_stack_instance] = nested_stack;
        ChainableStack::instance_ = nested_stack.get();
        recover_memory_nested();
      } catch (const std::exception& e) {
        local_queue.instance_stack_[local_queue.current_instance_].reset(
            new chainablestack_t(ChainableStack::queue().stack_id_));
        recover_memory_nested();
        throw;
      }
      */

      tls_scoped_stack_.local().execute([&] {
        if (child.sum_terms_.size() == 1) {
          sum_terms_.emplace_back(child.sum_terms_[0]);
        } else {
          sum_terms_.emplace_back(sum(child.sum_terms_));
        }
        child.sum_terms_.clear();
      });

      /*
      child_stack_ptr_.emplace_back(child.worker_stack_ptr_);
      child_stack_ptr_.insert(child_stack_ptr_.end(),
                              child.child_stack_ptr_.begin(),
                              child.child_stack_ptr_.end());
      */

      /*
         worker_stack_ptr_->var_stack_.insert(
             worker_stack_ptr_->var_stack_.end(),
             child.worker_stack_ptr_->var_stack_.begin(),
             child.worker_stack_ptr_->var_stack_.end());
         child.worker_stack_ptr_->var_stack_.clear();
         worker_stack_ptr_->var_nochain_stack_.insert(
             worker_stack_ptr_->var_nochain_stack_.end(),
             child.worker_stack_ptr_->var_nochain_stack_.begin(),
             child.worker_stack_ptr_->var_nochain_stack_.end());
         child.worker_stack_ptr_->var_nochain_stack_.clear();
      */
    }
  };

  T operator()(InputIt first, InputIt last, T init, BinaryFunction f) const {
    const std::size_t num_jobs = std::distance(first, last);

    // std::cout << "Running NEW var parallel_reduce_sum implementation ..."
    //          << std::endl;

    // All AD terms are written to thread-local AD tapes which are all
    // stored as part of the parent nochain stacks.
    chainablequeue_t& parent_queue = ChainableStack::queue();
    chainablestack_t& parent_stack = ChainableStack::instance();

    // we could tweak this for the
    // parent thread which can write
    // directly to it's own tape
    // tbb::enumerable_thread_specific<std::shared_ptr<chainablestack_t>>
    //    child_stacks(
    //        [&parent_stack]() { return parent_stack.get_child_stack(); });
    // todo: need get child stacks with the stack_id of the current tape!!!

    // partial_sum_buffer_t partial_sums;

    // tls_stack_ptr_t child_stacks(
    //    [&parent_stack]() { return parent_stack.get_child_stack();
    //    });

    tls_scoped_stack_t child_stacks(
        [&parent_stack]() { return ScopedChainableStack(parent_stack); });

    recursive_reducer worker(first, init, f, child_stacks);

    // static tbb::affinity_partitioner partitioner;
    tbb::auto_partitioner partitioner;
    // tbb::static_partitioner partitioner;
    // it seems that best performance is attained with the simple
    // partititioner and a reasonable grainsize
    // tbb::simple_partitioner partitioner;

    // TODO: add thread_local AD tapes

    // TODO: make grainsize a parameter??!!!
    tbb::parallel_reduce(tbb::blocked_range<std::size_t>(0, num_jobs, 1),
                         worker, partitioner);

    /*
    child_stacks.combine_each([&parent_stack](const stack_ptr_t& other_stack) {
      parent_stack.var_stack_.insert(parent_stack.var_stack_.end(),
                                     other_stack->var_stack_.begin(),
                                     other_stack->var_stack_.end());
      other_stack->var_stack_.clear();
      parent_stack.var_nochain_stack_.insert(
          parent_stack.var_nochain_stack_.end(),
          other_stack->var_nochain_stack_.begin(),
          other_stack->var_nochain_stack_.end());
      other_stack->var_nochain_stack_.clear();
    });
    */

    /*
    child_stacks.combine_each([&parent_stack](const ScopedChainableStack&
    other_scoped_stack) { const stack_ptr_t& other_stack =
    other_scoped_stack.local_stack_;
                                parent_stack.var_stack_.insert(parent_stack.var_stack_.end(),
                                                               other_stack->var_stack_.begin(),
                                                               other_stack->var_stack_.end());
                                other_stack->var_stack_.clear();
                                parent_stack.var_nochain_stack_.insert(
                                    parent_stack.var_nochain_stack_.end(),
                                    other_stack->var_nochain_stack_.begin(),
                                    other_stack->var_nochain_stack_.end());
                                other_stack->var_nochain_stack_.clear();
                              });
    */

    child_stacks.combine_each(
        [&parent_stack](ScopedChainableStack& child_scoped_stack) {
          child_scoped_stack.append_to_stack(parent_stack);
        });

    return sum(worker.sum_terms_);

    /*
    tbb::parallel_for(
        tbb::blocked_range<std::size_t>(0, num_jobs),
        [&](const tbb::blocked_range<size_t>& r) {
          // move current AD stack out
          // of the way in child thread
          chainablequeue_t& local_queue = ChainableStack::queue();

          try {
            start_nested();
            std::shared_ptr<chainablestack_t> managed_local_stack
                = child_stacks.local();
            const std::size_t nested_stack_instance
                = local_queue.current_instance_;
            std::shared_ptr<chainablestack_t> nested_stack
                = local_queue.instance_stack_[nested_stack_instance];
            local_queue.instance_stack_[nested_stack_instance]
                = managed_local_stack;
            ChainableStack::instance_ = managed_local_stack.get();
            auto elem = first;
            std::advance(elem, r.begin());
            for (std::size_t i = r.begin(); i != r.end(); ++elem, ++i) {
              f_eval[i] = f(*elem);
            }
            local_queue.instance_stack_[nested_stack_instance] = nested_stack;
            ChainableStack::instance_ = nested_stack.get();
            recover_memory_nested();
          } catch (const std::exception& e) {
            local_queue.instance_stack_[local_queue.current_instance_].reset(
                new chainablestack_t(ChainableStack::queue().stack_id_));
            recover_memory_nested();
            throw;
          }
        });

    child_stacks.combine_each(
        [&parent_stack](const std::shared_ptr<chainablestack_t>& other_stack) {
          parent_stack.var_stack_.insert(parent_stack.var_stack_.end(),
                                         other_stack->var_stack_.begin(),
                                         other_stack->var_stack_.end());
          other_stack->var_stack_.clear();
          parent_stack.var_nochain_stack_.insert(
              parent_stack.var_nochain_stack_.end(),
              other_stack->var_nochain_stack_.begin(),
              other_stack->var_nochain_stack_.end());
          other_stack->var_nochain_stack_.clear();
        });

    return std::move(f_eval);
    */
  }
};

}  // namespace internal
}  // namespace math
}  // namespace stan

#endif
