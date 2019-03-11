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

  struct recursive_reducer {
    InputIt first_;
    const BinaryFunction& f_;
    chainablestack_t& parent_stack_;
    stack_ptr_t worker_stack_ptr_;
    std::vector<var> sum_terms_;
    std::vector<stack_ptr_t> child_stack_ptr_;
    /*
    const std::size_t worker_id_;

    static std::size_t get_worker_count() {
      static std::atomic<std::size_t> worker_count{0};
      return worker_count.fetch_add(1);
    }
    */

    recursive_reducer(InputIt first, const T& init, const BinaryFunction& f,
                      chainablestack_t& parent_stack,
                      chainablequeue_t& parent_queue)
        : first_(first),
          f_(f),
          parent_stack_(parent_stack),
          worker_stack_ptr_(
              parent_queue.instance_stack_[parent_queue.current_instance_]),
          sum_terms_(1, init)  //,
                               // worker_id_(get_worker_count())
    {}

    recursive_reducer(recursive_reducer& other, tbb::split)
        : first_(other.first_),
          f_(other.f_),
          parent_stack_(other.parent_stack_),
          worker_stack_ptr_(parent_stack_.get_child_stack())  //,
    // worker_id_(get_worker_count())
    {
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
      try {
        start_nested();
        const std::size_t nested_stack_instance = local_queue.current_instance_;
        std::shared_ptr<chainablestack_t> nested_stack
            = local_queue.instance_stack_[nested_stack_instance];
        local_queue.instance_stack_[nested_stack_instance] = worker_stack_ptr_;
        ChainableStack::instance_ = worker_stack_ptr_.get();

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
    }

    void join(const recursive_reducer& child) {
      // std::cout << "Joining a child " << child.worker_id_ << " into worker "
      // << worker_id_ << std::endl;
      chainablequeue_t& local_queue = ChainableStack::queue();

      try {
        start_nested();
        const std::size_t nested_stack_instance = local_queue.current_instance_;
        std::shared_ptr<chainablestack_t> nested_stack
            = local_queue.instance_stack_[nested_stack_instance];
        local_queue.instance_stack_[nested_stack_instance] = worker_stack_ptr_;
        ChainableStack::instance_ = worker_stack_ptr_.get();

        sum_terms_.emplace_back(sum(child.sum_terms_));

        local_queue.instance_stack_[nested_stack_instance] = nested_stack;
        ChainableStack::instance_ = nested_stack.get();
        recover_memory_nested();
      } catch (const std::exception& e) {
        local_queue.instance_stack_[local_queue.current_instance_].reset(
            new chainablestack_t(ChainableStack::queue().stack_id_));
        recover_memory_nested();
        throw;
      }

      child_stack_ptr_.emplace_back(child.worker_stack_ptr_);
      child_stack_ptr_.insert(child_stack_ptr_.end(),
                              child.child_stack_ptr_.begin(),
                              child.child_stack_ptr_.end());

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

    recursive_reducer worker(first, init, f, parent_stack, parent_queue);

    // static tbb::affinity_partitioner partitioner;
    tbb::auto_partitioner partitioner;
    // tbb::static_partitioner partitioner;
    // it seems that best performance is attained with the simple
    // partititioner and a reasonable grainsize
    // tbb::simple_partitioner partitioner;

    // TODO: make grainsize a parameter??!!!
    tbb::parallel_reduce(tbb::blocked_range<std::size_t>(0, num_jobs, 1),
                         worker, partitioner);

    for (auto& child_stack_ptr : worker.child_stack_ptr_) {
      parent_stack.var_stack_.insert(parent_stack.var_stack_.end(),
                                     child_stack_ptr->var_stack_.begin(),
                                     child_stack_ptr->var_stack_.end());
      child_stack_ptr->var_stack_.clear();
      parent_stack.var_nochain_stack_.insert(
          parent_stack.var_nochain_stack_.end(),
          child_stack_ptr->var_nochain_stack_.begin(),
          child_stack_ptr->var_nochain_stack_.end());
      child_stack_ptr->var_nochain_stack_.clear();
    }

    // TODO: make this into a fast vectorized expression by
    // introducing a specialization to concurrent_vector of sum
    /*
    T sum(init);
    for (std::size_t i=0; i != partial_sums.size(); ++i) {
      sum += partial_sums[i];
    }
    */

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
