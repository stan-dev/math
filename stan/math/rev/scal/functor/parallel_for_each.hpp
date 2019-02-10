#ifndef STAN_MATH_REV_SCAL_FUNCTOR_PARALLEL_FOR_EACH_HPP
#define STAN_MATH_REV_SCAL_FUNCTOR_PARALLEL_FOR_EACH_HPP

#include <stan/math/prim/scal/functor/parallel_for_each.hpp>
#include <stan/math/prim/mat/fun/typedefs.hpp>

#include <stan/math/rev/core/nest_chainablestack.hpp>

#include <boost/iterator/counting_iterator.hpp>

#include <tbb/task_arena.h>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>

#include <iostream>
#include <iterator>
#include <map>
#include <tuple>
#include <vector>

namespace stan {
namespace math {
namespace internal {

template <class InputIt, class UnaryFunction>
struct parallel_map_impl<InputIt, UnaryFunction, var> {
  auto operator()(InputIt first, InputIt last, UnaryFunction f) const {
    typedef decltype(f(*first)) T_return_elem;
    typedef std::vector<decltype(f(*first))> T_return;
    typedef boost::counting_iterator<std::size_t> count_iter;

    const std::size_t num_jobs = std::distance(first, last);

    std::cout
        << "Running NEW var parallel_for_each implementation (non-nestable)..."
        << std::endl;

    typedef ChainableStack::AutodiffStackStorage chainablestack_t;
    typedef ChainableStack::AutodiffStackQueue chainablequeue_t;

    T_return f_eval(num_jobs);

    // All AD terms are written to thread-local AD tapes which are all
    // stored as part of the parent nochain stacks.
    chainablequeue_t& parent_queue = ChainableStack::queue();

    // we could tweak this for the
    // parent thread which can write
    // directly to it's own tape
    tbb::enumerable_thread_specific<std::shared_ptr<chainablestack_t>>
        child_stacks(
            [&parent_queue]() { return parent_queue.get_nochain_stack(); });

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
                new chainablestack_t());
            recover_memory_nested();
            throw;
          }
        });

    chainablestack_t& parent_tape = ChainableStack::instance();

    child_stacks.combine_each(
        [&parent_tape](const std::shared_ptr<chainablestack_t>& other_stack) {
          parent_tape.var_stack_.insert(parent_tape.var_stack_.end(),
                                        other_stack->var_stack_.begin(),
                                        other_stack->var_stack_.end());
          other_stack->var_stack_.clear();
          parent_tape.var_nochain_stack_.insert(
              parent_tape.var_nochain_stack_.end(),
              other_stack->var_nochain_stack_.begin(),
              other_stack->var_nochain_stack_.end());
          other_stack->var_nochain_stack_.clear();
        });

    return std::move(f_eval);
  }
};

}  // namespace internal
}  // namespace math
}  // namespace stan

#endif
