#ifndef STAN_MATH_REV_SCAL_FUNCTOR_PARALLEL_FOR_EACH_HPP
#define STAN_MATH_REV_SCAL_FUNCTOR_PARALLEL_FOR_EACH_HPP

#include <stan/math/prim/scal/functor/parallel_for_each.hpp>
#include <stan/math/prim/mat/fun/typedefs.hpp>

#include <stan/math/rev/core/nest_chainablestack.hpp>

#include <boost/iterator/counting_iterator.hpp>

#include <tbb/task_arena.h>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/combinable.h>

#include <iostream>
#include <iterator>
#include <map>
#include <vector>

namespace stan {
namespace math {
namespace internal {

template <class InputIt, class UnaryFunction>
struct parallel_map_impl<InputIt, UnaryFunction, var> {
  auto operator()(InputIt first, InputIt last, UnaryFunction f) const {
    typedef decltype(f(*first)) T_return_elem;
    typedef std::vector<decltype(f(*first))>
        T_return;
    typedef boost::counting_iterator<std::size_t> count_iter;

    const std::size_t num_jobs = std::distance(first, last);

    std::cout
        << "Running NEW var parallel_for_each implementation (non-nestable)..."
        << std::endl;

    typedef ChainableStack::AutodiffStackStorage chainablestack_t;

    T_return f_eval(num_jobs);

    tbb::this_task_arena::isolate([&]() {
                                    // stack_id => stack size map
                                    std::map<std::size_t, std::size_t> stack_starts;

                                    std::for_each(ChainableStack::thread_tapes_.begin(),
                                                  ChainableStack::thread_tapes_.end(),
                                                  [&](chainablestack_t& thread_stack) {
                                                    stack_starts.insert(std::make_pair(
                                                        thread_stack.id_, thread_stack.var_stack_.size()));
                                                  });

                                    const std::size_t parent_stack_id = ChainableStack::instance().id_;

                                    tbb::parallel_for( tbb::blocked_range<std::size_t>( 0, num_jobs ),
                                                       [&](const tbb::blocked_range<size_t>& r) {
                                                         auto elem = first;
                                                         std::advance(elem, r.begin());
                                                         for (std::size_t i = r.begin(); i != r.end(); ++elem, ++i) {
                                                           f_eval[i] = f(*elem);
                                                         }
                                                       });

                                    std::for_each(
                                        ChainableStack::thread_tapes_.begin(), ChainableStack::thread_tapes_.end(),
                                        [&](chainablestack_t& thread_stack) {
                                          if (thread_stack.id_ == parent_stack_id)
                                            return;
                                          const std::size_t thread_stack_size = thread_stack.var_stack_.size();
                                          if (thread_stack_size == 0)
                                            return;
                                          auto known = stack_starts.find(thread_stack.id_);
                                          if (known != stack_starts.end()) {
                                            if (known->second == thread_stack_size)
                                              return;
                                            std::cout << "registering remote AD tape (id = " << thread_stack.id_
                                                      << ") for block " << known->second << " to "
                                                      << thread_stack_size << std::endl;
                                            register_nested_chainablestack(
                                                thread_stack, known->second, thread_stack_size);
                                          } else {
                                            // the AD stack got created, so the starting
                                            // position is 0
                                            std::cout << "registering CREATED remote AD tape (id = "
                                                      << thread_stack.id_ << ") for block " << 0 << " to "
                                                      << thread_stack_size << std::endl;
                                            register_nested_chainablestack(thread_stack, 0,
                                                                           thread_stack_size);
                                          }
                                        });
                                  });

    return std::move(f_eval);
  }
};

}  // namespace internal
}  // namespace math
}  // namespace stan

#endif
