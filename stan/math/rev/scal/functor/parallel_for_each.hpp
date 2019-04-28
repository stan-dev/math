#ifndef STAN_MATH_REV_SCAL_FUNCTOR_PARALLEL_FOR_EACH_HPP
#define STAN_MATH_REV_SCAL_FUNCTOR_PARALLEL_FOR_EACH_HPP

#include <stan/math/prim/scal/functor/parallel_for_each.hpp>
#include <stan/math/prim/mat/fun/typedefs.hpp>

#include <stan/math/rev/core/scoped_chainablestack.hpp>

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
    typedef tbb::enumerable_thread_specific<ScopedChainableStack>
        tls_scoped_stack_t;

    T_return f_eval(num_jobs);

    // All AD terms are written to thread-local AD tapes which are all
    // stored as part of the parent nochain stacks.
    chainablestack_t& parent_stack = ChainableStack::instance();

    tls_scoped_stack_t tls_scoped_stacks(
        [&parent_stack]() { return ScopedChainableStack(parent_stack); });

    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, num_jobs),
                      [&](const tbb::blocked_range<size_t>& r) {
                        tls_scoped_stacks.local().execute([&] {
                          auto elem = first;
                          std::advance(elem, r.begin());
                          for (std::size_t i = r.begin(); i != r.end();
                               ++elem, ++i) {
                            f_eval[i] = f(*elem);
                          }
                        });
                      });

    tls_scoped_stacks.combine_each(
        [&parent_stack](ScopedChainableStack& child_scoped_stack) {
          child_scoped_stack.append_to_parent();
        });

    return std::move(f_eval);
  }
};

}  // namespace internal
}  // namespace math
}  // namespace stan

#endif
