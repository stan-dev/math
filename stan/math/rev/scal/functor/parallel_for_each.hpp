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
#include <tuple>
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

    // AD tape used / start / end
    typedef std::tuple<chainablestack_t&, std::size_t, std::size_t> nested_chunk_t;

    tbb::this_task_arena::isolate([&]() {
                                    // record the chunks executed on
                                    // which AD tapes
                                    tbb::combinable< std::vector<nested_chunk_t> > child_chunks;
                                    
                                    const std::size_t parent_ad_tape_idx = tbb::this_task_arena::current_thread_index();

                                    tbb::parallel_for( tbb::blocked_range<std::size_t>( 0, num_jobs ),
                                                       [&](const tbb::blocked_range<size_t>& r) {
                                                         chainablestack_t& local_tape = ChainableStack::instance();
                                                         const std::size_t start_stack_size = local_tape.var_stack_.size();
                                                         auto elem = first;
                                                         std::advance(elem, r.begin());
                                                         for (std::size_t i = r.begin(); i != r.end(); ++elem, ++i) {
                                                           f_eval[i] = f(*elem);
                                                         }
                                                         const std::size_t end_stack_size = local_tape.var_stack_.size();
                                                         if(parent_ad_tape_idx != tbb::this_task_arena::current_thread_index())
                                                           child_chunks.local().emplace_back(nested_chunk_t(local_tape, start_stack_size, end_stack_size));
                                                       });

                                    child_chunks.combine_each([](const std::vector<nested_chunk_t>& chunks) {
                                                                std::for_each(chunks.begin(), chunks.end(),
                                                                              [](const nested_chunk_t& nested_chunk) {
                                                                                register_nested_chainablestack(std::get<0>(nested_chunk),
                                                                                                               std::get<1>(nested_chunk),
                                                                                                               std::get<2>(nested_chunk));
                                                                              });
                                                              });
                                  });

    return std::move(f_eval);
  }
};

}  // namespace internal
}  // namespace math
}  // namespace stan

#endif
