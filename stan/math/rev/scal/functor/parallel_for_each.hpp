#ifndef STAN_MATH_REV_SCAL_FUNCTOR_PARALLEL_FOR_EACH_HPP
#define STAN_MATH_REV_SCAL_FUNCTOR_PARALLEL_FOR_EACH_HPP

#include <stan/math/prim/scal/functor/parallel_for_each.hpp>
#include <stan/math/prim/mat/fun/typedefs.hpp>

#include <stan/math/rev/core/nest_chainablestack.hpp>

#include <boost/iterator/counting_iterator.hpp>

#include <tbb/task_arena.h>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
//#include <tbb/enumerable_thread_specific.h>
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
    typedef std::vector<decltype(f(*first))> T_return;
    typedef boost::counting_iterator<std::size_t> count_iter;

    const std::size_t num_jobs = std::distance(first, last);

    std::cout
        << "Running NEW var parallel_for_each implementation (non-nestable)..."
        << std::endl;

    typedef ChainableStack::AutodiffStackStorage chainablestack_t;

    T_return f_eval(num_jobs);

    // AD tape used / start / end
    typedef std::tuple<chainablestack_t&, std::size_t, std::size_t>
        nested_chunk_t;

    tbb::this_task_arena::isolate([&]() {
      // record the chunks executed on
      // which AD tapes
      tbb::combinable<std::vector<nested_chunk_t> > child_chunks;

      // auto
      // tbb::enumerable_thread_specific<std::shared_ptr<AutodiffStackStorage>>
      // child_stacks;

      //          instance_nochain_stack_;

      const std::size_t parent_ad_tape_idx
          = tbb::this_task_arena::current_thread_index();
      // const std::size_t
      // parent_ad_tape_idx =
      // ChainableStack::instance().id_;

      chainablestack_t& parent_tape = ChainableStack::instance();

      tbb::parallel_for(
          tbb::blocked_range<std::size_t>(0, num_jobs),
          [&](const tbb::blocked_range<size_t>& r) {
            chainablestack_t& local_tape = ChainableStack::instance();
            const std::size_t start_stack_size = local_tape.var_stack_.size();
            auto elem = first;
            std::advance(elem, r.begin());
            for (std::size_t i = r.begin(); i != r.end(); ++elem, ++i) {
              f_eval[i] = f(*elem);
            }
            const std::size_t local_ad_tape_idx
                = tbb::this_task_arena::current_thread_index();
            const std::size_t end_stack_size = local_tape.var_stack_.size();
            if (parent_ad_tape_idx != local_ad_tape_idx)
              child_chunks.local().emplace_back(
                  nested_chunk_t(local_tape, start_stack_size, end_stack_size));
          });

      child_chunks.combine_each(
          [&parent_tape](const std::vector<nested_chunk_t>& chunks) {
            std::for_each(chunks.begin(), chunks.end(),
                          [&parent_tape](const nested_chunk_t& nested_chunk) {
                            chainablestack_t& other_tape
                                = std::get<0>(nested_chunk);
                            std::size_t start = std::get<1>(nested_chunk);
                            std::size_t end = std::get<2>(nested_chunk);
                            parent_tape.var_stack_.insert(
                                parent_tape.var_stack_.end(),
                                other_tape.var_stack_.begin() + start,
                                other_tape.var_stack_.begin() + end);

                            // register_nested_chainablestack(std::get<0>(nested_chunk),
                            //                               std::get<1>(nested_chunk),
                            //                               std::get<2>(nested_chunk));
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
