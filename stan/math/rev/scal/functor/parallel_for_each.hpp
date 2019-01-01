#ifndef STAN_MATH_REV_SCAL_FUNCTOR_PARALLEL_FOR_EACH_HPP
#define STAN_MATH_REV_SCAL_FUNCTOR_PARALLEL_FOR_EACH_HPP

#include <stan/math/parallel/for_each.hpp>
#include <stan/math/prim/mat/fun/typedefs.hpp>

#include <stan/math/rev/core/nest_chainablestack.hpp>

#include <boost/iterator/counting_iterator.hpp>

#include <vector>
#include <thread>

namespace stan {
namespace math {
namespace internal {

template <class InputIt, class UnaryFunction>
struct parallel_for_each_impl<InputIt, UnaryFunction, var> {
  auto operator()(InputIt first, InputIt last, UnaryFunction f) const {
    typedef decltype(f(*first)) T_return_elem;
    typedef Eigen::Matrix<var, Eigen::Dynamic, 1> T_return;
    typedef boost::counting_iterator<int> count_iter;

#ifdef STAN_THREADS
    constexpr std_par::execution::parallel_unsequenced_policy exec_policy
        = std_par::execution::par_unseq;
#else
    constexpr std_par::execution::sequenced_policy exec_policy
        = std_par::execution::seq;
#endif

    const int num_jobs = std::distance(first, last);

    typedef ChainableStack::AutodiffStackStorage chainablestack_t;

    std::vector<bool> stack_is_local(num_jobs, false);
    std::vector<chainablestack_t*> stack_used(num_jobs, nullptr);
    std::vector<std::size_t> stack_starts(num_jobs);
    std::vector<std::size_t> stack_ends(num_jobs);

    std::thread::id parent_thread = std::this_thread::get_id();

    std::vector<int> f_sizes(num_jobs);
    std::vector<T_return_elem> f_eval(num_jobs);

    std_par::for_each(
        exec_policy, count_iter(0), count_iter(num_jobs), [&](int i) -> void {
          InputIt elem = first;
          std::advance(elem, i);
          auto& elem_ref = *elem;
          stack_is_local[i] = parent_thread == std::this_thread::get_id();

          chainablestack_t* thread_stack
              = stack_is_local[i] ? nullptr : &ChainableStack::instance();
          stack_used[i] = thread_stack;

          if (!stack_is_local[i])
            stack_starts[i] = thread_stack->var_stack_.size();

          f_eval[i] = f(elem_ref);
          f_sizes[i] = num_elements(f_eval[i]);

          if (!stack_is_local[i])
            stack_ends[i] = thread_stack->var_stack_.size();
        });

    const int num_outputs = std::accumulate(f_sizes.begin(), f_sizes.end(), 0);
    T_return results(num_outputs);
    for (int i = 0, offset = 0, cur_stack_start = 0; i < num_jobs;
         offset += f_sizes[i], ++i) {
      //std::cout << "chunk i = " << i << std::endl;
      if (!stack_is_local[i]) {
        // if the current end == next start => then we can lump these
        // together
        if (i + 1 != num_jobs && stack_used[i + 1] == stack_used[i]
            && stack_starts[i + 1] == stack_ends[i]) {
          // do nothing as we merge the current and the next job results
          //std::cout << "merging block " << i << " and " << i + 1 << std::endl;
        } else {
          std::cout << "registering remote AD tape for blocks "
                    << cur_stack_start << " - " << i << std::endl;
          stan::math::register_nested_chainablestack(
              *stack_used[i], stack_starts[cur_stack_start], stack_ends[i]);
          cur_stack_start = i + 1;
        }
      } else {
        cur_stack_start = i;
      }
      results.block(offset, 0, f_sizes[i], 1).swap(f_eval[i]);
    }

    return results;
  }
};

}  // namespace internal
}  // namespace math
}  // namespace stan

#endif
