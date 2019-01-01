#ifndef STAN_MATH_REV_SCAL_FUNCTOR_PARALLEL_FOR_EACH_HPP
#define STAN_MATH_REV_SCAL_FUNCTOR_PARALLEL_FOR_EACH_HPP

#include <stan/math/parallel/for_each.hpp>
#include <stan/math/prim/mat/fun/typedefs.hpp>

#include <stan/math/rev/core/nest_chainablestack.hpp>

#include <boost/iterator/counting_iterator.hpp>

#include <map>
#include <utility>
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

    std::cout
        << "Running NEW var parallel_for_each implementation (non-nestable)..."
        << std::endl;

    typedef ChainableStack::AutodiffStackStorage chainablestack_t;

    // stack_id => stack size map
    std::map<std::size_t, std::size_t> stack_starts;

    std::for_each(ChainableStack::instance_.begin(),
                  ChainableStack::instance_.end(),
                  [&](chainablestack_t& thread_stack) {
                    stack_starts.insert(std::make_pair(
                        thread_stack.id_, thread_stack.var_stack_.size()));
                  });

    const std::size_t parent_stack_id = ChainableStack::instance().id_;

    std::vector<T_return_elem> f_eval(num_jobs);

    std_par::for_each(exec_policy, count_iter(0), count_iter(num_jobs),
                      [&](int i) -> void {
                        InputIt elem = first;
                        std::advance(elem, i);
                        f_eval[i] = f(*elem);
                      });

    std::for_each(
        ChainableStack::instance_.begin(), ChainableStack::instance_.end(),
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
            stan::math::register_nested_chainablestack(
                thread_stack, known->second, thread_stack_size);
          } else {
            // the AD stack got created, so the starting
            // position is 0
            std::cout << "registering CREATED remote AD tape (id = "
                      << thread_stack.id_ << ") for block " << 0 << " to "
                      << thread_stack_size << std::endl;
            stan::math::register_nested_chainablestack(thread_stack, 0,
                                                       thread_stack_size);
          }
        });

    return concatenate_row(f_eval);
  }
};

}  // namespace internal
}  // namespace math
}  // namespace stan

#endif
