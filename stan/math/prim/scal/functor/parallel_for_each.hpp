#ifndef STAN_MATH_PRIM_SCAL_FUNCTOR_PARALLEL_FOR_EACH_HPP
#define STAN_MATH_PRIM_SCAL_FUNCTOR_PARALLEL_FOR_EACH_HPP

#include <stan/math/parallel/for_each.hpp>
#include <stan/math/prim/scal/meta/return_type.hpp>
#include <stan/math/prim/mat/fun/typedefs.hpp>
#include <boost/iterator/counting_iterator.hpp>

#include <iostream>
#include <iterator>
#include <vector>

namespace stan {
namespace math {
namespace internal {

// base definition => compile error
template <class InputIt, class UnaryFunction, class T_return_type>
struct parallel_map_impl {};

template <class InputIt, class UnaryFunction>
struct parallel_map_impl<InputIt, UnaryFunction, double> {
  auto operator()(InputIt first, InputIt last, UnaryFunction f) const {
    typedef decltype(f(*first)) T_return_elem;
    typedef std::vector<decltype(f(*first))>
        T_return;
    typedef boost::counting_iterator<std::size_t> count_iter;

#ifdef STAN_THREADS
    constexpr std_par::execution::parallel_unsequenced_policy exec_policy
        = std_par::execution::par_unseq;
#else
    constexpr std_par::execution::sequenced_policy exec_policy
        = std_par::execution::seq;
#endif

    const std::size_t num_jobs = std::distance(first, last);

    std::cout << "Running base parallel_for_each implementation..."
              << std::endl;

    T_return f_eval(num_jobs);

    std_par::for_each(exec_policy, count_iter(0), count_iter(num_jobs),
                      [&](std::size_t i) -> void {
                        auto elem = first;
                        std::advance(elem, i);
                        f_eval[i] = f(*elem);
                      });

    return std::move(f_eval);
  }
};

}  // namespace internal

template <class InputIt, class UnaryFunction>
constexpr auto parallel_map(InputIt first, InputIt last, UnaryFunction f) {
  typedef typename return_type<decltype(f(*first))>::type return_base_t;
  return internal::parallel_map_impl<InputIt, UnaryFunction,
                                     return_base_t>()(first, last, f);
}

}  // namespace math
}  // namespace stan

#endif
