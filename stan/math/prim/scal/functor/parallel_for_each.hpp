#ifndef STAN_MATH_PRIM_SCAL_FUNCTOR_PARALLEL_FOR_EACH_HPP
#define STAN_MATH_PRIM_SCAL_FUNCTOR_PARALLEL_FOR_EACH_HPP

#include <stan/math/parallel/for_each.hpp>
#include <stan/math/prim/mat/fun/typedefs.hpp>
#include <stan/math/prim/mat/fun/concatenate_row.hpp>
#include <boost/iterator/counting_iterator.hpp>

#include <vector>

namespace stan {
namespace math {
namespace internal {

// base definition => compile error
template <class InputIt, class UnaryFunction, class T_return_type>
struct parallel_for_each_impl {};

template <class InputIt, class UnaryFunction>
struct parallel_for_each_impl<InputIt, UnaryFunction, double> {
  auto operator()(InputIt first, InputIt last, UnaryFunction f) const {
    typedef decltype(f(*first)) T_return_elem;
    typedef Eigen::Matrix<typename return_type<T_return_elem>::type,
                          Eigen::Dynamic, 1>
        T_return;
    typedef boost::counting_iterator<int> count_iter;

#ifdef STAN_THREADS
    constexpr std_par::execution::parallel_unsequenced_policy exec_policy
        = std_par::execution::par_unseq;
#else
    constexpr std_par::execution::sequenced_policy exec_policy
        = std_par::execution::seq;
#endif

    const int num_jobs = std::distance(first, last);

    std::cout << "Running base parallel_for_each implementation..."
              << std::endl;

    std::vector<T_return_elem> f_eval(num_jobs);

    std_par::for_each(exec_policy, count_iter(0), count_iter(num_jobs),
                      [&](int i) -> void {
                        InputIt elem = first;
                        std::advance(elem, i);
                        auto& elem_ref = *elem;
                        f_eval[i] = f(elem_ref);
                      });

    return concatenate_row(f_eval);
  }
};

}  // namespace internal

template <class InputIt, class UnaryFunction>
constexpr auto parallel_for_each(InputIt first, InputIt last, UnaryFunction f) {
  typedef typename return_type<decltype(f(*first))>::type return_base_t;
  return internal::parallel_for_each_impl<InputIt, UnaryFunction,
                                          return_base_t>()(first, last, f);
}

}  // namespace math
}  // namespace stan

#endif
