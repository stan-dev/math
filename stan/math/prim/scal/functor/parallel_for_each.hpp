#ifndef STAN_MATH_PRIM_SCAL_FUNCTOR_PARALLEL_FOR_EACH_HPP
#define STAN_MATH_PRIM_SCAL_FUNCTOR_PARALLEL_FOR_EACH_HPP

#include <stan/math/prim/scal/meta/return_type.hpp>
#include <stan/math/prim/mat/fun/typedefs.hpp>

#include <boost/iterator/counting_iterator.hpp>

#include <tbb/task_arena.h>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>

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

    const std::size_t num_jobs = std::distance(first, last);

    std::cout << "Running base parallel_for_each implementation..."
              << std::endl;

    T_return f_eval(num_jobs);

    tbb::parallel_for( tbb::blocked_range<std::size_t>( 0, num_jobs ),
                       [&](const tbb::blocked_range<size_t>& r) {
                         auto elem = first;
                         std::advance(elem, r.begin());
                         for (std::size_t i = r.begin(); i != r.end(); ++elem, ++i) {
                           f_eval[i] = f(*elem);
                         }
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
