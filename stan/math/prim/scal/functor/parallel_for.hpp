#ifndef STAN_MATH_PRIM_SCAL_FUNCTOR_PARALLEL_FOR_HPP
#define STAN_MATH_PRIM_SCAL_FUNCTOR_PARALLEL_FOR_HPP

#include <stan/math/parallel/for_each.hpp>
#include <stan/math/prim/mat/fun/typedefs.hpp>

#include <stan/math/prim/mat/functor/map_rect_reduce.hpp>
#include <stan/math/prim/mat/functor/map_rect_combine.hpp>
#include <boost/iterator/counting_iterator.hpp>

#include <vector>

namespace stan {
namespace math {

template <class Function, typename T1>
Eigen::Matrix<typename return_type<T1>::type, Eigen::Dynamic, 1> parallel_for(
    const Function& f, int start, int end,
    const Eigen::Matrix<T1, Eigen::Dynamic, 1>& arg1) {
  typedef typename return_type<T1>::type T_return_elem;
  typedef Eigen::Matrix<T_return_elem, Eigen::Dynamic, 1> T_return_job;
  typedef boost::counting_iterator<int> count_iter;

#ifdef STAN_THREADS
  constexpr std_par::execution::parallel_unsequenced_policy exec_policy
      = std_par::execution::par_unseq;
#else
  constexpr std_par::execution::sequenced_policy exec_policy
      = std_par::execution::seq;
#endif

  const int num_jobs = end - start;

  std::vector<int> f_sizes(num_jobs);
  std::vector<T_return_job> f_eval(num_jobs);

  std_par::for_each(exec_policy, count_iter(0), count_iter(num_jobs),
                    [&](int i) -> void {
                      f_eval[i] = f(i, arg1);
                      f_sizes[i] = f_eval[i].rows();
                    });

  const int num_outputs = std::accumulate(f_sizes.begin(), f_sizes.end(), 0);
  T_return_job results(num_outputs);
  for (int i = 0, offset = 0; i < num_jobs; offset += f_sizes[i], ++i) {
    results.block(offset, 0, f_sizes[i], 1) = f_eval[i];
  }

  return results;
}

}  // namespace math
}  // namespace stan

#endif
