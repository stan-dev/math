#ifdef STAN_THREADS
#include <stan/math/rev.hpp>
#include <test/unit/math/prim/functor/utils_threads.hpp>
#include <test/unit/math/prim/functor/reduce_sum_util.hpp>
#include <gtest/gtest.h>

namespace profiling_test {
stan::math::profile_map profiles_threading;
template <typename T>
struct grouped_count_lpdf {
  grouped_count_lpdf() {}

  // does the reduction in the sub-slice start to end
  template <typename VecInt1, typename VecT, typename VecInt2>
  inline T operator()(VecInt1&& sub_slice, std::size_t start, std::size_t end,
                      std::ostream* msgs, VecT&& lambda, VecInt2&& gidx) const {
    stan::math::profile<stan::math::var> p1("p1", profiles_threading);
    const std::size_t num_terms = end - start + 1;

    std::decay_t<VecT> lambda_slice(num_terms);
    for (std::size_t i = 0; i != num_terms; ++i)
      lambda_slice[i] = lambda[gidx[start + i]];

    return stan::math::poisson_lpmf(sub_slice, lambda_slice);
  }
};
}  // namespace profiling_test

TEST(Profiling, profile_threading) {
  using stan::math::var;
  using stan::math::test::get_new_msg;
  double lambda_d = 10.0;
  const std::size_t groups = 10;
  const std::size_t elems_per_group = 1000;
  const std::size_t elems = groups * elems_per_group;

  std::vector<int> data(elems);
  std::vector<int> gidx(elems);

  for (std::size_t i = 0; i != elems; ++i) {
    data[i] = i;
    gidx[i] = i / elems_per_group;
  }

  Eigen::Matrix<var, -1, 1> vlambda_v(groups);

  for (std::size_t i = 0; i != groups; ++i)
    vlambda_v[i] = i + 0.2;
  var lambda_v = vlambda_v[0];

  var poisson_lpdf
      = stan::math::reduce_sum<profiling_test::grouped_count_lpdf<var>>(
          data, 5, get_new_msg(), vlambda_v, gidx);
  poisson_lpdf.grad();
  stan::math::set_zero_all_adjoints();
  stan::math::recover_memory();
  EXPECT_GT(profiling_test::profiles_threading.size(), 0);
}
#endif
