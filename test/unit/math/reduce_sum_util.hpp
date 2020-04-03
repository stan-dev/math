#ifndef TEST_UNIT_MATH_REDUCE_SUM_UTIL
#define TEST_UNIT_MATH_REDUCE_SUM_UTIL

#include <stan/math.hpp>
#include <test/unit/math/test_ad.hpp>
#include <gtest/gtest.h>
#include <algorithm>
#include <sstream>
#include <tuple>
#include <vector>

namespace stan {
namespace math {
namespace test {
std::ostream* get_new_msg() {
  std::ostream* msgs = nullptr;
  return msgs;
}
// reduce functor which is the BinaryFunction
// here we use iterators which represent integer indices
template <typename T>
struct count_lpdf {
  count_lpdf() {}

  // does the reduction in the sub-slice start to end
  inline T operator()(std::size_t start, std::size_t end,
                      const std::vector<int>& sub_slice, std::ostream* msgs,
                      const std::vector<T>& lambda,
                      const std::vector<int>& idata) const {
    return stan::math::poisson_lpmf(sub_slice, lambda[0]);
  }
};

template <typename T>
struct nesting_count_lpdf {
  nesting_count_lpdf() {}

  // does the reduction in the sub-slice start to end
  inline T operator()(std::size_t start, std::size_t end,
                      const std::vector<int>& sub_slice, std::ostream* msgs,
                      const std::vector<T>& lambda,
                      const std::vector<int>& idata) const {
    return stan::math::reduce_sum<count_lpdf<T>>(sub_slice, 5, msgs, lambda,
                                                 idata);
  }
};

struct sum_lpdf {
  template <typename T, stan::require_stan_scalar_t<T>* = nullptr>
  T sum_(T arg) const {
    return arg;
  }

  template <typename EigMat, stan::require_eigen_t<EigMat>* = nullptr>
  auto sum_(EigMat&& arg) const {
    return stan::math::sum(arg);
  }

  template <typename Vec, stan::require_std_vector_t<Vec>* = nullptr>
  auto sum_(Vec&& arg) const {
    stan::scalar_type_t<Vec> sum = 0;
    for (size_t i = 0; i < arg.size(); ++i) {
      sum += sum_(arg[i]);
    }
    return sum;
  }

  template <typename T, typename... Args>
  inline auto operator()(std::size_t start, std::size_t end, T&& sub_slice,
                         std::ostream* msgs, Args&&... args) const {
    using return_type = stan::return_type_t<T, Args...>;

    return sum_(sub_slice)
           + sub_slice.size()
                 * stan::math::sum(std::vector<return_type>{
                       return_type(sum_(std::forward<Args>(args)))...});
  }
};

struct start_end_lpdf {
  template <typename T1, typename T2>
  inline auto operator()(std::size_t start, std::size_t end, T1&&,
                         std::ostream* msgs, T2&& data) const {
    stan::return_type_t<T1, T2> sum = 0;
    EXPECT_GE(start, 0);
    EXPECT_LE(end, data.size() - 1);
    for (size_t i = start; i <= end; i++) {
      sum += data[i];
    }
    return sum;
  }
};

/**
 * slice over the grouping variable which is a var
 */
template <typename T>
struct slice_group_count_lpdf {
  slice_group_count_lpdf() {}

  // does the reduction in the sub-slice start to end
  inline T operator()(std::size_t start, std::size_t end,
                      const std::vector<T>& lambda_slice, std::ostream* msgs,
                      const std::vector<int>& y,
                      const std::vector<int>& gsidx) const {
    const std::size_t num_groups = end - start + 1;
    T result = 0.0;
    for (std::size_t i = 0; i != num_groups; ++i) {
      std::vector<int> y_group(y.begin() + gsidx[start + i],
                               y.begin() + gsidx[start + i + 1]);
      result += stan::math::poisson_lpmf(y_group, lambda_slice[i]);
    }
    return result;
  }
};

/**
 * basic performance test for a hierarchical model
 */
template <typename T>
struct grouped_count_lpdf {
  grouped_count_lpdf() {}

  // does the reduction in the sub-slice start to end
  template <typename VecInt1, typename VecT, typename VecInt2>
  inline T operator()(std::size_t start, std::size_t end, VecInt1&& sub_slice,
                      std::ostream* msgs, VecT&& lambda, VecInt2&& gidx) const {
    const std::size_t num_terms = end - start + 1;
    // std::cout << "sub-slice " << start << " - " << end << "; num_terms = " <<
    // num_terms << "; size = " << sub_slice.size() << std::endl;
    std::decay_t<VecT> lambda_slice(num_terms);
    for (std::size_t i = 0; i != num_terms; ++i)
      lambda_slice[i] = lambda[gidx[start + i]];

    return stan::math::poisson_lpmf(sub_slice, lambda_slice);
  }
};

template <typename T1, typename T2, typename... Args>
void test_slices(T1 result, T2&& vec_value, Args&&... args) {
  tbb::task_scheduler_init default_scheduler;
  using stan::math::test::get_new_msg;
  using stan::math::test::sum_lpdf;

  std::vector<std::decay_t<T2>> data(5, vec_value);
  EXPECT_EQ(result, stan::math::reduce_sum_static<sum_lpdf>(
                        data, 1, get_new_msg(), args...))
      << "Failed for reduce_sum_static";
  EXPECT_EQ(result,
            stan::math::reduce_sum<sum_lpdf>(data, 1, get_new_msg(), args...))
      << "Failed for reduce_sum";
}

auto reduce_sum_static_int_sum_lpdf = [](auto&&... args) {
  return stan::math::reduce_sum_static<sum_lpdf>(std::vector<int>(2, 10.0), 1,
                                                 get_new_msg(), args...);
};

auto reduce_sum_static_sum_lpdf = [](auto&& data, auto&&... args) {
  return stan::math::reduce_sum_static<sum_lpdf>(data, 1, get_new_msg(),
                                                 args...);
};

auto reduce_sum_int_sum_lpdf = [](auto&&... args) {
  return stan::math::reduce_sum<sum_lpdf>(std::vector<int>(2, 10.0), 1,
                                          get_new_msg(), args...);
};

auto reduce_sum_sum_lpdf = [](auto&& data, auto&&... args) {
  return stan::math::reduce_sum<sum_lpdf>(data, 1, get_new_msg(), args...);
};

template <typename T1, typename T2>
void expect_ad_reduce_sum_lpdf(T1&& data, T2&& arg) {
  using stan::math::test::reduce_sum_int_sum_lpdf;
  using stan::math::test::reduce_sum_static_int_sum_lpdf;
  using stan::math::test::reduce_sum_static_sum_lpdf;
  using stan::math::test::reduce_sum_sum_lpdf;
  stan::test::expect_ad(reduce_sum_static_int_sum_lpdf, arg);
  stan::test::expect_ad(reduce_sum_static_sum_lpdf, data, arg);
  stan::test::expect_ad(reduce_sum_int_sum_lpdf, arg);
  stan::test::expect_ad(reduce_sum_sum_lpdf, data, arg);
}

template <int grainsize>
struct static_check_lpdf {
  template <typename T>
  inline auto operator()(std::size_t start, std::size_t end,
                         const std::vector<int>&, std::ostream* msgs,
                         const std::vector<T>& data) const {
    T sum = 0;
    // std::cout << "start: " << start << ", end: " << end << ", grainsize: " <<
    // grainsize << std::endl;
    EXPECT_LE(end - start + 1, grainsize);
    for (size_t i = start; i <= end; i++) {
      sum += data[i];
    }
    return sum;
  }
};

}  // namespace test
}  // namespace math
}  // namespace stan
#endif
