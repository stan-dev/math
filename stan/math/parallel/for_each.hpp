#ifndef STAN_MATH_PARALLEL_FOR_EACH_HPP
#define STAN_MATH_PARALLEL_FOR_EACH_HPP

#ifdef STAN_PSTL_CPP17

// Standard C++17 parallel facilities
#include <execution>
#include <numeric>
#include <algorithm>

#define std_par std

#elif STAN_PSTL_INTEL

// Intel Parallel STL
#include <pstl/execution>
#include <pstl/numeric>
#include <pstl/algorithm>

#define std_par std

#else

#include <stan/math/parallel/get_num_threads.hpp>

// fall-back solution which is based on C++11 features only (async) or
// the Intel TBB if available
#include <algorithm>
#include <type_traits>
#include <iterator>

#include <vector>
#include <thread>
#include <future>

#ifdef STAN_TBB
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#endif

#define std_par stan::math::internal

namespace stan {
namespace math {
namespace internal {

namespace execution {

// sequential execution policy:
struct sequenced_policy {
  typedef std::false_type parallel_mode;
  typedef std::false_type unsequenced_mode;
};

// parallel execution policy:
struct parallel_policy {
  typedef std::true_type parallel_mode;
  typedef std::false_type unsequenced_mode;
};
struct parallel_unsequenced_policy {
  typedef std::true_type parallel_mode;
  typedef std::true_type unsequenced_mode;
};

template <class T>
struct is_execution_policy : std::false_type {};
template <>
struct is_execution_policy<sequenced_policy> : std::true_type {};
template <>
struct is_execution_policy<parallel_policy> : std::true_type {};
template <>
struct is_execution_policy<parallel_unsequenced_policy> : std::true_type {};

// global execution policy objects:
static constexpr sequenced_policy seq{};
static constexpr parallel_policy par{};
static constexpr parallel_unsequenced_policy par_unseq{};

}  // namespace execution

template <class InputIt, class UnaryFunction>
UnaryFunction for_each_impl(InputIt first, InputIt last, UnaryFunction f,
                            std::true_type /* parallel mode */) {
#ifdef STAN_TBB
  // in case we have the TBB run the for_each loop using the Intel TBB
  // parallel_for_each
  tbb::parallel_for_each(first, last, f);
#else
  const int num_jobs = std::distance(first, last);
  const int num_threads = std::min(get_num_threads(), num_jobs);
  const int small_job_size = num_jobs / num_threads;
  const int num_big_jobs = num_jobs % num_threads;

  std::vector<int> job_chunk_size(num_threads, small_job_size);

  for (std::size_t i = 0; i < num_big_jobs; ++i)
    ++job_chunk_size[num_threads - i - 1];

  std::vector<std::future<void> > job_futures;

  InputIt cur_job_start = first;
  for (std::size_t i = 0; i < num_threads; ++i) {
    InputIt cur_job_end = cur_job_start;
    std::advance(cur_job_end, job_chunk_size[i]);
    job_futures.emplace_back(std::async(
        i == 0 ? std::launch::deferred : std::launch::async, [=, &f]() -> void {
          std::for_each<InputIt, UnaryFunction>(cur_job_start, cur_job_end, f);
        }));
    cur_job_start = cur_job_end;
  }

  for (std::size_t i = 0; i < num_threads; ++i)
    job_futures[i].wait();
#endif

  return std::move(f);
}

template <class InputIt, class UnaryFunction>
UnaryFunction for_each_impl(InputIt first, InputIt last, UnaryFunction f,
                            std::false_type /* parallel mode */) {
  return std::for_each(first, last, f);
}

/**
 * C++11 based implementation of the C++17 for_each which supports
 * sequenced and parallel execution policies, see
 * https://en.cppreference.com/w/cpp/algorithm/for_each for details.
 *
 * While the C++17 does not allow to control the number of threads,
 * this version creates at most STAN_NUM_THREADS concurrent threads
 * (which is defined in the environment of the running program).
 *
 * Note: The implemented parallel version should only be used with
 * random-access iterators.
 */
template <class ExecutionPolicy, class InputIt, class UnaryFunction>
UnaryFunction for_each(ExecutionPolicy& exec_policy, InputIt first,
                       InputIt last, UnaryFunction f) {
  return for_each_impl(first, last, f,
                       typename ExecutionPolicy::parallel_mode());
}

}  // namespace internal
}  // namespace math
}  // namespace stan

#endif

#endif
