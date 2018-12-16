#include <stan/math/parallel/for_each.hpp>
#include <boost/iterator/counting_iterator.hpp>

#include <gtest/gtest.h>

#include <string>
#include <vector>

struct for_each : public ::testing::Test {
  const int N = 11;
  std::vector<int> store;
  std::function<void(int)> fn = [&](int i) -> void { store[i] += i; };
  typedef boost::counting_iterator<int> count_iter;
  count_iter first;
  count_iter last;

  virtual void SetUp() {
    store = std::vector<int>(N, 0);

    first = count_iter(0);
    last = count_iter(N);
  }
};

// utility to set number of threads to use
void set_n_threads(std::size_t num_threads) {
  static char env_string[256];
  std::string num_threads_str = std::to_string(num_threads);
  snprintf(env_string, sizeof(env_string), "STAN_NUM_THREADS=%s",
           num_threads_str.c_str());
  putenv(env_string);
}

TEST_F(for_each, parallel_mode) {
  for (std::size_t i = 1; i < 2 * N + 1; ++i) {
    set_n_threads(i);

    std_par::for_each(std_par::execution::par, first, last, fn);

    for (std::size_t j = 0; j != N; ++j)
      EXPECT_EQ(store[j], i * j);
  }
}

TEST_F(for_each, parallel_unsequenced_mode) {
  for (std::size_t i = 1; i < 2 * N + 1; ++i) {
    set_n_threads(i);

    std_par::for_each(std_par::execution::par_unseq, first, last, fn);

    for (std::size_t j = 0; j != N; ++j)
      EXPECT_EQ(store[j], i * j);
  }
}

TEST_F(for_each, sequenced_mode) {
  std_par::for_each(std_par::execution::seq, first, last, fn);

  for (std::size_t i = 0; i != N; ++i)
    EXPECT_EQ(store[i], i);
}
