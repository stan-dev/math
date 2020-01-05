#ifndef TEST_UNIT_MATH_PRIM_FUNCTOR_UTILS_THREADS_HPP
#define TEST_UNIT_MATH_PRIM_FUNCTOR_UTILS_THREADS_HPP

#include <stdlib.h>
#include <string>

// utility to set number of threads to use
void set_n_threads(int num_threads) {
  static char env_string[256];
  std::string num_threads_str = std::to_string(num_threads);
  snprintf(env_string, sizeof(env_string), "STAN_NUM_THREADS=%s",
           num_threads_str.c_str());
  putenv(env_string);
}

// Can't easily use std::string as putenv require non-const char*
void set_n_threads(const char* value) {
  static char env_string[256];
  snprintf(env_string, sizeof(env_string), "STAN_NUM_THREADS=%s", value);
  putenv(env_string);
}

#endif
