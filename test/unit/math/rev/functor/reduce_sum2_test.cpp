#include <stan/math/prim/core.hpp>
#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <algorithm>
#include <sstream>
#include <tuple>
#include <vector>

std::ostream* msgs = nullptr;

template <typename T>
struct count_lpdf {
  count_lpdf() {}
  template <typename Vec>
  inline auto operator()(std::size_t start, std::size_t end,
                      Vec&& sub_slice, std::ostream* msgs,
                      const T& lambda,
                      int N) const {
    using stan::math::var;

    var sum = 0.0;
    for(int j = start; j < end; j++) {
      var lambda_mult = sub_slice[j - start] * lambda;
      for(int i = 0; i < N; i++) {
        sum += lambda_mult;
        lambda_mult *= lambda;
      }
    }

    return sum;
  }
};

TEST(v3_reduce_sum_benchmarks, reduce_sum_small) {
  using stan::math::var;

  stan::math::init_threadpool_tbb();

  std::vector<int> datasizes = { 1024, 4096, 16384 };
  std::vector<size_t> grainsizes = { 8, 16, 32, 64, 128, 256 };
  std::vector<int> worksizes = { 8, 16, 32, 64, 128, 256 };

  std::cout << "which_parallel, datasize, grainsize, worksize, time" << std::endl;
  for(auto datasize : datasizes) {
      for(auto worksize : worksizes) {
        for(auto grainsize : grainsizes) {
          std::vector<int> data(datasize, 1);

          var lambda_v = 0.5;

          double begin_norm_time = omp_get_wtime();
          var poisson_lpdf2 = 0.0;
          for(int i = 0; i < 100; i++) {
            poisson_lpdf2 += count_lpdf<var>()(0, data.size(), data, msgs, lambda_v, worksize);
          }
          double end_norm_time = omp_get_wtime();
          std::cout << "normie, " << datasize << ", " << grainsize <<
            ", " << worksize << ", " << end_norm_time - begin_norm_time << std::endl;

          double begin_par_time = omp_get_wtime();
          var poisson_lpdf = 0.0;
          for(int i = 0; i < 100; i++) {
            poisson_lpdf += stan::math::reduce_sum<count_lpdf<var>>(data, grainsize, msgs, lambda_v, worksize);
          }
          double end_par_time = omp_get_wtime();
          std::cout << "reduce_sum, " << datasize << ", " << grainsize <<
            ", " << worksize << ", " << end_par_time - begin_par_time << std::endl;
          stan::math::recover_memory();
        }
    }
  }

  stan::math::recover_memory();
}
