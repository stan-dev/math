

#include <stan/math.hpp>
#include <gtest/gtest.h>

#ifdef STAN_THREADS
TEST(MathFunctions, expInt) {
  using stan::math::exp;
  EXPECT_FLOAT_EQ(std::exp(3), exp(3));
  EXPECT_FLOAT_EQ(std::exp(3.1), exp(3.1));
  EXPECT_FLOAT_EQ(std::exp(3.0), exp(3.0));
}

TEST(MathFunctions, expVecN100) {
  using stan::math::exp;
  size_t N = 100;
  std::vector<double> vec(N);
  for (size_t i = 0; i < N; ++i) {
    vec[i] = i + 1;
  }
  std::vector<double> vec_test;
  EXPECT_NO_THROW(vec_test = stan::math::exp_test(vec));
  for (size_t i = 0; i < N; ++i) {
    EXPECT_FLOAT_EQ(std::exp(i + 1), vec_test[i]);
  }
}

TEST(MathFunctions, expVecBench) {
  // std timing includes
  using std::chrono::high_resolution_clock;
  using std::chrono::duration_cast;
  using std::chrono::duration;
  using std::chrono::milliseconds;

  // stan math includes
  using stan::math::exp;
  using stan::math::init_threadpool_tbb;
  size_t N = 10000; // we're computing exp 10000 times but  scaling number of threads
  // scaling Nthreads by squares N, N^2, N^3
  std::cout << "N,nThreads,msInt,msDouble\n";
  for (int i = 1; i < 10; ++i) { 
    size_t Nthreads = 2;
    Nthreads  = std::pow(Nthreads, i);
    stan::math::init_threadpool_tbb(Nthreads);
    std::vector<double> vec(N);
    for (size_t i = 0; i < N; ++i) {
      vec[i] = i + 1;
    }
    std::vector<double> vec_test;
    
    auto t1 = high_resolution_clock::now();
    EXPECT_NO_THROW(vec_test = stan::math::exp_test(vec));
    auto t2 = high_resolution_clock::now();

	/* Getting number of milliseconds as an integer. */
	auto ms_int = duration_cast<milliseconds>(t2 - t1);

	/* Getting number of milliseconds as a double. */
	duration<double, std::milli> ms_double = t2 - t1;

	std::cout << N << ",";
	std::cout << Nthreads << ",";
	std::cout << ms_int.count() << "ms,";
	std::cout << ms_double.count() << "ms\n";
  }
}

#else
TEST(MathFunctions, expInt) {
  using stan::math::exp;
  EXPECT_FLOAT_EQ(std::exp(3), exp(3));
  EXPECT_FLOAT_EQ(std::exp(3.1), exp(3.1));
  EXPECT_FLOAT_EQ(std::exp(3.0), exp(3.0));
}
#endif
