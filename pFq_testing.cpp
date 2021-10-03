#include <benchmark/benchmark.h>
#include <stan/math.hpp>

static void grad_pFq_bench(benchmark::State& state) {
  using namespace stan::math;
  vector_v a_v(2);
  a_v << 1, 1;
  vector_v b_v(1);
  b_v << 1;
  var z_v = 0.6;

  for (auto _ : state) {
    auto grad_tuple = grad_pFq(a_v, b_v, z_v);
  }
}
BENCHMARK(grad_pFq_bench);

static void grad_2F1_bench(benchmark::State& state) {
  using namespace stan::math;
  vector_v a_v(2);
  a_v << 1, 1;
  vector_v b_v(1);
  b_v << 1;
  var z_v = 0.6;

  double g_a1;
  double g_b1;

  for (auto _ : state) {
    grad_2F1(g_a1, g_b1, a_v[0].val(), a_v[1].val(), b_v[0].val(), z_v.val());
  }
}
BENCHMARK(grad_2F1_bench);

BENCHMARK_MAIN();