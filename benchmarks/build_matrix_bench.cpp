#include <benchmarks/build_matrix_benchmark.hpp>
#include <stan/math/rev/fun/build_matrix.hpp>

namespace stan {
namespace math {
namespace benchmark {
namespace build_matrix_callback {

struct builder {
  template <typename Vec, typename WeightMat, typename Functor>
  inline auto operator()(int m, const Vec& x, const WeightMat& weights,
                         const Functor& functor) const {
    return stan::math::build_matrix(m, m, functor, x, weights, m);
  }
};

static void BM_BuildMatrixSimpleF(::benchmark::State& state) {
  build_matrix::run(state, builder{}, build_matrix::simple_functor{});
}

static void BM_BuildMatrixComplicatedF(::benchmark::State& state) {
  build_matrix::run(state, builder{}, build_matrix::complicated_functor{});
}

BENCHMARK(BM_BuildMatrixSimpleF)->Apply(build_matrix::apply_cases);
BENCHMARK(BM_BuildMatrixComplicatedF)->Apply(build_matrix::apply_cases);

}  // namespace build_matrix_callback
}  // namespace benchmark
}  // namespace math
}  // namespace stan

int main(int argc, char** argv) {
  return stan::math::benchmark::build_matrix::benchmark_main(argc, argv);
}
