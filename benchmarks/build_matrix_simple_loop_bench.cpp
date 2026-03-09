#include <benchmarks/build_matrix_benchmark.hpp>

namespace stan {
namespace math {
namespace benchmark {
namespace build_matrix_simple_loop {

struct builder {
  template <typename Vec, typename WeightMat, typename Functor>
  inline auto operator()(int m, const Vec& x, const WeightMat& weights,
                         const Functor& functor) const {
    Eigen::Matrix<stan::math::var, -1, -1> out(m, m);
    for (int i = 0; i < m; ++i) {
      for (int j = 0; j < m; ++j) {
        out(i, j) = functor(i, j, x, weights, m);
      }
    }
    return out;
  }
};

static void BM_SimpleLoopSimpleF(::benchmark::State& state) {
  build_matrix::run(state, builder{}, build_matrix::simple_functor{});
}

static void BM_SimpleLoopComplicatedF(::benchmark::State& state) {
  build_matrix::run(state, builder{}, build_matrix::complicated_functor{});
}

BENCHMARK(BM_SimpleLoopSimpleF)->Apply(build_matrix::apply_cases);
BENCHMARK(BM_SimpleLoopComplicatedF)->Apply(build_matrix::apply_cases);

}  // namespace build_matrix_simple_loop
}  // namespace benchmark
}  // namespace math
}  // namespace stan

int main(int argc, char** argv) {
  return stan::math::benchmark::build_matrix::benchmark_main(argc, argv);
}
