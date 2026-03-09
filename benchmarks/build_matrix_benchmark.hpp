#ifndef STAN_BENCHMARKS_BUILD_MATRIX_BENCHMARK_HPP
#define STAN_BENCHMARKS_BUILD_MATRIX_BENCHMARK_HPP

#include <benchmark/benchmark.h>
#include <stan/math/rev.hpp>

#include <array>
#include <cstddef>

namespace stan {
namespace math {
namespace benchmark {
namespace build_matrix {

using var_vector = Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>;

constexpr std::array<int, 4> kInputSizes{{4, 16, 64, 256}};
constexpr std::array<int, 5> kMatrixSizes{{1, 2, 4, 8, 16}};

inline Eigen::VectorXd make_input_values(Eigen::Index p) {
  Eigen::VectorXd x(p);
  for (Eigen::Index i = 0; i < p; ++i) {
    x(i) = 0.1 + 0.01 * (i + 1);
  }
  return x;
}

inline Eigen::MatrixXd make_weights(Eigen::Index p, Eigen::Index m) {
  Eigen::MatrixXd weights(p, m * m);
  for (Eigen::Index col = 0; col < weights.cols(); ++col) {
    for (Eigen::Index row = 0; row < p; ++row) {
      weights(row, col)
          = 1.0 / static_cast<double>(row + col + 2)
            + 0.001 * static_cast<double>((row + 1) * ((col % m) + 1));
    }
  }
  return weights;
}

struct simple_functor {
  template <typename Vec, typename WeightMat>
  inline auto operator()(int i, int j, const Vec& x, const WeightMat& weights,
                         int m) const {
    using scalar_t = stan::return_type_t<Vec, WeightMat>;
    const Eigen::Index col = static_cast<Eigen::Index>(i * m + j);
    scalar_t sum = 0.0;
    for (Eigen::Index k = 0; k < x.size(); ++k) {
      sum += x.coeff(k) * weights.coeff(k, col);
    }
    return sum;
  }
};

struct complicated_functor {
  template <typename Vec, typename WeightMat>
  inline auto operator()(int i, int j, const Vec& x, const WeightMat& weights,
                         int m) const {
    using scalar_t = stan::return_type_t<Vec, WeightMat>;
    const Eigen::Index col = static_cast<Eigen::Index>(i * m + j);
    const double shift = 0.001 * static_cast<double>(i + j + 1);
    scalar_t sum = 0.0;
    for (Eigen::Index k = 0; k < x.size(); ++k) {
      const auto term = x.coeff(k) * weights.coeff(k, col) + shift;
      sum += stan::math::sin(term) + stan::math::cos(0.5 * term)
             + stan::math::log1p(stan::math::square(term));
    }
    return sum;
  }
};

inline void apply_cases(::benchmark::internal::Benchmark* bench) {
  for (const int p : kInputSizes) {
    for (const int m : kMatrixSizes) {
      bench->Args({p, m});
    }
  }
}

template <typename Builder, typename Functor>
inline void run(::benchmark::State& state, const Builder& builder,
                const Functor& functor) {
  try {
    const int p = static_cast<int>(state.range(0));
    const int m = static_cast<int>(state.range(1));
    const Eigen::VectorXd x_vals = make_input_values(p);
    const Eigen::MatrixXd weights = make_weights(p, m);

    state.counters["P"] = static_cast<double>(p);
    state.counters["M"] = static_cast<double>(m);
    state.counters["outputs"] = static_cast<double>(m * m);
    state.counters["work"] = static_cast<double>(p) * m * m;

    for (auto _ : state) {
      (void)_;
      var_vector x = x_vals;
      auto out = builder(m, x, weights, functor);
      stan::math::var total = stan::math::sum(out);
      total.grad();
      ::benchmark::DoNotOptimize(total.val());
      ::benchmark::ClobberMemory();
      stan::math::recover_memory();
    }
  } catch (const std::exception& e) {
    state.SkipWithError(e.what());
    stan::math::recover_memory();
  }
}

inline int benchmark_main(int argc, char** argv) {
  constexpr std::size_t kArenaBytes = 256 * 1024 * 1024;
  stan::math::ChainableStack::instance_->memalloc_.alloc(kArenaBytes);
  stan::math::recover_memory();

  ::benchmark::Initialize(&argc, argv);
  if (::benchmark::ReportUnrecognizedArguments(argc, argv)) {
    return 1;
  }
  ::benchmark::RunSpecifiedBenchmarks();
  return 0;
}

}  // namespace build_matrix
}  // namespace benchmark
}  // namespace math
}  // namespace stan

#endif
