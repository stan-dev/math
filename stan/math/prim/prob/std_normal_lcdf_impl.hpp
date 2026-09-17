#ifndef STAN_MATH_PRIM_PROB_STD_NORMAL_LCDF_IMPL_HPP
#define STAN_MATH_PRIM_PROB_STD_NORMAL_LCDF_IMPL_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/erfc.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1p.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <utility>
#include <vector>

namespace stan {
namespace math {
namespace internal {

// Cody (1969), Math. Comp. 23:631-637. Horner evaluation in r = 2/z^2
// avoids overflowing positive powers in the value and its derivatives.
template <typename T>
inline T std_normal_tail_correction(const T& r) {
  static constexpr double p[]
      = {0.000658749161529837803157, 0.0160837851487422766278,
         0.125781726111229246204,    0.360344899949804439429,
         0.305326634961232344035,    0.0163153871373020978498};
  static constexpr double q[]
      = {-0.00233520497626869185443, -0.0605183413124413191178,
         -0.527905102951428412248,   -1.87295284992346047209,
         -2.56852019228982242072,    -1.0};
  T numerator = p[5] * r + p[4];
  T denominator = q[5] * r + q[4];
  for (int i = 3; i >= 0; --i) {
    numerator = numerator * r + p[i];
    denominator = denominator * r + q[i];
  }
  return (numerator / denominator) / INV_SQRT_PI;
}

/** Scalar log Phi(z) and its slope, sharing the tail approximation or erfc.
 * The gradient is optional so primitive probability calls only pay for values.
 * The Cody crossover at |z| = sqrt(32) matches the erfc crossover in R pnorm.
 */
template <bool calc_grad, typename T, require_stan_scalar_t<T>* = nullptr>
inline std::pair<return_type_t<T>, return_type_t<T>> std_normal_lcdf_value_grad(
    const T& z_in) {
  using R = return_type_t<T>;
  const R z = z_in;
  if (z > 40) {
    return {0, 0};
  }
  if (z <= -4 * SQRT_TWO) {
    const R a = -z;
    const R log_a = log(a);
    const R inv_a = 1 / a;
    const R r = 2 * square(inv_a);
    const R correction = std_normal_tail_correction(r);
    const R value
        = -(0.5 * z) * z - HALF_LOG_TWO_PI - log_a + log1p(r * correction);
    if constexpr (calc_grad) {
      // Separate the leading a to preserve the 2/a^3 third derivative.
      return {value, a - (2 * correction * inv_a) / (1 + r * correction)};
    } else {
      return {value, 0};
    }
  }
  const R tail = 0.5 * erfc((z > 0 ? z : R(-z)) * INV_SQRT_TWO);
  const R value = z > 0 ? log1p(-tail) : log(tail);
  if constexpr (calc_grad) {
    return {value, INV_SQRT_TWO_PI * exp(-(0.5 * z) * z)
                       / (z > 0 ? R(1 - tail) : tail)};
  } else {
    return {value, 0};
  }
}

/** Elementwise values and slopes for normal-family vectorized expressions.
 * Keep primitive logarithms and exponentials in Eigen's vectorized path.
 * Extreme entries and higher-order autodiff use the scalar kernel.
 */
template <bool calc_grad, typename T, require_eigen_t<T>* = nullptr>
inline auto std_normal_lcdf_value_grad(const T& z) {
  using R = return_type_t<scalar_type_t<T>>;
  using Array = Eigen::Array<R, Eigen::Dynamic, 1>;
  Array values(z.size());
  Array slopes(calc_grad ? z.size() : 0);
  const auto scalar_evaluate = [&]() {
    for (Eigen::Index i = 0; i < z.size(); ++i) {
      const auto result = std_normal_lcdf_value_grad<calc_grad>(R(z.coeff(i)));
      values[i] = result.first;
      if constexpr (calc_grad) {
        slopes[i] = result.second;
      }
    }
  };
  if constexpr (std::is_arithmetic<R>::value) {
    // SIMD logarithms and exponentials provide no benefit for these inputs.
    if (z.size() < 2 || z.maxCoeff() <= -4 * SQRT_TWO || z.minCoeff() > 37) {
      scalar_evaluate();
      return std::make_pair(std::move(values), std::move(slopes));
    }
    std::vector<std::pair<Eigen::Index, std::pair<R, R>>> tails;
    for (Eigen::Index i = 0; i < z.size(); ++i) {
      const R zi = z.coeff(i);
      // The upper cutoff keeps subnormal probabilities out of SIMD math.
      const bool tail = zi <= -4 * SQRT_TWO || zi > 37;
      if (tail) {
        tails.emplace_back(i, std_normal_lcdf_value_grad<calc_grad>(zi));
      }
      values[i] = tail ? 0.5 : 0.5 * erfc(std::abs(zi) * INV_SQRT_TWO);
      if constexpr (calc_grad) {
        slopes[i] = tail ? 0 : zi;
      }
    }
    const Array upper = (-values).log1p();
    const Array lower = values.log();
    if constexpr (calc_grad) {
      slopes = INV_SQRT_TWO_PI * (-(0.5 * slopes) * slopes).exp();
    }
    for (Eigen::Index i = 0; i < z.size(); ++i) {
      const bool positive = z.coeff(i) > 0;
      if constexpr (calc_grad) {
        slopes[i] /= positive ? 1 - values[i] : values[i];
      }
      values[i] = positive ? upper[i] : lower[i];
    }
    for (const auto& entry : tails) {
      values[entry.first] = entry.second.first;
      if constexpr (calc_grad) {
        slopes[entry.first] = entry.second.second;
      }
    }
  } else {
    scalar_evaluate();
  }
  return std::make_pair(std::move(values), std::move(slopes));
}

}  // namespace internal
}  // namespace math
}  // namespace stan
#endif
