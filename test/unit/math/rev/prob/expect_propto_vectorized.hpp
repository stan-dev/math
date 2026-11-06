#ifndef TEST_UNIT_MATH_REV_PROB_EXPECT_PROPTO_VECTORIZED_HPP
#define TEST_UNIT_MATH_REV_PROB_EXPECT_PROPTO_VECTORIZED_HPP

#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <vector>

namespace stan {
namespace test {
namespace internal {

template <typename T, require_stan_scalar_t<T>* = nullptr>
inline const T& propto_element(const T& x, size_t) {
  return x;
}

template <typename T, require_not_stan_scalar_t<T>* = nullptr>
inline const auto& propto_element(const T& x, size_t i) {
  return x[i];
}

template <typename T, require_stan_scalar_t<T>* = nullptr>
inline void append_vars(std::vector<stan::math::var>& vars, const T& x) {
  if constexpr (is_var<T>::value) {
    vars.push_back(x);
  }
}

template <typename T, require_not_stan_scalar_t<T>* = nullptr>
inline void append_vars(std::vector<stan::math::var>& vars, const T& x) {
  for (size_t i = 0; i < static_cast<size_t>(x.size()); ++i) {
    append_vars(vars, x[i]);
  }
}

inline std::vector<double> adjoints(const std::vector<stan::math::var>& vars) {
  std::vector<double> adj;
  for (const auto& v : vars) {
    adj.push_back(v.adj());
  }
  return adj;
}

}  // namespace internal

/**
 * Expect a vectorized propto log density and its gradients to match the sum
 * of the elementwise scalar calls. Scalar arguments are shared by every scalar
 * call, so a scalar var's adjoint is summed over the elements, which checks
 * how the vectorized call broadcasts it.
 *
 * @param f Functor calling the propto log density.
 * @param args Scalars or std/Eigen vectors of double or var.
 */
template <typename F, typename... Args>
void expect_propto_vectorized(const F& f, const Args&... args) {
  std::vector<stan::math::var> vars;
  (internal::append_vars(vars, args), ...);
  ASSERT_FALSE(vars.empty()) << "propto needs an autodiff argument";
  const size_t N = stan::math::max_size(args...);

  stan::math::var lp = f(args...);
  stan::math::var lp_scalars = 0;
  for (size_t i = 0; i < N; ++i) {
    lp_scalars += f(internal::propto_element(args, i)...);
  }
  const double tol = 1e-12 * std::max(1.0, std::fabs(lp_scalars.val()));
  EXPECT_NEAR(lp_scalars.val(), lp.val(), tol);

  stan::math::set_zero_all_adjoints();
  lp.grad();
  const auto grad = internal::adjoints(vars);
  stan::math::set_zero_all_adjoints();
  lp_scalars.grad();
  const auto grad_scalars = internal::adjoints(vars);
  for (size_t i = 0; i < vars.size(); ++i) {
    EXPECT_NEAR(grad_scalars[i], grad[i],
                1e-12 * std::max(1.0, std::fabs(grad_scalars[i])))
        << "gradient " << i;
  }
}

}  // namespace test
}  // namespace stan
#endif
