#ifndef STAN_MATH_PRIM_PROB_EXPONENTIAL_QF_HPP
#define STAN_MATH_PRIM_PROB_EXPONENTIAL_QF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <vector>

namespace stan {
namespace math {

template <typename Tp, typename Tbeta,
          require_all_st_arithmetic<Tp, Tbeta>* = nullptr,
          require_any_eigen_vector_t<Tp, Tbeta>* = nullptr>
inline auto exponential_qf(const Tp& p, const Tbeta& beta) {
  return (-log1p(-as_array_or_scalar(p))
            / as_array_or_scalar(beta)).matrix();
}

template <typename Tp, typename Tbeta,
          require_all_arithmetic_t<Tp, Tbeta> * = nullptr>
inline auto exponential_qf(const Tp& p, const Tbeta& beta) {
  return -log1p(-p) / beta;
}

}  // namespace math
}  // namespace stan

#endif
