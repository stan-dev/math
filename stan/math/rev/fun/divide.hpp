#ifndef STAN_MATH_REV_FUN_DIVIDE_HPP
#define STAN_MATH_REV_FUN_DIVIDE_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/to_var.hpp>
#include <stan/math/rev/fun/typedefs.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Return matrix divided by scalar.
 *
 * @tparam Mat type of the matrix or expression
 * @tparam Scal type of the scalar
 * @param[in] m specified matrix or expression
 * @param[in] c specified scalar
 * @return matrix divided by the scalar
 */
template <typename Mat, typename Scal, require_eigen_t<Mat>* = nullptr,
          require_stan_scalar_t<Scal>* = nullptr,
          require_any_st_var<Mat, Scal>* = nullptr>
inline Eigen::Matrix<var, Mat::RowsAtCompileTime, Mat::ColsAtCompileTime>
divide(const Mat& m, const Scal& c) {
  const auto& m_ref = to_ref(m);
  const auto& c_ref = to_ref(c);

  auto arena_m = to_arena_if<!is_constant<Mat>::value>(m_ref);

  double invc = 1.0 / value_of(c_ref);

  using Mat_d = promote_scalar_t<double, Mat>;
  using Mat_v = promote_scalar_t<var, Mat>;

  arena_t<Mat_d> res_val = invc * value_of(m_ref);
  arena_t<Mat_v> res = res_val;

  reverse_pass_callback([arena_m, c, res, res_val, invc]() mutable {
    const auto& adj = to_ref(res.adj());

    if (!is_constant<Mat>::value) {
      forward_as<Mat_v>(arena_m).adj() += invc * adj;
    }
    if (!is_constant<Scal>::value) {
      forward_as<var>(c).adj() += -invc * (adj.array() * res_val.array()).sum();
    }
  });

  return res;
}

}  // namespace math
}  // namespace stan
#endif
