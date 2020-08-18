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
template <typename Mat, typename Scal,
	  require_eigen_t<Mat>* = nullptr,
	  require_stan_scalar_t<Scal>* = nullptr,
	  require_any_st_var<Mat, Scal>* = nullptr>
inline Eigen::Matrix<var, Mat::RowsAtCompileTime, Mat::ColsAtCompileTime>
divide(const Mat& m, Scal c) {
  auto arena_m = to_arena_if<!is_constant<Mat>::value>(m);

  double invc = 1.0 / value_of(c);

  using Mat_v = Eigen::Matrix<var,
			      Mat::RowsAtCompileTime,
			      Mat::ColsAtCompileTime>;
  
  arena_matrix<Mat_v> res = invc * value_of(m);

  reverse_pass_callback([=]() mutable {
    if (!is_constant<Mat>::value) {
      accumulate_adjoints(arena_m, invc * res.adj());
    }
    if (!is_constant<Scal>::value) {
      accumulate_adjoints(c, -invc *
			  (res.adj().array() * res.val().array()).sum());
    }
  });

  return res;
}

}  // namespace math
}  // namespace stan
#endif
