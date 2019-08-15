#ifndef STAN_MATH_PRIM_MAT_FUN_QUAD_FORM_DIAG_HPP
#define STAN_MATH_PRIM_MAT_FUN_QUAD_FORM_DIAG_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/mat/err/check_square.hpp>
#include <stan/math/prim/mat/err/check_vector.hpp>
#include <stan/math/prim/scal/err/check_size_match.hpp>

namespace stan {
namespace math {

template <typename T1, typename T2, typename = enable_if_all_eigen<T1, T2>>
inline auto quad_form_diag(const T1& mat, const T2& vec) {
  check_vector("quad_form_diag", "vec", vec);
  check_square("quad_form_diag", "mat", mat);
  check_size_match("quad_form_diag", "rows of mat", mat.rows(), "size of vec",
                   vec.size());
  return vec.asDiagonal() * mat * vec.asDiagonal();
}

}  // namespace math
}  // namespace stan
#endif
