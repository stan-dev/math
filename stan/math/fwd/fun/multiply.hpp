#ifndef STAN_MATH_FWD_FUN_MULTIPLY_HPP
#define STAN_MATH_FWD_FUN_MULTIPLY_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/fun/typedefs.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/fwd/fun/to_fvar.hpp>
#include <stan/math/prim/fun/dot_product.hpp>
#include <stan/math/prim/fun/multiply.hpp>

namespace stan {
namespace math {

template <typename Mat1, typename Mat2,
          require_all_eigen_vt<is_fvar, Mat1, Mat2>* = nullptr,
          require_vt_same<Mat1, Mat2>* = nullptr,
          require_not_eigen_row_and_col_t<Mat1, Mat2>* = nullptr>
inline auto multiply(Mat1&& m1, Mat2&& m2) {
  check_multiplicable("multiply", "m1", m1, "m2", m2);
  decltype(auto) m1_ref = to_ref(m1);
  decltype(auto) m2_ref = to_ref(m2);
  return to_fvar(multiply(m1_ref.val(), m2_ref.val()),
                 add(multiply(m1_ref.val(), m2_ref.d()),
                     multiply(m1_ref.d(), m2_ref.val())));
}

template <typename Mat1, typename Mat2,
          require_eigen_vt<is_fvar, Mat1>* = nullptr,
          require_eigen_vt<std::is_floating_point, Mat2>* = nullptr,
          require_not_eigen_row_and_col_t<Mat1, Mat2>* = nullptr>
inline auto multiply(Mat1&& m1, Mat2&& m2) {
  check_multiplicable("multiply", "m1", m1, "m2", m2);
  decltype(auto) m1_ref = to_ref(m1);
  decltype(auto) m2_ref = to_ref(m2);
  return to_fvar(multiply(m1_ref.val(), m2_ref), multiply(m1_ref.d(), m2_ref));
}

template <typename Mat1, typename Mat2,
          require_eigen_vt<std::is_floating_point, Mat1>* = nullptr,
          require_eigen_vt<is_fvar, Mat2>* = nullptr,
          require_not_eigen_row_and_col_t<Mat1, Mat2>* = nullptr>
inline auto multiply(Mat1&& m1, Mat2&& m2) {
  check_multiplicable("multiply", "m1", m1, "m2", m2);
  decltype(auto) m1_ref = to_ref(m1);
  decltype(auto) m2_ref = to_ref(m2);
  return to_fvar(multiply(m1_ref, m2_ref.val()), multiply(m1_ref, m2_ref.d()));
}

}  // namespace math
}  // namespace stan
#endif
