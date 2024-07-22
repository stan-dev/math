#ifndef STAN_MATH_REV_META_MODIFY_EIGEN_OPTIONS_HPP
#define STAN_MATH_REV_META_MODIFY_EIGEN_OPTIONS_HPP

#include <stan/math/rev/core/var.hpp>
#include <stan/math/prim/meta/modify_eigen_options.hpp>

namespace stan {
namespace math {
namespace internal {

template <typename Mat, int NewOptions>
struct change_eigen_options_impl<var_value<Mat>, NewOptions,
                                 require_eigen_matrix_base_t<Mat>> {
  using type = var_value<Eigen::Matrix<
      typename Mat::Scalar, Mat::RowsAtCompileTime, Mat::ColsAtCompileTime,
      NewOptions, Mat::MaxRowsAtCompileTime, Mat::MaxColsAtCompileTime>>;
};

template <typename Mat, int NewOptions>
struct change_eigen_options_impl<var_value<Mat>, NewOptions,
                                 require_eigen_array_t<Mat>> {
  using type = var_value<Eigen::Array<
      typename Mat::Scalar, Mat::RowsAtCompileTime, Mat::ColsAtCompileTime,
      NewOptions, Mat::MaxRowsAtCompileTime, Mat::MaxColsAtCompileTime>>;
};

}  // namespace internal
}  // namespace math
}  // namespace stan

#endif
