#ifndef STAN_MATH_PRIM_META_MODIFY_EIGEN_OPTIONS_HPP
#define STAN_MATH_PRIM_META_MODIFY_EIGEN_OPTIONS_HPP

#include <stan/math/prim/meta/is_eigen_matrix_base.hpp>
#include <stan/math/prim/meta/plain_type.hpp>

namespace stan {
namespace math {
namespace internal {
/**
 * Change the options of an Eigen matrix or array.
 * @tparam Mat type of the matrix or array
 * @tparam NewOptions new options for the matrix or array
 */
template <typename Mat, int NewOptions, typename = void>
struct change_eigen_options_impl {};

template <typename Mat, int NewOptions>
struct change_eigen_options_impl<Mat, NewOptions,
                                 require_eigen_matrix_base_t<Mat>> {
  using type
      = Eigen::Matrix<typename Mat::Scalar, Mat::RowsAtCompileTime,
                      Mat::ColsAtCompileTime, NewOptions,
                      Mat::MaxRowsAtCompileTime, Mat::MaxColsAtCompileTime>;
};

template <typename Mat, int NewOptions>
struct change_eigen_options_impl<Mat, NewOptions, require_eigen_array_t<Mat>> {
  using type
      = Eigen::Array<typename Mat::Scalar, Mat::RowsAtCompileTime,
                     Mat::ColsAtCompileTime, NewOptions,
                     Mat::MaxRowsAtCompileTime, Mat::MaxColsAtCompileTime>;
};
}  // namespace internal
/**
 * Change the options of an Eigen matrix or array.
 * @tparam Mat type of the matrix or array
 * @tparam NewOptions new options for the matrix or array
 */
template <typename Mat, int NewOptions>
using change_eigen_options_t = typename internal::change_eigen_options_impl<
    plain_type_t<std::decay_t<Mat>>, NewOptions>::type;

}  // namespace math
}  // namespace stan

#endif
