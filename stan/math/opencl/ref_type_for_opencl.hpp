#ifndef STAN_MATH_OPENCL_REF_TYPE_FOR_OPENCL_HPP
#define STAN_MATH_OPENCL_REF_TYPE_FOR_OPENCL_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/pinned_matrix.hpp>
#include <stan/math/prim/meta/value_type.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {
/**
 * Determines appropriate type for assigning expression of given type to, so
 * that the resulting type has directly accessible contiguous colum-major data,
 * which is needed to copy to OpenCL device for construction of matrix_cl.
 *
 * THis is similar to `ref_type` except this also copies expressions, which do
 * not have contiguous data (in both dimensions) or have row major storage
 * order.
 *
 * Warning: if a variable of this type could be assigned a rvalue, make sure
 * template parameter `T` is of correct reference type (rvalue).
 * @tparam T type to determine reference for
 */
template <typename T, typename = void>
struct ref_type_for_opencl {
  using T_val = std::remove_reference_t<T>;
  using T_plain_col_major = std::conditional_t<
      std::is_same<typename Eigen::internal::traits<T_val>::XprKind,
                   Eigen::MatrixXpr>::value,
      Eigen::Matrix<value_type_t<T>, T_val::RowsAtCompileTime,
                    T_val::ColsAtCompileTime>,
      Eigen::Array<value_type_t<T>, T_val::RowsAtCompileTime,
                   T_val::ColsAtCompileTime>>;
  using T_pinned = math::pinned_matrix<T_plain_col_major>;
  using T_optionally_ref
      = std::conditional_t<std::is_rvalue_reference<T>::value, T_val, const T&>;
  using T_val_derived = std::decay_t<decltype(std::declval<T_val>().derived())>;
  // Setting Outer stride of Ref to 0 (default) won't actually check that
  // expression has contiguous outer stride. Instead we need to check that
  // evaluator flags contain LinearAccessBit and PacketAccessBit.
  using type = std::conditional_t<
      Eigen::internal::traits<Eigen::Ref<std::decay_t<T_plain_col_major>>>::
              template match<T_val_derived>::MatchAtCompileTime
          && (Eigen::internal::evaluator<T_val_derived>::Flags
              & Eigen::LinearAccessBit)
          && (Eigen::internal::evaluator<T_val_derived>::Flags
              & Eigen::PacketAccessBit),
      T_optionally_ref, T_pinned>;
};

template <typename T>
struct ref_type_for_opencl<T, require_not_eigen_t<T>> {
  using type = std::conditional_t<std::is_rvalue_reference<T>::value,
                                  std::remove_reference_t<T>, const T&>;
};

template <typename T>
struct ref_type_for_opencl<T, require_arena_matrix_t<T>> {
  using type =
      typename ref_type_for_opencl<typename std::decay_t<T>::Base>::type;
};

template <typename T>
using ref_type_for_opencl_t = typename ref_type_for_opencl<T>::type;

}  // namespace math
}  // namespace stan

#endif
#endif
