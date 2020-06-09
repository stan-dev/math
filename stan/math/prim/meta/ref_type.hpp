#ifndef STAN_MATH_PRIM_META_REF_TYPE_HPP
#define STAN_MATH_PRIM_META_REF_TYPE_HPP

#include <stan/math/prim/meta/is_eigen.hpp>
#include <stan/math/prim/meta/is_vector.hpp>
#include <stan/math/prim/meta/plain_type.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

#include <type_traits>

namespace stan {

/**
 * Determines appropriate type for assigning expression of given type to,
 * to evaluate expensive expressions, but not make a copy if T involves no
 * calculations. This works similarly as `Eigen::Ref`. It also handles
 * rvalue references, so it can be used with perfect forwarding.
 *
 * Warning: if a variable of this type could be assigned a rvalue, make sure
 * template parameter `T` is of correct reference type (rvalue).
 * @tparam T type to determine reference for
 */
template <typename T, typename = void>
struct ref_type {
  using T_plain = plain_type_t<T>;
  using T_optionally_ref
      = std::conditional_t<std::is_rvalue_reference<T>::value,
                           std::remove_reference_t<T>, const T&>;
  using type = std::conditional_t<
      Eigen::internal::traits<Eigen::Ref<std::decay_t<T_plain>>>::
          template match<std::decay_t<T>>::MatchAtCompileTime,
      T_optionally_ref, T_plain>;
};

template <typename T>
struct ref_type<T, require_not_eigen_t<T>> {
  using type
      = std::conditional_t<std::is_rvalue_reference<T>::value, T, const T&>;
};

template <typename T>
using ref_type_t = typename ref_type<T>::type;

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
  using T_optionally_ref
      = std::conditional_t<std::is_rvalue_reference<T>::value, T_val, const T&>;
  // Setting Outer stride of Ref to 0 (default) won't actually check that
  // expression has contiguous outer stride. Instead we need to check that
  // evaluator flags contain LinearAccessBit.
  using type = std::conditional_t<
      Eigen::internal::traits<Eigen::Ref<std::decay_t<T_plain_col_major>>>::
              template match<T_val>::MatchAtCompileTime
          && (Eigen::internal::evaluator<T_val>::Flags
              & Eigen::LinearAccessBit),
      T_optionally_ref, T_plain_col_major>;
};

template <typename T>
struct ref_type_for_opencl<T, require_not_eigen_t<T>> {
  using type = std::conditional_t<std::is_rvalue_reference<T>::value,
                                  std::remove_reference_t<T>, const T&>;
};

template <typename T>
using ref_type_for_opencl_t = typename ref_type_for_opencl<T>::type;

}  // namespace stan

#endif
