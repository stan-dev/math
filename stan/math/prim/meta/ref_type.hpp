#ifndef STAN_MATH_PRIM_META_REF_TYPE_HPP
#define STAN_MATH_PRIM_META_REF_TYPE_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/is_constant.hpp>
#include <stan/math/prim/meta/is_eigen.hpp>
#include <stan/math/prim/meta/is_arena_matrix.hpp>
#include <stan/math/prim/meta/is_vector.hpp>
#include <stan/math/prim/meta/plain_type.hpp>
#include <type_traits>

namespace stan {

/**
 * If the condition is true determines appropriate type for assigning expression
 * of given type to, to evaluate expensive expressions, but not make a copy if T
 * involves no calculations. This works similarly as `Eigen::Ref`. It also
 * handles rvalue references, so it can be used with perfect forwarding. If the
 * condition is false the expression will never be evaluated.
 *
 * Warning: if a variable of this type could be assigned a rvalue, make sure
 * template parameter `T` is of correct reference type (rvalue).
 * @tparam T type to determine reference for
 */
template <bool Condition, typename T, typename = void>
struct ref_type_if {
  using type = std::conditional_t<std::is_rvalue_reference<T>::value,
                                  std::remove_reference_t<T>, const T&>;
};

template <bool Condition, typename T>
struct ref_type_if<
    Condition, T,
    require_all_t<is_eigen<T>, bool_constant<!is_arena_matrix<T>::value>>> {
  using T_plain = plain_type_t<T>;
  using T_optionally_ref
      = std::conditional_t<std::is_rvalue_reference<T>::value,
                           std::remove_reference_t<T>, const T&>;
  using T_dec = std::decay_t<decltype(std::declval<T>().derived())>;

  using type = std::conditional_t<
      Eigen::internal::traits<Eigen::Ref<std::decay_t<T_plain>>>::
              template match<T_dec>::MatchAtCompileTime
          || !Condition,
      T_optionally_ref, T_plain>;
};

template <bool Condition, typename T>
struct ref_type_if<Condition, T, require_arena_matrix_t<T>> {
  using type =
      typename ref_type_if<Condition, typename std::decay_t<T>::Base>::type;
};

template <typename T>
using ref_type_t = typename ref_type_if<true, T>::type;

template <bool Condition, typename T>
using ref_type_if_t = typename ref_type_if<Condition, T>::type;

template <typename T>
using ref_type_if_not_constant_t =
    typename ref_type_if<!is_constant<T>::value, T>::type;

}  // namespace stan

#endif
