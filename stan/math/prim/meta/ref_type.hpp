#ifndef STAN_MATH_PRIM_META_REF_TYPE_HPP
#define STAN_MATH_PRIM_META_REF_TYPE_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/is_eigen.hpp>
#include <stan/math/prim/meta/is_arena_matrix.hpp>
#include <stan/math/prim/meta/is_vector.hpp>
#include <stan/math/prim/meta/plain_type.hpp>
#include <type_traits>

namespace stan {

/** /ingroup type_trait
 * If the condition is true determines appropriate type for assigning expression
 * of given type to, to evaluate expensive expressions, but not make a copy if T
 * involves no calculations. This works similarly as [`Eigen::Ref`](https://eigen.tuxfamily.org/dox/classEigen_1_1Ref.html).
 * The deduction scheme used in `ref_type_if` also handles rvalue references,
 * so it can be used with perfect forwarding. If the condition is false the
 * expression will never be evaluated.
 *
 * Internally, this type trait compares `T` and `Eigen::Ref<T>` to see whether
 * the defines all of the following  trait flags as true for each type.
 * 1. HasDirectAccess
 *  - The underlying array of coefficients can be directly accessed as a plain strided array
 *    (This is true for objects with ownership of memory but not for expressions).
 * 2. StorageOrderMatch
 *  - Both are either RowMajor or ColumnMajor
 * 3. InnerStrideMatch && OuterStrideMatch
 *  - The inner and outer strides are exactly equal
 * 4. AlignmentMatch
 *  - The memory alignment of both are the same
 * 5.  ScalarTypeMatch
 *  - The underlying scalar types are the same
 *
 * If any of the above conditionals fail then the type returned is the same as
 * from `plain_type`. If all of the above conditions pass then this type trait's
 * `type` holds a `const T&` for lvalue inputs and `T` for rvalue references.
 *
 * Warning: if a variable of this type could be assigned a rvalue, make sure
 * template parameter `T` is of correct reference type (rvalue).
 * @tparam Condition If false, `type` will hold `T`.
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

/** /ingroup type_trait
 * See the docs for `ret_type_if`
 */
template <typename T>
using ref_type_t = typename ref_type_if<true, T>::type;

template <bool Condition, typename T>
using ref_type_if_t = typename ref_type_if<Condition, T>::type;

}  // namespace stan

#endif
