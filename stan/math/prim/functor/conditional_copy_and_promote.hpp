#ifndef STAN_MATH_PRIM_FUNCTOR_CONDITIONAL_COPY_AND_PROMOTE_HPP
#define STAN_MATH_PRIM_FUNCTOR_CONDITIONAL_COPY_AND_PROMOTE_HPP

#include <stan/math/prim/functor/apply.hpp>
#include <stan/math/prim/functor/map_if.hpp>
#include <stan/math/prim/functor/make_holder_tuple.hpp>
#include <stan/math/prim/fun/eval.hpp>
#include <stan/math/prim/fun/promote_scalar.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <stan/math/prim/meta/is_tuple.hpp>
#include <stan/math/prim/meta/is_var.hpp>
#include <stan/math/prim/meta/is_vector.hpp>
#include <stan/math/prim/meta/scalar_type.hpp>
#include <cstdint>
#include <utility>
#include <vector>

namespace stan::math::internal {

/**
 * Decide if object should be deep or shallow copied when
 * using @ref conditional_copy_and_promote .
 */
enum class COPY_TYPE : uint8_t { SHALLOW = 0, DEEP = 1 };

/**
 * Conditional copy and promote a type's scalar type to a `PromotedType`.
 * @tparam Filter type trait with a static constexpr bool member `value`
 *  that is true if the type should be promoted. Otherwise, the type is
 *  left unchanged.
 * @tparam PromotedType type to promote the scalar to.
 * @tparam CopyType type of copy to perform.
 * @tparam Args variadic arguments.
 * @param args variadic arguments to conditionally copy and promote.
 * @return a tuple where each element is either a reference to the original
 * argument or a promoted copy of the argument.
 */
template <template <typename...> class Filter, typename PromotedType = double,
          COPY_TYPE CopyType = COPY_TYPE::DEEP, typename... Args>
inline auto conditional_copy_and_promote(Args&&... args) {
  return map_if<Filter>(
      [](auto&& arg) {
        if constexpr (is_tuple_v<decltype(arg)>) {
          return stan::math::apply(
              [](auto&&... inner_args) {
                return make_holder_tuple(
                    conditional_copy_and_promote<Filter, PromotedType,
                                                 CopyType>(
                        std::forward<decltype(inner_args)>(inner_args))...);
              },
              std::forward<decltype(arg)>(arg));
        } else if constexpr (is_std_vector_v<decltype(arg)>) {
          std::vector<decltype(conditional_copy_and_promote<
                               Filter, PromotedType, CopyType>(arg[0]))>
              ret;
          for (std::size_t i = 0; i < arg.size(); ++i) {
            ret.push_back(
                conditional_copy_and_promote<Filter, PromotedType, CopyType>(
                    arg[i]));
          }
          return ret;
        } else {
          if constexpr (CopyType == COPY_TYPE::DEEP) {
            return stan::math::eval(promote_scalar<PromotedType>(
                value_of_rec(std::forward<decltype(arg)>(arg))));
          } else if (CopyType == COPY_TYPE::SHALLOW) {
            if constexpr (std::is_same_v<PromotedType,
                                         scalar_type_t<decltype(arg)>>) {
              return std::forward<decltype(arg)>(arg);
            } else {
              return stan::math::eval(promote_scalar<PromotedType>(
                  std::forward<decltype(arg)>(arg)));
            }
          }
        }
      },
      std::forward<Args>(args)...);
}

/**
 * Conditional deep copy types with a `var` scalar type to `PromotedType`.
 * @tparam PromotedType type to promote the scalar to.
 * @tparam Args variadic arguments.
 * @param args variadic arguments to conditionally copy and promote.
 * @return a tuple where each element is either a reference to the original
 * argument or a promoted copy of the argument.
 */
template <typename PromotedType, typename... Args>
inline auto deep_copy_vargs(Args&&... args) {
  return conditional_copy_and_promote<is_any_var_scalar, PromotedType,
                                      COPY_TYPE::DEEP>(
      std::forward<Args>(args)...);
}

/**
 * Conditional shallow copy types with a `var` scalar type to `PromotedType`.
 * @note This function is useful whenever you are inside of nested autodiff
 *  and want to allow the input arguments from an outer autodiff to be used
 *  in an inner autodiff without making a hard copy of the input arguments.
 * @tparam PromotedType type to promote the scalar to.
 * @tparam Args variadic arguments.
 * @param args variadic arguments to conditionally copy and promote.
 * @return a tuple where each element is either a reference to the original
 * argument or a promoted copy of the argument.
 */
template <typename PromotedType, typename... Args>
inline auto shallow_copy_vargs(Args&&... args) {
  return conditional_copy_and_promote<is_any_var_scalar, PromotedType,
                                      COPY_TYPE::SHALLOW>(
      std::forward<Args>(args)...);
}

}  // namespace stan::math::internal

#endif
