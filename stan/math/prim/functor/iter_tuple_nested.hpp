#ifndef STAN_MATH_PRIM_FUNCTOR_ITER_TUPLE_N_HPP
#define STAN_MATH_PRIM_FUNCTOR_ITER_TUPLE_N_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/functor/for_each.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <utility>
namespace stan::math {
/**
 * Iterate and nest into a tuple or std::vector to apply `f` to each matrix or
 * scalar type.
 * @tparam F a functor with `operator()(Arg&&)`
 * @tparam Types types of arguments to `f`
 * @param f functor to apply
 * @param args arguments to apply `f` to. If `args` is a tuple or std::vector,
 * this function will nest until it finds an Eigen type or scalar then apply `f`
 * to that value.
 * @return void, all arguments are passed by reference and this function will
 * only create side effects.
 */
template <typename F, typename... Types>
inline void iter_tuple_nested(F&& f, Types&&... args) {
  constexpr bool is_vec_container
      = (is_std_vector_v<Types> && ...)
        && (!is_stan_scalar<value_type_t<Types>>::value && ...);
  if constexpr ((is_tuple_v<Types> && ...)) {
    stan::math::for_each(
        [&f](auto&&... args_i) {
          iter_tuple_nested(f, std::forward<decltype(args_i)>(args_i)...);
        },
        std::forward<Types>(args)...);
  } else if constexpr (is_vec_container) {
    const auto vec_size = max_size(args...);
    for (Eigen::Index i = 0; i < vec_size; ++i) {
      iter_tuple_nested(f, args[i]...);
    }
  } else {
    f(std::forward<Types>(args)...);
  }
}

}  // namespace stan::math
#endif
