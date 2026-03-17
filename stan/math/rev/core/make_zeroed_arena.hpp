#ifndef STAN_MATH_REV_CORE_MAKE_ZEROED_ARENA_HPP
#define STAN_MATH_REV_CORE_MAKE_ZEROED_ARENA_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/functor/filter_map.hpp>

namespace stan::math::internal {
/**
 * Creates an arena type that is the same type as the input and initialized with
 * zeros
 * @tparam Input tuple, std::vector, Eigen type, or scalar
 * @param input The input to be converted to an arena type
 * @return An arena type with the same structure as the input and initialized to
 * zeros
 */
template <typename Input>
inline constexpr auto make_zeroed_arena(Input&& input) {
  if constexpr (is_tuple_v<Input>) {
    return stan::math::filter_map<is_any_var_scalar>(
        [](auto&& output_i) { return make_zeroed_arena(output_i); }, input);
  } else if constexpr (is_std_vector_v<Input>) {
    if constexpr (!is_var_v<value_type_t<Input>>) {
      const auto output_size = input.size();
      arena_t<std::vector<decltype(make_zeroed_arena(input[0]))>> ret;
      ret.reserve(output_size);
      for (Eigen::Index i = 0; i < output_size; ++i) {
        ret.push_back(make_zeroed_arena(input[i]));
      }
      return ret;
    } else {
      return arena_t<std::vector<double>>(input.size(), 0.0);
    }
  } else if constexpr (is_eigen_v<Input>) {
    return arena_t<promote_scalar_t<double, Input>>(
        plain_type_t<promote_scalar_t<double, Input>>::Zero(input.rows(),
                                                            input.cols()));
  } else if constexpr (is_var<Input>::value) {
    return static_cast<double>(0.0);
  }
}

}  // namespace stan::math::internal

#endif
