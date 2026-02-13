#ifndef STAN_MATH_REV_FUN_ADJOINT_OF_HPP
#define STAN_MATH_REV_FUN_ADJOINT_OF_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/prim/functor/map.hpp>

namespace stan::math {

/**
 * Returns a reference to a variable's adjoint.
 * @note In the case of a `std::vector`, this will return back an Eigen map of a vector of the adjoints.
 * @tparam T Any object contains a scalar type of `var`
 * @param x a var
 * @return reference to `x`'s adjoint
 */
template <typename T>
inline decltype(auto) adjoint_of(T&& x) noexcept {
  if constexpr (!is_all_var_scalar_v<T>) {
    static_assert(sizeof(std::decay_t<T>*) == 0, "Non-var types do not have an adjoint.");
  } else {
    return map([](auto&& x_) -> decltype(auto) {
      using arg_t = decltype(x_);
      using arg_decay_t = std::decay_t<decltype(x_)>;
      if constexpr (is_tuple_v<arg_decay_t>) {
        return stan::math::apply([](auto&&... args) -> decltype(auto) {
          return make_holder_tuple(adjoint_of(std::forward<decltype(args)>(args))...);
        }, std::forward<arg_t>(x_));
      } else if constexpr (is_var_v<arg_decay_t>) {
        return std::forward<arg_t>(x_).adj();
      } else if constexpr (is_eigen_v<arg_decay_t>) {
        return make_holder([](auto&& x_hold_) -> decltype(auto) {
          return std::forward<decltype(x_hold_)>(x_hold_).adj();
        }, std::forward<arg_t>(x_));
      } else if constexpr (is_std_vector_v<arg_decay_t>) {
        if constexpr (is_var_v<value_type_t<arg_decay_t>>) {
          return make_holder([](auto&& x_hold_) -> decltype(auto) {
            constexpr bool const_x = std::is_const_v<std::remove_reference_t<decltype(x_hold_)>>;
            using base_map_t = Eigen::Matrix<var, -1, 1>;
            using map_t = std::conditional_t<const_x, const base_map_t, base_map_t>;
            return Eigen::Map<map_t>(x_hold_.data(), x_hold_.size()).adj();
          }, std::forward<arg_t>(x_));
        } else {
          if constexpr (std::is_rvalue_reference_v<arg_t&&>) {
              std::vector<decltype(adjoint_of(std::move(x[0])))> adjoints;
              adjoints.reserve(x_.size());
              for (auto&& elem : x_) {
                adjoints.push_back(adjoint_of(std::move(elem)));
              }
              return adjoints;
          } else {
              std::vector<decltype(adjoint_of(x[0]))> adjoints;
              adjoints.reserve(x_.size());
              for (auto&& elem : x_) {
                adjoints.push_back(adjoint_of(elem));
              }
              return adjoints;
          }
        }
      }
    }, std::forward<T>(x));
  }
}
}  // namespace stan::math

#endif  // ADJOINT_OF_HPP
