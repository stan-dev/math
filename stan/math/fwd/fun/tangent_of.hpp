#ifndef STAN_MATH_REV_FUN_TANGENT_OF_HPP
#define STAN_MATH_REV_FUN_TANGENT_OF_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/prim/functor/map.hpp>

namespace stan::math {

/**
 * Returns a reference to a variable's tangent.
 * @note In the case of a `std::vector`, this will return back an Eigen map of a vector of the tangents.
 * @tparam T Any object contains a scalar type of `fvar`
 * @param x a fvar
 * @return reference to `x`'s tangent
 */
template <typename T>
inline decltype(auto) tangent_of(T&& x) noexcept {
  if constexpr (!is_all_fvar_scalar_v<T>) {
    static_assert(sizeof(std::decay_t<T>*) == 0, "Non-fvar types do not have an tangent.");
  } else {
    return map([](auto&& x_) -> decltype(auto) {
      using arg_t = decltype(x_);
      using arg_decay_t = std::decay_t<decltype(x_)>;
      if constexpr (is_tuple_v<arg_decay_t>) {
        return stan::math::apply([](auto&&... args) -> decltype(auto) {
          return make_holder_tuple(tangent_of(std::forward<decltype(args)>(args))...);
        }, std::forward<arg_t>(x_));
      } else if constexpr (is_fvar_v<arg_decay_t>) {
        return std::forward<arg_t>(x_).d();
      } else if constexpr (is_eigen_v<arg_decay_t>) {
        return make_holder([](auto&& x_hold_) -> decltype(auto) {
          return std::forward<decltype(x_hold_)>(x_hold_).d();
        }, std::forward<arg_t>(x_));
      } else if constexpr (is_std_vector_v<arg_decay_t>) {
        if constexpr (is_fvar_v<value_type_t<arg_decay_t>>) {
          return make_holder([](auto&& x_hold_) -> decltype(auto) {
            constexpr bool const_x = std::is_const_v<std::remove_reference_t<decltype(x_hold_)>>;
            using base_map_t = Eigen::Matrix<value_type_t<arg_decay_t>, -1, 1>;
            using map_t = std::conditional_t<const_x, const base_map_t, base_map_t>;
            return Eigen::Map<map_t>(x_hold_.data(), x_hold_.size()).d();
          }, std::forward<arg_t>(x_));
        } else {
          if constexpr (std::is_rvalue_reference_v<arg_t&&>) {
              std::vector<decltype(tangent_of(std::move(x[0])))> tangents;
              tangents.reserve(x_.size());
              for (auto&& elem : x_) {
                tangents.push_back(tangent_of(std::move(elem)));
              }
              return tangents;
          } else {
              std::vector<decltype(tangent_of(x[0]))> tangents;
              tangents.reserve(x_.size());
              for (auto&& elem : x_) {
                tangents.push_back(tangent_of(elem));
              }
              return tangents;
          }
        }
      }
    }, std::forward<T>(x));
  }
}
}  // namespace stan::math

#endif  // TANGENT_OF_HPP
