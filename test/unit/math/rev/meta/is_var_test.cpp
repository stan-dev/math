#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/meta/var_tuple_filter.hpp>
#include <gtest/gtest.h>

TEST(MetaTraitsRevScal, is_var) {
  using stan::is_var;
  EXPECT_TRUE(is_var<stan::math::var>::value);
}

TEST(MetaTraitsRevScal, is_var_value) {
  using stan::is_var_value;
  EXPECT_TRUE((is_var_value<stan::math::var>::value));
  EXPECT_TRUE((is_var_value<
               stan::math::var_value<Eigen::Matrix<double, -1, -1>>>::value));
}
using stan::math::var;
template <typename... Types>
using remaining_vars_
    = std::tuple_size<stan::math::var_to_vari_filter_t<Types...>>;

using x_vis_size_ = std::tuple_size<stan::math::var_to_vari_filter_t<
    double, var, double, var, double, Eigen::Matrix<var, -1, -1>, var, double>>;

template <typename... Types>
using var_position_ = std::integral_constant<
    size_t,
    x_vis_size_::value - remaining_vars_<std::decay_t<Types>...>::value>;

void blahh() {}

template <typename T, typename... Pargs>
void blahh(T&& x, Pargs&&... args) {
  std::cout << "\n" << var_position_<T, Pargs...>::value << "\n";
  blahh(args...);
}

TEST(MetaTraitsRevScal, is_var_valueeee) {
  std::cout << "\n" << x_vis_size_::value << "\n";
  blahh(2.0, var(1), 2, var(1), 1.0, Eigen::Matrix<var, -1, -1>(1, 1), var(1),
        1);
}
