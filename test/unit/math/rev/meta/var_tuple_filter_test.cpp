#include <stan/math/rev/meta/var_tuple_filter.hpp>
#include <test/unit/pretty_print_types.hpp>
#include <gtest/gtest.h>

TEST(MetaTraitsRevScal, var_tuple_filter) {
  using stan::math::var_value;
  using stan::math::var;
  using stan::math::var_to_vari_filter_t;
  using stan::math::test::type_name;
  using checker = var_to_vari_filter_t<var, double, double,
   Eigen::Matrix<var, -1, -1>, std::vector<var>>;
  std::cout << "\n" << type_name<checker>() << "\n";
}
