#include <stan/math/rev/meta.hpp>
#include <test/unit/pretty_print_types.hpp>
#include <gtest/gtest.h>

TEST(MetaTraitsRevScal, var_tuple_filter) {
  using stan::math::var;
  using stan::math::var_to_vari_filter_t;
  using stan::math::var_value;
  using stan::math::test::type_name;

  std::cout << type_name<stan::value_type_t<Eigen::Matrix<var, -1, -1>>>() << std::endl << std::endl;
  
  std::cout << type_name<stan::get_var_vari_value_t<Eigen::Matrix<var, -1, -1>>>() << std::endl;
  
  using checker
    = var_to_vari_filter_t<var, double, double, Eigen::Matrix<var, -1, -1>, var_value<Eigen::MatrixXd>,
                             std::vector<var>>;

  std::cout << "\n" << type_name<checker>() << "\n";
}
