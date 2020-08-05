#include <stan/math/rev/functor/adj_arg.hpp>
#include <test/unit/util.hpp>
#include <test/unit/pretty_print_types.hpp>
#include <gtest/gtest.h>
#include <algorithm>
#include <sstream>
#include <tuple>
#include <vector>

TEST(AgradRev, adj_arg_types) {
  using stan::math::adj_arg_t;
  using stan::math::var;
  using std::is_same;
  using eig_mat = Eigen::Matrix<double, -1, -1>;
  using eig_vec = Eigen::Matrix<double, -1, 1>;
  using eig_row_vec = Eigen::Matrix<double, 1, -1>;
  using std_vec = std::vector<double>;
  using map_eig_mat = Eigen::Map<eig_mat>;
  using map_eig_vec = Eigen::Map<eig_vec>;;
  using map_eig_row_vec = Eigen::Map<eig_row_vec>;;
  EXPECT_TRUE((is_same<map_eig_mat, adj_arg_t<eig_mat>>::value));
  EXPECT_TRUE((is_same<map_eig_vec, adj_arg_t<eig_vec>>::value));
  EXPECT_TRUE((is_same<map_eig_vec, adj_arg_t<std_vec>>::value));
  EXPECT_TRUE((is_same<map_eig_row_vec, adj_arg_t<eig_row_vec>>::value));
}

TEST(AgradRev, setup_adj_arg_types_fill) {
  using stan::math::adj_arg_t;
  using stan::math::var;
  using stan::math::setup_adj_arg;
  using stan::math::value_of;
  using eig_mat = Eigen::Matrix<var, -1, -1>;
  using eig_vec = Eigen::Matrix<var, -1, 1>;
  using eig_row_vec = Eigen::Matrix<var, 1, -1>;
  using std_vec = std::vector<var>;
  using arg_eig_mat = adj_arg_t<eig_mat>;
  using arg_eig_vec = adj_arg_t<eig_vec>;
  using arg_std_vec = adj_arg_t<std_vec>;
  using arg_eig_row_vec = adj_arg_t<eig_row_vec>;

  var A(1);
  adj_arg_t<var> A_arg = setup_adj_arg<var>(value_of(A));
  A_arg = A.val();
  EXPECT_EQ(A_arg, A.val());

  eig_mat eig_mat_v = eig_mat::Random(10, 10);
  arg_eig_mat arg_eig_mat_d = setup_adj_arg<eig_mat>(eig_mat_v.rows(), eig_mat_v.cols());
  arg_eig_mat_d = value_of(eig_mat_v);
  EXPECT_MATRIX_EQ(arg_eig_mat_d, value_of(eig_mat_v));

  eig_vec eig_vec_v = eig_vec::Random(10);
  arg_eig_vec arg_eig_vec_d = setup_adj_arg<eig_vec>(eig_vec_v.size());
  arg_eig_vec_d = value_of(eig_vec_v);
  EXPECT_MATRIX_EQ(arg_eig_vec_d, value_of(eig_vec_v));


  std_vec std_vec_v;
  for (int i = 0; i < 10; i++) {
    std_vec_v.push_back(var(i));
  }
  arg_std_vec arg_std_vec_d = setup_adj_arg<std_vec>(std_vec_v.size());
  arg_std_vec_d = adj_arg_t<std_vec>(value_of(std_vec_v).data(), std_vec_v.size());
  for (int i = 0; i < 10; i++) {
    EXPECT_FLOAT_EQ(arg_std_vec_d(i), std_vec_v[i].val());
  }

  eig_row_vec eig_row_vec_v = eig_row_vec::Random(10);
  arg_eig_row_vec arg_eig_row_vec_d = setup_adj_arg<eig_row_vec>(eig_row_vec_v.size());
  arg_eig_row_vec_d = value_of(eig_row_vec_v);
  EXPECT_MATRIX_EQ(arg_eig_row_vec_d, value_of(eig_row_vec_v));

}

TEST(AgradRev, setup_adj_arg_types_fill_override) {
  using stan::math::adj_arg_t;
  using stan::math::var;
  using stan::math::setup_adj_arg;
  using stan::math::value_of;
  using eig_mat = Eigen::Matrix<double, -1, -1>;
  using eig_vec = Eigen::Matrix<double, -1, 1>;
  using eig_row_vec = Eigen::Matrix<double, 1, -1>;
  using std_vec = std::vector<double>;
  using arg_eig_mat = adj_arg_t<eig_mat>;
  using arg_eig_vec = adj_arg_t<eig_vec>;
  using arg_std_vec = adj_arg_t<std_vec>;
  using arg_eig_row_vec = adj_arg_t<eig_row_vec>;

  double A(1);
  adj_arg_t<double> A_arg = setup_adj_arg<double, true>(value_of(A));
  A_arg = A;
  EXPECT_EQ(A_arg, A);

  eig_mat eig_mat_v = eig_mat::Random(10, 10);
  arg_eig_mat arg_eig_mat_d = setup_adj_arg<eig_mat, true>(eig_mat_v.rows(), eig_mat_v.cols());
  arg_eig_mat_d = value_of(eig_mat_v);
  EXPECT_MATRIX_EQ(arg_eig_mat_d, value_of(eig_mat_v));

  eig_vec eig_vec_v = eig_vec::Random(10);
  arg_eig_vec arg_eig_vec_d = setup_adj_arg<eig_vec, true>(eig_vec_v.size());
  arg_eig_vec_d = value_of(eig_vec_v);
  EXPECT_MATRIX_EQ(arg_eig_vec_d, value_of(eig_vec_v));

  std_vec std_vec_v;
  for (int i = 0; i < 10; i++) {
    std_vec_v.push_back(i);
  }
  arg_std_vec arg_std_vec_d = setup_adj_arg<std_vec, true>(std_vec_v.size());
  arg_std_vec_d = adj_arg_t<std_vec>(value_of(std_vec_v).data(), std_vec_v.size());
  for (int i = 0; i < 10; i++) {
    EXPECT_FLOAT_EQ(arg_std_vec_d(i), std_vec_v[i]);
  }


  eig_row_vec eig_row_vec_v = eig_row_vec::Random(10);
  arg_eig_row_vec arg_eig_row_vec_d = setup_adj_arg<eig_row_vec, true>(eig_row_vec_v.size());
  arg_eig_row_vec_d = value_of(eig_row_vec_v);
  EXPECT_MATRIX_EQ(arg_eig_row_vec_d, value_of(eig_row_vec_v));

}

TEST(AgradRev, setup_adj_arg_types_empty) {
  using stan::math::adj_arg_t;
  using stan::math::var;
  using stan::math::setup_adj_arg;
  using stan::math::value_of;
  using eig_mat = Eigen::Matrix<var, -1, -1>;
  using eig_vec = Eigen::Matrix<var, -1, 1>;
  using eig_row_vec = Eigen::Matrix<var, 1, -1>;
  using std_vec = std::vector<var>;
  using arg_eig_mat = adj_arg_t<eig_mat>;
  using arg_eig_vec = adj_arg_t<eig_vec>;
  using arg_std_vec = adj_arg_t<std_vec>;
  using arg_eig_row_vec = adj_arg_t<eig_row_vec>;

  double A(1);
  adj_arg_t<double> A_arg = setup_adj_arg<double, false>(value_of(A));
  EXPECT_EQ(A_arg, 0);

  eig_mat eig_mat_v = eig_mat::Random(10, 10);
  arg_eig_mat arg_eig_mat_d = setup_adj_arg<eig_mat, false>(eig_mat_v.rows(), eig_mat_v.cols());
  EXPECT_EQ(nullptr, arg_eig_mat_d.data());

  eig_vec eig_vec_v = eig_vec::Random(10);
  arg_eig_vec arg_eig_vec_d = setup_adj_arg<eig_vec, false>(eig_vec_v.size());
  EXPECT_EQ(nullptr, arg_eig_vec_d.data());


  std_vec std_vec_v;
  for (int i = 0; i < 10; i++) {
    std_vec_v.push_back(var(i));
  }
  arg_std_vec arg_std_vec_d = setup_adj_arg<std_vec, false>(std_vec_v.size());
  EXPECT_EQ(nullptr, arg_std_vec_d.data());


  eig_row_vec eig_row_vec_v = eig_row_vec::Random(10);
  arg_eig_row_vec arg_eig_row_vec_d = setup_adj_arg<eig_row_vec, false>(eig_row_vec_v.size());
  EXPECT_EQ(nullptr, arg_eig_row_vec_d.data());

}
