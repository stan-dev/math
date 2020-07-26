#include <stan/math/rev/functor/adj_op.hpp>
#include <test/unit/pretty_print_types.hpp>
#include <gtest/gtest.h>
#include <algorithm>
#include <sstream>
#include <tuple>
#include <vector>

namespace stan {
namespace math {
namespace test {

// Test we can do assignment to the map in the var op for matrices
template <typename T>
void test_adj_op_assign_accessors_v(T&& x) {
  for (int i = 0; i < x.size(); ++i) {
    x.map()(i) = i;
  }
  for (int i = 0; i < x.size(); ++i) {
    EXPECT_EQ(x.map()(i), static_cast<double>(i))
     << "Failed for type:" << type_name<T>();
  }
  for (int i = 0; i < x.size(); ++i) {
    x(i) = x.size() - i;
  }
  for (int i = 0; i < x.size(); ++i) {
    EXPECT_EQ(x.map()(i), static_cast<double>(x.size() - i))
     << "Failed for type:" << type_name<T>();
  }
}

// Test for correct sizes in var op for vectors
template <typename T>
void test_adj_op_sizes_v(size_t n) {
  stan::math::adj_op<T> x_v(n);
  if (stan::is_eigen_row_vector<T>::value) {
    EXPECT_EQ(x_v.rows(), 1) << "Failed for type:" << type_name<T>();
    EXPECT_EQ(x_v.cols(), n) << "Failed for type:" << type_name<T>();
    EXPECT_EQ(x_v.size(), n) << "Failed for type:" << type_name<T>();
  } else {
    EXPECT_EQ(x_v.rows(), n) << "Failed for type:" << type_name<T>();
    EXPECT_EQ(x_v.cols(), 1) << "Failed for type:" << type_name<T>();
    EXPECT_EQ(x_v.size(), n) << "Failed for type:" << type_name<T>();
  }
  test_adj_op_assign_accessors_v(x_v);
}

// Test sizes in var op
template <typename T>
void test_adj_op_sizes_v(size_t n, size_t m) {
  stan::math::adj_op<T> x_v(n, m);
  EXPECT_EQ(x_v.rows(), n);
  EXPECT_EQ(x_v.cols(), m);
  EXPECT_EQ(x_v.size(), n * m);
  test_adj_op_assign_accessors_v(x_v);
}

// Test we throw for non-var op assignment.
template <typename T>
void test_adj_op_assign_accessors_d(T&& x) {
  for (int i = 0; i < x.size(); ++i) {
    EXPECT_THROW((x.map()(i) = i), std::domain_error);
  }
  for (int i = 0; i < x.size(); ++i) {
    EXPECT_THROW((x(i) = x.size() - i), std::domain_error);
  }
}

// Test for 0 sizes in non-var op matrices
template <typename T>
void test_adj_op_sizes_d(size_t n) {
  stan::math::adj_op<T> x_v(n);
  EXPECT_EQ(x_v.rows(), 0) << "Failed for type:" << type_name<T>();
  EXPECT_EQ(x_v.cols(), 0) << "Failed for type:" << type_name<T>();
  EXPECT_EQ(x_v.size(), 0) << "Failed for type:" << type_name<T>();
  test_adj_op_assign_accessors_d(x_v);
}

// Test for 0 sizes in non-var op vectors
template <typename T>
void test_adj_op_sizes_d(size_t n, size_t m) {
  stan::math::adj_op<T> x_v(n, m);
  EXPECT_EQ(x_v.rows(), 0) << "Failed for type:" << type_name<T>();
  EXPECT_EQ(x_v.cols(), 0) << "Failed for type:" << type_name<T>();
  EXPECT_EQ(x_v.size(), 0) << "Failed for type:" << type_name<T>();
  test_adj_op_assign_accessors_d(x_v);
}

}
}
}
TEST(AgradRev, adj_op_assign_sizes) {
  using stan::math::test::test_adj_op_sizes_v;
  using stan::math::test::test_adj_op_sizes_d;
  using stan::math::adj_op;
  using stan::math::var;
  using eig_mat_v = Eigen::Matrix<var, -1, -1>;
  test_adj_op_sizes_v<eig_mat_v>(5, 10);

  using eig_mat_d = Eigen::Matrix<double, -1, -1>;
  adj_op<eig_mat_d> adj_mat_d(5, 10);
  test_adj_op_sizes_d<eig_mat_d>(5, 10);

  using eig_vec_v = Eigen::Matrix<var, -1, 1>;
  test_adj_op_sizes_v<eig_vec_v>(5);
  using eig_vec_d = Eigen::Matrix<double, -1, 1>;
  test_adj_op_sizes_d<eig_vec_d>(5);


  using eig_rowvec_v = Eigen::Matrix<var, 1, -1>;
  test_adj_op_sizes_v<eig_rowvec_v>(5);
  using eig_rowvec_d = Eigen::Matrix<double, 1, -1>;
  test_adj_op_sizes_d<eig_rowvec_d>(5);

  using std_vec_v = std::vector<var>;
  test_adj_op_sizes_v<std_vec_v>(5);
  using std_vec_d = std::vector<double>;
  test_adj_op_sizes_d<std_vec_d>(5);

  adj_op<var> adj_v(5);
  EXPECT_EQ(adj_v.rows(), 0);
  EXPECT_EQ(adj_v.cols(), 0);
  EXPECT_EQ(adj_v.size(), 1);
  adj_v.map() = 2.0;
  EXPECT_EQ(adj_v.map(), 2.0);
  adj_op<double> adj_d(5);
  EXPECT_EQ(adj_d.rows(), 0);
  EXPECT_EQ(adj_d.cols(), 0);
  EXPECT_EQ(adj_d.size(), 0);
  EXPECT_THROW((adj_d.map() = 1), std::domain_error);

}

namespace stan {
namespace math {
namespace test {
template <typename T1, typename T2>
void test_adj_op_assign_containers(size_t n) {
  T1 x_d = T1(n);
  for (int i = 0; i < x_d.size(); ++i) {
    x_d[i] = i;
  }
  T2 x_var(n);
  for (int i = 0; i < x_d.size(); ++i) {
    x_var[i] = x_d[i];
  }
  stan::math::adj_op<T2> adj_v(x_var);
  stan::math::adj_op<T1, true> adj_d(x_d);
  for (int i = 0; i < x_d.size(); ++i) {
    EXPECT_EQ(adj_v(i), static_cast<double>(i));
    EXPECT_EQ(adj_d(i), static_cast<double>(i));
  }
}

template <typename T1, typename T2>
void test_adj_op_assign_containers(size_t n, size_t m) {
  T1 x_d = T1(n, m);
  for (int i = 0; i < x_d.size(); ++i) {
    x_d(i) = i;
  }
  T2 x_var(x_d);
  stan::math::adj_op<T2> adj_v(x_var);
  stan::math::adj_op<T1, true> adj_d(x_d);
  for (int i = 0; i < x_d.size(); ++i) {
    EXPECT_EQ(adj_v(i), static_cast<double>(i));
    EXPECT_EQ(adj_d(i), static_cast<double>(i));
  }
}
}
}
}
TEST(AgradRev, adj_op_assign_containers) {
  using stan::math::var;
  using stan::math::adj_op;
  using stan::math::test::test_adj_op_assign_containers;
  using eig_mat_v = Eigen::Matrix<var, -1, -1>;
  using eig_mat_d = Eigen::Matrix<double, -1, -1>;
  test_adj_op_assign_containers<eig_mat_d, eig_mat_v>(5, 5);
  using eig_vec_v = Eigen::Matrix<var, -1, 1>;
  using eig_vec_d = Eigen::Matrix<double, -1, 1>;
  test_adj_op_assign_containers<eig_vec_d, eig_vec_v>(5);
  using eig_rowvec_v = Eigen::Matrix<var, 1, -1>;
  using eig_rowvec_d = Eigen::Matrix<double, 1, -1>;
  test_adj_op_assign_containers<eig_rowvec_d, eig_rowvec_v>(5);
  using std_vec_v = std::vector<var>;
  using std_vec_d = std::vector<double>;
  test_adj_op_assign_containers<std_vec_d, std_vec_v>(5);
}

TEST(AgradRev, adj_op_assign_scalar) {
  using stan::math::adj_op;
  using stan::math::var;
  adj_op<var> adj_vd(3.0);
  EXPECT_EQ(adj_vd.map(), 3.0);
  adj_op<var> adj_v(var(3.0));
  EXPECT_EQ(adj_v.map(), 3.0);
  adj_op<double, true> adj_d(3.0);
  EXPECT_EQ(adj_d.map(), 3.0);
}
