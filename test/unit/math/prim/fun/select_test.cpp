#include <stan/math/prim.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

TEST(MathFunctions, select_scalar_scalar) {
  using stan::math::select;

  EXPECT_FLOAT_EQ(2.0, select(true, 2.0, 5.5));
  EXPECT_FLOAT_EQ(5.5, select(false, 2.0, 5.5));

  double dbl_a = -5.5;
  int int_b = 6;

  int int_res = select(false, dbl_a, int_b);
  double dbl_res = select(true, dbl_a, int_b);

  EXPECT_EQ(int_res, int_b);
  EXPECT_FLOAT_EQ(dbl_res, dbl_a);
}

TEST(MathFunctions, select_std_std) {
  using stan::math::select;

  std::vector<double> std_a{1.2, 63, 1.25};
  std::vector<double> std_b{-1.5, 2.1, -47.20};

  EXPECT_STD_VECTOR_EQ(std_a, select(true, std_a, std_b));
  EXPECT_STD_VECTOR_EQ(std_b, select(false, std_a, std_b));

  std::vector<std::vector<double>> std_a_nested{std_a, std_a, std_a};
  std::vector<std::vector<double>> std_b_nested{std_b, std_b, std_b};

  EXPECT_STD_VECTOR_EQ(std_a_nested, select(true, std_a_nested, std_b_nested));
  EXPECT_STD_VECTOR_EQ(std_b_nested, select(false, std_a_nested, std_b_nested));

  std::vector<double> std_a_short{1.2, 63};
  std::vector<double> std_b_short{-1.5, 2.1};

  EXPECT_THROW(select(true, std_a, std_b_short), std::invalid_argument);
  EXPECT_THROW(select(true, std_a_short, std_b), std::invalid_argument);
}

TEST(MathFunctions, select_eigen_eigen) {
  using stan::math::select;

  Eigen::VectorXd eig_a(3);
  eig_a << 1.2, 63, 1.25;
  Eigen::VectorXd eig_b(3);
  eig_b << -1.5, 2.1, -47.20;

  EXPECT_STD_VECTOR_EQ(eig_a, select(true, eig_a, eig_b));
  EXPECT_STD_VECTOR_EQ(eig_b, select(false, eig_a, eig_b));

  std::vector<Eigen::VectorXd> eig_a_nested{eig_a, eig_a, eig_a};
  std::vector<Eigen::VectorXd> eig_b_nested{eig_b, eig_b, eig_b};

  EXPECT_STD_VECTOR_EQ(eig_a_nested, select(true, eig_a_nested, eig_b_nested));
  EXPECT_STD_VECTOR_EQ(eig_b_nested, select(false, eig_a_nested, eig_b_nested));

  Eigen::VectorXd eig_a_short(2);
  eig_a_short << 1.2, 63;
  Eigen::VectorXd eig_b_short(2);
  eig_b_short << -1.5, 2.1;

  EXPECT_THROW(select(true, eig_a, eig_b_short), std::invalid_argument);
  EXPECT_THROW(select(true, eig_a_short, eig_b), std::invalid_argument);
}

TEST(MathFunctions, select_scalar_container) {
  using stan::math::select;

  std::vector<double> std_a{1.2, 63, 1.25};
  std::vector<double> std_b{-1.5, 2.1, -47.20};

  std::vector<double> promoted_a{0.5, 0.5, 0.5};
  std::vector<double> promoted_b{10.55, 10.55, 10.55};

  EXPECT_STD_VECTOR_EQ(std_a, select(true, std_a, 0.5));
  EXPECT_STD_VECTOR_EQ(std_b, select(false, 10.55, std_b));

  // Scalar value returned as a container of the same type and size as argument
  EXPECT_STD_VECTOR_EQ(promoted_a, select(false, std_a, 0.5));
  EXPECT_STD_VECTOR_EQ(promoted_b, select(true, 10.55, std_b));

  std::vector<std::vector<double>> std_a_nested{std_a, std_a, std_a};
  std::vector<std::vector<double>> std_b_nested{std_b, std_b, std_b};

  std::vector<std::vector<double>> promoted_a_nested{promoted_a, promoted_a,
                                                     promoted_a};
  std::vector<std::vector<double>> promoted_b_nested{promoted_b, promoted_b,
                                                     promoted_b};

  EXPECT_STD_VECTOR_EQ(std_a_nested, select(true, std_a_nested, 0.5));
  EXPECT_STD_VECTOR_EQ(std_b_nested, select(false, 10.55, std_b_nested));

  EXPECT_STD_VECTOR_EQ(promoted_a_nested, select(false, std_a_nested, 0.5));
  EXPECT_STD_VECTOR_EQ(promoted_b_nested, select(true, 10.55, std_b_nested));
}

TEST(MathFunctions, select_scalar_nested_eig) {
  using stan::math::select;

  Eigen::MatrixXd eigmat_a = Eigen::MatrixXd::Random(2, 2);
  std::vector<std::vector<Eigen::MatrixXd>> nested_eigmat_a
      = {{eigmat_a}, {eigmat_a}};

  double dbl_a = -1.25;
  Eigen::MatrixXd dbl_a_promoted
      = Eigen::MatrixXd::Constant(eigmat_a.rows(), eigmat_a.cols(), dbl_a);
  std::vector<std::vector<Eigen::MatrixXd>> nested_dbl_a_promoted
      = {{dbl_a_promoted}, {dbl_a_promoted}};

  EXPECT_STD_VECTOR_EQ(nested_dbl_a_promoted,
                       select(true, dbl_a, nested_eigmat_a));
  EXPECT_STD_VECTOR_EQ(nested_dbl_a_promoted,
                       select(false, nested_eigmat_a, dbl_a));

  EXPECT_STD_VECTOR_EQ(nested_eigmat_a, select(false, dbl_a, nested_eigmat_a));
  EXPECT_STD_VECTOR_EQ(nested_eigmat_a, select(true, nested_eigmat_a, dbl_a));
}

TEST(MathFunctions, select_array_bool_scalar_scalar) {
  using stan::math::select;

  Eigen::ArrayXd val(3);
  val << 1, -1, 1;

  Eigen::ArrayXd a_array(3);
  a_array << 0.1, 0.1, 0.1;

  Eigen::ArrayXd b_array(3);
  b_array << 2.5, 2.5, 2.5;

  Eigen::ArrayXd a_b_a(3);
  a_b_a << 0.1, 2.5, 0.1;

  Eigen::ArrayXd b_a_b(3);
  b_a_b << 2.5, 0.1, 2.5;

  EXPECT_STD_VECTOR_EQ(a_array, select(val > -5, 0.1, 2.5));
  EXPECT_STD_VECTOR_EQ(b_array, select(val < -5, 0.1, 2.5));
  EXPECT_STD_VECTOR_EQ(a_b_a, select(val > 0, 0.1, 2.5));
  EXPECT_STD_VECTOR_EQ(b_a_b, select(val < 0, 0.1, 2.5));
}

TEST(MathFunctions, select_array_bool) {
  using stan::math::select;

  Eigen::ArrayXd val(3);
  val << 1, -1, 1;

  double a = 0.1;
  double b = 2.5;

  Eigen::ArrayXd a_array(3);
  a_array << a, a, a;

  Eigen::ArrayXd b_array(3);
  b_array << b, b, b;

  Eigen::ArrayXd a_b_a(3);
  a_b_a << a, b, a;

  Eigen::ArrayXd b_a_b(3);
  b_a_b << b, a, b;

  EXPECT_STD_VECTOR_EQ(a_array, select(val > -5, a_array, b_array));
  EXPECT_STD_VECTOR_EQ(a_array, select(val > -5, a, b_array));

  EXPECT_STD_VECTOR_EQ(b_array, select(val < -5, a_array, b_array));
  EXPECT_STD_VECTOR_EQ(b_array, select(val < -5, a_array, b));

  EXPECT_STD_VECTOR_EQ(a_b_a, select(val > 0, a_array, b_array));
  EXPECT_STD_VECTOR_EQ(a_b_a, select(val > 0, a, b_array));
  EXPECT_STD_VECTOR_EQ(a_b_a, select(val > 0, a_array, b));

  EXPECT_STD_VECTOR_EQ(b_a_b, select(val < 0, a_array, b_array));
  EXPECT_STD_VECTOR_EQ(b_a_b, select(val < 0, a, b_array));
  EXPECT_STD_VECTOR_EQ(b_a_b, select(val < 0, a_array, b));

  Eigen::ArrayXd val_short(2);
  val_short << 1, -1;

  Eigen::ArrayXd a_array_short(2);
  a_array_short << a, a;

  Eigen::ArrayXd b_array_short(2);
  b_array_short << b, b;

  EXPECT_THROW(select(val_short < 0, a_array, b_array), std::invalid_argument);
  EXPECT_THROW(select(val_short < 0, a, b_array), std::invalid_argument);
  EXPECT_THROW(select(val_short < 0, a_array, b), std::invalid_argument);

  EXPECT_THROW(select(val < 0, a_array_short, b_array), std::invalid_argument);
  EXPECT_THROW(select(val < 0, a_array_short, b), std::invalid_argument);

  EXPECT_THROW(select(val < 0, a_array, b_array_short), std::invalid_argument);
  EXPECT_THROW(select(val < 0, a, b_array_short), std::invalid_argument);

  EXPECT_THROW(select(val < 0, b_array_short, b_array_short),
               std::invalid_argument);
}
