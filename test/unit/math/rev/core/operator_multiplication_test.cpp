#include <stan/math/rev/core.hpp>
#include <test/unit/math/rev/fun/util.hpp>
#include <gtest/gtest.h>
#include <vector>
#include <limits>

TEST(MathRev, multiplication_scalar_scalar) {
  using stan::math::var;
  var x = 2.0;
  var y = 3.5;
  var lp = x * y;
  lp.grad();

  EXPECT_FLOAT_EQ(lp.val(), 7.0);
  EXPECT_FLOAT_EQ(x.adj(), 3.5);
  EXPECT_FLOAT_EQ(y.adj(), 2.0);
  stan::math::recover_memory();
}

TEST(MathRev, multiplication_scalar_vector) {
  using stan::math::sum;
  using stan::math::var;
  using stan::math::var_value;
  int N = 3;
  var x = 2.0;
  Eigen::VectorXd y_val(N);
  for(size_t i = 0; i < y_val.size(); ++i)
    y_val(i) = i + 1;
  var_value<Eigen::VectorXd> y = y_val;
  var lp = sum(x * y);
  lp.grad();

  EXPECT_FLOAT_EQ(lp.val(), 12.0);
  EXPECT_FLOAT_EQ(x.adj(), 6.0);
  EXPECT_FLOAT_EQ(y.adj()(0), 2.0);
  EXPECT_FLOAT_EQ(y.adj()(1), 2.0);
  EXPECT_FLOAT_EQ(y.adj()(2), 2.0);
  stan::math::set_zero_all_adjoints();
  var lp2 = sum(y * x);
  lp2.grad();
  EXPECT_FLOAT_EQ(lp.val(), 12.0);
  EXPECT_FLOAT_EQ(x.adj(), 6.0);
  EXPECT_FLOAT_EQ(y.adj()(0), 2.0);
  EXPECT_FLOAT_EQ(y.adj()(1), 2.0);
  EXPECT_FLOAT_EQ(y.adj()(2), 2.0);
  stan::math::recover_memory();
}

/*TEST(MathRev, TestVarEigenVectors) {
  using stan::math::sum;
  using stan::math::var;
  using stan::math::var_value;
  Eigen::Matrix<double, -1, 1> x_vals(9);
  Eigen::Matrix<double, 1, -1> y_vals(9);
  x_vals << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  y_vals << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  var_value<Eigen::Matrix<double, -1, 1>> x = x_vals;
  var_value<Eigen::Matrix<double, 1, -1>> y = y_vals;
  var lp = sum(x * 5.0);
  lp += sum(y * lp);
  stan::math::recover_memory();
}

TEST(MathRev, TestVarEigen) {
  using stan::math::sum;
  using stan::math::var;
  using stan::math::var_value;
  Eigen::Matrix<double, -1, -1> x_vals(3, 3);
  Eigen::Matrix<double, -1, -1> y_vals(3, 3);
  x_vals << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  y_vals << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  var_value<Eigen::Matrix<double, -1, -1>> x = x_vals;
  var_value<Eigen::Matrix<double, -1, -1>> y = y_vals;
  var lp = 0;
  var_value<Eigen::Matrix<double, -1, -1>> mul_xy = x * y;
  var sum_mul_xy = sum(mul_xy);

  lp -= sum_mul_xy;
  lp.grad();
  puts("-------------");
  std::cout << "lp static val: \n" << lp.val() << "\n";
  std::cout << "lp static adj: \n" << lp.adj() << "\n";
  puts("---------");
  std::cout << "sum_mul_xy static val: \n" << sum_mul_xy.val() << "\n";
  std::cout << "sum_mul_xy static adj: \n" << sum_mul_xy.adj() << "\n";
  puts("---------");

  std::cout << "mul_xy stat val: \n" << mul_xy.val() << "\n";
  std::cout << "mul_xy stat adj: \n" << mul_xy.adj() << "\n";
  puts("---------");

  std::cout << "x stat val: \n" << x.val() << "\n";
  std::cout << "x stat adj: \n" << x.adj() << "\n";
  puts("---------");

  std::cout << "y stat val: \n" << y.val() << "\n";
  std::cout << "y stat adj: \n" << y.adj() << "\n";

  stan::math::set_zero_all_adjoints();
  puts("");
  puts("---------");
  puts("Dynamic Matrix");
  puts("---------");
  puts("");
  Eigen::Matrix<var, -1, -1> x_dyn = x_vals;
  Eigen::Matrix<var, -1, -1> y_dyn = y_vals;
  Eigen::Matrix<var, -1, -1> mul_xy_dyn = multiply(x_dyn, y_dyn);
  var lp_dyn = 0;
  var sum_mul_xy_dyn = sum(mul_xy_dyn);
  lp_dyn -= sum_mul_xy_dyn;
  lp_dyn.grad();
  std::cout << "lp dyn val: \n" << lp_dyn.val() << "\n";
  std::cout << "lp dyn adj: \n" << lp_dyn.adj() << "\n";
  puts("---------");
  std::cout << "sum_mul_xy dyn val: \n" << sum_mul_xy_dyn.val() << "\n";
  std::cout << "sum_mul_xy dyn adj: \n" << sum_mul_xy_dyn.adj() << "\n";
  puts("---------");

  std::cout << "mul_xy dyn val: \n" << mul_xy_dyn.val() << "\n";
  std::cout << "mul_xy dyn adj: \n" << mul_xy_dyn.adj() << "\n";
  puts("---------");

  std::cout << "x dyn val: \n" << x_dyn.val() << "\n";
  std::cout << "x dyn adj: \n" << x_dyn.adj() << "\n";
  puts("---------");

  std::cout << "y dyn val: \n" << y_dyn.val() << "\n";
  std::cout << "y dyn adj: \n" << y_dyn.adj() << "\n";
  stan::math::recover_memory();
}

TEST(MathRev, TestVarEigenMatColVec) {
  using stan::math::sum;
  using stan::math::var;
  using stan::math::var_value;
  Eigen::Matrix<double, -1, -1> x_vals(3, 3);
  Eigen::Matrix<double, -1, 1> y_vals(3);
  x_vals << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  y_vals << 1, 2, 3;
  var_value<Eigen::Matrix<double, -1, -1>> x = x_vals;
  var_value<Eigen::Matrix<double, -1, 1>> y = y_vals;
  var lp = 0;
  auto mul_xy = x * y;
  auto sum_mul_xy = sum(mul_xy);
  lp -= sum_mul_xy;
  lp.grad();
  puts("-------------");
  std::cout << "lp static val: \n" << lp.val() << "\n";
  std::cout << "lp static adj: \n" << lp.adj() << "\n";
  puts("---------");
  std::cout << "sum_mul_xy static val: \n" << sum_mul_xy.val() << "\n";
  std::cout << "sum_mul_xy static adj: \n" << sum_mul_xy.adj() << "\n";
  puts("---------");

  std::cout << "mul_xy stat val: \n" << mul_xy.val() << "\n";
  std::cout << "mul_xy stat adj: \n" << mul_xy.adj() << "\n";
  puts("---------");

  std::cout << "x stat val: \n" << x.val() << "\n";
  std::cout << "x stat adj: \n" << x.adj() << "\n";
  puts("---------");

  std::cout << "y stat val: \n" << y.val() << "\n";
  std::cout << "y stat adj: \n" << y.adj() << "\n";

  stan::math::set_zero_all_adjoints();
  puts("");
  puts("---------");
  puts("Dynamic Matrix");
  puts("---------");
  puts("");
  Eigen::Matrix<var, -1, -1> x_dyn = x_vals;
  Eigen::Matrix<var, -1, 1> y_dyn = y_vals;
  Eigen::Matrix<var, -1, 1> mul_xy_dyn = multiply(x_dyn, y_dyn);
  var lp_dyn = 0;
  var sum_mul_xy_dyn = sum(mul_xy_dyn);
  lp_dyn -= sum_mul_xy_dyn;
  lp_dyn.grad();
  std::cout << "lp dyn val: \n" << lp_dyn.val() << "\n";
  std::cout << "lp dyn adj: \n" << lp_dyn.adj() << "\n";
  puts("---------");
  std::cout << "sum_mul_xy dyn val: \n" << sum_mul_xy_dyn.val() << "\n";
  std::cout << "sum_mul_xy dyn adj: \n" << sum_mul_xy_dyn.adj() << "\n";
  puts("---------");

  std::cout << "mul_xy dyn val: \n" << mul_xy_dyn.val() << "\n";
  std::cout << "mul_xy dyn adj: \n" << mul_xy_dyn.adj() << "\n";
  puts("---------");

  std::cout << "x dyn val: \n" << x_dyn.val() << "\n";
  std::cout << "x dyn adj: \n" << x_dyn.adj() << "\n";
  puts("---------");

  std::cout << "y dyn val: \n" << y_dyn.val() << "\n";
  std::cout << "y dyn adj: \n" << y_dyn.adj() << "\n";
  stan::math::recover_memory();
}

TEST(MathRev, TestVarEigenRowVecMat) {
  using stan::math::sum;
  using stan::math::var;
  using stan::math::var_value;
  Eigen::Matrix<double, 1, -1> x_vals(3);
  Eigen::Matrix<double, -1, -1> y_vals(3, 3);
  x_vals << 1, 2, 3;
  y_vals << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  var_value<Eigen::Matrix<double, 1, -1>> x = x_vals;
  var_value<Eigen::Matrix<double, -1, -1>> y = y_vals;
  var lp = 0;
  auto mul_xy = x * y;
  auto sum_mul_xy = sum(mul_xy);
  lp -= sum_mul_xy;
  lp.grad();
  puts("-------------");
  std::cout << "lp static val: \n" << lp.val() << "\n";
  std::cout << "lp static adj: \n" << lp.adj() << "\n";
  puts("---------");
  std::cout << "sum_mul_xy static val: \n" << sum_mul_xy.val() << "\n";
  std::cout << "sum_mul_xy static adj: \n" << sum_mul_xy.adj() << "\n";
  puts("---------");

  std::cout << "mul_xy stat val: \n" << mul_xy.val() << "\n";
  std::cout << "mul_xy stat adj: \n" << mul_xy.adj() << "\n";
  puts("---------");

  std::cout << "x stat val: \n" << x.val() << "\n";
  std::cout << "x stat adj: \n" << x.adj() << "\n";
  puts("---------");

  std::cout << "y stat val: \n" << y.val() << "\n";
  std::cout << "y stat adj: \n" << y.adj() << "\n";

  stan::math::set_zero_all_adjoints();
  puts("");
  puts("---------");
  puts("Dynamic Matrix");
  puts("---------");
  puts("");
  Eigen::Matrix<var, 1, -1> x_dyn = x_vals;
  Eigen::Matrix<var, -1, -1> y_dyn = y_vals;
  Eigen::Matrix<var, 1, -1> mul_xy_dyn = multiply(x_dyn, y_dyn);
  var lp_dyn = 0;
  var sum_mul_xy_dyn = sum(mul_xy_dyn);
  lp_dyn -= sum_mul_xy_dyn;
  lp_dyn.grad();
  std::cout << "lp dyn val: \n" << lp_dyn.val() << "\n";
  std::cout << "lp dyn adj: \n" << lp_dyn.adj() << "\n";
  puts("---------");
  std::cout << "sum_mul_xy dyn val: \n" << sum_mul_xy_dyn.val() << "\n";
  std::cout << "sum_mul_xy dyn adj: \n" << sum_mul_xy_dyn.adj() << "\n";
  puts("---------");

  std::cout << "mul_xy dyn val: \n" << mul_xy_dyn.val() << "\n";
  std::cout << "mul_xy dyn adj: \n" << mul_xy_dyn.adj() << "\n";
  puts("---------");

  std::cout << "x dyn val: \n" << x_dyn.val() << "\n";
  std::cout << "x dyn adj: \n" << x_dyn.adj() << "\n";
  puts("---------");

  std::cout << "y dyn val: \n" << y_dyn.val() << "\n";
  std::cout << "y dyn adj: \n" << y_dyn.adj() << "\n";
  stan::math::recover_memory();
}

TEST(MathRev, TestVarEigenColVecRowVec) {
  using stan::math::sum;
  using stan::math::var;
  using stan::math::var_value;
  Eigen::Matrix<double, 1, -1> x_vals(3);
  Eigen::Matrix<double, -1, 1> y_vals(3);
  x_vals << 1, 2, 3;
  y_vals << 1, 2, 3;
  var_value<Eigen::Matrix<double, 1, -1>> x = x_vals;
  var_value<Eigen::Matrix<double, -1, 1>> y = y_vals;
  var lp = 0;
  auto mul_xy = x * y;
  auto sum_mul_xy = sum(mul_xy);
  lp -= sum_mul_xy;
  lp.grad();
  puts("-------------");
  std::cout << "lp static val: \n" << lp.val() << "\n";
  std::cout << "lp static adj: \n" << lp.adj() << "\n";
  puts("---------");
  std::cout << "sum_mul_xy static val: \n" << sum_mul_xy.val() << "\n";
  std::cout << "sum_mul_xy static adj: \n" << sum_mul_xy.adj() << "\n";
  puts("---------");

  std::cout << "mul_xy stat val: \n" << mul_xy.val() << "\n";
  std::cout << "mul_xy stat adj: \n" << mul_xy.adj() << "\n";
  puts("---------");

  std::cout << "x stat val: \n" << x.val() << "\n";
  std::cout << "x stat adj: \n" << x.adj() << "\n";
  puts("---------");

  std::cout << "y stat val: \n" << y.val() << "\n";
  std::cout << "y stat adj: \n" << y.adj() << "\n";

  stan::math::set_zero_all_adjoints();
  puts("");
  puts("---------");
  puts("Dynamic Matrix");
  puts("---------");
  puts("");
  Eigen::Matrix<var, 1, -1> x_dyn = x_vals;
  Eigen::Matrix<var, -1, 1> y_dyn = y_vals;
  auto mul_xy_dyn = multiply(x_dyn, y_dyn);
  var lp_dyn = 0;
  var sum_mul_xy_dyn = sum(mul_xy_dyn);
  lp_dyn -= sum_mul_xy_dyn;
  lp_dyn.grad();
  std::cout << "lp dyn val: \n" << lp_dyn.val() << "\n";
  std::cout << "lp dyn adj: \n" << lp_dyn.adj() << "\n";
  puts("---------");
  std::cout << "sum_mul_xy dyn val: \n" << sum_mul_xy_dyn.val() << "\n";
  std::cout << "sum_mul_xy dyn adj: \n" << sum_mul_xy_dyn.adj() << "\n";
  puts("---------");

  std::cout << "mul_xy dyn val: \n" << mul_xy_dyn.val() << "\n";
  std::cout << "mul_xy dyn adj: \n" << mul_xy_dyn.adj() << "\n";
  puts("---------");

  std::cout << "x dyn val: \n" << x_dyn.val() << "\n";
  std::cout << "x dyn adj: \n" << x_dyn.adj() << "\n";
  puts("---------");

  std::cout << "y dyn val: \n" << y_dyn.val() << "\n";
  std::cout << "y dyn adj: \n" << y_dyn.adj() << "\n";
  stan::math::recover_memory();
}

TEST(MathRev, TestVarEigenRowVecColVec) {
  using stan::math::sum;
  using stan::math::var;
  using stan::math::var_value;
  Eigen::Matrix<double, -1, 1> x_vals(3);
  Eigen::Matrix<double, 1, -1> y_vals(3);
  x_vals << 1, 2, 3;
  y_vals << 1, 2, 3;
  var_value<Eigen::Matrix<double, -1, 1>> x = x_vals;
  var_value<Eigen::Matrix<double, 1, -1>> y = y_vals;
  var lp = 0;
  auto mul_xy = x * y;
  auto sum_mul_xy = sum(mul_xy);
  lp -= sum_mul_xy;
  lp.grad();
  puts("-------------");
  std::cout << "lp static val: \n" << lp.val() << "\n";
  std::cout << "lp static adj: \n" << lp.adj() << "\n";
  puts("---------");
  std::cout << "sum_mul_xy static val: \n" << sum_mul_xy.val() << "\n";
  std::cout << "sum_mul_xy static adj: \n" << sum_mul_xy.adj() << "\n";
  puts("---------");

  std::cout << "mul_xy stat val: \n" << mul_xy.val() << "\n";
  std::cout << "mul_xy stat adj: \n" << mul_xy.adj() << "\n";
  puts("---------");

  std::cout << "x stat val: \n" << x.val() << "\n";
  std::cout << "x stat adj: \n" << x.adj() << "\n";
  puts("---------");

  std::cout << "y stat val: \n" << y.val() << "\n";
  std::cout << "y stat adj: \n" << y.adj() << "\n";

  stan::math::set_zero_all_adjoints();
  puts("");
  puts("---------");
  puts("Dynamic Matrix");
  puts("---------");
  puts("");
  Eigen::Matrix<var, -1, 1> x_dyn = x_vals;
  Eigen::Matrix<var, 1, -1> y_dyn = y_vals;
  auto mul_xy_dyn = multiply(x_dyn, y_dyn);
  var lp_dyn = 0;
  var sum_mul_xy_dyn = sum(mul_xy_dyn);
  lp_dyn -= sum_mul_xy_dyn;
  lp_dyn.grad();
  std::cout << "lp dyn val: \n" << lp_dyn.val() << "\n";
  std::cout << "lp dyn adj: \n" << lp_dyn.adj() << "\n";
  puts("---------");
  std::cout << "sum_mul_xy dyn val: \n" << sum_mul_xy_dyn.val() << "\n";
  std::cout << "sum_mul_xy dyn adj: \n" << sum_mul_xy_dyn.adj() << "\n";
  puts("---------");

  std::cout << "mul_xy dyn val: \n" << mul_xy_dyn.val() << "\n";
  std::cout << "mul_xy dyn adj: \n" << mul_xy_dyn.adj() << "\n";
  puts("---------");

  std::cout << "x dyn val: \n" << x_dyn.val() << "\n";
  std::cout << "x dyn adj: \n" << x_dyn.adj() << "\n";
  puts("---------");

  std::cout << "y dyn val: \n" << y_dyn.val() << "\n";
  std::cout << "y dyn adj: \n" << y_dyn.adj() << "\n";
  stan::math::recover_memory();
}

TEST(MathRev, TestVarEigenVarDbl) {
  using stan::math::sum;
  using stan::math::var;
  using stan::math::var_value;
  Eigen::Matrix<double, -1, -1> x_vals(3, 3);
  Eigen::Matrix<double, -1, -1> y_vals(3, 3);
  x_vals << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  y_vals << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  var_value<Eigen::Matrix<double, -1, -1>> x = x_vals;
  var lp = 0;
  auto mul_xy = x * y_vals;
  auto sum_mul_xy = sum(mul_xy);
  lp -= sum_mul_xy;
  lp.grad();
  puts("-------------");
  std::cout << "lp static val: \n" << lp.val() << "\n";
  std::cout << "lp static adj: \n" << lp.adj() << "\n";
  puts("---------");
  std::cout << "sum_mul_xy static val: \n" << sum_mul_xy.val() << "\n";
  std::cout << "sum_mul_xy static adj: \n" << sum_mul_xy.adj() << "\n";
  puts("---------");

  std::cout << "mul_xy stat val: \n" << mul_xy.val() << "\n";
  std::cout << "mul_xy stat adj: \n" << mul_xy.adj() << "\n";
  puts("---------");

  std::cout << "x stat val: \n" << x.val() << "\n";
  std::cout << "x stat adj: \n" << x.adj() << "\n";

  stan::math::set_zero_all_adjoints();
  puts("");
  puts("---------");
  puts("Dynamic Matrix");
  puts("---------");
  puts("");
  Eigen::Matrix<var, -1, -1> x_dyn = x_vals;
  Eigen::Matrix<var, -1, -1> mul_xy_dyn = multiply(x_dyn, y_vals);
  var lp_dyn = 0;
  var sum_mul_xy_dyn = sum(mul_xy_dyn);
  lp_dyn -= sum_mul_xy_dyn;
  lp_dyn.grad();
  std::cout << "lp dyn val: \n" << lp_dyn.val() << "\n";
  std::cout << "lp dyn adj: \n" << lp_dyn.adj() << "\n";
  puts("---------");
  std::cout << "sum_mul_xy dyn val: \n" << sum_mul_xy_dyn.val() << "\n";
  std::cout << "sum_mul_xy dyn adj: \n" << sum_mul_xy_dyn.adj() << "\n";
  puts("---------");

  std::cout << "mul_xy dyn val: \n" << mul_xy_dyn.val() << "\n";
  std::cout << "mul_xy dyn adj: \n" << mul_xy_dyn.adj() << "\n";
  puts("---------");

  std::cout << "x dyn val: \n" << x_dyn.val() << "\n";
  std::cout << "x dyn adj: \n" << x_dyn.adj() << "\n";
  puts("---------");
  stan::math::recover_memory();
}

TEST(MathRev, TestVarEigenDblVal) {
  using stan::math::sum;
  using stan::math::var;
  using stan::math::var_value;
  Eigen::Matrix<double, -1, -1> x_vals(3, 3);
  Eigen::Matrix<double, -1, -1> y_vals(3, 3);
  x_vals << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  y_vals << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  var_value<Eigen::Matrix<double, -1, -1>> y = y_vals;
  var lp = 0;
  auto mul_xy = x_vals * y;
  auto sum_mul_xy = sum(mul_xy);
  lp -= sum_mul_xy;
  lp.grad();
  puts("-------------");
  std::cout << "lp static val: \n" << lp.val() << "\n";
  std::cout << "lp static adj: \n" << lp.adj() << "\n";
  puts("---------");
  std::cout << "sum_mul_xy static val: \n" << sum_mul_xy.val() << "\n";
  std::cout << "sum_mul_xy static adj: \n" << sum_mul_xy.adj() << "\n";
  puts("---------");

  std::cout << "mul_xy stat val: \n" << mul_xy.val() << "\n";
  std::cout << "mul_xy stat adj: \n" << mul_xy.adj() << "\n";
  puts("---------");

  std::cout << "y stat val: \n" << y.val() << "\n";
  std::cout << "y stat adj: \n" << y.adj() << "\n";

  stan::math::set_zero_all_adjoints();
  puts("");
  puts("---------");
  puts("Dynamic Matrix");
  puts("---------");
  puts("");
  Eigen::Matrix<var, -1, -1> y_dyn = y_vals;
  Eigen::Matrix<var, -1, -1> mul_xy_dyn = multiply(x_vals, y_dyn);
  var lp_dyn = 0;
  var sum_mul_xy_dyn = sum(mul_xy_dyn);
  lp_dyn -= sum_mul_xy_dyn;
  lp_dyn.grad();
  std::cout << "lp dyn val: \n" << lp_dyn.val() << "\n";
  std::cout << "lp dyn adj: \n" << lp_dyn.adj() << "\n";
  puts("---------");
  std::cout << "sum_mul_xy dyn val: \n" << sum_mul_xy_dyn.val() << "\n";
  std::cout << "sum_mul_xy dyn adj: \n" << sum_mul_xy_dyn.adj() << "\n";
  puts("---------");

  std::cout << "mul_xy dyn val: \n" << mul_xy_dyn.val() << "\n";
  std::cout << "mul_xy dyn adj: \n" << mul_xy_dyn.adj() << "\n";
  puts("---------");

  std::cout << "y dyn val: \n" << y_dyn.val() << "\n";
  std::cout << "y dyn adj: \n" << y_dyn.adj() << "\n";
  stan::math::recover_memory();
}
*/
