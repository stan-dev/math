#include <stan/math/rev/core.hpp>
#include <test/unit/math/expect_near_rel.hpp>
#include <test/unit/math/rev/fun/util.hpp>
#include <gtest/gtest.h>
#include <vector>
#include <limits>

using var_matrix = stan::math::var_value<Eigen::MatrixXd>;
using var_vector = stan::math::var_value<Eigen::VectorXd>;
using var_row_vector  = stan::math::var_value<Eigen::RowVectorXd>;
using stan::math::var;
using stan::math::var_value;
using stan::math::sum;

auto make_var_vector(size_t N = 3) {
  Eigen::Matrix<var, Eigen::Dynamic, 1> out(N);
  for(size_t n = 0; n < N; ++n) {
    out(n) = n + 1;
  }
  return out;
}

auto make_var_row_vector(size_t N = 3) {
  Eigen::Matrix<var, 1, Eigen::Dynamic> out(N);
  for(size_t n = 0; n < N; ++n) {
    out(n) = n + 1;
  }
  return out;
}

auto make_var_matrix(size_t N = 3, size_t M = 3) {
  Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> out(N, M);
  for(size_t n = 0; n < N * M; ++n) {
    out(n) = n + 1;
  }
  return out;
}

template <typename F, int R, int C>
void expect_ad2(const F& f, const Eigen::Matrix<double, R, C>& x) {
  Eigen::Matrix<var, R, C> xv = x;
  var_value<Eigen::Matrix<double, R, C>> vx = x;

  auto o = f(x).eval();
  auto ov = f(xv);
  auto vo = f(vx);

  stan::math::set_zero_all_adjoints();
  sum(ov).grad();
  auto xv_adj = xv.adj().eval();
  stan::math::set_zero_all_adjoints();
  sum(vo).grad();
  auto vx_adj = vx.adj().eval();

  auto g = [&](const Eigen::VectorXd& flat_x) {
    Eigen::Matrix<double, R, C> x_shaped_variable(x.rows(), x.cols());
    for(size_t i = 0; i < flat_x.size(); ++i)
      x_shaped_variable(i) = flat_x(i);

    auto o = f(x_shaped_variable);

    return sum(o);
  };

  double tmp;
  Eigen::VectorXd flat_x(x.size());
  Eigen::VectorXd grad_x;
  for(size_t i = 0; i < flat_x.size(); ++i)
    flat_x(i) = x(i);
  stan::math::finite_diff_gradient_auto(g, flat_x, tmp, grad_x);
  Eigen::Matrix<double, R, C> x_shaped_grad(x.rows(), x.cols());
  for(size_t i = 0; i < grad_x.size(); ++i)
    x_shaped_grad(i) = grad_x(i);

  stan::test::expect_near_rel("matrix of vars function output", o, value_of(ov));
  stan::test::expect_near_rel("var matrix function output", o, value_of(vo));
  stan::test::expect_near_rel("matrix of vars argument adjoints", x_shaped_grad, xv_adj);
  stan::test::expect_near_rel("var matrix argument adjoints", x_shaped_grad, vx_adj);
}

/*TEST(MathRev, cast_vector) {
  Eigen::VectorXd vec = Eigen::VectorXd::Random(3);

  auto f = [&](const auto& arg) {
    return arg;
  };

  expect_ad2(f, vec);
  //expect_ad2(f2, vec);
  }*/

TEST(MathRev, scalar_vector) {
  double s = 1.5;
  Eigen::VectorXd vec = Eigen::VectorXd::Random(3);

  auto f1 = [&](const auto& arg) {
    return s * arg;
  };

  auto f2 = [&](const auto& arg) {
    return arg * s;
  };

  expect_ad2(f1, vec);
  expect_ad2(f2, vec);
}

TEST(MathRev, scalar_row_vector) {
  double s = 1.5;
  Eigen::RowVectorXd rvec = Eigen::RowVectorXd::Random(3);

  auto f1 = [&](const auto& arg) {
    return s * arg;
  };

  auto f2 = [&](const auto& arg) {
    return arg * s;
  };

  expect_ad2(f1, rvec);
  expect_ad2(f2, rvec);
}

TEST(MathRev, scalar_matrix) {
  double s = 1.5;
  Eigen::MatrixXd mat = Eigen::MatrixXd::Random(3, 2);

  auto f1 = [&](const auto& arg) {
    return s * arg;
  };

  auto f2 = [&](const auto& arg) {
    return arg * s;
  };

  expect_ad2(f1, mat);
  expect_ad2(f2, mat);
}

TEST(MathRev, vector_row_vector) {
  Eigen::VectorXd vec = Eigen::VectorXd::Random(3);
  Eigen::RowVectorXd rvec = Eigen::RowVectorXd::Random(3);

  /*auto f1 = [&](const auto& arg) {
    return arg * vec;
    };*/

  auto f2 = [&](const auto& arg) {
    return arg * rvec;
  };

  //expect_ad2(f1, rvec);
  expect_ad2(f2, vec);
}

TEST(MathRev, vector_matrix) {
  Eigen::VectorXd vec = Eigen::VectorXd::Random(3);
  Eigen::MatrixXd mat = Eigen::MatrixXd::Random(3, 3);

  auto f1 = [&](const auto& arg) {
    return arg * vec;
  };

  auto f2 = [&](const auto& arg) {
    return mat * arg;
  };

  expect_ad2(f1, mat);
  //expect_ad2(f2, vec);
}

TEST(MathRev, row_vector_matrix) {
  Eigen::RowVectorXd rvec = Eigen::RowVectorXd::Random(3);
  Eigen::MatrixXd mat = Eigen::MatrixXd::Random(3, 3);

  auto f1 = [&](const auto& arg) {
    return arg * mat;
  };

  auto f2 = [&](const auto& arg) {
    return rvec * arg;
  };

  expect_ad2(f1, rvec);
  expect_ad2(f2, mat);
}

/*TEST(MathRev, matrix_matrix) {
  var_matrix mat1 = make_var_matrix(3);
  var_matrix mat2 = make_var_matrix(3);

  var_matrix tmp = mat1 * mat2;

  sum(tmp).grad();
  }*/

/*TEST(MathRev, TestVarEigen) {
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
