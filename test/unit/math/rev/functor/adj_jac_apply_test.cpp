#include <stan/math/rev/core.hpp>
#include <test/unit/math/rev/util.hpp>
#include <gtest/gtest.h>
#include <algorithm>
#include <sstream>
#include <tuple>
#include <vector>

/*
 * Check scalar return type
 */
template <typename T>
struct ScalarSinFunctor {
  stan::math::adj_op<T> x_;
  explicit ScalarSinFunctor(const T& x) : x_(x) {}

  double operator()(const double& x) {
    return sin(x_.map());
  }

  auto multiply_adjoint_jacobian(const double& adj) {
    return std::make_tuple(cos(x_.map()) * adj);
  }
};

TEST(AgradRev, test_scalar_sin_stack) {
  using stan::math::var;
  stan::math::var x1, y1;
  x1 = 1.0;

  y1 = stan::math::adj_jac_apply<ScalarSinFunctor<var>>(x1);

  test::check_varis_on_stack(y1);
}

TEST(AgradRev, test_scalar_sin_values) {
  using stan::math::var;
  stan::math::var x1, y1;
  x1 = 1.0;

  y1 = stan::math::adj_jac_apply<ScalarSinFunctor<var>>(x1);

  EXPECT_NEAR(y1.val(), 0.841470984807897, 1e-10);
}

TEST(AgradRev, test_scalar_sin_jac) {
  using stan::math::var;
  stan::math::var x1, y1;
  x1 = 1.0;

  y1 = stan::math::adj_jac_apply<ScalarSinFunctor<var>>(x1);

  y1.grad();
  EXPECT_NEAR(x1.adj(), 0.5403023058681398, 1e-10);
}

/*
 * Check std::vector return type
 */
template <typename T>
struct StdVectorSinFunctor {
  stan::math::adj_op<T> x_;
  // double* x_;
  // int N_;
  explicit StdVectorSinFunctor(const T& x) : x_(x) {}

  std::vector<double> operator()(const std::vector<double>& x) {
    std::vector<double> out(x_.size());
    for (int i = 0; i < x_.size(); ++i) {
      out[i] = sin(x[i]);
    }
    return out;
  }

  auto multiply_adjoint_jacobian(const std::vector<double>& adj) {
    std::vector<double> adj_jac(x_.size());
    for (int i = 0; i < x_.size(); ++i) {
      adj_jac[i] = cos(x_(i)) * adj[i];
    }
    return std::make_tuple(adj_jac);
  }
};

TEST(AgradRev, test_std_vector_sin_stack) {
  std::vector<stan::math::var> x1(2), y1;
  x1[0] = 1.0;
  x1[1] = 2.0;

  y1 = stan::math::adj_jac_apply<
      StdVectorSinFunctor<std::vector<stan::math::var>>>(x1);

  test::check_varis_on_stack(y1);
}

TEST(AgradRev, test_std_vector_sin_values) {
  std::vector<stan::math::var> x1(2), y1;
  x1[0] = 1.0;
  x1[1] = 2.0;

  y1 = stan::math::adj_jac_apply<
      StdVectorSinFunctor<std::vector<stan::math::var>>>(x1);

  EXPECT_NEAR(y1[0].val(), 0.841470984807897, 1e-10);
  EXPECT_NEAR(y1[1].val(), 0.909297426825682, 1e-10);
}

TEST(AgradRev, test_std_vector_sin_jac) {
  std::vector<stan::math::var> x1(2), y1;
  x1[0] = 1.0;
  x1[1] = 2.0;

  y1 = stan::math::adj_jac_apply<
      StdVectorSinFunctor<std::vector<stan::math::var>>>(x1);

  y1[0].grad();
  EXPECT_NEAR(x1[0].adj(), 0.5403023058681398, 1e-10);
  EXPECT_NEAR(x1[1].adj(), 0.0, 1e-10);

  stan::math::set_zero_all_adjoints();
  y1[1].grad();
  EXPECT_FLOAT_EQ(x1[0].adj(), 0.0);
  EXPECT_FLOAT_EQ(x1[1].adj(), -0.4161468365471424);
}

/**
 * Test Eigen::VectorXd return types
 */
template <typename T>
struct SinFunctor {
  stan::math::adj_op<T> x_;
  explicit SinFunctor(const T& x) : x_(x) {}
  Eigen::VectorXd operator()(const Eigen::VectorXd& x) {
    Eigen::VectorXd out(x_.size());
    for (int n = 0; n < x_.size(); ++n) {
      out(n) = sin(x(n));
    }
    return out;
  }

  auto multiply_adjoint_jacobian(const Eigen::VectorXd& adj) {
    Eigen::VectorXd out(x_.size());
    for (int n = 0; n < x_.size(); ++n) {
      out(n) = cos(x_(n)) * adj(n);
    }
    return std::make_tuple(out);
  }
};

TEST(AgradRev, test_vector_sin_stack) {
  using eig_vec = Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>;
  eig_vec x1(1), x2(2), y1(1), y2(2);
  x1 << 1.0;
  x2 << 2.0, 1.0;

  y1 = stan::math::adj_jac_apply<SinFunctor<eig_vec>>(x1);
  y2 = stan::math::adj_jac_apply<SinFunctor<eig_vec>>(x2);

  test::check_varis_on_stack(y1);
  test::check_varis_on_stack(y2);
}

TEST(AgradRev, test_vector_sin_values) {
  using eig_vec = Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>;
  eig_vec x1(1), x2(2), y1(1), y2(2);
  x1 << 1.0;
  x2 << 2.0, 1.0;

  y1 = stan::math::adj_jac_apply<SinFunctor<eig_vec>>(x1);
  y2 = stan::math::adj_jac_apply<SinFunctor<eig_vec>>(x2);

  EXPECT_FLOAT_EQ(y1(0).val(), 0.841470984807897);
  EXPECT_FLOAT_EQ(y2(0).val(), 0.909297426825682);
  EXPECT_FLOAT_EQ(y2(1).val(), 0.841470984807897);
}

TEST(AgradRev, test_vector_sin_multiple_jac) {
  using eig_vec = Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>;
  eig_vec x1(1), x2(2), y1(1), y2(2);
  x1 << 1.0;
  x2 << 2.0, 1.0;

  y1 = stan::math::adj_jac_apply<SinFunctor<eig_vec>>(x1);
  y2 = stan::math::adj_jac_apply<SinFunctor<eig_vec>>(x2);

  y1(0).grad();
  EXPECT_FLOAT_EQ(x1(0).adj(), 0.5403023058681398);
  EXPECT_FLOAT_EQ(x2(0).adj(), 0.0);
  EXPECT_FLOAT_EQ(x2(1).adj(), 0.0);

  stan::math::set_zero_all_adjoints();

  y2(0).grad();
  EXPECT_FLOAT_EQ(x1(0).adj(), 0.0);
  EXPECT_FLOAT_EQ(x2(0).adj(), -0.4161468365471424);
  EXPECT_FLOAT_EQ(x2(1).adj(), 0.0);

  stan::math::set_zero_all_adjoints();

  y2(1).grad();
  EXPECT_FLOAT_EQ(x1(0).adj(), 0.0);
  EXPECT_FLOAT_EQ(x2(0).adj(), 0.0);
  EXPECT_FLOAT_EQ(x2(1).adj(), 0.5403023058681398);

  stan::math::set_zero_all_adjoints();

  stan::math::var sum_y2 = (1.73 * y2(0) + 1.57 * y2(1));
  sum_y2.grad();
  EXPECT_FLOAT_EQ(x1(0).adj(), 0.0);
  EXPECT_FLOAT_EQ(x2(0).adj(), 1.73 * -0.4161468365471424);
  EXPECT_FLOAT_EQ(x2(1).adj(), 1.57 * 0.5403023058681398);
}

/**
 * Test Eigen::RowVectorXd return types
 */
template <typename T>
struct RowVectorSinFunctor {
  stan::math::adj_op<T> x_;
  //  int N_;
  //  double* x_mem_;
  explicit RowVectorSinFunctor(const T& x) : x_(x) {}
  Eigen::RowVectorXd operator()(const Eigen::RowVectorXd& x) {
    Eigen::RowVectorXd out(x_.size());
    for (int n = 0; n < x_.size(); ++n) {
      out(n) = sin(x(n));
    }

    return out;
  }

  auto multiply_adjoint_jacobian(const Eigen::RowVectorXd& adj) {
    Eigen::RowVectorXd out(x_.size());
    for (int n = 0; n < x_.size(); ++n) {
      out(n) = cos(x_(n)) * adj(n);
    }
    return std::make_tuple(out);
  }
};

TEST(AgradRev, test_row_vector_sin_stack) {
  using eig_row_vec = Eigen::Matrix<stan::math::var, 1, Eigen::Dynamic>;
  eig_row_vec x1(1), x2(2), y1(1), y2(2);
  x1 << 1.0;
  x2 << 2.0, 1.0;

  y1 = stan::math::adj_jac_apply<RowVectorSinFunctor<eig_row_vec>>(x1);
  y2 = stan::math::adj_jac_apply<RowVectorSinFunctor<eig_row_vec>>(x2);

  test::check_varis_on_stack(y1);
  test::check_varis_on_stack(y2);
}

TEST(AgradRev, test_row_vector_sin_values) {
  using eig_row_vec = Eigen::Matrix<stan::math::var, 1, Eigen::Dynamic>;
  eig_row_vec x1(1), x2(2), y1(1), y2(2);
  x1 << 1.0;
  x2 << 2.0, 1.0;

  y1 = stan::math::adj_jac_apply<RowVectorSinFunctor<eig_row_vec>>(x1);
  y2 = stan::math::adj_jac_apply<RowVectorSinFunctor<eig_row_vec>>(x2);

  EXPECT_FLOAT_EQ(y1(0).val(), 0.841470984807897);
  EXPECT_FLOAT_EQ(y2(0).val(), 0.909297426825682);
  EXPECT_FLOAT_EQ(y2(1).val(), 0.841470984807897);
}

TEST(AgradRev, test_row_vector_sin_multiple_jac) {
  using eig_row_vec = Eigen::Matrix<stan::math::var, 1, Eigen::Dynamic>;
  eig_row_vec x1(1), x2(2), y1(1), y2(2);
  x1 << 1.0;
  x2 << 2.0, 1.0;

  y1 = stan::math::adj_jac_apply<RowVectorSinFunctor<eig_row_vec>>(x1);
  y2 = stan::math::adj_jac_apply<RowVectorSinFunctor<eig_row_vec>>(x2);

  y1(0).grad();
  EXPECT_FLOAT_EQ(x1(0).adj(), 0.5403023058681398);
  EXPECT_FLOAT_EQ(x2(0).adj(), 0.0);
  EXPECT_FLOAT_EQ(x2(1).adj(), 0.0);

  stan::math::set_zero_all_adjoints();

  y2(0).grad();
  EXPECT_FLOAT_EQ(x1(0).adj(), 0.0);
  EXPECT_FLOAT_EQ(x2(0).adj(), -0.4161468365471424);
  EXPECT_FLOAT_EQ(x2(1).adj(), 0.0);

  stan::math::set_zero_all_adjoints();

  y2(1).grad();
  EXPECT_FLOAT_EQ(x1(0).adj(), 0.0);
  EXPECT_FLOAT_EQ(x2(0).adj(), 0.0);
  EXPECT_FLOAT_EQ(x2(1).adj(), 0.5403023058681398);

  stan::math::set_zero_all_adjoints();

  stan::math::var sum_y2 = (1.73 * y2(0) + 1.57 * y2(1));
  sum_y2.grad();
  EXPECT_FLOAT_EQ(x1(0).adj(), 0.0);
  EXPECT_FLOAT_EQ(x2(0).adj(), 1.73 * -0.4161468365471424);
  EXPECT_FLOAT_EQ(x2(1).adj(), 1.57 * 0.5403023058681398);
}

/**
 * Test Eigen::MatrixXd return types
 */
template <typename T>
struct MatrixSinFunctor {
  stan::math::adj_op<T> x_;
  explicit MatrixSinFunctor(const T& x) : x_(x) {}
  Eigen::MatrixXd operator()(const Eigen::MatrixXd& x) {
    Eigen::MatrixXd out(x_.rows(), x_.cols());
    for (int n = 0; n < x_.size(); ++n) {
      out(n) = sin(x(n));
    }

    return out;
  }

  auto multiply_adjoint_jacobian(const Eigen::MatrixXd& adj) {
    Eigen::MatrixXd out(x_.rows(), x_.cols());
    for (int n = 0; n < x_.size(); ++n) {
      out(n) = cos(x_(n)) * adj(n);
    }
    return std::make_tuple(out);
  }
};

TEST(AgradRev, test_matrix_sin_stack) {
  using eig_mat
      = Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic>;
  eig_mat x(2, 2), y(2, 2);
  x << 2.0, 1.0, 0.0, -1.0;
  y = stan::math::adj_jac_apply<MatrixSinFunctor<eig_mat>>(x);
  test::check_varis_on_stack(y);
}

TEST(AgradRev, test_matrix_sin_values) {
  using eig_mat
      = Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic>;
  eig_mat x(2, 2), y(2, 2);
  x << 2.0, 1.0, 0.0, -1.0;
  y = stan::math::adj_jac_apply<MatrixSinFunctor<eig_mat>>(x);
  EXPECT_FLOAT_EQ(y(0, 0).val(), 0.909297426825682);
  EXPECT_FLOAT_EQ(y(0, 1).val(), 0.841470984807897);
  EXPECT_FLOAT_EQ(y(1, 0).val(), 0.0);
  EXPECT_FLOAT_EQ(y(1, 1).val(), -0.841470984807897);
}

TEST(AgradRev, test_matrix_sin_multiple_jac) {
  using eig_mat
      = Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic>;
  eig_mat x(2, 2), y(2, 2);
  x << 2.0, 1.0, 0.0, -1.0;
  y = stan::math::adj_jac_apply<MatrixSinFunctor<eig_mat>>(x);
  y(0, 0).grad();
  EXPECT_FLOAT_EQ(x(0, 0).adj(), -0.4161468365471424);
  EXPECT_FLOAT_EQ(x(0, 1).adj(), 0.0);
  EXPECT_FLOAT_EQ(x(1, 0).adj(), 0.0);
  EXPECT_FLOAT_EQ(x(1, 1).adj(), 0.0);

  stan::math::set_zero_all_adjoints();
  y(0, 1).grad();
  EXPECT_FLOAT_EQ(x(0, 0).adj(), 0.0);
  EXPECT_FLOAT_EQ(x(0, 1).adj(), 0.5403023058681398);
  EXPECT_FLOAT_EQ(x(1, 0).adj(), 0.0);
  EXPECT_FLOAT_EQ(x(1, 1).adj(), 0.0);

  stan::math::set_zero_all_adjoints();
  y(1, 0).grad();
  EXPECT_FLOAT_EQ(x(0, 0).adj(), 0.0);
  EXPECT_FLOAT_EQ(x(0, 1).adj(), 0.0);
  EXPECT_FLOAT_EQ(x(1, 0).adj(), 1.0);
  EXPECT_FLOAT_EQ(x(1, 1).adj(), 0.0);

  stan::math::set_zero_all_adjoints();
  y(1, 1).grad();
  EXPECT_FLOAT_EQ(x(0, 0).adj(), 0.0);
  EXPECT_FLOAT_EQ(x(0, 1).adj(), 0.0);
  EXPECT_FLOAT_EQ(x(1, 0).adj(), 0.0);
  EXPECT_FLOAT_EQ(x(1, 1).adj(), 0.5403023058681398);
}

/**
 * Test that all types that are documented to be supported can actually be
 * included in functor. Check initialization of is_vars_ and offsets_ are to
 * specification.
 */
struct WeirdArgumentListFunctor1 {
  template <typename... Args>
  WeirdArgumentListFunctor1(const Args&... args) {}
  Eigen::VectorXd operator()(double, int, const double&, const int&,
      std::vector<double>, std::vector<int>, const std::vector<double>&,
      const std::vector<int>&, Eigen::Matrix<double, Eigen::Dynamic, 1>,
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>,
      Eigen::Matrix<double, 2, Eigen::Dynamic>, Eigen::Matrix<double, 5, 1>,
      const Eigen::Matrix<double, Eigen::Dynamic, 1>&,
      const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>&,
      const Eigen::Matrix<double, 2, Eigen::Dynamic>&,
      const Eigen::Matrix<double, 5, 1>&) {
    return Eigen::VectorXd(1);
  }

  auto multiply_adjoint_jacobian(const Eigen::VectorXd& y_adj) {
    return std::make_tuple(
        double(), int(), double(), int(), std::vector<double>(),
        std::vector<int>(), std::vector<double>(), std::vector<int>(),
        Eigen::Matrix<double, Eigen::Dynamic, 1>(),
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(),
        Eigen::Matrix<double, 2, Eigen::Dynamic>(),
        Eigen::Matrix<double, 5, 1>(),
        Eigen::Matrix<double, Eigen::Dynamic, 1>(),
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(),
        Eigen::Matrix<double, 2, Eigen::Dynamic>(),
        Eigen::Matrix<double, 5, 1>());
  }
};

template <typename F, typename... Targs>
auto make_vari_for_test(const Targs&... args) {
  auto vi = new stan::math::adj_jac_vari<F, Targs...>(args...);
  return vi;
}

TEST(AgradRev,
     test_weird_argument_list_functor_compiles_and_sets_is_var_and_offsets_) {
  int i;
  double d;
  stan::math::var v(5.0);
  std::vector<int> vi(2);
  std::vector<double> vd(2);
  std::vector<stan::math::var> vv(2, 0);
  Eigen::Matrix<double, Eigen::Dynamic, 1> ed1(3);
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> ed2(2, 4);
  Eigen::Matrix<double, 2, Eigen::Dynamic> ed3(2, 5);
  Eigen::Matrix<double, 5, 1> ed4;
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> ev1(3);
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> ev2(2, 4);
  Eigen::Matrix<stan::math::var, 2, Eigen::Dynamic> ev3(2, 5);
  Eigen::Matrix<stan::math::var, 5, 1> ev4;
  ev1.setZero();
  ev2.setZero();
  ev3.setZero();
  ev4.setZero();
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> y1
      = stan::math::adj_jac_apply<WeirdArgumentListFunctor1>(
          d, i, d, i, vd, vi, vd, vi, ed1, ed2, ed3, ed4, ed1, ed2, ed3, ed4);

  y1(0).grad();

  auto vi1 = make_vari_for_test<WeirdArgumentListFunctor1>(
      d, i, d, i, vd, vi, vd, vi, ed1, ed2, ed3, ed4, ed1, ed2, ed3, ed4);

  EXPECT_EQ(vi1->is_var_,
            (std::array<bool, 16>(
                {{false, false, false, false, false, false, false, false, false,
                  false, false, false, false, false, false, false}})));

  stan::math::var(vi1).grad();

  auto vi2 = make_vari_for_test<WeirdArgumentListFunctor1>(
      v, i, d, i, vv, vi, vd, vi, ev1, ed2, ev3, ed4, ev1, ed2, ev3, ed4);

  EXPECT_EQ(vi2->is_var_,
            (std::array<bool, 16>(
                {{true, false, false, false, true, false, false, false, true,
                  false, true, false, true, false, true, false}})));

  stan::math::var(vi2).grad();

  auto vi3 = make_vari_for_test<WeirdArgumentListFunctor1>(
      d, i, v, i, vd, vi, vv, vi, ed1, ev2, ed3, ev4, ed1, ev2, ed3, ev4);

  EXPECT_EQ(vi3->is_var_,
            (std::array<bool, 16>(
                {{false, false, true, false, false, false, true, false, false,
                  true, false, true, false, true, false, true}})));

  stan::math::var(vi3).grad();

  auto vi4 = make_vari_for_test<WeirdArgumentListFunctor1>(
      v, i, d, i, vd, vi, vv, vi, ev1, ed2, ed3, ev4, ev1, ed2, ed3, ev4);

  EXPECT_EQ(vi4->is_var_,
            (std::array<bool, 16>(
                {{true, false, false, false, false, false, true, false, true,
                  false, false, true, true, false, false, true}})));

  stan::math::var(vi4).grad();
}

/**
 * Test to make sure variable values get passed forward and adjoint values get
 * passed back for all var types. Mix in some integer types for good measure
 *
 * Repeat this test while also individually making some of the could-be-autodiff
 * variables doubles instead
 */
struct CheckAdjointsPassingThrough {
  int size_vd;
  int rows_ed1;
  int rows_ed2;
  int cols_ed2;
  int cols_ed3;
  template <typename... Args>
  CheckAdjointsPassingThrough(const Args&... args) {}
  Eigen::VectorXd operator()(const double& d,
      const std::vector<double>& vd, const int&,
      const Eigen::Matrix<double, Eigen::Dynamic, 1>& ed1,
      const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& ed2,
      const std::vector<int>&,
      const Eigen::Matrix<double, 1, Eigen::Dynamic>& ed3,
      const Eigen::Matrix<double, 1, 1>& ed4) {
    size_vd = vd.size();
    rows_ed1 = ed1.rows();
    rows_ed2 = ed2.rows();
    cols_ed2 = ed2.cols();
    cols_ed3 = ed3.cols();
    Eigen::VectorXd out(1 + size_vd + rows_ed1 + rows_ed2 * cols_ed2 + cols_ed3
                        + 1);

    out(0) = d;
    for (int i = 0; i < size_vd; i++)
      out(1 + i) = vd[i];
    for (int i = 0; i < rows_ed1; i++)
      out(1 + size_vd + i) = ed1(i);
    for (int i = 0; i < rows_ed2 * cols_ed2; i++)
      out(1 + size_vd + rows_ed1 + i) = ed2(i);
    for (int i = 0; i < cols_ed3; i++)
      out(1 + size_vd + rows_ed1 + rows_ed2 * cols_ed2 + i) = ed3(i);
    out(1 + size_vd + rows_ed1 + rows_ed2 * cols_ed2 + cols_ed3) = ed4(0);

    return out;
  }

  auto multiply_adjoint_jacobian(const Eigen::VectorXd& y_adj) {
    double d;
    std::vector<double> vd(size_vd);
    Eigen::Matrix<double, Eigen::Dynamic, 1> ed1(rows_ed1);
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> ed2(rows_ed2,
                                                              cols_ed2);
    Eigen::Matrix<double, 1, Eigen::Dynamic> ed3(cols_ed3);
    Eigen::Matrix<double, 1, 1> ed4;
    d = y_adj(0);
    for (int i = 0; i < size_vd; i++)
      vd[i] = y_adj(1 + i);
    for (int i = 0; i < rows_ed1; i++)
      ed1(i) = y_adj(1 + size_vd + i);
    for (int i = 0; i < rows_ed2 * cols_ed2; i++)
      ed2(i) = y_adj(1 + size_vd + rows_ed1 + i);
    for (int i = 0; i < cols_ed3; i++)
      ed3(i) = y_adj(1 + size_vd + rows_ed1 + rows_ed2 * cols_ed2 + i);
    ed4(0) = y_adj(1 + size_vd + rows_ed1 + rows_ed2 * cols_ed2 + cols_ed3);
    return std::make_tuple(d, vd, int(), ed1, ed2, std::vector<int>(), ed3,
                           ed4);
  }
};

TEST(AgradRev, test_pass_through_working_all_var_types) {
  int size_vd = 5, rows_ed1 = 3, rows_ed2 = 2, cols_ed2 = 3, cols_ed3 = 4;
  stan::math::var d = 1.0;
  std::vector<stan::math::var> vd(size_vd, 1.0);
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> ed1(rows_ed1);
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> ed2(rows_ed2,
                                                                     cols_ed2);
  std::vector<int> vi(5, 5);
  Eigen::Matrix<stan::math::var, 1, Eigen::Dynamic> ed3(cols_ed3);
  Eigen::Matrix<stan::math::var, 1, 1> ed4;

  for (int i = 0; i < size_vd; i++)
    vd[i] = 1.0;

  for (int i = 0; i < rows_ed1; i++)
    ed1(i) = 1.0;

  for (int i = 0; i < rows_ed2 * cols_ed2; i++)
    ed2(i) = 1.0;

  for (int i = 0; i < cols_ed3; i++)
    ed3(i) = 1.0;

  ed4(0) = 1.0;

  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> y
      = stan::math::adj_jac_apply<CheckAdjointsPassingThrough>(
          d, vd, 5, ed1, ed2, vi, ed3, ed4);

  y(0).grad();
  EXPECT_FLOAT_EQ(y(0).val(), d.val());
  EXPECT_FLOAT_EQ(d.adj(), 1.0);
  for (int i = 0; i < size_vd; i++)
    EXPECT_FLOAT_EQ(vd[i].adj(), 0);
  for (int i = 0; i < rows_ed1; i++)
    EXPECT_FLOAT_EQ(ed1(i).adj(), 0);
  for (int i = 0; i < rows_ed2 * cols_ed2; i++)
    EXPECT_FLOAT_EQ(ed2(i).adj(), 0);
  for (int i = 0; i < cols_ed3; i++)
    EXPECT_FLOAT_EQ(ed3(i).adj(), 0);
  EXPECT_FLOAT_EQ(ed4(0).adj(), 0);
  stan::math::set_zero_all_adjoints();

  for (int j = 0; j < size_vd; j++) {
    y(1 + j).grad();
    EXPECT_FLOAT_EQ(y(1 + j).val(), vd[j].val());
    EXPECT_FLOAT_EQ(vd[j].adj(), 1.0);
    EXPECT_FLOAT_EQ(d.adj(), 0.0);
    for (int i = 0; i < size_vd; i++)
      if (i != j) {
        EXPECT_FLOAT_EQ(vd[i].adj(), 0);
      }
    for (int i = 0; i < rows_ed1; i++)
      EXPECT_FLOAT_EQ(ed1(i).adj(), 0);
    for (int i = 0; i < rows_ed2 * cols_ed2; i++)
      EXPECT_FLOAT_EQ(ed2(i).adj(), 0);
    for (int i = 0; i < cols_ed3; i++)
      EXPECT_FLOAT_EQ(ed3(i).adj(), 0);
    EXPECT_FLOAT_EQ(ed4(0).adj(), 0);
    stan::math::set_zero_all_adjoints();
  }

  for (int j = 0; j < rows_ed1; j++) {
    y(1 + size_vd + j).grad();
    EXPECT_FLOAT_EQ(y(1 + size_vd + j).val(), ed1(j).val());
    EXPECT_FLOAT_EQ(ed1(j).adj(), 1.0);
    EXPECT_FLOAT_EQ(d.adj(), 0.0);
    for (int i = 0; i < size_vd; i++)
      EXPECT_FLOAT_EQ(vd[i].adj(), 0);
    for (int i = 0; i < rows_ed1; i++)
      if (i != j) {
        EXPECT_FLOAT_EQ(ed1(i).adj(), 0);
      }
    for (int i = 0; i < rows_ed2 * cols_ed2; i++)
      EXPECT_FLOAT_EQ(ed2(i).adj(), 0);
    for (int i = 0; i < cols_ed3; i++)
      EXPECT_FLOAT_EQ(ed3(i).adj(), 0);
    EXPECT_FLOAT_EQ(ed4(0).adj(), 0);
    stan::math::set_zero_all_adjoints();
  }

  for (int j = 0; j < rows_ed2 * cols_ed2; j++) {
    y(1 + size_vd + rows_ed1 + j).grad();
    EXPECT_FLOAT_EQ(y(1 + size_vd + rows_ed1 + j).val(), ed2(j).val());
    EXPECT_FLOAT_EQ(ed2(j).adj(), 1.0);
    EXPECT_FLOAT_EQ(d.adj(), 0.0);
    for (int i = 0; i < size_vd; i++)
      EXPECT_FLOAT_EQ(vd[i].adj(), 0);
    for (int i = 0; i < rows_ed1; i++)
      EXPECT_FLOAT_EQ(ed1(i).adj(), 0);
    for (int i = 0; i < rows_ed2 * cols_ed2; i++)
      if (i != j) {
        EXPECT_FLOAT_EQ(ed2(i).adj(), 0);
      }
    for (int i = 0; i < cols_ed3; i++)
      EXPECT_FLOAT_EQ(ed3(i).adj(), 0);
    EXPECT_FLOAT_EQ(ed4(0).adj(), 0);
    stan::math::set_zero_all_adjoints();
  }

  for (int j = 0; j < cols_ed3; j++) {
    y(1 + size_vd + rows_ed1 + rows_ed2 * cols_ed2 + j).grad();
    EXPECT_FLOAT_EQ(y(1 + size_vd + rows_ed1 + rows_ed2 * cols_ed2 + j).val(),
                    ed3(j).val());
    EXPECT_FLOAT_EQ(ed3(j).adj(), 1.0);
    EXPECT_FLOAT_EQ(d.adj(), 0.0);
    for (int i = 0; i < size_vd; i++)
      EXPECT_FLOAT_EQ(vd[i].adj(), 0);
    for (int i = 0; i < rows_ed1; i++)
      EXPECT_FLOAT_EQ(ed1(i).adj(), 0);
    for (int i = 0; i < rows_ed2 * cols_ed2; i++)
      EXPECT_FLOAT_EQ(ed2(i).adj(), 0);
    for (int i = 0; i < cols_ed3; i++)
      if (i != j) {
        EXPECT_FLOAT_EQ(ed3(i).adj(), 0);
      }
    EXPECT_FLOAT_EQ(ed4(0).adj(), 0);
    stan::math::set_zero_all_adjoints();
  }

  y(1 + size_vd + rows_ed1 + rows_ed2 * cols_ed2 + cols_ed3).grad();
  EXPECT_FLOAT_EQ(
      y(1 + size_vd + rows_ed1 + rows_ed2 * cols_ed2 + cols_ed3).val(),
      ed4(0).val());
  EXPECT_FLOAT_EQ(d.adj(), 0.0);
  for (int i = 0; i < size_vd; i++)
    EXPECT_FLOAT_EQ(vd[i].adj(), 0);
  for (int i = 0; i < rows_ed1; i++)
    EXPECT_FLOAT_EQ(ed1(i).adj(), 0);
  for (int i = 0; i < rows_ed2 * cols_ed2; i++)
    EXPECT_FLOAT_EQ(ed2(i).adj(), 0);
  for (int i = 0; i < cols_ed3; i++)
    EXPECT_FLOAT_EQ(ed3(i).adj(), 0);
  EXPECT_FLOAT_EQ(ed4(0).adj(), 1.0);
  stan::math::set_zero_all_adjoints();
}

TEST(AgradRev, test_pass_through_working_all_var_types_different_shapes) {
  int size_vd = 3, rows_ed1 = 7, rows_ed2 = 3, cols_ed2 = 5, cols_ed3 = 1;
  stan::math::var d = 1.0;
  std::vector<stan::math::var> vd(size_vd, 1.0);
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> ed1(rows_ed1);
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> ed2(rows_ed2,
                                                                     cols_ed2);
  std::vector<int> vi(5, 5);
  Eigen::Matrix<stan::math::var, 1, Eigen::Dynamic> ed3(cols_ed3);
  Eigen::Matrix<stan::math::var, 1, 1> ed4;

  for (int i = 0; i < size_vd; i++)
    vd[i] = 1.0;

  for (int i = 0; i < rows_ed1; i++)
    ed1(i) = 1.0;

  for (int i = 0; i < rows_ed2 * cols_ed2; i++)
    ed2(i) = 1.0;

  for (int i = 0; i < cols_ed3; i++)
    ed3(i) = 1.0;

  ed4(0) = 1.0;

  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> y
      = stan::math::adj_jac_apply<CheckAdjointsPassingThrough>(
          d, vd, 3, ed1, ed2, vi, ed3, ed4);

  y(0).grad();
  EXPECT_FLOAT_EQ(y(0).val(), d.val());
  EXPECT_FLOAT_EQ(d.adj(), 1.0);
  for (int i = 0; i < size_vd; i++)
    EXPECT_FLOAT_EQ(vd[i].adj(), 0);
  for (int i = 0; i < rows_ed1; i++)
    EXPECT_FLOAT_EQ(ed1(i).adj(), 0);
  for (int i = 0; i < rows_ed2 * cols_ed2; i++)
    EXPECT_FLOAT_EQ(ed2(i).adj(), 0);
  for (int i = 0; i < cols_ed3; i++)
    EXPECT_FLOAT_EQ(ed3(i).adj(), 0);
  EXPECT_FLOAT_EQ(ed4(0).adj(), 0);
  stan::math::set_zero_all_adjoints();

  for (int j = 0; j < size_vd; j++) {
    y(1 + j).grad();
    EXPECT_FLOAT_EQ(y(1 + j).val(), vd[j].val());
    EXPECT_FLOAT_EQ(vd[j].adj(), 1.0);
    EXPECT_FLOAT_EQ(d.adj(), 0.0);
    for (int i = 0; i < size_vd; i++)
      if (i != j) {
        EXPECT_FLOAT_EQ(vd[i].adj(), 0);
      }
    for (int i = 0; i < rows_ed1; i++)
      EXPECT_FLOAT_EQ(ed1(i).adj(), 0);
    for (int i = 0; i < rows_ed2 * cols_ed2; i++)
      EXPECT_FLOAT_EQ(ed2(i).adj(), 0);
    for (int i = 0; i < cols_ed3; i++)
      EXPECT_FLOAT_EQ(ed3(i).adj(), 0);
    EXPECT_FLOAT_EQ(ed4(0).adj(), 0);
    stan::math::set_zero_all_adjoints();
  }

  for (int j = 0; j < rows_ed1; j++) {
    y(1 + size_vd + j).grad();
    EXPECT_FLOAT_EQ(y(1 + size_vd + j).val(), ed1(j).val());
    EXPECT_FLOAT_EQ(ed1(j).adj(), 1.0);
    EXPECT_FLOAT_EQ(d.adj(), 0.0);
    for (int i = 0; i < size_vd; i++)
      EXPECT_FLOAT_EQ(vd[i].adj(), 0);
    for (int i = 0; i < rows_ed1; i++)
      if (i != j) {
        EXPECT_FLOAT_EQ(ed1(i).adj(), 0);
      }
    for (int i = 0; i < rows_ed2 * cols_ed2; i++)
      EXPECT_FLOAT_EQ(ed2(i).adj(), 0);
    for (int i = 0; i < cols_ed3; i++)
      EXPECT_FLOAT_EQ(ed3(i).adj(), 0);
    EXPECT_FLOAT_EQ(ed4(0).adj(), 0);
    stan::math::set_zero_all_adjoints();
  }

  for (int j = 0; j < rows_ed2 * cols_ed2; j++) {
    y(1 + size_vd + rows_ed1 + j).grad();
    EXPECT_FLOAT_EQ(y(1 + size_vd + rows_ed1 + j).val(), ed2(j).val());
    EXPECT_FLOAT_EQ(ed2(j).adj(), 1.0);
    EXPECT_FLOAT_EQ(d.adj(), 0.0);
    for (int i = 0; i < size_vd; i++)
      EXPECT_FLOAT_EQ(vd[i].adj(), 0);
    for (int i = 0; i < rows_ed1; i++)
      EXPECT_FLOAT_EQ(ed1(i).adj(), 0);
    for (int i = 0; i < rows_ed2 * cols_ed2; i++)
      if (i != j) {
        EXPECT_FLOAT_EQ(ed2(i).adj(), 0);
      }
    for (int i = 0; i < cols_ed3; i++)
      EXPECT_FLOAT_EQ(ed3(i).adj(), 0);
    EXPECT_FLOAT_EQ(ed4(0).adj(), 0);
    stan::math::set_zero_all_adjoints();
  }

  for (int j = 0; j < cols_ed3; j++) {
    y(1 + size_vd + rows_ed1 + rows_ed2 * cols_ed2 + j).grad();
    EXPECT_FLOAT_EQ(y(1 + size_vd + rows_ed1 + rows_ed2 * cols_ed2 + j).val(),
                    ed3(j).val());
    EXPECT_FLOAT_EQ(ed3(j).adj(), 1.0);
    EXPECT_FLOAT_EQ(d.adj(), 0.0);
    for (int i = 0; i < size_vd; i++)
      EXPECT_FLOAT_EQ(vd[i].adj(), 0);
    for (int i = 0; i < rows_ed1; i++)
      EXPECT_FLOAT_EQ(ed1(i).adj(), 0);
    for (int i = 0; i < rows_ed2 * cols_ed2; i++)
      EXPECT_FLOAT_EQ(ed2(i).adj(), 0);
    for (int i = 0; i < cols_ed3; i++)
      if (i != j) {
        EXPECT_FLOAT_EQ(ed3(i).adj(), 0);
      }
    EXPECT_FLOAT_EQ(ed4(0).adj(), 0);
    stan::math::set_zero_all_adjoints();
  }

  y(1 + size_vd + rows_ed1 + rows_ed2 * cols_ed2 + cols_ed3).grad();
  EXPECT_FLOAT_EQ(
      y(1 + size_vd + rows_ed1 + rows_ed2 * cols_ed2 + cols_ed3).val(),
      ed4(0).val());
  EXPECT_FLOAT_EQ(d.adj(), 0.0);
  for (int i = 0; i < size_vd; i++)
    EXPECT_FLOAT_EQ(vd[i].adj(), 0);
  for (int i = 0; i < rows_ed1; i++)
    EXPECT_FLOAT_EQ(ed1(i).adj(), 0);
  for (int i = 0; i < rows_ed2 * cols_ed2; i++)
    EXPECT_FLOAT_EQ(ed2(i).adj(), 0);
  for (int i = 0; i < cols_ed3; i++)
    EXPECT_FLOAT_EQ(ed3(i).adj(), 0);
  EXPECT_FLOAT_EQ(ed4(0).adj(), 1.0);
  stan::math::set_zero_all_adjoints();
}

TEST(AgradRev, test_pass_through_working_all_var_types_double_test_1) {
  int size_vd = 3, rows_ed1 = 7, rows_ed2 = 3, cols_ed2 = 5, cols_ed3 = 1;
  double d = 1.0;
  std::vector<stan::math::var> vd(size_vd, 1.0);
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> ed1(rows_ed1);
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> ed2(rows_ed2,
                                                                     cols_ed2);
  std::vector<int> vi(5, 5);
  Eigen::Matrix<stan::math::var, 1, Eigen::Dynamic> ed3(cols_ed3);
  Eigen::Matrix<stan::math::var, 1, 1> ed4;

  for (int i = 0; i < size_vd; i++)
    vd[i] = 1.0;

  for (int i = 0; i < rows_ed1; i++)
    ed1(i) = 1.0;

  for (int i = 0; i < rows_ed2 * cols_ed2; i++)
    ed2(i) = 1.0;

  for (int i = 0; i < cols_ed3; i++)
    ed3(i) = 1.0;

  ed4(0) = 1.0;

  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> y
      = stan::math::adj_jac_apply<CheckAdjointsPassingThrough>(
          d, vd, 3, ed1, ed2, vi, ed3, ed4);

  y(0).grad();
  EXPECT_FLOAT_EQ(y(0).val(), d);
  for (int i = 0; i < size_vd; i++)
    EXPECT_FLOAT_EQ(vd[i].adj(), 0);
  for (int i = 0; i < rows_ed1; i++)
    EXPECT_FLOAT_EQ(ed1(i).adj(), 0);
  for (int i = 0; i < rows_ed2 * cols_ed2; i++)
    EXPECT_FLOAT_EQ(ed2(i).adj(), 0);
  for (int i = 0; i < cols_ed3; i++)
    EXPECT_FLOAT_EQ(ed3(i).adj(), 0);
  EXPECT_FLOAT_EQ(ed4(0).adj(), 0);
  stan::math::set_zero_all_adjoints();

  for (int j = 0; j < size_vd; j++) {
    y(1 + j).grad();
    EXPECT_FLOAT_EQ(y(1 + j).val(), vd[j].val());
    EXPECT_FLOAT_EQ(vd[j].adj(), 1.0);
    for (int i = 0; i < size_vd; i++)
      if (i != j) {
        EXPECT_FLOAT_EQ(vd[i].adj(), 0);
      }
    for (int i = 0; i < rows_ed1; i++)
      EXPECT_FLOAT_EQ(ed1(i).adj(), 0);
    for (int i = 0; i < rows_ed2 * cols_ed2; i++)
      EXPECT_FLOAT_EQ(ed2(i).adj(), 0);
    for (int i = 0; i < cols_ed3; i++)
      EXPECT_FLOAT_EQ(ed3(i).adj(), 0);
    EXPECT_FLOAT_EQ(ed4(0).adj(), 0);
    stan::math::set_zero_all_adjoints();
  }

  for (int j = 0; j < rows_ed1; j++) {
    y(1 + size_vd + j).grad();
    EXPECT_FLOAT_EQ(y(1 + size_vd + j).val(), ed1(j).val());
    EXPECT_FLOAT_EQ(ed1(j).adj(), 1.0);
    for (int i = 0; i < size_vd; i++)
      EXPECT_FLOAT_EQ(vd[i].adj(), 0);
    for (int i = 0; i < rows_ed1; i++)
      if (i != j) {
        EXPECT_FLOAT_EQ(ed1(i).adj(), 0);
      }
    for (int i = 0; i < rows_ed2 * cols_ed2; i++)
      EXPECT_FLOAT_EQ(ed2(i).adj(), 0);
    for (int i = 0; i < cols_ed3; i++)
      EXPECT_FLOAT_EQ(ed3(i).adj(), 0);
    EXPECT_FLOAT_EQ(ed4(0).adj(), 0);
    stan::math::set_zero_all_adjoints();
  }

  for (int j = 0; j < rows_ed2 * cols_ed2; j++) {
    y(1 + size_vd + rows_ed1 + j).grad();
    EXPECT_FLOAT_EQ(y(1 + size_vd + rows_ed1 + j).val(), ed2(j).val());
    EXPECT_FLOAT_EQ(ed2(j).adj(), 1.0);
    for (int i = 0; i < size_vd; i++)
      EXPECT_FLOAT_EQ(vd[i].adj(), 0);
    for (int i = 0; i < rows_ed1; i++)
      EXPECT_FLOAT_EQ(ed1(i).adj(), 0);
    for (int i = 0; i < rows_ed2 * cols_ed2; i++)
      if (i != j) {
        EXPECT_FLOAT_EQ(ed2(i).adj(), 0);
      }
    for (int i = 0; i < cols_ed3; i++)
      EXPECT_FLOAT_EQ(ed3(i).adj(), 0);
    EXPECT_FLOAT_EQ(ed4(0).adj(), 0);
    stan::math::set_zero_all_adjoints();
  }

  for (int j = 0; j < cols_ed3; j++) {
    y(1 + size_vd + rows_ed1 + rows_ed2 * cols_ed2 + j).grad();
    EXPECT_FLOAT_EQ(y(1 + size_vd + rows_ed1 + rows_ed2 * cols_ed2 + j).val(),
                    ed3(j).val());
    EXPECT_FLOAT_EQ(ed3(j).adj(), 1.0);
    for (int i = 0; i < size_vd; i++)
      EXPECT_FLOAT_EQ(vd[i].adj(), 0);
    for (int i = 0; i < rows_ed1; i++)
      EXPECT_FLOAT_EQ(ed1(i).adj(), 0);
    for (int i = 0; i < rows_ed2 * cols_ed2; i++)
      EXPECT_FLOAT_EQ(ed2(i).adj(), 0);
    for (int i = 0; i < cols_ed3; i++)
      if (i != j) {
        EXPECT_FLOAT_EQ(ed3(i).adj(), 0);
      }
    EXPECT_FLOAT_EQ(ed4(0).adj(), 0);
    stan::math::set_zero_all_adjoints();
  }

  y(1 + size_vd + rows_ed1 + rows_ed2 * cols_ed2 + cols_ed3).grad();
  EXPECT_FLOAT_EQ(
      y(1 + size_vd + rows_ed1 + rows_ed2 * cols_ed2 + cols_ed3).val(),
      ed4(0).val());
  for (int i = 0; i < size_vd; i++)
    EXPECT_FLOAT_EQ(vd[i].adj(), 0);
  for (int i = 0; i < rows_ed1; i++)
    EXPECT_FLOAT_EQ(ed1(i).adj(), 0);
  for (int i = 0; i < rows_ed2 * cols_ed2; i++)
    EXPECT_FLOAT_EQ(ed2(i).adj(), 0);
  for (int i = 0; i < cols_ed3; i++)
    EXPECT_FLOAT_EQ(ed3(i).adj(), 0);
  EXPECT_FLOAT_EQ(ed4(0).adj(), 1.0);
  stan::math::set_zero_all_adjoints();
}

TEST(AgradRev, test_pass_through_working_all_var_types_double_test_2) {
  int size_vd = 3, rows_ed1 = 7, rows_ed2 = 3, cols_ed2 = 5, cols_ed3 = 1;
  stan::math::var d = 1.0;
  std::vector<double> vd(size_vd, 1.0);
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> ed1(rows_ed1);
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> ed2(rows_ed2,
                                                                     cols_ed2);
  std::vector<int> vi(5, 5);
  Eigen::Matrix<stan::math::var, 1, Eigen::Dynamic> ed3(cols_ed3);
  Eigen::Matrix<stan::math::var, 1, 1> ed4;

  for (int i = 0; i < size_vd; i++)
    vd[i] = 1.0;

  for (int i = 0; i < rows_ed1; i++)
    ed1(i) = 1.0;

  for (int i = 0; i < rows_ed2 * cols_ed2; i++)
    ed2(i) = 1.0;

  for (int i = 0; i < cols_ed3; i++)
    ed3(i) = 1.0;

  ed4(0) = 1.0;

  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> y
      = stan::math::adj_jac_apply<CheckAdjointsPassingThrough>(
          d, vd, 3, ed1, ed2, vi, ed3, ed4);

  y(0).grad();
  EXPECT_FLOAT_EQ(y(0).val(), d.val());
  EXPECT_FLOAT_EQ(d.adj(), 1.0);
  for (int i = 0; i < rows_ed1; i++)
    EXPECT_FLOAT_EQ(ed1(i).adj(), 0);
  for (int i = 0; i < rows_ed2 * cols_ed2; i++)
    EXPECT_FLOAT_EQ(ed2(i).adj(), 0);
  for (int i = 0; i < cols_ed3; i++)
    EXPECT_FLOAT_EQ(ed3(i).adj(), 0);
  EXPECT_FLOAT_EQ(ed4(0).adj(), 0);
  stan::math::set_zero_all_adjoints();

  for (int j = 0; j < size_vd; j++) {
    y(1 + j).grad();
    EXPECT_FLOAT_EQ(y(1 + j).val(), vd[j]);
    for (int i = 0; i < rows_ed1; i++)
      EXPECT_FLOAT_EQ(ed1(i).adj(), 0);
    for (int i = 0; i < rows_ed2 * cols_ed2; i++)
      EXPECT_FLOAT_EQ(ed2(i).adj(), 0);
    for (int i = 0; i < cols_ed3; i++)
      EXPECT_FLOAT_EQ(ed3(i).adj(), 0);
    EXPECT_FLOAT_EQ(ed4(0).adj(), 0);
    stan::math::set_zero_all_adjoints();
  }

  for (int j = 0; j < rows_ed1; j++) {
    y(1 + size_vd + j).grad();
    EXPECT_FLOAT_EQ(y(1 + size_vd + j).val(), ed1(j).val());
    EXPECT_FLOAT_EQ(ed1(j).adj(), 1.0);
    for (int i = 0; i < rows_ed1; i++)
      if (i != j) {
        EXPECT_FLOAT_EQ(ed1(i).adj(), 0);
      }
    for (int i = 0; i < rows_ed2 * cols_ed2; i++)
      EXPECT_FLOAT_EQ(ed2(i).adj(), 0);
    for (int i = 0; i < cols_ed3; i++)
      EXPECT_FLOAT_EQ(ed3(i).adj(), 0);
    EXPECT_FLOAT_EQ(ed4(0).adj(), 0);
    stan::math::set_zero_all_adjoints();
  }

  for (int j = 0; j < rows_ed2 * cols_ed2; j++) {
    y(1 + size_vd + rows_ed1 + j).grad();
    EXPECT_FLOAT_EQ(y(1 + size_vd + rows_ed1 + j).val(), ed2(j).val());
    EXPECT_FLOAT_EQ(ed2(j).adj(), 1.0);
    for (int i = 0; i < rows_ed1; i++)
      EXPECT_FLOAT_EQ(ed1(i).adj(), 0);
    for (int i = 0; i < rows_ed2 * cols_ed2; i++)
      if (i != j) {
        EXPECT_FLOAT_EQ(ed2(i).adj(), 0);
      }
    for (int i = 0; i < cols_ed3; i++)
      EXPECT_FLOAT_EQ(ed3(i).adj(), 0);
    EXPECT_FLOAT_EQ(ed4(0).adj(), 0);
    stan::math::set_zero_all_adjoints();
  }

  for (int j = 0; j < cols_ed3; j++) {
    y(1 + size_vd + rows_ed1 + rows_ed2 * cols_ed2 + j).grad();
    EXPECT_FLOAT_EQ(y(1 + size_vd + rows_ed1 + rows_ed2 * cols_ed2 + j).val(),
                    ed3(j).val());
    EXPECT_FLOAT_EQ(ed3(j).adj(), 1.0);
    for (int i = 0; i < rows_ed1; i++)
      EXPECT_FLOAT_EQ(ed1(i).adj(), 0);
    for (int i = 0; i < rows_ed2 * cols_ed2; i++)
      EXPECT_FLOAT_EQ(ed2(i).adj(), 0);
    for (int i = 0; i < cols_ed3; i++)
      if (i != j) {
        EXPECT_FLOAT_EQ(ed3(i).adj(), 0);
      }
    EXPECT_FLOAT_EQ(ed4(0).adj(), 0);
    stan::math::set_zero_all_adjoints();
  }

  y(1 + size_vd + rows_ed1 + rows_ed2 * cols_ed2 + cols_ed3).grad();
  EXPECT_FLOAT_EQ(
      y(1 + size_vd + rows_ed1 + rows_ed2 * cols_ed2 + cols_ed3).val(),
      ed4(0).val());
  for (int i = 0; i < rows_ed1; i++)
    EXPECT_FLOAT_EQ(ed1(i).adj(), 0);
  for (int i = 0; i < rows_ed2 * cols_ed2; i++)
    EXPECT_FLOAT_EQ(ed2(i).adj(), 0);
  for (int i = 0; i < cols_ed3; i++)
    EXPECT_FLOAT_EQ(ed3(i).adj(), 0);
  EXPECT_FLOAT_EQ(ed4(0).adj(), 1.0);
  stan::math::set_zero_all_adjoints();
}

TEST(AgradRev, test_pass_through_working_all_var_types_double_test_3) {
  int size_vd = 3, rows_ed1 = 7, rows_ed2 = 3, cols_ed2 = 5, cols_ed3 = 1;
  stan::math::var d = 1.0;
  std::vector<stan::math::var> vd(size_vd, 1.0);
  Eigen::Matrix<double, Eigen::Dynamic, 1> ed1(rows_ed1);
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> ed2(rows_ed2,
                                                                     cols_ed2);
  std::vector<int> vi(5, 5);
  Eigen::Matrix<stan::math::var, 1, Eigen::Dynamic> ed3(cols_ed3);
  Eigen::Matrix<stan::math::var, 1, 1> ed4;

  for (int i = 0; i < size_vd; i++)
    vd[i] = 1.0;

  for (int i = 0; i < rows_ed1; i++)
    ed1(i) = 1.0;

  for (int i = 0; i < rows_ed2 * cols_ed2; i++)
    ed2(i) = 1.0;

  for (int i = 0; i < cols_ed3; i++)
    ed3(i) = 1.0;

  ed4(0) = 1.0;

  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> y
      = stan::math::adj_jac_apply<CheckAdjointsPassingThrough>(
          d, vd, 3, ed1, ed2, vi, ed3, ed4);

  y(0).grad();
  EXPECT_FLOAT_EQ(y(0).val(), d.val());
  EXPECT_FLOAT_EQ(d.adj(), 1.0);
  for (int i = 0; i < size_vd; i++)
    EXPECT_FLOAT_EQ(vd[i].adj(), 0);
  for (int i = 0; i < rows_ed2 * cols_ed2; i++)
    EXPECT_FLOAT_EQ(ed2(i).adj(), 0);
  for (int i = 0; i < cols_ed3; i++)
    EXPECT_FLOAT_EQ(ed3(i).adj(), 0);
  EXPECT_FLOAT_EQ(ed4(0).adj(), 0);
  stan::math::set_zero_all_adjoints();

  for (int j = 0; j < size_vd; j++) {
    y(1 + j).grad();
    EXPECT_FLOAT_EQ(y(1 + j).val(), vd[j].val());
    EXPECT_FLOAT_EQ(vd[j].adj(), 1.0);
    EXPECT_FLOAT_EQ(d.adj(), 0.0);
    for (int i = 0; i < size_vd; i++)
      if (i != j) {
        EXPECT_FLOAT_EQ(vd[i].adj(), 0);
      }
    for (int i = 0; i < rows_ed2 * cols_ed2; i++)
      EXPECT_FLOAT_EQ(ed2(i).adj(), 0);
    for (int i = 0; i < cols_ed3; i++)
      EXPECT_FLOAT_EQ(ed3(i).adj(), 0);
    EXPECT_FLOAT_EQ(ed4(0).adj(), 0);
    stan::math::set_zero_all_adjoints();
  }

  for (int j = 0; j < rows_ed1; j++) {
    y(1 + size_vd + j).grad();
    EXPECT_FLOAT_EQ(y(1 + size_vd + j).val(), ed1(j));
    EXPECT_FLOAT_EQ(d.adj(), 0.0);
    for (int i = 0; i < size_vd; i++)
      EXPECT_FLOAT_EQ(vd[i].adj(), 0);
    for (int i = 0; i < rows_ed2 * cols_ed2; i++)
      EXPECT_FLOAT_EQ(ed2(i).adj(), 0);
    for (int i = 0; i < cols_ed3; i++)
      EXPECT_FLOAT_EQ(ed3(i).adj(), 0);
    EXPECT_FLOAT_EQ(ed4(0).adj(), 0);
    stan::math::set_zero_all_adjoints();
  }

  for (int j = 0; j < rows_ed2 * cols_ed2; j++) {
    y(1 + size_vd + rows_ed1 + j).grad();
    EXPECT_FLOAT_EQ(y(1 + size_vd + rows_ed1 + j).val(), ed2(j).val());
    EXPECT_FLOAT_EQ(ed2(j).adj(), 1.0);
    EXPECT_FLOAT_EQ(d.adj(), 0.0);
    for (int i = 0; i < size_vd; i++)
      EXPECT_FLOAT_EQ(vd[i].adj(), 0);
    for (int i = 0; i < rows_ed2 * cols_ed2; i++)
      if (i != j) {
        EXPECT_FLOAT_EQ(ed2(i).adj(), 0);
      }
    for (int i = 0; i < cols_ed3; i++)
      EXPECT_FLOAT_EQ(ed3(i).adj(), 0);
    EXPECT_FLOAT_EQ(ed4(0).adj(), 0);
    stan::math::set_zero_all_adjoints();
  }

  for (int j = 0; j < cols_ed3; j++) {
    y(1 + size_vd + rows_ed1 + rows_ed2 * cols_ed2 + j).grad();
    EXPECT_FLOAT_EQ(y(1 + size_vd + rows_ed1 + rows_ed2 * cols_ed2 + j).val(),
                    ed3(j).val());
    EXPECT_FLOAT_EQ(ed3(j).adj(), 1.0);
    EXPECT_FLOAT_EQ(d.adj(), 0.0);
    for (int i = 0; i < size_vd; i++)
      EXPECT_FLOAT_EQ(vd[i].adj(), 0);
    for (int i = 0; i < rows_ed2 * cols_ed2; i++)
      EXPECT_FLOAT_EQ(ed2(i).adj(), 0);
    for (int i = 0; i < cols_ed3; i++)
      if (i != j) {
        EXPECT_FLOAT_EQ(ed3(i).adj(), 0);
      }
    EXPECT_FLOAT_EQ(ed4(0).adj(), 0);
    stan::math::set_zero_all_adjoints();
  }

  y(1 + size_vd + rows_ed1 + rows_ed2 * cols_ed2 + cols_ed3).grad();
  EXPECT_FLOAT_EQ(
      y(1 + size_vd + rows_ed1 + rows_ed2 * cols_ed2 + cols_ed3).val(),
      ed4(0).val());
  EXPECT_FLOAT_EQ(d.adj(), 0.0);
  for (int i = 0; i < size_vd; i++)
    EXPECT_FLOAT_EQ(vd[i].adj(), 0);
  for (int i = 0; i < rows_ed2 * cols_ed2; i++)
    EXPECT_FLOAT_EQ(ed2(i).adj(), 0);
  for (int i = 0; i < cols_ed3; i++)
    EXPECT_FLOAT_EQ(ed3(i).adj(), 0);
  EXPECT_FLOAT_EQ(ed4(0).adj(), 1.0);
  stan::math::set_zero_all_adjoints();
}

TEST(AgradRev, test_pass_through_working_all_var_types_double_test_4) {
  int size_vd = 3, rows_ed1 = 7, rows_ed2 = 3, cols_ed2 = 5, cols_ed3 = 1;
  stan::math::var d = 1.0;
  std::vector<stan::math::var> vd(size_vd, 1.0);
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> ed1(rows_ed1);
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> ed2(rows_ed2, cols_ed2);
  std::vector<int> vi(5, 5);
  Eigen::Matrix<stan::math::var, 1, Eigen::Dynamic> ed3(cols_ed3);
  Eigen::Matrix<stan::math::var, 1, 1> ed4;

  for (int i = 0; i < size_vd; i++)
    vd[i] = 1.0;

  for (int i = 0; i < rows_ed1; i++)
    ed1(i) = 1.0;

  for (int i = 0; i < rows_ed2 * cols_ed2; i++)
    ed2(i) = 1.0;

  for (int i = 0; i < cols_ed3; i++)
    ed3(i) = 1.0;

  ed4(0) = 1.0;

  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> y
      = stan::math::adj_jac_apply<CheckAdjointsPassingThrough>(
          d, vd, 3, ed1, ed2, vi, ed3, ed4);

  y(0).grad();
  EXPECT_FLOAT_EQ(y(0).val(), d.val());
  EXPECT_FLOAT_EQ(d.adj(), 1.0);
  for (int i = 0; i < size_vd; i++)
    EXPECT_FLOAT_EQ(vd[i].adj(), 0);
  for (int i = 0; i < rows_ed1; i++)
    EXPECT_FLOAT_EQ(ed1(i).adj(), 0);
  for (int i = 0; i < cols_ed3; i++)
    EXPECT_FLOAT_EQ(ed3(i).adj(), 0);
  EXPECT_FLOAT_EQ(ed4(0).adj(), 0);
  stan::math::set_zero_all_adjoints();

  for (int j = 0; j < size_vd; j++) {
    y(1 + j).grad();
    EXPECT_FLOAT_EQ(y(1 + j).val(), vd[j].val());
    EXPECT_FLOAT_EQ(vd[j].adj(), 1.0);
    EXPECT_FLOAT_EQ(d.adj(), 0.0);
    for (int i = 0; i < size_vd; i++)
      if (i != j) {
        EXPECT_FLOAT_EQ(vd[i].adj(), 0);
      }
    for (int i = 0; i < rows_ed1; i++)
      EXPECT_FLOAT_EQ(ed1(i).adj(), 0);
    for (int i = 0; i < cols_ed3; i++)
      EXPECT_FLOAT_EQ(ed3(i).adj(), 0);
    EXPECT_FLOAT_EQ(ed4(0).adj(), 0);
    stan::math::set_zero_all_adjoints();
  }

  for (int j = 0; j < rows_ed1; j++) {
    y(1 + size_vd + j).grad();
    EXPECT_FLOAT_EQ(y(1 + size_vd + j).val(), ed1(j).val());
    EXPECT_FLOAT_EQ(ed1(j).adj(), 1.0);
    EXPECT_FLOAT_EQ(d.adj(), 0.0);
    for (int i = 0; i < size_vd; i++)
      EXPECT_FLOAT_EQ(vd[i].adj(), 0);
    for (int i = 0; i < rows_ed1; i++)
      if (i != j) {
        EXPECT_FLOAT_EQ(ed1(i).adj(), 0);
      }
    for (int i = 0; i < cols_ed3; i++)
      EXPECT_FLOAT_EQ(ed3(i).adj(), 0);
    EXPECT_FLOAT_EQ(ed4(0).adj(), 0);
    stan::math::set_zero_all_adjoints();
  }

  for (int j = 0; j < rows_ed2 * cols_ed2; j++) {
    y(1 + size_vd + rows_ed1 + j).grad();
    EXPECT_FLOAT_EQ(y(1 + size_vd + rows_ed1 + j).val(), ed2(j));
    EXPECT_FLOAT_EQ(d.adj(), 0.0);
    for (int i = 0; i < size_vd; i++)
      EXPECT_FLOAT_EQ(vd[i].adj(), 0);
    for (int i = 0; i < rows_ed1; i++)
      EXPECT_FLOAT_EQ(ed1(i).adj(), 0);
    for (int i = 0; i < cols_ed3; i++)
      EXPECT_FLOAT_EQ(ed3(i).adj(), 0);
    EXPECT_FLOAT_EQ(ed4(0).adj(), 0);
    stan::math::set_zero_all_adjoints();
  }

  for (int j = 0; j < cols_ed3; j++) {
    y(1 + size_vd + rows_ed1 + rows_ed2 * cols_ed2 + j).grad();
    EXPECT_FLOAT_EQ(y(1 + size_vd + rows_ed1 + rows_ed2 * cols_ed2 + j).val(),
                    ed3(j).val());
    EXPECT_FLOAT_EQ(ed3(j).adj(), 1.0);
    EXPECT_FLOAT_EQ(d.adj(), 0.0);
    for (int i = 0; i < size_vd; i++)
      EXPECT_FLOAT_EQ(vd[i].adj(), 0);
    for (int i = 0; i < rows_ed1; i++)
      EXPECT_FLOAT_EQ(ed1(i).adj(), 0);
    for (int i = 0; i < cols_ed3; i++)
      if (i != j) {
        EXPECT_FLOAT_EQ(ed3(i).adj(), 0);
      }
    EXPECT_FLOAT_EQ(ed4(0).adj(), 0);
    stan::math::set_zero_all_adjoints();
  }

  y(1 + size_vd + rows_ed1 + rows_ed2 * cols_ed2 + cols_ed3).grad();
  EXPECT_FLOAT_EQ(
      y(1 + size_vd + rows_ed1 + rows_ed2 * cols_ed2 + cols_ed3).val(),
      ed4(0).val());
  EXPECT_FLOAT_EQ(d.adj(), 0.0);
  for (int i = 0; i < size_vd; i++)
    EXPECT_FLOAT_EQ(vd[i].adj(), 0);
  for (int i = 0; i < rows_ed1; i++)
    EXPECT_FLOAT_EQ(ed1(i).adj(), 0);
  for (int i = 0; i < cols_ed3; i++)
    EXPECT_FLOAT_EQ(ed3(i).adj(), 0);
  EXPECT_FLOAT_EQ(ed4(0).adj(), 1.0);
  stan::math::set_zero_all_adjoints();
}

TEST(AgradRev, test_pass_through_working_all_var_types_double_test_5) {
  int size_vd = 3, rows_ed1 = 7, rows_ed2 = 3, cols_ed2 = 5, cols_ed3 = 1;
  stan::math::var d = 1.0;
  std::vector<stan::math::var> vd(size_vd, 1.0);
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> ed1(rows_ed1);
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> ed2(rows_ed2,
                                                                     cols_ed2);
  std::vector<int> vi(5, 5);
  Eigen::Matrix<double, 1, Eigen::Dynamic> ed3(cols_ed3);
  Eigen::Matrix<stan::math::var, 1, 1> ed4;

  for (int i = 0; i < size_vd; i++)
    vd[i] = 1.0;

  for (int i = 0; i < rows_ed1; i++)
    ed1(i) = 1.0;

  for (int i = 0; i < rows_ed2 * cols_ed2; i++)
    ed2(i) = 1.0;

  for (int i = 0; i < cols_ed3; i++)
    ed3(i) = 1.0;

  ed4(0) = 1.0;

  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> y
      = stan::math::adj_jac_apply<CheckAdjointsPassingThrough>(
          d, vd, 3, ed1, ed2, vi, ed3, ed4);

  y(0).grad();
  EXPECT_FLOAT_EQ(y(0).val(), d.val());
  EXPECT_FLOAT_EQ(d.adj(), 1.0);
  for (int i = 0; i < size_vd; i++)
    EXPECT_FLOAT_EQ(vd[i].adj(), 0);
  for (int i = 0; i < rows_ed1; i++)
    EXPECT_FLOAT_EQ(ed1(i).adj(), 0);
  for (int i = 0; i < rows_ed2 * cols_ed2; i++)
    EXPECT_FLOAT_EQ(ed2(i).adj(), 0);
  EXPECT_FLOAT_EQ(ed4(0).adj(), 0);
  stan::math::set_zero_all_adjoints();

  for (int j = 0; j < size_vd; j++) {
    y(1 + j).grad();
    EXPECT_FLOAT_EQ(y(1 + j).val(), vd[j].val());
    EXPECT_FLOAT_EQ(vd[j].adj(), 1.0);
    EXPECT_FLOAT_EQ(d.adj(), 0.0);
    for (int i = 0; i < size_vd; i++)
      if (i != j) {
        EXPECT_FLOAT_EQ(vd[i].adj(), 0);
      }
    for (int i = 0; i < rows_ed1; i++)
      EXPECT_FLOAT_EQ(ed1(i).adj(), 0);
    for (int i = 0; i < rows_ed2 * cols_ed2; i++)
      EXPECT_FLOAT_EQ(ed2(i).adj(), 0);
    EXPECT_FLOAT_EQ(ed4(0).adj(), 0);
    stan::math::set_zero_all_adjoints();
  }

  for (int j = 0; j < rows_ed1; j++) {
    y(1 + size_vd + j).grad();
    EXPECT_FLOAT_EQ(y(1 + size_vd + j).val(), ed1(j).val());
    EXPECT_FLOAT_EQ(ed1(j).adj(), 1.0);
    EXPECT_FLOAT_EQ(d.adj(), 0.0);
    for (int i = 0; i < size_vd; i++)
      EXPECT_FLOAT_EQ(vd[i].adj(), 0);
    for (int i = 0; i < rows_ed1; i++)
      if (i != j) {
        EXPECT_FLOAT_EQ(ed1(i).adj(), 0);
      }
    for (int i = 0; i < rows_ed2 * cols_ed2; i++)
      EXPECT_FLOAT_EQ(ed2(i).adj(), 0);
    EXPECT_FLOAT_EQ(ed4(0).adj(), 0);
    stan::math::set_zero_all_adjoints();
  }

  for (int j = 0; j < rows_ed2 * cols_ed2; j++) {
    y(1 + size_vd + rows_ed1 + j).grad();
    EXPECT_FLOAT_EQ(y(1 + size_vd + rows_ed1 + j).val(), ed2(j).val());
    EXPECT_FLOAT_EQ(ed2(j).adj(), 1.0);
    EXPECT_FLOAT_EQ(d.adj(), 0.0);
    for (int i = 0; i < size_vd; i++)
      EXPECT_FLOAT_EQ(vd[i].adj(), 0);
    for (int i = 0; i < rows_ed1; i++)
      EXPECT_FLOAT_EQ(ed1(i).adj(), 0);
    for (int i = 0; i < rows_ed2 * cols_ed2; i++)
      if (i != j) {
        EXPECT_FLOAT_EQ(ed2(i).adj(), 0);
      }
    EXPECT_FLOAT_EQ(ed4(0).adj(), 0);
    stan::math::set_zero_all_adjoints();
  }

  for (int j = 0; j < cols_ed3; j++) {
    y(1 + size_vd + rows_ed1 + rows_ed2 * cols_ed2 + j).grad();
    EXPECT_FLOAT_EQ(y(1 + size_vd + rows_ed1 + rows_ed2 * cols_ed2 + j).val(),
                    ed3(j));
    EXPECT_FLOAT_EQ(d.adj(), 0.0);
    for (int i = 0; i < size_vd; i++)
      EXPECT_FLOAT_EQ(vd[i].adj(), 0);
    for (int i = 0; i < rows_ed1; i++)
      EXPECT_FLOAT_EQ(ed1(i).adj(), 0);
    for (int i = 0; i < rows_ed2 * cols_ed2; i++)
      EXPECT_FLOAT_EQ(ed2(i).adj(), 0);
    EXPECT_FLOAT_EQ(ed4(0).adj(), 0);
    stan::math::set_zero_all_adjoints();
  }

  y(1 + size_vd + rows_ed1 + rows_ed2 * cols_ed2 + cols_ed3).grad();
  EXPECT_FLOAT_EQ(
      y(1 + size_vd + rows_ed1 + rows_ed2 * cols_ed2 + cols_ed3).val(),
      ed4(0).val());
  EXPECT_FLOAT_EQ(d.adj(), 0.0);
  for (int i = 0; i < size_vd; i++)
    EXPECT_FLOAT_EQ(vd[i].adj(), 0);
  for (int i = 0; i < rows_ed1; i++)
    EXPECT_FLOAT_EQ(ed1(i).adj(), 0);
  for (int i = 0; i < rows_ed2 * cols_ed2; i++)
    EXPECT_FLOAT_EQ(ed2(i).adj(), 0);
  EXPECT_FLOAT_EQ(ed4(0).adj(), 1.0);
  stan::math::set_zero_all_adjoints();
}

TEST(AgradRev, test_pass_through_working_all_var_types_double_test_6) {
  int size_vd = 3, rows_ed1 = 7, rows_ed2 = 3, cols_ed2 = 5, cols_ed3 = 1;
  stan::math::var d = 1.0;
  std::vector<stan::math::var> vd(size_vd, 1.0);
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> ed1(rows_ed1);
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> ed2(rows_ed2,
                                                                     cols_ed2);
  std::vector<int> vi(5, 5);
  Eigen::Matrix<stan::math::var, 1, Eigen::Dynamic> ed3(cols_ed3);
  Eigen::Matrix<double, 1, 1> ed4;

  for (int i = 0; i < size_vd; i++)
    vd[i] = 1.0;

  for (int i = 0; i < rows_ed1; i++)
    ed1(i) = 1.0;

  for (int i = 0; i < rows_ed2 * cols_ed2; i++)
    ed2(i) = 1.0;

  for (int i = 0; i < cols_ed3; i++)
    ed3(i) = 1.0;

  ed4(0) = 1.0;

  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> y
      = stan::math::adj_jac_apply<CheckAdjointsPassingThrough>(
          d, vd, 3, ed1, ed2, vi, ed3, ed4);

  y(0).grad();
  EXPECT_FLOAT_EQ(y(0).val(), d.val());
  EXPECT_FLOAT_EQ(d.adj(), 1.0);
  for (int i = 0; i < size_vd; i++)
    EXPECT_FLOAT_EQ(vd[i].adj(), 0);
  for (int i = 0; i < rows_ed1; i++)
    EXPECT_FLOAT_EQ(ed1(i).adj(), 0);
  for (int i = 0; i < rows_ed2 * cols_ed2; i++)
    EXPECT_FLOAT_EQ(ed2(i).adj(), 0);
  for (int i = 0; i < cols_ed3; i++)
    EXPECT_FLOAT_EQ(ed3(i).adj(), 0);
  stan::math::set_zero_all_adjoints();

  for (int j = 0; j < size_vd; j++) {
    y(1 + j).grad();
    EXPECT_FLOAT_EQ(y(1 + j).val(), vd[j].val());
    EXPECT_FLOAT_EQ(vd[j].adj(), 1.0);
    EXPECT_FLOAT_EQ(d.adj(), 0.0);
    for (int i = 0; i < size_vd; i++)
      if (i != j) {
        EXPECT_FLOAT_EQ(vd[i].adj(), 0);
      }
    for (int i = 0; i < rows_ed1; i++)
      EXPECT_FLOAT_EQ(ed1(i).adj(), 0);
    for (int i = 0; i < rows_ed2 * cols_ed2; i++)
      EXPECT_FLOAT_EQ(ed2(i).adj(), 0);
    for (int i = 0; i < cols_ed3; i++)
      EXPECT_FLOAT_EQ(ed3(i).adj(), 0);
    stan::math::set_zero_all_adjoints();
  }

  for (int j = 0; j < rows_ed1; j++) {
    y(1 + size_vd + j).grad();
    EXPECT_FLOAT_EQ(y(1 + size_vd + j).val(), ed1(j).val());
    EXPECT_FLOAT_EQ(ed1(j).adj(), 1.0);
    EXPECT_FLOAT_EQ(d.adj(), 0.0);
    for (int i = 0; i < size_vd; i++)
      EXPECT_FLOAT_EQ(vd[i].adj(), 0);
    for (int i = 0; i < rows_ed1; i++)
      if (i != j) {
        EXPECT_FLOAT_EQ(ed1(i).adj(), 0);
      }
    for (int i = 0; i < rows_ed2 * cols_ed2; i++)
      EXPECT_FLOAT_EQ(ed2(i).adj(), 0);
    for (int i = 0; i < cols_ed3; i++)
      EXPECT_FLOAT_EQ(ed3(i).adj(), 0);
    stan::math::set_zero_all_adjoints();
  }

  for (int j = 0; j < rows_ed2 * cols_ed2; j++) {
    y(1 + size_vd + rows_ed1 + j).grad();
    EXPECT_FLOAT_EQ(y(1 + size_vd + rows_ed1 + j).val(), ed2(j).val());
    EXPECT_FLOAT_EQ(ed2(j).adj(), 1.0);
    EXPECT_FLOAT_EQ(d.adj(), 0.0);
    for (int i = 0; i < size_vd; i++)
      EXPECT_FLOAT_EQ(vd[i].adj(), 0);
    for (int i = 0; i < rows_ed1; i++)
      EXPECT_FLOAT_EQ(ed1(i).adj(), 0);
    for (int i = 0; i < rows_ed2 * cols_ed2; i++)
      if (i != j) {
        EXPECT_FLOAT_EQ(ed2(i).adj(), 0);
      }
    for (int i = 0; i < cols_ed3; i++)
      EXPECT_FLOAT_EQ(ed3(i).adj(), 0);
    stan::math::set_zero_all_adjoints();
  }

  for (int j = 0; j < cols_ed3; j++) {
    y(1 + size_vd + rows_ed1 + rows_ed2 * cols_ed2 + j).grad();
    EXPECT_FLOAT_EQ(y(1 + size_vd + rows_ed1 + rows_ed2 * cols_ed2 + j).val(),
                    ed3(j).val());
    EXPECT_FLOAT_EQ(ed3(j).adj(), 1.0);
    EXPECT_FLOAT_EQ(d.adj(), 0.0);
    for (int i = 0; i < size_vd; i++)
      EXPECT_FLOAT_EQ(vd[i].adj(), 0);
    for (int i = 0; i < rows_ed1; i++)
      EXPECT_FLOAT_EQ(ed1(i).adj(), 0);
    for (int i = 0; i < rows_ed2 * cols_ed2; i++)
      EXPECT_FLOAT_EQ(ed2(i).adj(), 0);
    for (int i = 0; i < cols_ed3; i++)
      if (i != j) {
        EXPECT_FLOAT_EQ(ed3(i).adj(), 0);
      }
    stan::math::set_zero_all_adjoints();
  }

  y(1 + size_vd + rows_ed1 + rows_ed2 * cols_ed2 + cols_ed3).grad();
  EXPECT_FLOAT_EQ(
      y(1 + size_vd + rows_ed1 + rows_ed2 * cols_ed2 + cols_ed3).val(), ed4(0));
  EXPECT_FLOAT_EQ(d.adj(), 0.0);
  for (int i = 0; i < size_vd; i++)
    EXPECT_FLOAT_EQ(vd[i].adj(), 0);
  for (int i = 0; i < rows_ed1; i++)
    EXPECT_FLOAT_EQ(ed1(i).adj(), 0);
  for (int i = 0; i < rows_ed2 * cols_ed2; i++)
    EXPECT_FLOAT_EQ(ed2(i).adj(), 0);
  for (int i = 0; i < cols_ed3; i++)
    EXPECT_FLOAT_EQ(ed3(i).adj(), 0);
  stan::math::set_zero_all_adjoints();
}
/**
 * Test a functor with multiple input types that takes advantage of needs_adj_
 * functionality
 */

template <typename T1, typename T2, typename T3, typename T4>
struct SinCosFunctor {
  stan::math::adj_op<T1> x1_;
  stan::math::adj_op<T2> x2_;
  stan::math::adj_op<T3> x3_;
  stan::math::adj_op<T4> x4_;
  SinCosFunctor(const T1& x1, const T2& x2, const T3& x3, const T4& x4)
      : x1_(x1), x2_(x2), x3_(x3), x4_(x4) {}
  Eigen::VectorXd operator()(const Eigen::VectorXd& x1, const int& x2,
                             const std::vector<int>& x3,
                             const std::vector<double>& x4) {
    stan::math::check_matching_sizes("SinCosFunctor", "x1", x1, "x4", x4);
    Eigen::VectorXd out(x1.size());
    for (int n = 0; n < x1.size(); ++n) {
      out(n) = sin(x1(n)) + cos(x4[n]);
    }
    return out;
  }

  auto multiply_adjoint_jacobian(const Eigen::VectorXd& adj) {
    Eigen::VectorXd out1;
    std::vector<double> out4;
    if (stan::is_var<stan::scalar_type_t<T1>>::value) {
      out1.resize(x1_.size());
      for (int n = 0; n < x1_.size(); ++n) {
        out1(n) = cos(x1_(n)) * adj(n);
      }
    }
    if (stan::is_var<stan::scalar_type_t<T4>>::value) {
      out4.resize(x4_.size());
      for (int n = 0; n < x4_.size(); ++n) {
        out4[n] = -sin(x4_(n)) * adj(n);
      }
    }
    return std::make_tuple(out1, 0, std::vector<int>(), out4);
  }
};

TEST(AgradRev, test_sincos_stack) {
  using eig_vec = Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>;
  using std_vec = std::vector<stan::math::var>;
  eig_vec x11(1), x21(2), y1(1), y2(2);
  std_vec x12(1), x22(2);
  x11 << 1.0;
  x12[0] = 5.0;
  x21 << 2.0, 1.0;
  x22[0] = 5.0;
  x22[1] = 3.0;
  y1 = stan::math::adj_jac_apply<
      SinCosFunctor<eig_vec, int, std::vector<int>, std_vec>>(
      x11, 0, std::vector<int>(5, 0), x12);
  y2 = stan::math::adj_jac_apply<
      SinCosFunctor<eig_vec, int, std::vector<int>, std_vec>>(
      x21, 0, std::vector<int>(5, 0), x22);

  test::check_varis_on_stack(y1);
  test::check_varis_on_stack(y2);
}

TEST(AgradRev, test_sincos_values) {
  using eig_vec = Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>;
  using std_vec = std::vector<stan::math::var>;
  eig_vec x11(1), x21(2), y1(1), y2(2);
  std_vec x12(1), x22(2);
  x11 << 1.0;
  x12[0] = 5.0;
  x21 << 2.0, 1.0;
  x22[0] = 5.0;
  x22[1] = 3.0;
  using sin_cos_func = SinCosFunctor<eig_vec, int, std::vector<int>, std_vec>;
  y1 = stan::math::adj_jac_apply<sin_cos_func>(x11, 0, std::vector<int>(5, 0),
                                               x12);
  y2 = stan::math::adj_jac_apply<sin_cos_func>(x21, 0, std::vector<int>(5, 0),
                                               x22);

  EXPECT_FLOAT_EQ(y1(0).val(), 1.125133170271123);
  EXPECT_FLOAT_EQ(y2(0).val(), 1.192959612288908);
  EXPECT_FLOAT_EQ(y2(1).val(), -0.1485215117925489);
}

TEST(AgradRev, test_sincos_multiple_jac_vv) {
  using eig_vec = Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>;
  using std_vec = std::vector<stan::math::var>;
  eig_vec x11(1), x21(2), y1(1), y2(2);
  std_vec x12(1), x22(2);
  x11 << 1.0;
  x12[0] = 5.0;
  x21 << 2.0, 1.0;
  x22[0] = 5.0;
  x22[1] = 3.0;
  using sin_cos_func = SinCosFunctor<eig_vec, int, std::vector<int>, std_vec>;
  y1 = stan::math::adj_jac_apply<sin_cos_func>(x11, 0, std::vector<int>(5, 0),
                                               x12);
  y2 = stan::math::adj_jac_apply<sin_cos_func>(x21, 0, std::vector<int>(5, 0),
                                               x22);

  y1(0).grad();
  EXPECT_FLOAT_EQ(x11(0).adj(), 0.5403023058681398);
  EXPECT_FLOAT_EQ(x12[0].adj(), 0.958924274663139);
  EXPECT_FLOAT_EQ(x21(0).adj(), 0.0);
  EXPECT_FLOAT_EQ(x21(1).adj(), 0.0);
  EXPECT_FLOAT_EQ(x22[0].adj(), 0.0);
  EXPECT_FLOAT_EQ(x22[1].adj(), 0.0);

  stan::math::set_zero_all_adjoints();

  y2(0).grad();
  EXPECT_FLOAT_EQ(x11(0).adj(), 0.0);
  EXPECT_FLOAT_EQ(x12[0].adj(), 0.0);
  EXPECT_FLOAT_EQ(x21(0).adj(), -0.4161468365471424);
  EXPECT_FLOAT_EQ(x21(1).adj(), 0.0);
  EXPECT_FLOAT_EQ(x22[0].adj(), 0.958924274663139);
  EXPECT_FLOAT_EQ(x22[1].adj(), 0.0);

  stan::math::set_zero_all_adjoints();

  y2(1).grad();
  EXPECT_FLOAT_EQ(x11(0).adj(), 0.0);
  EXPECT_FLOAT_EQ(x12[0].adj(), 0.0);
  EXPECT_FLOAT_EQ(x21(0).adj(), 0.0);
  EXPECT_FLOAT_EQ(x21(1).adj(), 0.5403023058681398);
  EXPECT_FLOAT_EQ(x22[0].adj(), 0.0);
  EXPECT_FLOAT_EQ(x22[1].adj(), -0.1411200080598672);

  stan::math::set_zero_all_adjoints();

  stan::math::var sum_y2 = (1.73 * y2(0) + 1.57 * y2(1));
  sum_y2.grad();
  EXPECT_FLOAT_EQ(x11(0).adj(), 0.0);
  EXPECT_FLOAT_EQ(x12[0].adj(), 0.0);
  EXPECT_FLOAT_EQ(x21(0).adj(), 1.73 * -0.4161468365471424);
  EXPECT_FLOAT_EQ(x21(1).adj(), 1.57 * 0.5403023058681398);
  EXPECT_FLOAT_EQ(x22[0].adj(), 1.73 * 0.958924274663139);
  EXPECT_FLOAT_EQ(x22[1].adj(), 1.57 * -0.1411200080598672);
}

TEST(AgradRev, test_sincos_multiple_jac_dv) {
  using eig_vec_v = Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>;
  using eig_vec_d = Eigen::Matrix<double, Eigen::Dynamic, 1>;
  using std_vec_v = std::vector<stan::math::var>;
  eig_vec_v y1(1), y2(2);
  eig_vec_d x11(1), x21(2);
  std_vec_v x12(1), x22(2);
  x11 << 1.0;
  x12[0] = 5.0;
  x21 << 2.0, 1.0;
  x22[0] = 5.0;
  x22[1] = 3.0;
  using sin_cos_func
      = SinCosFunctor<eig_vec_d, int, std::vector<int>, std_vec_v>;
  y1 = stan::math::adj_jac_apply<sin_cos_func>(x11, 0, std::vector<int>(5, 0),
                                               x12);
  y2 = stan::math::adj_jac_apply<sin_cos_func>(x21, 0, std::vector<int>(5, 0),
                                               x22);

  y1(0).grad();
  EXPECT_FLOAT_EQ(x12[0].adj(), 0.958924274663139);
  EXPECT_FLOAT_EQ(x22[0].adj(), 0.0);
  EXPECT_FLOAT_EQ(x22[1].adj(), 0.0);

  stan::math::set_zero_all_adjoints();

  y2(0).grad();
  EXPECT_FLOAT_EQ(x12[0].adj(), 0.0);
  EXPECT_FLOAT_EQ(x22[0].adj(), 0.958924274663139);
  EXPECT_FLOAT_EQ(x22[1].adj(), 0.0);

  stan::math::set_zero_all_adjoints();

  y2(1).grad();
  EXPECT_FLOAT_EQ(x12[0].adj(), 0.0);
  EXPECT_FLOAT_EQ(x22[0].adj(), 0.0);
  EXPECT_FLOAT_EQ(x22[1].adj(), -0.1411200080598672);

  stan::math::set_zero_all_adjoints();

  stan::math::var sum_y2 = (1.73 * y2(0) + 1.57 * y2(1));
  sum_y2.grad();
  EXPECT_FLOAT_EQ(x12[0].adj(), 0.0);
  EXPECT_FLOAT_EQ(x22[0].adj(), 1.73 * 0.958924274663139);
  EXPECT_FLOAT_EQ(x22[1].adj(), 1.57 * -0.1411200080598672);
}

TEST(AgradRev, test_sincos_multiple_jac_vd) {
  using eig_vec_v = Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>;
  using std_vec_d = std::vector<double>;
  eig_vec_v x11(1), x21(2), y1(1), y2(2);
  std_vec_d x12(1), x22(2);
  x11 << 1.0;
  x12[0] = 5.0;
  x21 << 2.0, 1.0;
  x22[0] = 5.0;
  x22[1] = 3.0;
  using sin_cos_func
      = SinCosFunctor<eig_vec_v, int, std::vector<int>, std_vec_d>;
  y1 = stan::math::adj_jac_apply<sin_cos_func>(x11, 0, std::vector<int>(5, 0),
                                               x12);
  y2 = stan::math::adj_jac_apply<sin_cos_func>(x21, 0, std::vector<int>(5, 0),
                                               x22);

  y1(0).grad();
  EXPECT_FLOAT_EQ(x11(0).adj(), 0.5403023058681398);
  EXPECT_FLOAT_EQ(x21(0).adj(), 0.0);
  EXPECT_FLOAT_EQ(x21(1).adj(), 0.0);

  stan::math::set_zero_all_adjoints();

  y2(0).grad();
  EXPECT_FLOAT_EQ(x11(0).adj(), 0.0);
  EXPECT_FLOAT_EQ(x21(0).adj(), -0.4161468365471424);
  EXPECT_FLOAT_EQ(x21(1).adj(), 0.0);

  stan::math::set_zero_all_adjoints();

  y2(1).grad();
  EXPECT_FLOAT_EQ(x11(0).adj(), 0.0);
  EXPECT_FLOAT_EQ(x21(0).adj(), 0.0);
  EXPECT_FLOAT_EQ(x21(1).adj(), 0.5403023058681398);

  stan::math::set_zero_all_adjoints();

  stan::math::var sum_y2 = (1.73 * y2(0) + 1.57 * y2(1));
  sum_y2.grad();
  EXPECT_FLOAT_EQ(x11(0).adj(), 0.0);
  EXPECT_FLOAT_EQ(x21(0).adj(), 1.73 * -0.4161468365471424);
  EXPECT_FLOAT_EQ(x21(1).adj(), 1.57 * 0.5403023058681398);
}

/**
 * Test a multi-argument functor with an Eigen VectorXd input and a double input
 */
template <typename T1, typename T2>
struct SinCosFunctor2 {
  stan::math::adj_op<T1> x1_;
  stan::math::adj_op<T2> x2_;
  SinCosFunctor2(const T1& x1, const T2& x2) : x1_(x1), x2_(x2) {}

  Eigen::VectorXd operator()(const Eigen::VectorXd& x1, const double& x2) {
    Eigen::VectorXd out(x1.size());
    for (int n = 0; n < x1.size(); ++n) {
      out(n) = sin(x1(n)) + cos(x2);
    }
    return out;
  }

  auto multiply_adjoint_jacobian(const Eigen::VectorXd& adj) {
    Eigen::VectorXd out1;
    double out2 = 0.0;
    if (x1_.needs_adj) {
      out1.resize(x1_.size());
      for (int n = 0; n < x1_.size(); ++n) {
        out1(n) = cos(x1_(n)) * adj(n);
      }
    }
    if (x2_.needs_adj) {
      for (int n = 0; n < x1_.size(); ++n) {
        out2 += -sin(x2_.map()) * adj(n);
      }
    }
    return std::make_tuple(out1, out2);
  }
};

TEST(AgradRev, test_eigen_vector_scalar_stack) {
  using eig_vec = Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>;
  using stan::math::var;
  eig_vec x11(1), x21(2), y1(1), y2(2);
  var x12, x22;
  x11 << 1.0;
  x12 = 5.0;
  x21 << 2.0, 1.0;
  x22 = 3.0;

  y1 = stan::math::adj_jac_apply<SinCosFunctor2<eig_vec, var>>(x11, x12);
  y2 = stan::math::adj_jac_apply<SinCosFunctor2<eig_vec, var>>(x21, x22);

  test::check_varis_on_stack(y1);
  test::check_varis_on_stack(y2);
}

TEST(AgradRev, test_eigen_vector_scalar_values) {
  using eig_vec = Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>;
  using stan::math::var;
  eig_vec x11(1), x21(2), y1(1), y2(2);
  var x12, x22;
  x11 << 1.0;
  x12 = 5.0;
  x21 << 2.0, 1.0;
  x22 = 3.0;

  y1 = stan::math::adj_jac_apply<SinCosFunctor2<eig_vec, var>>(x11, x12);
  y2 = stan::math::adj_jac_apply<SinCosFunctor2<eig_vec, var>>(x21, x22);

  EXPECT_FLOAT_EQ(y1(0).val(), 1.125133170271123);
  EXPECT_FLOAT_EQ(y2(0).val(), -0.0806950697747637);
  EXPECT_FLOAT_EQ(y2(1).val(), -0.1485215117925489);
}

TEST(AgradRev, test_eigen_vector_scalar_multiple_jac_vv) {
  using eig_vec = Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>;
  using stan::math::var;
  eig_vec x11(1), x21(2), y1(1), y2(2);
  var x12, x22;
  x11 << 1.0;
  x12 = 5.0;
  x21 << 2.0, 1.0;
  x22 = 3.0;

  y1 = stan::math::adj_jac_apply<SinCosFunctor2<eig_vec, var>>(x11, x12);
  y2 = stan::math::adj_jac_apply<SinCosFunctor2<eig_vec, var>>(x21, x22);

  y1(0).grad();
  EXPECT_FLOAT_EQ(x11(0).adj(), 0.5403023058681398);
  EXPECT_FLOAT_EQ(x12.adj(), 0.958924274663139);
  EXPECT_FLOAT_EQ(x21(0).adj(), 0.0);
  EXPECT_FLOAT_EQ(x21(1).adj(), 0.0);
  EXPECT_FLOAT_EQ(x22.adj(), 0.0);

  stan::math::set_zero_all_adjoints();

  y2(0).grad();
  EXPECT_FLOAT_EQ(x11(0).adj(), 0.0);
  EXPECT_FLOAT_EQ(x12.adj(), 0.0);
  EXPECT_FLOAT_EQ(x21(0).adj(), -0.4161468365471424);
  EXPECT_FLOAT_EQ(x21(1).adj(), 0.0);
  EXPECT_FLOAT_EQ(x22.adj(), -0.1411200080598672);

  stan::math::set_zero_all_adjoints();

  y2(1).grad();
  EXPECT_FLOAT_EQ(x11(0).adj(), 0.0);
  EXPECT_FLOAT_EQ(x12.adj(), 0.0);
  EXPECT_FLOAT_EQ(x21(0).adj(), 0.0);
  EXPECT_FLOAT_EQ(x21(1).adj(), 0.5403023058681398);
  EXPECT_FLOAT_EQ(x22.adj(), -0.1411200080598672);

  stan::math::set_zero_all_adjoints();

  stan::math::var sum_y2 = (1.73 * y2(0) + 1.57 * y2(1));
  sum_y2.grad();
  EXPECT_FLOAT_EQ(x11(0).adj(), 0.0);
  EXPECT_FLOAT_EQ(x12.adj(), 0.0);
  EXPECT_FLOAT_EQ(x21(0).adj(), 1.73 * -0.4161468365471424);
  EXPECT_FLOAT_EQ(x21(1).adj(), 1.57 * 0.5403023058681398);
  EXPECT_FLOAT_EQ(x22.adj(), (1.73 + 1.57) * -0.1411200080598672);
}

TEST(AgradRev, test_eigen_vector_scalar_multiple_jac_dv) {
  using eig_vec_v = Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>;
  using eig_vec_d = Eigen::Matrix<double, Eigen::Dynamic, 1>;
  using stan::math::var;
  eig_vec_v y1(1), y2(2);
  eig_vec_d x11(1), x21(2);
  var x12, x22;
  x11 << 1.0;
  x12 = 5.0;
  x21 << 2.0, 1.0;
  x22 = 3.0;
  using sin_cos_fun = SinCosFunctor2<eig_vec_v, var>;
  y1 = stan::math::adj_jac_apply<sin_cos_fun>(x11, x12);
  y2 = stan::math::adj_jac_apply<sin_cos_fun>(x21, x22);

  y1(0).grad();
  EXPECT_FLOAT_EQ(x12.adj(), 0.958924274663139);
  EXPECT_FLOAT_EQ(x22.adj(), 0.0);

  stan::math::set_zero_all_adjoints();

  y2(0).grad();
  EXPECT_FLOAT_EQ(x12.adj(), 0.0);
  EXPECT_FLOAT_EQ(x22.adj(), -0.1411200080598672);

  stan::math::set_zero_all_adjoints();

  y2(1).grad();
  EXPECT_FLOAT_EQ(x12.adj(), 0.0);
  EXPECT_FLOAT_EQ(x22.adj(), -0.1411200080598672);

  stan::math::set_zero_all_adjoints();

  stan::math::var sum_y2 = (1.73 * y2(0) + 1.57 * y2(1));
  sum_y2.grad();
  EXPECT_FLOAT_EQ(x12.adj(), 0.0);
  EXPECT_FLOAT_EQ(x22.adj(), (1.73 + 1.57) * -0.1411200080598672);
}

TEST(AgradRev, test_eigen_vector_scalar_multiple_jac_vd) {
  using eig_vec_v = Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>;
  eig_vec_v x11(1), x21(2), y1(1), y2(2);
  double x12, x22;
  x11 << 1.0;
  x12 = 5.0;
  x21 << 2.0, 1.0;
  x22 = 3.0;
  using sin_cos_func = SinCosFunctor2<eig_vec_v, double>;
  y1 = stan::math::adj_jac_apply<sin_cos_func>(x11, x12);
  y2 = stan::math::adj_jac_apply<sin_cos_func>(x21, x22);

  y1(0).grad();
  EXPECT_FLOAT_EQ(x11(0).adj(), 0.5403023058681398);
  EXPECT_FLOAT_EQ(x21(0).adj(), 0.0);
  EXPECT_FLOAT_EQ(x21(1).adj(), 0.0);

  stan::math::set_zero_all_adjoints();

  y2(0).grad();
  EXPECT_FLOAT_EQ(x11(0).adj(), 0.0);
  EXPECT_FLOAT_EQ(x21(0).adj(), -0.4161468365471424);
  EXPECT_FLOAT_EQ(x21(1).adj(), 0.0);

  stan::math::set_zero_all_adjoints();

  y2(1).grad();
  EXPECT_FLOAT_EQ(x11(0).adj(), 0.0);
  EXPECT_FLOAT_EQ(x21(0).adj(), 0.0);
  EXPECT_FLOAT_EQ(x21(1).adj(), 0.5403023058681398);

  stan::math::set_zero_all_adjoints();

  stan::math::var sum_y2 = (1.73 * y2(0) + 1.57 * y2(1));
  sum_y2.grad();
  EXPECT_FLOAT_EQ(x11(0).adj(), 0.0);
  EXPECT_FLOAT_EQ(x21(0).adj(), 1.73 * -0.4161468365471424);
  EXPECT_FLOAT_EQ(x21(1).adj(), 1.57 * 0.5403023058681398);
}

/**
 * Test a functor with a double then Eigen::VectorXd input
 * idk why these numbers are swapped
 */
template <typename T1, typename T2>
struct SinCosFunctor3 {
  int N_;
  stan::math::adj_op<T1> x1_;
  stan::math::adj_op<T2> x2_;
  template <typename S1, typename S2>
  SinCosFunctor3(const S1& x1, const S2& x2) : x1_(x1), x2_(x2) {}
  Eigen::VectorXd operator()(const double& x1, const Eigen::VectorXd& x2) {
    N_ = x2.size();
    Eigen::VectorXd out(N_);
    for (int n = 0; n < N_; ++n) {
      out(n) = sin(x2(n)) + cos(x1);
    }

    return out;
  }

  auto multiply_adjoint_jacobian(const Eigen::VectorXd& adj) {
    Eigen::VectorXd out2;
    double out1 = 0.0;
    if (x1_.needs_adj) {
      for (int n = 0; n < N_; ++n) {
        out1 += -sin(x1_.map()) * adj(n);
      }
    }
    if (x2_.needs_adj) {
      out2.resize(N_);
      for (int n = 0; n < N_; ++n) {
        out2(n) = cos(x2_(n)) * adj(n);
      }
    }
    return std::make_tuple(out1, out2);
  }
};

TEST(AgradRev, test_sincos_scalar_eigen_vector_stack) {
  using eig_vec = Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>;
  using stan::math::var;
  eig_vec x11(1), x21(2), y1(1), y2(2);
  var x12, x22;
  x11 << 1.0;
  x12 = 5.0;
  x21 << 2.0, 1.0;
  x22 = 3.0;
  using sin_cos_fun = SinCosFunctor3<stan::math::var, eig_vec>;
  y1 = stan::math::adj_jac_apply<sin_cos_fun>(x12, x11);
  y2 = stan::math::adj_jac_apply<sin_cos_fun>(x22, x21);
  test::check_varis_on_stack(y1);
  test::check_varis_on_stack(y2);
}

TEST(AgradRev, test_sincos_scalar_eigen_vector_values) {
  using eig_vec = Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>;
  using stan::math::var;
  eig_vec x11(1), x21(2), y1(1), y2(2);
  stan::math::var x12, x22;
  x11 << 1.0;
  x12 = 5.0;
  x21 << 2.0, 1.0;
  x22 = 3.0;
  using sin_cos_fun = SinCosFunctor3<stan::math::var, eig_vec>;
  y1 = stan::math::adj_jac_apply<sin_cos_fun>(x12, x11);
  y2 = stan::math::adj_jac_apply<sin_cos_fun>(x22, x21);

  EXPECT_FLOAT_EQ(y1(0).val(), 1.125133170271123);
  EXPECT_FLOAT_EQ(y2(0).val(), -0.0806950697747637);
  EXPECT_FLOAT_EQ(y2(1).val(), -0.1485215117925489);
}

TEST(AgradRev, test_sincos_scalar_eigen_vector_multiple_jac_vv) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> x11(1), x21(2), y1(1),
      y2(2);
  stan::math::var x12, x22;
  x11 << 1.0;
  x12 = 5.0;
  x21 << 2.0, 1.0;
  x22 = 3.0;
  using sin_cos_func
      = SinCosFunctor3<stan::math::var,
                       Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>>;
  y1 = stan::math::adj_jac_apply<sin_cos_func>(x12, x11);
  y2 = stan::math::adj_jac_apply<sin_cos_func>(x22, x21);
  /*
    y1(0).grad();
    EXPECT_FLOAT_EQ(x11(0).adj(), 0.5403023058681398);
    EXPECT_FLOAT_EQ(x12.adj(), 0.958924274663139);
    EXPECT_FLOAT_EQ(x21(0).adj(), 0.0);
    EXPECT_FLOAT_EQ(x21(1).adj(), 0.0);
    EXPECT_FLOAT_EQ(x22.adj(), 0.0);

    stan::math::set_zero_all_adjoints();

    y2(0).grad();
    EXPECT_FLOAT_EQ(x11(0).adj(), 0.0);
    EXPECT_FLOAT_EQ(x12.adj(), 0.0);
    EXPECT_FLOAT_EQ(x21(0).adj(), -0.4161468365471424);
    EXPECT_FLOAT_EQ(x21(1).adj(), 0.0);
    EXPECT_FLOAT_EQ(x22.adj(), -0.1411200080598672);

    stan::math::set_zero_all_adjoints();

    y2(1).grad();
    EXPECT_FLOAT_EQ(x11(0).adj(), 0.0);
    EXPECT_FLOAT_EQ(x12.adj(), 0.0);
    EXPECT_FLOAT_EQ(x21(0).adj(), 0.0);
    EXPECT_FLOAT_EQ(x21(1).adj(), 0.5403023058681398);
    EXPECT_FLOAT_EQ(x22.adj(), -0.1411200080598672);

    stan::math::set_zero_all_adjoints();
  */
  stan::math::var sum_y2 = (1.73 * y2(0) + 1.57 * y2(1));
  sum_y2.grad();
  EXPECT_FLOAT_EQ(x11(0).adj(), 0.0);
  EXPECT_FLOAT_EQ(x12.adj(), 0.0);
  EXPECT_FLOAT_EQ(x21(0).adj(), 1.73 * -0.4161468365471424);
  EXPECT_FLOAT_EQ(x21(1).adj(), 1.57 * 0.5403023058681398);
  EXPECT_FLOAT_EQ(x22.adj(), (1.73 + 1.57) * -0.1411200080598672);
}

TEST(AgradRev, test_sincos_scalar_eigen_vector_multiple_jac_vd) {
  using eig_vec_v = Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>;
  using eig_vec_d = Eigen::Matrix<double, Eigen::Dynamic, 1>;
  using stan::math::var;
  eig_vec_v y1(1), y2(2);
  eig_vec_d x11(1), x21(2);
  var x12, x22;
  x11 << 1.0;
  x12 = 5.0;
  x21 << 2.0, 1.0;
  x22 = 3.0;
  using sin_cos_func = SinCosFunctor3<var, eig_vec_d>;
  y1 = stan::math::adj_jac_apply<sin_cos_func>(x12, x11);
  y2 = stan::math::adj_jac_apply<sin_cos_func>(x22, x21);

  y1(0).grad();
  EXPECT_FLOAT_EQ(x12.adj(), 0.958924274663139);
  EXPECT_FLOAT_EQ(x22.adj(), 0.0);

  stan::math::set_zero_all_adjoints();

  y2(0).grad();
  EXPECT_FLOAT_EQ(x12.adj(), 0.0);
  EXPECT_FLOAT_EQ(x22.adj(), -0.1411200080598672);

  stan::math::set_zero_all_adjoints();

  y2(1).grad();
  EXPECT_FLOAT_EQ(x12.adj(), 0.0);
  EXPECT_FLOAT_EQ(x22.adj(), -0.1411200080598672);

  stan::math::set_zero_all_adjoints();

  stan::math::var sum_y2 = (1.73 * y2(0) + 1.57 * y2(1));
  sum_y2.grad();
  EXPECT_FLOAT_EQ(x12.adj(), 0.0);
  EXPECT_FLOAT_EQ(x22.adj(), (1.73 + 1.57) * -0.1411200080598672);
}

TEST(AgradRev, test_sincos_scalar_eigen_vector_multiple_jac_dv) {
  using eig_vec_v = Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>;
  eig_vec_v x11(1), x21(2), y1(1), y2(2);
  double x12, x22;
  x11 << 1.0;
  x12 = 5.0;
  x21 << 2.0, 1.0;
  x22 = 3.0;
  using sin_cos_func = SinCosFunctor3<double, eig_vec_v>;
  y1 = stan::math::adj_jac_apply<sin_cos_func>(x12, x11);
  y2 = stan::math::adj_jac_apply<sin_cos_func>(x22, x21);

  y1(0).grad();
  EXPECT_FLOAT_EQ(x11(0).adj(), 0.5403023058681398);
  EXPECT_FLOAT_EQ(x21(0).adj(), 0.0);
  EXPECT_FLOAT_EQ(x21(1).adj(), 0.0);

  stan::math::set_zero_all_adjoints();

  y2(0).grad();
  EXPECT_FLOAT_EQ(x11(0).adj(), 0.0);
  EXPECT_FLOAT_EQ(x21(0).adj(), -0.4161468365471424);
  EXPECT_FLOAT_EQ(x21(1).adj(), 0.0);

  stan::math::set_zero_all_adjoints();

  y2(1).grad();
  EXPECT_FLOAT_EQ(x11(0).adj(), 0.0);
  EXPECT_FLOAT_EQ(x21(0).adj(), 0.0);
  EXPECT_FLOAT_EQ(x21(1).adj(), 0.5403023058681398);

  stan::math::set_zero_all_adjoints();

  stan::math::var sum_y2 = (1.73 * y2(0) + 1.57 * y2(1));
  sum_y2.grad();
  EXPECT_FLOAT_EQ(x11(0).adj(), 0.0);
  EXPECT_FLOAT_EQ(x21(0).adj(), 1.73 * -0.4161468365471424);
  EXPECT_FLOAT_EQ(x21(1).adj(), 1.57 * 0.5403023058681398);
}
