#include <stan/math/rev/core.hpp>
#include <test/unit/math/rev/mat/util.hpp>
#include <gtest/gtest.h>
#include <algorithm>
#include <sstream>
#include <tuple>
#include <vector>

/*
 * Check scalar return type
 */
struct ScalarSinFunctor {
  double x_;
  template <std::size_t size>
  double operator()(const std::array<bool, size>& needs_adj, const double& x) {
    x_ = x;

    return sin(x_);
  }

  template <std::size_t size>
  auto multiply_adjoint_jacobian(const std::array<bool, size>& needs_adj,
                                 const double& adj) {
    return std::make_tuple(cos(x_) * adj);
  }
};

TEST(AgradRev, test_scalar_sin_stack) {
  stan::math::var x1, y1;
  x1 = 1.0;

  y1 = stan::math::adj_jac_apply<ScalarSinFunctor>(x1);

  test::check_varis_on_stack(y1);
}

TEST(AgradRev, test_scalar_sin_values) {
  stan::math::var x1, y1;
  x1 = 1.0;

  y1 = stan::math::adj_jac_apply<ScalarSinFunctor>(x1);

  EXPECT_NEAR(y1.val(), 0.841470984807897, 1e-10);
}

TEST(AgradRev, test_scalar_sin_jac) {
  stan::math::var x1, y1;
  x1 = 1.0;

  y1 = stan::math::adj_jac_apply<ScalarSinFunctor>(x1);

  y1.grad();
  EXPECT_NEAR(x1.adj(), 0.5403023058681398, 1e-10);
}

/*
 * Check std::vector return type
 */
struct StdVectorSinFunctor {
  double* x_;
  int N_;

  template <std::size_t size>
  std::vector<double> operator()(const std::array<bool, size>& needs_adj,
                                 const std::vector<double>& x) {
    N_ = x.size();
    std::vector<double> out(N_);

    x_ = stan::math::ChainableStack::instance().memalloc_.alloc_array<double>(
        N_);

    for (int i = 0; i < N_; ++i) {
      out[i] = sin(x[i]);
      x_[i] = x[i];
    }

    return out;
  }

  template <std::size_t size>
  auto multiply_adjoint_jacobian(const std::array<bool, size>& needs_adj,
                                 const std::vector<double>& adj) {
    std::vector<double> adj_jac(N_);
    for (int i = 0; i < N_; ++i) {
      adj_jac[i] = cos(x_[i]) * adj[i];
    }
    return std::make_tuple(adj_jac);
  }
};

TEST(AgradRev, test_std_vector_sin_stack) {
  std::vector<stan::math::var> x1(2), y1;
  x1[0] = 1.0;
  x1[1] = 2.0;

  y1 = stan::math::adj_jac_apply<StdVectorSinFunctor>(x1);

  test::check_varis_on_stack(y1);
}

TEST(AgradRev, test_std_vector_sin_values) {
  std::vector<stan::math::var> x1(2), y1;
  x1[0] = 1.0;
  x1[1] = 2.0;

  y1 = stan::math::adj_jac_apply<StdVectorSinFunctor>(x1);

  EXPECT_NEAR(y1[0].val(), 0.841470984807897, 1e-10);
  EXPECT_NEAR(y1[1].val(), 0.909297426825682, 1e-10);
}

TEST(AgradRev, test_std_vector_sin_jac) {
  std::vector<stan::math::var> x1(2), y1;
  x1[0] = 1.0;
  x1[1] = 2.0;

  y1 = stan::math::adj_jac_apply<StdVectorSinFunctor>(x1);

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
struct SinFunctor {
  int N_;
  double* x_mem_;
  template <std::size_t size>
  Eigen::VectorXd operator()(const std::array<bool, size>& needs_adj,
                             const Eigen::VectorXd& x) {
    N_ = x.size();
    Eigen::VectorXd out(N_);
    x_mem_
        = stan::math::ChainableStack::instance().memalloc_.alloc_array<double>(
            N_);

    for (int n = 0; n < N_; ++n) {
      x_mem_[n] = x(n);
      out(n) = sin(x(n));
    }

    return out;
  }

  template <std::size_t size>
  auto multiply_adjoint_jacobian(const std::array<bool, size>& needs_adj,
                                 const Eigen::VectorXd& adj) {
    Eigen::VectorXd out(N_);

    for (int n = 0; n < N_; ++n) {
      out(n) = cos(x_mem_[n]) * adj(n);
    }

    return std::make_tuple(out);
  }
};

TEST(AgradRev, test_vector_sin_stack) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> x1(1), x2(2), y1(1), y2(2);
  x1 << 1.0;
  x2 << 2.0, 1.0;

  y1 = stan::math::adj_jac_apply<SinFunctor>(x1);
  y2 = stan::math::adj_jac_apply<SinFunctor>(x2);

  test::check_varis_on_stack(y1);
  test::check_varis_on_stack(y2);
}

TEST(AgradRev, test_vector_sin_values) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> x1(1), x2(2), y1(1), y2(2);
  x1 << 1.0;
  x2 << 2.0, 1.0;

  y1 = stan::math::adj_jac_apply<SinFunctor>(x1);
  y2 = stan::math::adj_jac_apply<SinFunctor>(x2);

  EXPECT_FLOAT_EQ(y1(0).val(), 0.841470984807897);
  EXPECT_FLOAT_EQ(y2(0).val(), 0.909297426825682);
  EXPECT_FLOAT_EQ(y2(1).val(), 0.841470984807897);
}

TEST(AgradRev, test_vector_sin_multiple_jac) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> x1(1), x2(2), y1(1), y2(2);
  x1 << 1.0;
  x2 << 2.0, 1.0;

  y1 = stan::math::adj_jac_apply<SinFunctor>(x1);
  y2 = stan::math::adj_jac_apply<SinFunctor>(x2);

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
struct RowVectorSinFunctor {
  int N_;
  double* x_mem_;
  template <std::size_t size>
  Eigen::RowVectorXd operator()(const std::array<bool, size>& needs_adj,
                                const Eigen::RowVectorXd& x) {
    N_ = x.size();
    Eigen::RowVectorXd out(N_);
    x_mem_
        = stan::math::ChainableStack::instance().memalloc_.alloc_array<double>(
            N_);

    for (int n = 0; n < N_; ++n) {
      x_mem_[n] = x(n);
      out(n) = sin(x(n));
    }

    return out;
  }

  template <std::size_t size>
  auto multiply_adjoint_jacobian(const std::array<bool, size>& needs_adj,
                                 const Eigen::RowVectorXd& adj) {
    Eigen::RowVectorXd out(N_);

    for (int n = 0; n < N_; ++n) {
      out(n) = cos(x_mem_[n]) * adj(n);
    }

    return std::make_tuple(out);
  }
};

TEST(AgradRev, test_row_vector_sin_stack) {
  Eigen::Matrix<stan::math::var, 1, Eigen::Dynamic> x1(1), x2(2), y1(1), y2(2);
  x1 << 1.0;
  x2 << 2.0, 1.0;

  y1 = stan::math::adj_jac_apply<RowVectorSinFunctor>(x1);
  y2 = stan::math::adj_jac_apply<RowVectorSinFunctor>(x2);

  test::check_varis_on_stack(y1);
  test::check_varis_on_stack(y2);
}

TEST(AgradRev, test_row_vector_sin_values) {
  Eigen::Matrix<stan::math::var, 1, Eigen::Dynamic> x1(1), x2(2), y1(1), y2(2);
  x1 << 1.0;
  x2 << 2.0, 1.0;

  y1 = stan::math::adj_jac_apply<RowVectorSinFunctor>(x1);
  y2 = stan::math::adj_jac_apply<RowVectorSinFunctor>(x2);

  EXPECT_FLOAT_EQ(y1(0).val(), 0.841470984807897);
  EXPECT_FLOAT_EQ(y2(0).val(), 0.909297426825682);
  EXPECT_FLOAT_EQ(y2(1).val(), 0.841470984807897);
}

TEST(AgradRev, test_row_vector_sin_multiple_jac) {
  Eigen::Matrix<stan::math::var, 1, Eigen::Dynamic> x1(1), x2(2), y1(1), y2(2);
  x1 << 1.0;
  x2 << 2.0, 1.0;

  y1 = stan::math::adj_jac_apply<RowVectorSinFunctor>(x1);
  y2 = stan::math::adj_jac_apply<RowVectorSinFunctor>(x2);

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
struct MatrixSinFunctor {
  int N_;
  int M_;
  double* x_mem_;
  template <std::size_t size>
  Eigen::MatrixXd operator()(const std::array<bool, size>& needs_adj,
                             const Eigen::MatrixXd& x) {
    N_ = x.rows();
    M_ = x.cols();
    Eigen::MatrixXd out(N_, M_);
    x_mem_
        = stan::math::ChainableStack::instance().memalloc_.alloc_array<double>(
            N_ * M_);

    for (int n = 0; n < N_ * M_; ++n) {
      x_mem_[n] = x(n);
      out(n) = sin(x(n));
    }

    return out;
  }

  template <std::size_t size>
  auto multiply_adjoint_jacobian(const std::array<bool, size>& needs_adj,
                                 const Eigen::MatrixXd& adj) {
    Eigen::MatrixXd out(N_, M_);

    for (int n = 0; n < N_ * M_; ++n) {
      out(n) = cos(x_mem_[n]) * adj(n);
    }

    return std::make_tuple(out);
  }
};

TEST(AgradRev, test_matrix_sin_stack) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> x(2, 2),
      y(2, 2);
  x << 2.0, 1.0, 0.0, -1.0;

  y = stan::math::adj_jac_apply<MatrixSinFunctor>(x);

  test::check_varis_on_stack(y);
}

TEST(AgradRev, test_matrix_sin_values) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> x(2, 2),
      y(2, 2);
  x << 2.0, 1.0, 0.0, -1.0;

  y = stan::math::adj_jac_apply<MatrixSinFunctor>(x);

  EXPECT_FLOAT_EQ(y(0, 0).val(), 0.909297426825682);
  EXPECT_FLOAT_EQ(y(0, 1).val(), 0.841470984807897);
  EXPECT_FLOAT_EQ(y(1, 0).val(), 0.0);
  EXPECT_FLOAT_EQ(y(1, 1).val(), -0.841470984807897);
}

TEST(AgradRev, test_matrix_sin_multiple_jac) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> x(2, 2),
      y(2, 2);
  x << 2.0, 1.0, 0.0, -1.0;

  y = stan::math::adj_jac_apply<MatrixSinFunctor>(x);

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
  template <size_t size>
  Eigen::VectorXd operator()(
      std::array<bool, size> needs_adj, double, int, const double&, const int&,
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

  template <size_t size>
  auto multiply_adjoint_jacobian(const std::array<bool, size>& needs_adj,
                                 const Eigen::VectorXd& y_adj) {
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
  auto vi = new stan::math::adj_jac_vari<F, Targs...>();

  (*vi)(args...);

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

  EXPECT_EQ(vi1->offsets_, (std::array<int, 16>({{0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                  0, 0, 0, 0, 0, 0}})));

  stan::math::var(vi1).grad();

  auto vi2 = make_vari_for_test<WeirdArgumentListFunctor1>(
      v, i, d, i, vv, vi, vd, vi, ev1, ed2, ev3, ed4, ev1, ed2, ev3, ed4);

  EXPECT_EQ(vi2->is_var_,
            (std::array<bool, 16>(
                {{true, false, false, false, true, false, false, false, true,
                  false, true, false, true, false, true, false}})));

  EXPECT_EQ(vi2->offsets_, (std::array<int, 16>({{0, 1, 1, 1, 1, 3, 3, 3, 3, 6,
                                                  6, 16, 16, 19, 19, 29}})));

  stan::math::var(vi2).grad();

  auto vi3 = make_vari_for_test<WeirdArgumentListFunctor1>(
      d, i, v, i, vd, vi, vv, vi, ed1, ev2, ed3, ev4, ed1, ev2, ed3, ev4);

  EXPECT_EQ(vi3->is_var_,
            (std::array<bool, 16>(
                {{false, false, true, false, false, false, true, false, false,
                  true, false, true, false, true, false, true}})));

  EXPECT_EQ(vi3->offsets_, (std::array<int, 16>({{0, 0, 0, 1, 1, 1, 1, 3, 3, 3,
                                                  11, 11, 16, 16, 24, 24}})));

  stan::math::var(vi3).grad();

  auto vi4 = make_vari_for_test<WeirdArgumentListFunctor1>(
      v, i, d, i, vd, vi, vv, vi, ev1, ed2, ed3, ev4, ev1, ed2, ed3, ev4);

  EXPECT_EQ(vi4->is_var_,
            (std::array<bool, 16>(
                {{true, false, false, false, false, false, true, false, true,
                  false, false, true, true, false, false, true}})));

  EXPECT_EQ(vi4->offsets_, (std::array<int, 16>({{0, 1, 1, 1, 1, 1, 1, 3, 3, 6,
                                                  6, 6, 11, 14, 14, 14}})));

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
  template <size_t size>
  Eigen::VectorXd operator()(
      std::array<bool, size> needs_adj, const double& d,
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

  template <size_t size>
  auto multiply_adjoint_jacobian(const std::array<bool, size>& needs_adj,
                                 const Eigen::VectorXd& y_adj) {
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
struct SinCosFunctor {
  int N_;
  double* x1_mem_;
  double* x4_mem_;

  template <std::size_t size>
  Eigen::VectorXd operator()(const std::array<bool, size>& needs_adj,
                             const Eigen::VectorXd& x1, const int& x2,
                             const std::vector<int>& x3,
                             const std::vector<double>& x4) {
    stan::math::check_matching_sizes("SinCosFunctor", "x1", x1, "x4", x4);
    N_ = x1.size();
    Eigen::VectorXd out(N_);

    if (needs_adj[0]) {
      x1_mem_ = stan::math::ChainableStack::instance()
                    .memalloc_.alloc_array<double>(N_);
      std::copy(x1.data(), x1.data() + N_, x1_mem_);
    }

    EXPECT_FALSE(needs_adj[1]);
    EXPECT_FALSE(needs_adj[2]);

    if (needs_adj[3]) {
      x4_mem_ = stan::math::ChainableStack::instance()
                    .memalloc_.alloc_array<double>(N_);
      std::copy(x4.data(), x4.data() + N_, x4_mem_);
    }

    for (int n = 0; n < N_; ++n) {
      out(n) = sin(x1(n)) + cos(x4[n]);
    }

    return out;
  }

  template <std::size_t size>
  auto multiply_adjoint_jacobian(const std::array<bool, size>& needs_adj,
                                 const Eigen::VectorXd& adj) {
    Eigen::VectorXd out1;
    std::vector<double> out4;

    if (needs_adj[0]) {
      out1.resize(N_);
      for (int n = 0; n < N_; ++n) {
        out1(n) = cos(x1_mem_[n]) * adj(n);
      }
    }

    EXPECT_FALSE(needs_adj[1]);
    EXPECT_FALSE(needs_adj[2]);

    if (needs_adj[3]) {
      out4.resize(N_);
      for (int n = 0; n < N_; ++n) {
        out4[n] = -sin(x4_mem_[n]) * adj(n);
      }
    }

    return std::make_tuple(out1, 0, std::vector<int>(), out4);
  }
};

TEST(AgradRev, test_sincos_stack) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> x11(1), x21(2), y1(1),
      y2(2);
  std::vector<stan::math::var> x12(1), x22(2);
  x11 << 1.0;
  x12[0] = 5.0;
  x21 << 2.0, 1.0;
  x22[0] = 5.0;
  x22[1] = 3.0;

  y1 = stan::math::adj_jac_apply<SinCosFunctor>(x11, 0, std::vector<int>(5, 0),
                                                x12);
  y2 = stan::math::adj_jac_apply<SinCosFunctor>(x21, 0, std::vector<int>(5, 0),
                                                x22);

  test::check_varis_on_stack(y1);
  test::check_varis_on_stack(y2);
}

TEST(AgradRev, test_sincos_values) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> x11(1), x21(2), y1(1),
      y2(2);
  std::vector<stan::math::var> x12(1), x22(2);
  x11 << 1.0;
  x12[0] = 5.0;
  x21 << 2.0, 1.0;
  x22[0] = 5.0;
  x22[1] = 3.0;

  y1 = stan::math::adj_jac_apply<SinCosFunctor>(x11, 0, std::vector<int>(5, 0),
                                                x12);
  y2 = stan::math::adj_jac_apply<SinCosFunctor>(x21, 0, std::vector<int>(5, 0),
                                                x22);

  EXPECT_FLOAT_EQ(y1(0).val(), 1.125133170271123);
  EXPECT_FLOAT_EQ(y2(0).val(), 1.192959612288908);
  EXPECT_FLOAT_EQ(y2(1).val(), -0.1485215117925489);
}

TEST(AgradRev, test_sincos_multiple_jac_vv) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> x11(1), x21(2), y1(1),
      y2(2);
  std::vector<stan::math::var> x12(1), x22(2);
  x11 << 1.0;
  x12[0] = 5.0;
  x21 << 2.0, 1.0;
  x22[0] = 5.0;
  x22[1] = 3.0;

  y1 = stan::math::adj_jac_apply<SinCosFunctor>(x11, 0, std::vector<int>(5, 0),
                                                x12);
  y2 = stan::math::adj_jac_apply<SinCosFunctor>(x21, 0, std::vector<int>(5, 0),
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
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> y1(1), y2(2);
  std::vector<stan::math::var> x12(1), x22(2);
  Eigen::Matrix<double, Eigen::Dynamic, 1> x11(1), x21(2);
  x11 << 1.0;
  x12[0] = 5.0;
  x21 << 2.0, 1.0;
  x22[0] = 5.0;
  x22[1] = 3.0;

  y1 = stan::math::adj_jac_apply<SinCosFunctor>(x11, 0, std::vector<int>(5, 0),
                                                x12);
  y2 = stan::math::adj_jac_apply<SinCosFunctor>(x21, 0, std::vector<int>(5, 0),
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
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> x11(1), x21(2), y1(1),
      y2(2);
  std::vector<double> x12(1), x22(2);
  x11 << 1.0;
  x12[0] = 5.0;
  x21 << 2.0, 1.0;
  x22[0] = 5.0;
  x22[1] = 3.0;

  y1 = stan::math::adj_jac_apply<SinCosFunctor>(x11, 0, std::vector<int>(5, 0),
                                                x12);
  y2 = stan::math::adj_jac_apply<SinCosFunctor>(x21, 0, std::vector<int>(5, 0),
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
struct SinCosFunctor2 {
  int N_;
  double* x1_mem_;
  double x2_;

  template <std::size_t size>
  Eigen::VectorXd operator()(const std::array<bool, size>& needs_adj,
                             const Eigen::VectorXd& x1, const double& x2) {
    N_ = x1.size();
    Eigen::VectorXd out(N_);

    if (needs_adj[0]) {
      x1_mem_ = stan::math::ChainableStack::instance()
                    .memalloc_.alloc_array<double>(N_);
      std::copy(x1.data(), x1.data() + N_, x1_mem_);
    }

    if (needs_adj[1]) {
      x2_ = x2;
    }

    for (int n = 0; n < N_; ++n) {
      out(n) = sin(x1(n)) + cos(x2);
    }

    return out;
  }

  template <std::size_t size>
  auto multiply_adjoint_jacobian(const std::array<bool, size>& needs_adj,
                                 const Eigen::VectorXd& adj) {
    Eigen::VectorXd out1;
    double out2 = 0.0;

    if (needs_adj[0]) {
      out1.resize(N_);
      for (int n = 0; n < N_; ++n) {
        out1(n) = cos(x1_mem_[n]) * adj(n);
      }
    }

    if (needs_adj[1]) {
      for (int n = 0; n < N_; ++n) {
        out2 += -sin(x2_) * adj(n);
      }
    }

    return std::make_tuple(out1, out2);
  }
};

TEST(AgradRev, test_eigen_vector_scalar_stack) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> x11(1), x21(2), y1(1),
      y2(2);
  stan::math::var x12, x22;
  x11 << 1.0;
  x12 = 5.0;
  x21 << 2.0, 1.0;
  x22 = 3.0;

  y1 = stan::math::adj_jac_apply<SinCosFunctor2>(x11, x12);
  y2 = stan::math::adj_jac_apply<SinCosFunctor2>(x21, x22);

  test::check_varis_on_stack(y1);
  test::check_varis_on_stack(y2);
}

TEST(AgradRev, test_eigen_vector_scalar_values) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> x11(1), x21(2), y1(1),
      y2(2);
  stan::math::var x12, x22;
  x11 << 1.0;
  x12 = 5.0;
  x21 << 2.0, 1.0;
  x22 = 3.0;

  y1 = stan::math::adj_jac_apply<SinCosFunctor2>(x11, x12);
  y2 = stan::math::adj_jac_apply<SinCosFunctor2>(x21, x22);

  EXPECT_FLOAT_EQ(y1(0).val(), 1.125133170271123);
  EXPECT_FLOAT_EQ(y2(0).val(), -0.0806950697747637);
  EXPECT_FLOAT_EQ(y2(1).val(), -0.1485215117925489);
}

TEST(AgradRev, test_eigen_vector_scalar_multiple_jac_vv) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> x11(1), x21(2), y1(1),
      y2(2);
  stan::math::var x12, x22;
  x11 << 1.0;
  x12 = 5.0;
  x21 << 2.0, 1.0;
  x22 = 3.0;

  y1 = stan::math::adj_jac_apply<SinCosFunctor2>(x11, x12);
  y2 = stan::math::adj_jac_apply<SinCosFunctor2>(x21, x22);

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
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> y1(1), y2(2);
  Eigen::Matrix<double, Eigen::Dynamic, 1> x11(1), x21(2);
  stan::math::var x12, x22;
  x11 << 1.0;
  x12 = 5.0;
  x21 << 2.0, 1.0;
  x22 = 3.0;

  y1 = stan::math::adj_jac_apply<SinCosFunctor2>(x11, x12);
  y2 = stan::math::adj_jac_apply<SinCosFunctor2>(x21, x22);

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
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> x11(1), x21(2), y1(1),
      y2(2);
  double x12, x22;
  x11 << 1.0;
  x12 = 5.0;
  x21 << 2.0, 1.0;
  x22 = 3.0;

  y1 = stan::math::adj_jac_apply<SinCosFunctor2>(x11, x12);
  y2 = stan::math::adj_jac_apply<SinCosFunctor2>(x21, x22);

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
 */
struct SinCosFunctor3 {
  int N_;
  double* x1_mem_;
  double x2_;

  template <std::size_t size>
  Eigen::VectorXd operator()(const std::array<bool, size>& needs_adj,
                             const double& x2, const Eigen::VectorXd& x1) {
    N_ = x1.size();
    Eigen::VectorXd out(N_);

    if (needs_adj[1]) {
      x1_mem_ = stan::math::ChainableStack::instance()
                    .memalloc_.alloc_array<double>(N_);
      std::copy(x1.data(), x1.data() + N_, x1_mem_);
    }

    if (needs_adj[0]) {
      x2_ = x2;
    }

    for (int n = 0; n < N_; ++n) {
      out(n) = sin(x1(n)) + cos(x2);
    }

    return out;
  }

  template <std::size_t size>
  auto multiply_adjoint_jacobian(const std::array<bool, size>& needs_adj,
                                 const Eigen::VectorXd& adj) {
    Eigen::VectorXd out1;
    double out2 = 0.0;

    if (needs_adj[1]) {
      out1.resize(N_);
      for (int n = 0; n < N_; ++n) {
        out1(n) = cos(x1_mem_[n]) * adj(n);
      }
    }

    if (needs_adj[0]) {
      for (int n = 0; n < N_; ++n) {
        out2 += -sin(x2_) * adj(n);
      }
    }

    return std::make_tuple(out2, out1);
  }
};

TEST(AgradRev, test_sincos_scalar_eigen_vector_stack) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> x11(1), x21(2), y1(1),
      y2(2);
  stan::math::var x12, x22;
  x11 << 1.0;
  x12 = 5.0;
  x21 << 2.0, 1.0;
  x22 = 3.0;

  y1 = stan::math::adj_jac_apply<SinCosFunctor3>(x12, x11);
  y2 = stan::math::adj_jac_apply<SinCosFunctor3>(x22, x21);

  test::check_varis_on_stack(y1);
  test::check_varis_on_stack(y2);
}

TEST(AgradRev, test_sincos_scalar_eigen_vector_values) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> x11(1), x21(2), y1(1),
      y2(2);
  stan::math::var x12, x22;
  x11 << 1.0;
  x12 = 5.0;
  x21 << 2.0, 1.0;
  x22 = 3.0;

  y1 = stan::math::adj_jac_apply<SinCosFunctor3>(x12, x11);
  y2 = stan::math::adj_jac_apply<SinCosFunctor3>(x22, x21);

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

  y1 = stan::math::adj_jac_apply<SinCosFunctor3>(x12, x11);
  y2 = stan::math::adj_jac_apply<SinCosFunctor3>(x22, x21);

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

TEST(AgradRev, test_sincos_scalar_eigen_vector_multiple_jac_vd) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> y1(1), y2(2);
  Eigen::Matrix<double, Eigen::Dynamic, 1> x11(1), x21(2);
  stan::math::var x12, x22;
  x11 << 1.0;
  x12 = 5.0;
  x21 << 2.0, 1.0;
  x22 = 3.0;

  y1 = stan::math::adj_jac_apply<SinCosFunctor3>(x12, x11);
  y2 = stan::math::adj_jac_apply<SinCosFunctor3>(x22, x21);

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
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> x11(1), x21(2), y1(1),
      y2(2);
  double x12, x22;
  x11 << 1.0;
  x12 = 5.0;
  x21 << 2.0, 1.0;
  x22 = 3.0;

  y1 = stan::math::adj_jac_apply<SinCosFunctor3>(x12, x11);
  y2 = stan::math::adj_jac_apply<SinCosFunctor3>(x22, x21);

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
