#include <stan/math/rev/core.hpp>
#include <test/unit/math/rev/mat/util.hpp>
#include <gtest/gtest.h>
#include <sstream>

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

TEST(AgradRev, test_sin_stack) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> x1(1), x2(2), y1(1), y2(2);
  x1 << 1.0;
  x2 << 2.0, 1.0;

  y1 = stan::math::adj_jac_apply<SinFunctor>(x1);
  y2 = stan::math::adj_jac_apply<SinFunctor>(x2);

  test::check_varis_on_stack(y1);
  test::check_varis_on_stack(y2);
}

TEST(AgradRev, test_sin_values) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> x1(1), x2(2), y1(1), y2(2);
  x1 << 1.0;
  x2 << 2.0, 1.0;

  y1 = stan::math::adj_jac_apply<SinFunctor>(x1);
  y2 = stan::math::adj_jac_apply<SinFunctor>(x2);

  EXPECT_NEAR(y1(0).val(), 0.841470984807897, 1e-10);
  EXPECT_NEAR(y2(0).val(), 0.909297426825682, 1e-10);
  EXPECT_NEAR(y2(1).val(), 0.841470984807897, 1e-10);
}

TEST(AgradRev, test_sin_multiple_jac) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> x1(1), x2(2), y1(1), y2(2);
  x1 << 1.0;
  x2 << 2.0, 1.0;

  y1 = stan::math::adj_jac_apply<SinFunctor>(x1);
  y2 = stan::math::adj_jac_apply<SinFunctor>(x2);

  y1(0).grad();
  EXPECT_NEAR(x1(0).adj(), 0.5403023058681398, 1e-10);
  EXPECT_NEAR(x2(0).adj(), 0.0, 1e-10);
  EXPECT_NEAR(x2(1).adj(), 0.0, 1e-10);

  stan::math::set_zero_all_adjoints();

  y2(0).grad();
  EXPECT_NEAR(x1(0).adj(), 0.0, 1e-10);
  EXPECT_NEAR(x2(0).adj(), -0.4161468365471424, 1e-10);
  EXPECT_NEAR(x2(1).adj(), 0.0, 1e-10);

  stan::math::set_zero_all_adjoints();

  y2(1).grad();
  EXPECT_NEAR(x1(0).adj(), 0.0, 1e-10);
  EXPECT_NEAR(x2(0).adj(), 0.0, 1e-10);
  EXPECT_NEAR(x2(1).adj(), 0.5403023058681398, 1e-10);

  stan::math::set_zero_all_adjoints();

  stan::math::var sum_y2 = (1.73 * y2(0) + 1.57 * y2(1));
  sum_y2.grad();
  EXPECT_NEAR(x1(0).adj(), 0.0, 1e-10);
  EXPECT_NEAR(x2(0).adj(), 1.73 * -0.4161468365471424, 1e-10);
  EXPECT_NEAR(x2(1).adj(), 1.57 * 0.5403023058681398, 1e-10);
}

struct SinCosFunctor {
  int N_;
  double* x1_mem_;
  double* x2_mem_;

  template <std::size_t size>
  Eigen::VectorXd operator()(const std::array<bool, size>& needs_adj,
                             const Eigen::VectorXd& x1,
                             const Eigen::VectorXd& x2) {
    stan::math::check_matching_sizes("SinCosFunctor", "x1", x1, "x2", x2);
    N_ = x1.size();
    Eigen::VectorXd out(N_);

    if (needs_adj[0]) {
      x1_mem_ = stan::math::ChainableStack::instance()
                    .memalloc_.alloc_array<double>(N_);
      std::copy(x1.data(), x1.data() + N_, x1_mem_);
    }

    if (needs_adj[1]) {
      x2_mem_ = stan::math::ChainableStack::instance()
                    .memalloc_.alloc_array<double>(N_);
      std::copy(x2.data(), x2.data() + N_, x2_mem_);
    }

    for (int n = 0; n < N_; ++n) {
      out(n) = sin(x1(n)) + cos(x2(n));
    }

    return out;
  }

  template <std::size_t size>
  auto multiply_adjoint_jacobian(const std::array<bool, size>& needs_adj,
                                 const Eigen::VectorXd& adj) {
    Eigen::VectorXd out1;
    Eigen::VectorXd out2;

    if (needs_adj[0]) {
      out1.resize(N_);
      for (int n = 0; n < N_; ++n) {
        out1(n) = cos(x1_mem_[n]) * adj(n);
      }
    }

    if (needs_adj[1]) {
      out2.resize(N_);
      for (int n = 0; n < N_; ++n) {
        out2(n) = -sin(x2_mem_[n]) * adj(n);
      }
    }

    return std::make_tuple(out1, out2);
  }
};

TEST(AgradRev, test_sincos_stack) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> x11(1), x12(1), x21(2),
      x22(2), y1(1), y2(2);
  x11 << 1.0;
  x12 << 5.0;
  x21 << 2.0, 1.0;
  x22 << 5.0, 3.0;

  y1 = stan::math::adj_jac_apply<SinCosFunctor>(x11, x12);
  y2 = stan::math::adj_jac_apply<SinCosFunctor>(x21, x22);

  test::check_varis_on_stack(y1);
  test::check_varis_on_stack(y2);
}

TEST(AgradRev, test_sincos_values) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> x11(1), x12(1), x21(2),
      x22(2), y1(1), y2(2);
  x11 << 1.0;
  x12 << 5.0;
  x21 << 2.0, 1.0;
  x22 << 5.0, 3.0;

  y1 = stan::math::adj_jac_apply<SinCosFunctor>(x11, x12);
  y2 = stan::math::adj_jac_apply<SinCosFunctor>(x21, x22);

  EXPECT_NEAR(y1(0).val(), 1.125133170271123, 1e-10);
  EXPECT_NEAR(y2(0).val(), 1.192959612288908, 1e-10);
  EXPECT_NEAR(y2(1).val(), -0.1485215117925489, 1e-10);
}

TEST(AgradRev, test_sincos_multiple_jac_vv) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> x11(1), x12(1), x21(2),
      x22(2), y1(1), y2(2);
  x11 << 1.0;
  x12 << 5.0;
  x21 << 2.0, 1.0;
  x22 << 5.0, 3.0;

  y1 = stan::math::adj_jac_apply<SinCosFunctor>(x11, x12);
  y2 = stan::math::adj_jac_apply<SinCosFunctor>(x21, x22);

  y1(0).grad();
  EXPECT_NEAR(x11(0).adj(), 0.5403023058681398, 1e-10);
  EXPECT_NEAR(x12(0).adj(), 0.958924274663139, 1e-10);
  EXPECT_NEAR(x21(0).adj(), 0.0, 1e-10);
  EXPECT_NEAR(x21(1).adj(), 0.0, 1e-10);
  EXPECT_NEAR(x22(0).adj(), 0.0, 1e-10);
  EXPECT_NEAR(x22(1).adj(), 0.0, 1e-10);

  stan::math::set_zero_all_adjoints();

  y2(0).grad();
  EXPECT_NEAR(x11(0).adj(), 0.0, 1e-10);
  EXPECT_NEAR(x12(0).adj(), 0.0, 1e-10);
  EXPECT_NEAR(x21(0).adj(), -0.4161468365471424, 1e-10);
  EXPECT_NEAR(x21(1).adj(), 0.0, 1e-10);
  EXPECT_NEAR(x22(0).adj(), 0.958924274663139, 1e-10);
  EXPECT_NEAR(x22(1).adj(), 0.0, 1e-10);

  stan::math::set_zero_all_adjoints();

  y2(1).grad();
  EXPECT_NEAR(x11(0).adj(), 0.0, 1e-10);
  EXPECT_NEAR(x12(0).adj(), 0.0, 1e-10);
  EXPECT_NEAR(x21(0).adj(), 0.0, 1e-10);
  EXPECT_NEAR(x21(1).adj(), 0.5403023058681398, 1e-10);
  EXPECT_NEAR(x22(0).adj(), 0.0, 1e-10);
  EXPECT_NEAR(x22(1).adj(), -0.1411200080598672, 1e-10);

  stan::math::set_zero_all_adjoints();

  stan::math::var sum_y2 = (1.73 * y2(0) + 1.57 * y2(1));
  sum_y2.grad();
  EXPECT_NEAR(x11(0).adj(), 0.0, 1e-10);
  EXPECT_NEAR(x12(0).adj(), 0.0, 1e-10);
  EXPECT_NEAR(x21(0).adj(), 1.73 * -0.4161468365471424, 1e-10);
  EXPECT_NEAR(x21(1).adj(), 1.57 * 0.5403023058681398, 1e-10);
  EXPECT_NEAR(x22(0).adj(), 1.73 * 0.958924274663139, 1e-10);
  EXPECT_NEAR(x22(1).adj(), 1.57 * -0.1411200080598672, 1e-10);
}

TEST(AgradRev, test_sincos_multiple_jac_dv) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> x12(1), x22(2), y1(1),
      y2(2);
  Eigen::Matrix<double, Eigen::Dynamic, 1> x11(1), x21(2);
  x11 << 1.0;
  x12 << 5.0;
  x21 << 2.0, 1.0;
  x22 << 5.0, 3.0;

  y1 = stan::math::adj_jac_apply<SinCosFunctor>(x11, x12);
  y2 = stan::math::adj_jac_apply<SinCosFunctor>(x21, x22);

  y1(0).grad();
  EXPECT_NEAR(x12(0).adj(), 0.958924274663139, 1e-10);
  EXPECT_NEAR(x22(0).adj(), 0.0, 1e-10);
  EXPECT_NEAR(x22(1).adj(), 0.0, 1e-10);

  stan::math::set_zero_all_adjoints();

  y2(0).grad();
  EXPECT_NEAR(x12(0).adj(), 0.0, 1e-10);
  EXPECT_NEAR(x22(0).adj(), 0.958924274663139, 1e-10);
  EXPECT_NEAR(x22(1).adj(), 0.0, 1e-10);

  stan::math::set_zero_all_adjoints();

  y2(1).grad();
  EXPECT_NEAR(x12(0).adj(), 0.0, 1e-10);
  EXPECT_NEAR(x22(0).adj(), 0.0, 1e-10);
  EXPECT_NEAR(x22(1).adj(), -0.1411200080598672, 1e-10);

  stan::math::set_zero_all_adjoints();

  stan::math::var sum_y2 = (1.73 * y2(0) + 1.57 * y2(1));
  sum_y2.grad();
  EXPECT_NEAR(x12(0).adj(), 0.0, 1e-10);
  EXPECT_NEAR(x22(0).adj(), 1.73 * 0.958924274663139, 1e-10);
  EXPECT_NEAR(x22(1).adj(), 1.57 * -0.1411200080598672, 1e-10);
}

TEST(AgradRev, test_sincos_multiple_jac_vd) {
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> x11(1), x21(2), y1(1),
      y2(2);
  Eigen::Matrix<double, Eigen::Dynamic, 1> x12(1), x22(2);
  x11 << 1.0;
  x12 << 5.0;
  x21 << 2.0, 1.0;
  x22 << 5.0, 3.0;

  y1 = stan::math::adj_jac_apply<SinCosFunctor>(x11, x12);
  y2 = stan::math::adj_jac_apply<SinCosFunctor>(x21, x22);

  y1(0).grad();
  EXPECT_NEAR(x11(0).adj(), 0.5403023058681398, 1e-10);
  EXPECT_NEAR(x21(0).adj(), 0.0, 1e-10);
  EXPECT_NEAR(x21(1).adj(), 0.0, 1e-10);

  stan::math::set_zero_all_adjoints();

  y2(0).grad();
  EXPECT_NEAR(x11(0).adj(), 0.0, 1e-10);
  EXPECT_NEAR(x21(0).adj(), -0.4161468365471424, 1e-10);
  EXPECT_NEAR(x21(1).adj(), 0.0, 1e-10);

  stan::math::set_zero_all_adjoints();

  y2(1).grad();
  EXPECT_NEAR(x11(0).adj(), 0.0, 1e-10);
  EXPECT_NEAR(x21(0).adj(), 0.0, 1e-10);
  EXPECT_NEAR(x21(1).adj(), 0.5403023058681398, 1e-10);

  stan::math::set_zero_all_adjoints();

  stan::math::var sum_y2 = (1.73 * y2(0) + 1.57 * y2(1));
  sum_y2.grad();
  EXPECT_NEAR(x11(0).adj(), 0.0, 1e-10);
  EXPECT_NEAR(x21(0).adj(), 1.73 * -0.4161468365471424, 1e-10);
  EXPECT_NEAR(x21(1).adj(), 1.57 * 0.5403023058681398, 1e-10);
}

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

  EXPECT_NEAR(y1(0).val(), 1.125133170271123, 1e-10);
  EXPECT_NEAR(y2(0).val(), -0.0806950697747637, 1e-10);
  EXPECT_NEAR(y2(1).val(), -0.1485215117925489, 1e-10);
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
  EXPECT_NEAR(x11(0).adj(), 0.5403023058681398, 1e-10);
  EXPECT_NEAR(x12.adj(), 0.958924274663139, 1e-10);
  EXPECT_NEAR(x21(0).adj(), 0.0, 1e-10);
  EXPECT_NEAR(x21(1).adj(), 0.0, 1e-10);
  EXPECT_NEAR(x22.adj(), 0.0, 1e-10);

  stan::math::set_zero_all_adjoints();

  y2(0).grad();
  EXPECT_NEAR(x11(0).adj(), 0.0, 1e-10);
  EXPECT_NEAR(x12.adj(), 0.0, 1e-10);
  EXPECT_NEAR(x21(0).adj(), -0.4161468365471424, 1e-10);
  EXPECT_NEAR(x21(1).adj(), 0.0, 1e-10);
  EXPECT_NEAR(x22.adj(), -0.1411200080598672, 1e-10);

  stan::math::set_zero_all_adjoints();

  y2(1).grad();
  EXPECT_NEAR(x11(0).adj(), 0.0, 1e-10);
  EXPECT_NEAR(x12.adj(), 0.0, 1e-10);
  EXPECT_NEAR(x21(0).adj(), 0.0, 1e-10);
  EXPECT_NEAR(x21(1).adj(), 0.5403023058681398, 1e-10);
  EXPECT_NEAR(x22.adj(), -0.1411200080598672, 1e-10);

  stan::math::set_zero_all_adjoints();

  stan::math::var sum_y2 = (1.73 * y2(0) + 1.57 * y2(1));
  sum_y2.grad();
  EXPECT_NEAR(x11(0).adj(), 0.0, 1e-10);
  EXPECT_NEAR(x12.adj(), 0.0, 1e-10);
  EXPECT_NEAR(x21(0).adj(), 1.73 * -0.4161468365471424, 1e-10);
  EXPECT_NEAR(x21(1).adj(), 1.57 * 0.5403023058681398, 1e-10);
  EXPECT_NEAR(x22.adj(), (1.73 + 1.57) * -0.1411200080598672, 1e-10);
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
  EXPECT_NEAR(x12.adj(), 0.958924274663139, 1e-10);
  EXPECT_NEAR(x22.adj(), 0.0, 1e-10);

  stan::math::set_zero_all_adjoints();

  y2(0).grad();
  EXPECT_NEAR(x12.adj(), 0.0, 1e-10);
  EXPECT_NEAR(x22.adj(), -0.1411200080598672, 1e-10);

  stan::math::set_zero_all_adjoints();

  y2(1).grad();
  EXPECT_NEAR(x12.adj(), 0.0, 1e-10);
  EXPECT_NEAR(x22.adj(), -0.1411200080598672, 1e-10);

  stan::math::set_zero_all_adjoints();

  stan::math::var sum_y2 = (1.73 * y2(0) + 1.57 * y2(1));
  sum_y2.grad();
  EXPECT_NEAR(x12.adj(), 0.0, 1e-10);
  EXPECT_NEAR(x22.adj(), (1.73 + 1.57) * -0.1411200080598672, 1e-10);
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
  EXPECT_NEAR(x11(0).adj(), 0.5403023058681398, 1e-10);
  EXPECT_NEAR(x21(0).adj(), 0.0, 1e-10);
  EXPECT_NEAR(x21(1).adj(), 0.0, 1e-10);

  stan::math::set_zero_all_adjoints();

  y2(0).grad();
  EXPECT_NEAR(x11(0).adj(), 0.0, 1e-10);
  EXPECT_NEAR(x21(0).adj(), -0.4161468365471424, 1e-10);
  EXPECT_NEAR(x21(1).adj(), 0.0, 1e-10);

  stan::math::set_zero_all_adjoints();

  y2(1).grad();
  EXPECT_NEAR(x11(0).adj(), 0.0, 1e-10);
  EXPECT_NEAR(x21(0).adj(), 0.0, 1e-10);
  EXPECT_NEAR(x21(1).adj(), 0.5403023058681398, 1e-10);

  stan::math::set_zero_all_adjoints();

  stan::math::var sum_y2 = (1.73 * y2(0) + 1.57 * y2(1));
  sum_y2.grad();
  EXPECT_NEAR(x11(0).adj(), 0.0, 1e-10);
  EXPECT_NEAR(x21(0).adj(), 1.73 * -0.4161468365471424, 1e-10);
  EXPECT_NEAR(x21(1).adj(), 1.57 * 0.5403023058681398, 1e-10);
}

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

  EXPECT_NEAR(y1(0).val(), 1.125133170271123, 1e-10);
  EXPECT_NEAR(y2(0).val(), -0.0806950697747637, 1e-10);
  EXPECT_NEAR(y2(1).val(), -0.1485215117925489, 1e-10);
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
  EXPECT_NEAR(x11(0).adj(), 0.5403023058681398, 1e-10);
  EXPECT_NEAR(x12.adj(), 0.958924274663139, 1e-10);
  EXPECT_NEAR(x21(0).adj(), 0.0, 1e-10);
  EXPECT_NEAR(x21(1).adj(), 0.0, 1e-10);
  EXPECT_NEAR(x22.adj(), 0.0, 1e-10);

  stan::math::set_zero_all_adjoints();

  y2(0).grad();
  EXPECT_NEAR(x11(0).adj(), 0.0, 1e-10);
  EXPECT_NEAR(x12.adj(), 0.0, 1e-10);
  EXPECT_NEAR(x21(0).adj(), -0.4161468365471424, 1e-10);
  EXPECT_NEAR(x21(1).adj(), 0.0, 1e-10);
  EXPECT_NEAR(x22.adj(), -0.1411200080598672, 1e-10);

  stan::math::set_zero_all_adjoints();

  y2(1).grad();
  EXPECT_NEAR(x11(0).adj(), 0.0, 1e-10);
  EXPECT_NEAR(x12.adj(), 0.0, 1e-10);
  EXPECT_NEAR(x21(0).adj(), 0.0, 1e-10);
  EXPECT_NEAR(x21(1).adj(), 0.5403023058681398, 1e-10);
  EXPECT_NEAR(x22.adj(), -0.1411200080598672, 1e-10);

  stan::math::set_zero_all_adjoints();

  stan::math::var sum_y2 = (1.73 * y2(0) + 1.57 * y2(1));
  sum_y2.grad();
  EXPECT_NEAR(x11(0).adj(), 0.0, 1e-10);
  EXPECT_NEAR(x12.adj(), 0.0, 1e-10);
  EXPECT_NEAR(x21(0).adj(), 1.73 * -0.4161468365471424, 1e-10);
  EXPECT_NEAR(x21(1).adj(), 1.57 * 0.5403023058681398, 1e-10);
  EXPECT_NEAR(x22.adj(), (1.73 + 1.57) * -0.1411200080598672, 1e-10);
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
  EXPECT_NEAR(x12.adj(), 0.958924274663139, 1e-10);
  EXPECT_NEAR(x22.adj(), 0.0, 1e-10);

  stan::math::set_zero_all_adjoints();

  y2(0).grad();
  EXPECT_NEAR(x12.adj(), 0.0, 1e-10);
  EXPECT_NEAR(x22.adj(), -0.1411200080598672, 1e-10);

  stan::math::set_zero_all_adjoints();

  y2(1).grad();
  EXPECT_NEAR(x12.adj(), 0.0, 1e-10);
  EXPECT_NEAR(x22.adj(), -0.1411200080598672, 1e-10);

  stan::math::set_zero_all_adjoints();

  stan::math::var sum_y2 = (1.73 * y2(0) + 1.57 * y2(1));
  sum_y2.grad();
  EXPECT_NEAR(x12.adj(), 0.0, 1e-10);
  EXPECT_NEAR(x22.adj(), (1.73 + 1.57) * -0.1411200080598672, 1e-10);
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
  EXPECT_NEAR(x11(0).adj(), 0.5403023058681398, 1e-10);
  EXPECT_NEAR(x21(0).adj(), 0.0, 1e-10);
  EXPECT_NEAR(x21(1).adj(), 0.0, 1e-10);

  stan::math::set_zero_all_adjoints();

  y2(0).grad();
  EXPECT_NEAR(x11(0).adj(), 0.0, 1e-10);
  EXPECT_NEAR(x21(0).adj(), -0.4161468365471424, 1e-10);
  EXPECT_NEAR(x21(1).adj(), 0.0, 1e-10);

  stan::math::set_zero_all_adjoints();

  y2(1).grad();
  EXPECT_NEAR(x11(0).adj(), 0.0, 1e-10);
  EXPECT_NEAR(x21(0).adj(), 0.0, 1e-10);
  EXPECT_NEAR(x21(1).adj(), 0.5403023058681398, 1e-10);

  stan::math::set_zero_all_adjoints();

  stan::math::var sum_y2 = (1.73 * y2(0) + 1.57 * y2(1));
  sum_y2.grad();
  EXPECT_NEAR(x11(0).adj(), 0.0, 1e-10);
  EXPECT_NEAR(x21(0).adj(), 1.73 * -0.4161468365471424, 1e-10);
  EXPECT_NEAR(x21(1).adj(), 1.57 * 0.5403023058681398, 1e-10);
}
