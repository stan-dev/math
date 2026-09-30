#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/mix/fun/multiply_operator_util.hpp>

// Eigen `operator*` matrix products of autodiff matrices are evaluated by
// Stan's product hooks. Each test checks that the hook is selected for every
// autodiff scalar and that values and derivatives are correct.

using stan::test::multiply_operator::expect_assignment_hooked;
using stan::test::multiply_operator::expect_hooked;
using stan::test::multiply_operator::test_matrix;

namespace {
auto product = [](const auto& a, const auto& b) { return a * b; };
auto product_eval = [](const auto& a, const auto& b) { return (a * b).eval(); };
}  // namespace

TEST(mathMixMultiplyOperator, matrixMatrix) {
  Eigen::MatrixXd a = test_matrix(3, 4, 0.1);
  Eigen::MatrixXd b = test_matrix(4, 2, -0.2);
  expect_hooked(product, a, b);
  stan::test::expect_ad(product_eval, a, b);
  Eigen::MatrixXd empty(0, 4);
  stan::test::expect_ad(product_eval, empty, b);
}

TEST(mathMixMultiplyOperator, matrixVector) {
  Eigen::MatrixXd a = test_matrix(3, 4, 0.1);
  Eigen::VectorXd v = test_matrix(4, 1, 0.3);
  expect_hooked(product, a, v);
  stan::test::expect_ad(product_eval, a, v);
}

TEST(mathMixMultiplyOperator, rowVectorMatrix) {
  Eigen::RowVectorXd rv = test_matrix(1, 3, 0.2);
  Eigen::MatrixXd a = test_matrix(3, 4, 0.1);
  expect_hooked(product, rv, a);
  stan::test::expect_ad(product_eval, rv, a);
}

TEST(mathMixMultiplyOperator, innerAndOuterProduct) {
  Eigen::RowVectorXd rv = test_matrix(1, 4, -0.1);
  Eigen::VectorXd v = test_matrix(4, 1, 0.3);
  expect_hooked(product, rv, v);
  expect_hooked(product, v, rv);
  stan::test::expect_ad(product_eval, rv, v);
  stan::test::expect_ad(product_eval, v, rv);
}

TEST(mathMixMultiplyOperator, fixedSize) {
  // 2x2 uses Eigen's lazy coefficient-based product, 4x4 the
  // coefficient-based one
  expect_hooked(product, Eigen::Matrix2d(test_matrix(2, 2, 0.1)),
                Eigen::Matrix2d(test_matrix(2, 2, -0.3)));
  expect_hooked(product, Eigen::Matrix4d(test_matrix(4, 4, 0.1)),
                Eigen::Matrix4d(test_matrix(4, 4, -0.3)));
  auto f_fixed = [](const auto& a, const auto& b) {
    using a_scalar = typename std::decay_t<decltype(a)>::Scalar;
    using b_scalar = typename std::decay_t<decltype(b)>::Scalar;
    using ret_scalar = stan::return_type_t<a_scalar, b_scalar>;
    const Eigen::Matrix<a_scalar, 2, 2> a2 = a;
    const Eigen::Matrix<b_scalar, 2, 2> b2 = b;
    const Eigen::Matrix<ret_scalar, 2, 2> c = a2 * b2;
    return Eigen::Matrix<ret_scalar, -1, -1>(c);
  };
  stan::test::expect_ad(f_fixed, test_matrix(2, 2, 0.1),
                        test_matrix(2, 2, -0.3));
}

TEST(mathMixMultiplyOperator, lazyProduct) {
  auto lazy = [](const auto& a, const auto& b) { return a.lazyProduct(b); };
  auto lazy_eval
      = [](const auto& a, const auto& b) { return a.lazyProduct(b).eval(); };
  Eigen::MatrixXd a = test_matrix(2, 3, 0.1);
  Eigen::MatrixXd b = test_matrix(3, 2, -0.2);
  expect_hooked(lazy, a, b);
  stan::test::expect_ad(lazy_eval, a, b);
}

TEST(mathMixMultiplyOperator, expressionOperands) {
  auto f_t = [](const auto& a, const auto& b) {
    return a.transpose() * b.block(0, 0, 3, 2);
  };
  auto f_t_eval
      = [&f_t](const auto& a, const auto& b) { return f_t(a, b).eval(); };
  Eigen::MatrixXd a = test_matrix(3, 4, 0.1);
  Eigen::MatrixXd b = test_matrix(3, 3, -0.2);
  expect_hooked(f_t, a, b);
  stan::test::expect_ad(f_t_eval, a, b);
  auto f_chain = [](const auto& a, const auto& b) {
    return (a * b * b.transpose()).eval();
  };
  stan::test::expect_ad(f_chain, a, test_matrix(4, 2, 0.4));
}

TEST(mathMixMultiplyOperator, assignment) {
  expect_assignment_hooked();
  auto f_add_sub = [](const auto& a, const auto& b) {
    auto c = (a * b).eval();
    c += a * b;
    c -= a * b.transpose();
    return c;
  };
  auto f_noalias = [](const auto& a, const auto& b) {
    Eigen::Matrix<stan::return_type_t<decltype(a), decltype(b)>, -1, -1> c(
        a.rows(), b.cols());
    c.noalias() = a * b;
    return c;
  };
  Eigen::MatrixXd a = test_matrix(3, 3, 0.1);
  Eigen::MatrixXd b = test_matrix(3, 3, -0.2);
  stan::test::expect_ad(f_add_sub, a, b);
  stan::test::expect_ad(f_noalias, a, b);
}

TEST(mathMixMultiplyOperator, triangularTimesVector) {
  // Eigen's triangular kernels call the matrix-vector kernel directly. They
  // do not support mixing autodiff and double scalars, so both operands are
  // promoted to the common scalar.
  auto f = [](const auto& a, const auto& v) {
    using scalar = stan::return_type_t<decltype(a), decltype(v)>;
    const Eigen::Matrix<scalar, -1, -1> a_s = a;
    const Eigen::Matrix<scalar, -1, 1> v_s = v;
    return (a_s.template triangularView<Eigen::Lower>() * v_s).eval();
  };
  Eigen::MatrixXd a = test_matrix(9, 9, 0.1);
  Eigen::VectorXd v = test_matrix(9, 1, 0.3);
  stan::test::expect_ad(f, a, v);
}
