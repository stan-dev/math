#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/mix/fun/multiply_operator_util.hpp>

// Values and derivatives of Eigen `operator*` matrix products of autodiff
// matrices. Which Stan function computes each form is tested in
// test/unit/math/rev/core/eigen_product_hook_test.cpp.

using stan::test::multiply_operator::test_matrix;

namespace {
// Eigen cannot multiply fvar and double matrices until fvar products are
// routed to Stan's functions, so these tests run in reverse mode only.
constexpr bool reverse_only = true;
auto product_eval = [](const auto& a, const auto& b) { return (a * b).eval(); };
}  // namespace

TEST(mathMixMultiplyOperator, matrixMatrix) {
  Eigen::MatrixXd a = test_matrix(3, 4, 0.1);
  Eigen::MatrixXd b = test_matrix(4, 2, -0.2);
  stan::test::expect_ad<reverse_only>(product_eval, a, b);
  Eigen::MatrixXd empty(0, 4);
  stan::test::expect_ad<reverse_only>(product_eval, empty, b);
}

TEST(mathMixMultiplyOperator, matrixVector) {
  Eigen::MatrixXd a = test_matrix(3, 4, 0.1);
  Eigen::VectorXd v = test_matrix(4, 1, 0.3);
  stan::test::expect_ad<reverse_only>(product_eval, a, v);
}

TEST(mathMixMultiplyOperator, rowVectorMatrix) {
  Eigen::RowVectorXd rv = test_matrix(1, 3, 0.2);
  Eigen::MatrixXd a = test_matrix(3, 4, 0.1);
  stan::test::expect_ad<reverse_only>(product_eval, rv, a);
}

TEST(mathMixMultiplyOperator, innerAndOuterProduct) {
  Eigen::RowVectorXd rv = test_matrix(1, 4, -0.1);
  Eigen::VectorXd v = test_matrix(4, 1, 0.3);
  stan::test::expect_ad<reverse_only>(product_eval, rv, v);
  stan::test::expect_ad<reverse_only>(product_eval, v, rv);
}

TEST(mathMixMultiplyOperator, fixedSize) {
  // 2x2 uses Eigen's lazy coefficient-based product
  auto f_fixed = [](const auto& a, const auto& b) {
    using a_scalar = typename std::decay_t<decltype(a)>::Scalar;
    using b_scalar = typename std::decay_t<decltype(b)>::Scalar;
    using ret_scalar = stan::return_type_t<a_scalar, b_scalar>;
    const Eigen::Matrix<a_scalar, 2, 2> a2 = a;
    const Eigen::Matrix<b_scalar, 2, 2> b2 = b;
    const Eigen::Matrix<ret_scalar, 2, 2> c = a2 * b2;
    return Eigen::Matrix<ret_scalar, -1, -1>(c);
  };
  stan::test::expect_ad<reverse_only>(f_fixed, test_matrix(2, 2, 0.1),
                                      test_matrix(2, 2, -0.3));
}

TEST(mathMixMultiplyOperator, lazyProduct) {
  auto lazy_eval
      = [](const auto& a, const auto& b) { return a.lazyProduct(b).eval(); };
  Eigen::MatrixXd a = test_matrix(2, 3, 0.1);
  Eigen::MatrixXd b = test_matrix(3, 2, -0.2);
  stan::test::expect_ad<reverse_only>(lazy_eval, a, b);
}

TEST(mathMixMultiplyOperator, expressionOperands) {
  auto f_t = [](const auto& a, const auto& b) {
    return (a.transpose() * b.block(0, 0, 3, 2)).eval();
  };
  Eigen::MatrixXd a = test_matrix(3, 4, 0.1);
  Eigen::MatrixXd b = test_matrix(3, 3, -0.2);
  stan::test::expect_ad<reverse_only>(f_t, a, b);
  auto f_chain = [](const auto& a, const auto& b) {
    return (a * b * b.transpose()).eval();
  };
  stan::test::expect_ad<reverse_only>(f_chain, a, test_matrix(4, 2, 0.4));
}

TEST(mathMixMultiplyOperator, assignment) {
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
    c.noalias() += a * b;
    return c;
  };
  Eigen::MatrixXd a = test_matrix(3, 3, 0.1);
  Eigen::MatrixXd b = test_matrix(3, 3, -0.2);
  stan::test::expect_ad<reverse_only>(f_add_sub, a, b);
  stan::test::expect_ad<reverse_only>(f_noalias, a, b);
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
