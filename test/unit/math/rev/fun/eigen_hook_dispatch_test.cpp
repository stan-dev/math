#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <complex>
#include <type_traits>
#include <utility>

// Checks which Eigen struct each product expression instantiates. Stan's
// specializations of Eigen's evaluators and assignments derive from
// `stan::math::internal::eigen_ad_hook`, so a type that derives from it is
// handled by Stan's functions. Each type is built as Eigen builds it: the
// expression types come from `decltype` of the expression, the assignment
// source is that type without `const` (Eigen deduces it from a `const Src&`
// argument), and the assignment functor is `Op<Dst::Scalar, Src::Scalar>`.

namespace {

using var = stan::math::var;
template <typename T>
using mat = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

template <typename T>
constexpr bool is_hooked
    = std::is_base_of<stan::math::internal::eigen_ad_hook, T>::value;

template <typename Xpr>
constexpr bool evaluator_hooked = is_hooked<Eigen::internal::evaluator<Xpr>>;

template <template <typename, typename> class Op, typename Dst, typename Src>
constexpr bool assignment_hooked = is_hooked<Eigen::internal::Assignment<
    Dst, std::decay_t<Src>, Op<typename Dst::Scalar, typename Src::Scalar>>>;

template <template <typename, typename> class Op, typename Dst, typename Src>
void expect_assignment_hooked() {
  EXPECT_TRUE((assignment_hooked<Op, Dst, Src>));
  EXPECT_TRUE((assignment_hooked<Op, Eigen::Map<Dst>, Src>));
  EXPECT_TRUE((assignment_hooked<Op, Eigen::Block<Dst>, Src>));
}

template <typename L, typename R>
struct operands {
  using lhs = L;
  using rhs = R;
};

template <typename T>
class EigenHookDispatch : public ::testing::Test {};

using ad_operands = ::testing::Types<
    operands<mat<var>, mat<var>>, operands<mat<var>, mat<double>>,
    operands<mat<double>, mat<var>>, operands<mat<var>, mat<int>>,
    operands<mat<int>, mat<var>>>;
TYPED_TEST_SUITE(EigenHookDispatch, ad_operands);

}  // namespace

TYPED_TEST(EigenHookDispatch, matrix_product_evaluator) {
  using L = typename TypeParam::lhs;
  using R = typename TypeParam::rhs;
  using product_t
      = decltype(std::declval<const L&>() * std::declval<const R&>());
  EXPECT_TRUE(evaluator_hooked<product_t>);
}

TYPED_TEST(EigenHookDispatch, matrix_product_assignment) {
  using L = typename TypeParam::lhs;
  using R = typename TypeParam::rhs;
  using product_t
      = decltype(std::declval<const L&>() * std::declval<const R&>());
  using alias_free_t = Eigen::Product<L, R, Eigen::AliasFreeProduct>;
  using Eigen::internal::add_assign_op;
  using Eigen::internal::assign_op;
  using Eigen::internal::sub_assign_op;
  expect_assignment_hooked<assign_op, mat<var>, product_t>();
  expect_assignment_hooked<add_assign_op, mat<var>, product_t>();
  expect_assignment_hooked<sub_assign_op, mat<var>, product_t>();
  expect_assignment_hooked<assign_op, mat<var>, alias_free_t>();
  expect_assignment_hooked<add_assign_op, mat<var>, alias_free_t>();
  expect_assignment_hooked<sub_assign_op, mat<var>, alias_free_t>();
}

TYPED_TEST(EigenHookDispatch, cwise_product_evaluator) {
  using L = typename TypeParam::lhs;
  using R = typename TypeParam::rhs;
  using array_product_t = decltype(std::declval<const L&>().array()
                                   * std::declval<const R&>().array());
  using scalar_times_matrix_t
      = decltype(std::declval<typename L::Scalar>() * std::declval<const R&>());
  using matrix_times_scalar_t
      = decltype(std::declval<const L&>() * std::declval<typename R::Scalar>());
  EXPECT_TRUE(evaluator_hooked<array_product_t>);
  EXPECT_TRUE(evaluator_hooked<scalar_times_matrix_t>);
  EXPECT_TRUE(evaluator_hooked<matrix_times_scalar_t>);
}

TEST(EigenHookDispatch, non_ad_products_use_eigen) {
  using Eigen::internal::assign_op;
  using double_product_t = decltype(std::declval<const mat<double>&>()
                                    * std::declval<const mat<double>&>());
  EXPECT_FALSE(evaluator_hooked<double_product_t>);
  EXPECT_FALSE((assignment_hooked<assign_op, mat<double>, double_product_t>));
  // a double product assigned to var has no autodiff to record
  EXPECT_FALSE((assignment_hooked<assign_op, mat<var>, double_product_t>));
  EXPECT_FALSE(
      evaluator_hooked<decltype(std::declval<const mat<double>&>().array()
                                * std::declval<const mat<double>&>().array())>);
}

TEST(EigenHookDispatch, lazy_and_complex_products_use_eigen) {
  using lazy_t = decltype(std::declval<const mat<var>&>().lazyProduct(
      std::declval<const mat<var>&>()));
  using complex_t = decltype(std::declval<const mat<std::complex<var>>&>()
                             * std::declval<const mat<std::complex<var>>&>());
  EXPECT_FALSE(evaluator_hooked<lazy_t>);
  EXPECT_FALSE(evaluator_hooked<complex_t>);
}

template <typename T>
constexpr bool is_cwise_product
    = stan::math::internal::is_eigen_cwise_product<T>::value;

TEST(EigenHookDispatch, is_eigen_cwise_product) {
  // Eigen's operators return `const` expressions, which the trait decays
  const mat<var> a(2, 2);
  const mat<double> b(2, 2);
  const var s = 2;
  EXPECT_TRUE(is_cwise_product<decltype(a.cwiseProduct(b))>);
  EXPECT_TRUE(is_cwise_product<decltype(a.array() * b.array())>);
  EXPECT_TRUE(is_cwise_product<decltype(s * a)>);
  EXPECT_TRUE(is_cwise_product<decltype(b * s)>);
  EXPECT_FALSE(is_cwise_product<decltype(a * b)>);
  EXPECT_FALSE(is_cwise_product<decltype(a + b)>);
  EXPECT_FALSE(is_cwise_product<mat<var>>);
}

TEST(EigenHookDispatch, inner_product_assignment) {
  using Eigen::internal::add_assign_op;
  using Eigen::internal::assign_op;
  using Eigen::internal::sub_assign_op;
  using row_vec = Eigen::Matrix<var, 1, Eigen::Dynamic>;
  using col_vec = Eigen::Matrix<var, Eigen::Dynamic, 1>;
  using inner_t = decltype(std::declval<const row_vec&>()
                           * std::declval<const col_vec&>());
  using dst_t = Eigen::Matrix<var, 1, 1>;
  EXPECT_TRUE((assignment_hooked<assign_op, dst_t, inner_t>));
  EXPECT_TRUE((assignment_hooked<add_assign_op, dst_t, inner_t>));
  EXPECT_TRUE((assignment_hooked<sub_assign_op, dst_t, inner_t>));
}

TEST(EigenHookDispatch, triangular_product_assignment_uses_eigen) {
  using Eigen::internal::add_assign_op;
  using Eigen::internal::assign_op;
  using Eigen::internal::sub_assign_op;
  using triangular_t
      = decltype(std::declval<const mat<var>&>().triangularView<Eigen::Lower>()
                 * std::declval<const mat<var>&>());
  EXPECT_FALSE((assignment_hooked<assign_op, mat<var>, triangular_t>));
  EXPECT_FALSE((assignment_hooked<add_assign_op, mat<var>, triangular_t>));
  EXPECT_FALSE((assignment_hooked<sub_assign_op, mat<var>, triangular_t>));
}
