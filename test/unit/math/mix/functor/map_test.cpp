#include <stan/math/mix.hpp>
#include <test/unit/math/prim/util.hpp>
#include <test/unit/pretty_print_types.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <numeric>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace {

struct mixed_output_functor {
  int operator()(double x) const { return static_cast<int>(x * 2.0); }

  double operator()(const stan::math::var& x) const { return x.val(); }

  Eigen::Index operator()(const Eigen::VectorXd& x) const { return x.size(); }

  double operator()(const Eigen::MatrixXd& x) const { return x.sum(); }

  std::size_t operator()(const std::vector<int>& x) const { return x.size(); }

  std::vector<int> operator()(const std::vector<std::vector<int>>& x) const {
    std::vector<int> out;
    out.reserve(x.size());
    for (const auto& inner : x) {
      out.push_back(std::accumulate(inner.begin(), inner.end(), 0));
    }
    return out;
  }
};

TEST(MathMixFunctor, map_scalar_support_and_return_type_changes) {
  int i = 4;
  auto&& i_ref = stan::math::map(
      [](auto& x) -> auto& { return x; }, i);
  EXPECT_SAME_TYPE(int&, decltype(i_ref));
  i_ref = 11;
  EXPECT_EQ(i, 11);

  const int ci = 5;
  auto&& ci_ref = stan::math::map(
      [](const auto& x) -> const auto& { return x; }, ci);
  EXPECT_SAME_TYPE(const int&, decltype(ci_ref));

  auto as_bool = stan::math::map([](double x) { return x > 0.0; }, -1.5);
  EXPECT_TYPE(bool, as_bool);
  EXPECT_FALSE(as_bool);

  stan::math::var ad = 2.0;
  auto ad_out = stan::math::map([](const auto& x) { return x * 3.0; }, ad);
  EXPECT_TYPE(stan::math::var, ad_out);
  EXPECT_FLOAT_EQ(ad_out.val(), 6.0);

  stan::math::fvar<stan::math::var> fvad(2.0, 1.0);
  auto fvad_out
      = stan::math::map([](const auto& x) { return x * 4.0; }, fvad);
  EXPECT_TYPE(stan::math::fvar<stan::math::var>, fvad_out);
  EXPECT_FLOAT_EQ(fvad_out.val_.val(), 8.0);
  EXPECT_FLOAT_EQ(fvad_out.d_.val(), 4.0);
}

TEST(MathMixFunctor, map_non_tuple_accepts_rvalue_scalar) {
  auto out = stan::math::map([](double&& x) { return x + 1.0; }, 2.5);
  EXPECT_TYPE(double, out);
  EXPECT_FLOAT_EQ(out, 3.5);
}

TEST(MathMixFunctor, map_non_tuple_eigen_vector) {
  Eigen::VectorXd eigen_vec(3);
  eigen_vec << 1.0, 2.0, 3.0;
  auto vec_out = stan::math::map(
                      [](const Eigen::VectorXd& x) {
                        return (2.0 * x.array()).matrix();
                      },
                      eigen_vec)
                     .eval();
  EXPECT_SAME_TYPE(Eigen::VectorXd, std::decay_t<decltype(vec_out)>);
  Eigen::VectorXd expected_vec(3);
  expected_vec << 2.0, 4.0, 6.0;
  EXPECT_MATRIX_EQ(vec_out, expected_vec);
}

TEST(MathMixFunctor, map_non_tuple_eigen_matrix) {
  Eigen::MatrixXd eigen_mat(2, 2);
  eigen_mat << 1.0, 2.0, 3.0, 4.0;
  auto mat_out
      = stan::math::map([](const Eigen::MatrixXd& x) { return (x * 2.0); },
                        eigen_mat)
            .eval();
  EXPECT_SAME_TYPE(Eigen::MatrixXd, std::decay_t<decltype(mat_out)>);
  EXPECT_MATRIX_EQ(mat_out, (eigen_mat * 2.0).eval());
}

TEST(MathMixFunctor, map_non_tuple_accepts_eigen_vector_expression) {
  Eigen::VectorXd lhs(3);
  lhs << 1.0, 2.0, 3.0;
  Eigen::VectorXd rhs(3);
  rhs << 4.0, 5.0, 6.0;

  auto vec_out
      = stan::math::map([](const auto& x) { return x.eval(); }, lhs + rhs)
            .eval();
  EXPECT_SAME_TYPE(Eigen::VectorXd, std::decay_t<decltype(vec_out)>);
  EXPECT_MATRIX_EQ(vec_out, (lhs + rhs).eval());
}

TEST(MathMixFunctor, map_non_tuple_accepts_eigen_matrix_expression) {
  Eigen::MatrixXd lhs(2, 2);
  lhs << 1.0, 2.0, 3.0, 4.0;
  Eigen::MatrixXd rhs(2, 2);
  rhs << -1.0, 0.5, 1.5, -2.0;

  auto mat_out
      = stan::math::map([](const auto& x) { return x.eval(); }, lhs * rhs)
            .eval();
  EXPECT_SAME_TYPE(Eigen::MatrixXd, std::decay_t<decltype(mat_out)>);
  EXPECT_MATRIX_EQ(mat_out, (lhs * rhs).eval());
}

TEST(MathMixFunctor, map_non_tuple_std_vector) {
  std::vector<int> std_vec{1, 2, 3};
  auto std_vec_out = stan::math::map(
      [](const std::vector<int>& x) {
        std::vector<double> out;
        out.reserve(x.size());
        for (int v : x) {
          out.push_back(v + 0.5);
        }
        return out;
      },
      std_vec);
  EXPECT_SAME_TYPE(std::vector<double>, std::decay_t<decltype(std_vec_out)>);
  ASSERT_EQ(std_vec_out.size(), std_vec.size());
  EXPECT_FLOAT_EQ(std_vec_out[0], 1.5);
  EXPECT_FLOAT_EQ(std_vec_out[1], 2.5);
  EXPECT_FLOAT_EQ(std_vec_out[2], 3.5);
}

TEST(MathMixFunctor, map_non_tuple_accepts_rvalue_std_vector) {
  auto std_vec_out = stan::math::map(
      [](std::vector<int>&& x) {
        x.push_back(4);
        return std::move(x);
      },
      std::vector<int>{1, 2, 3});
  EXPECT_SAME_TYPE(std::vector<int>, std::decay_t<decltype(std_vec_out)>);
  std::vector<int> expected_std_vec{1, 2, 3, 4};
  EXPECT_EQ(std_vec_out, expected_std_vec);
}

TEST(MathMixFunctor, map_non_tuple_nested_std_vector) {
  std::vector<std::vector<int>> nested_vec{{1, 2, 3}, {4}, {}};
  auto nested_out = stan::math::map(
      [](const std::vector<std::vector<int>>& x) {
        std::vector<int> out;
        out.reserve(x.size());
        for (const auto& inner : x) {
          out.push_back(std::accumulate(inner.begin(), inner.end(), 0));
        }
        return out;
      },
      nested_vec);
  EXPECT_SAME_TYPE(std::vector<int>, std::decay_t<decltype(nested_out)>);
  std::vector<int> expected_nested{6, 4, 0};
  EXPECT_EQ(nested_out, expected_nested);
}

TEST(MathMixFunctor, map_non_tuple_accepts_rvalue_nested_std_vector) {
  auto nested_out = stan::math::map(
      [](std::vector<std::vector<int>>&& x) {
        x.push_back({7, 8});
        return std::move(x);
      },
      std::vector<std::vector<int>>{{1, 2}, {3, 4}});
  EXPECT_SAME_TYPE(std::vector<std::vector<int>>,
                   std::decay_t<decltype(nested_out)>);
  ASSERT_EQ(nested_out.size(), 3);
  EXPECT_EQ(nested_out.back(), (std::vector<int>{7, 8}));
}

TEST(MathMixFunctor, map_tuple_mixed_inputs_respect_output_types) {
  using stan::test::is_ref_element_v;
  using stan::test::is_same_tuple_element_v;

  double scalar = 1.5;
  stan::math::var ad_scalar = 2.5;
  Eigen::VectorXd eigen_vec(2);
  eigen_vec << 3.0, 4.0;
  Eigen::MatrixXd eigen_mat(2, 2);
  eigen_mat << 1.0, 2.0, 3.0, 4.0;
  std::vector<int> std_vec{5, 6, 7};
  std::vector<std::vector<int>> nested_vec{{1, 2}, {3, 4, 5}};

  auto out = stan::math::map(
      mixed_output_functor{},
      std::forward_as_tuple(scalar, ad_scalar, eigen_vec, eigen_mat, std_vec,
                            nested_vec));
  using out_t = decltype(out);
  EXPECT_EQ(std::tuple_size_v<out_t>, 6);
  EXPECT_TRUE((is_same_tuple_element_v<0, out_t, int>));
  EXPECT_TRUE((is_same_tuple_element_v<1, out_t, double>));
  EXPECT_TRUE((is_same_tuple_element_v<2, out_t, Eigen::Index>));
  EXPECT_TRUE((is_same_tuple_element_v<3, out_t, double>));
  EXPECT_TRUE((is_same_tuple_element_v<4, out_t, std::size_t>));
  EXPECT_TRUE((is_same_tuple_element_v<5, out_t, std::vector<int>>));
  EXPECT_FALSE((is_ref_element_v<0, out_t>));
  EXPECT_FALSE((is_ref_element_v<1, out_t>));
  EXPECT_FALSE((is_ref_element_v<2, out_t>));
  EXPECT_FALSE((is_ref_element_v<3, out_t>));
  EXPECT_FALSE((is_ref_element_v<4, out_t>));
  EXPECT_FALSE((is_ref_element_v<5, out_t>));

  EXPECT_EQ(std::get<0>(out), 3);
  EXPECT_FLOAT_EQ(std::get<1>(out), 2.5);
  EXPECT_EQ(std::get<2>(out), 2);
  EXPECT_FLOAT_EQ(std::get<3>(out), 10.0);
  EXPECT_EQ(std::get<4>(out), 3);
  std::vector<int> expected_nested{3, 12};
  EXPECT_EQ(std::get<5>(out), expected_nested);
}

TEST(MathMixFunctor, map_tuple_preserves_lvalue_references_when_returned) {
  using stan::test::is_lvalue_ref_element_v;

  int scalar = 2;
  stan::math::var ad_scalar = 3.0;
  Eigen::VectorXd eigen_vec(2);
  eigen_vec << 1.0, 2.0;
  std::vector<int> std_vec{5, 6};
  std::vector<std::vector<int>> nested_vec{{1}, {2, 3}};

  auto out = stan::math::map(
      [](auto&& x) -> decltype(auto) { return std::forward<decltype(x)>(x); },
      std::forward_as_tuple(scalar, ad_scalar, eigen_vec, std_vec, nested_vec));
  using out_t = decltype(out);
  EXPECT_TRUE((is_lvalue_ref_element_v<0, out_t>));
  EXPECT_TRUE((is_lvalue_ref_element_v<1, out_t>));
  EXPECT_TRUE((is_lvalue_ref_element_v<2, out_t>));
  EXPECT_TRUE((is_lvalue_ref_element_v<3, out_t>));
  EXPECT_TRUE((is_lvalue_ref_element_v<4, out_t>));

  std::get<0>(out) = 11;
  std::get<1>(out) = 7.0;
  std::get<2>(out)(0) = -1.0;
  std::get<3>(out)[1] = 42;
  std::get<4>(out)[1].push_back(9);

  EXPECT_EQ(scalar, 11);
  EXPECT_FLOAT_EQ(ad_scalar.val(), 7.0);
  EXPECT_FLOAT_EQ(eigen_vec(0), -1.0);
  EXPECT_EQ(std_vec[1], 42);
  EXPECT_EQ(nested_vec[1].back(), 9);
}

TEST(MathMixFunctor, map_tuple_rvalues_are_stored_by_value) {
  using stan::test::is_ref_element_v;

  auto out = stan::math::map(
      [](auto&& x) -> decltype(auto) { return std::forward<decltype(x)>(x); },
      std::make_tuple(4.0, std::vector<int>{1, 2, 3},
                      Eigen::VectorXd::Constant(2, 5.0)));
  using out_t = decltype(out);
  EXPECT_FALSE((is_ref_element_v<0, out_t>));
  EXPECT_FALSE((is_ref_element_v<1, out_t>));
  EXPECT_FALSE((is_ref_element_v<2, out_t>));

  EXPECT_FLOAT_EQ(std::get<0>(out), 4.0);
  std::vector<int> expected_std_vec{1, 2, 3};
  EXPECT_EQ(std::get<1>(out), expected_std_vec);
  EXPECT_FLOAT_EQ(std::get<2>(out)(0), 5.0);
  EXPECT_FLOAT_EQ(std::get<2>(out)(1), 5.0);
}

TEST(MathMixFunctor, map_empty_tuple_returns_empty_tuple) {
  auto out = stan::math::map([](auto&& x) { return x; }, std::tuple<>{});
  EXPECT_EQ(std::tuple_size_v<decltype(out)>, 0);
}

}  // namespace
