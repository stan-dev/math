#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/mix/util.hpp>
#include <test/unit/math/rev/fun/util.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(AgradMix, array_value_of) {
  using stan::math::fvar;
  using stan::math::value_of;
  using stan::math::var;
  using std::vector;

  vector<fvar<fvar<var>>> a(10);
  for (size_t i = 0; i < 10; ++i)
    a[i] = fvar<fvar<var>>(i);
  vector<fvar<fvar<var>>> b(5);
  for (size_t i = 0; i < 5; ++i)
    b[i] = fvar<fvar<var>>(10 + i);

  vector<fvar<var>> d_a = value_of(a);
  vector<fvar<var>> d_b = value_of(b);

  for (int i = 0; i < 5; ++i)
    EXPECT_FLOAT_EQ(b[i].val_.val_.val(), d_b[i].val_.val());

  for (int i = 0; i < 10; ++i)
    EXPECT_FLOAT_EQ(a[i].val_.val_.val(), d_a[i].val_.val());
}

TEST(AgradMix, matrix_value_of) {
  using stan::math::fvar;
  using stan::math::value_of;
  using stan::math::var;
  using std::vector;

  vector<double> a_vals;

  for (size_t i = 0; i < 10; ++i)
    a_vals.push_back(i + 1);

  vector<double> b_vals;

  for (size_t i = 10; i < 15; ++i)
    b_vals.push_back(i + 1);

  Eigen::Matrix<double, 2, 5> a;
  ::fill(a_vals, a);
  Eigen::Matrix<double, 5, 1> b;
  ::fill(b_vals, b);

  Eigen::Matrix<var, 2, 5> v_a;
  ::fill(a_vals, v_a);
  Eigen::Matrix<var, 5, 1> v_b;
  ::fill(b_vals, v_b);

  Eigen::Matrix<fvar<var>, 2, 5> fv_a;
  ::fill(a_vals, fv_a);
  Eigen::Matrix<fvar<var>, 5, 1> fv_b;
  ::fill(b_vals, fv_b);

  Eigen::Matrix<fvar<fvar<var>>, 2, 5> ffv_a;
  ::fill(a_vals, ffv_a);
  Eigen::Matrix<fvar<fvar<var>>, 5, 1> ffv_b;
  ::fill(b_vals, ffv_b);

  Eigen::MatrixXd d_a = value_of(a);
  Eigen::VectorXd d_b = value_of(b);
  Eigen::MatrixXd d_v_a = value_of(v_a);
  Eigen::MatrixXd d_v_b = value_of(v_b);
  Eigen::Matrix<var, -1, -1> d_fv_a = value_of(fv_a);
  Eigen::Matrix<var, -1, -1> d_fv_b = value_of(fv_b);
  Eigen::Matrix<fvar<var>, -1, -1> d_ffv_a = value_of(ffv_a);
  Eigen::Matrix<fvar<var>, -1, -1> d_ffv_b = value_of(ffv_b);

  for (Eigen::Index i = 0; i < 5; ++i) {
    EXPECT_FLOAT_EQ(b(i), d_b(i));
    EXPECT_FLOAT_EQ(b(i), d_v_b(i));
    EXPECT_FLOAT_EQ(b(i), d_fv_b(i).val());
    EXPECT_FLOAT_EQ(b(i), d_ffv_b(i).val_.val());
  }

  for (Eigen::Index i = 0; i < 2; ++i)
    for (Eigen::Index j = 0; j < 5; ++j) {
      EXPECT_FLOAT_EQ(a(i, j), d_a(i, j));
      EXPECT_FLOAT_EQ(a(i, j), d_v_a(i, j));
      EXPECT_FLOAT_EQ(a(i, j), d_fv_a(i, j).val());
      EXPECT_FLOAT_EQ(a(i, j), d_ffv_a(i, j).val_.val());
    }
}

TEST(AgradMix, tuple_value_of) {
  using stan::math::fvar;
  using stan::math::value_of;
  using stan::math::var;
  using std::vector;

  vector<double> a_vals;

  for (size_t i = 0; i < 10; ++i)
    a_vals.push_back(i + 1);

  vector<double> b_vals;

  for (size_t i = 10; i < 15; ++i)
    b_vals.push_back(i + 1);

  Eigen::Matrix<double, 2, 5> a;
  ::fill(a_vals, a);
  Eigen::Matrix<var, 2, 5> v_a;
  ::fill(a_vals, v_a);
  Eigen::Matrix<fvar<var>, 2, 5> fv_a;
  ::fill(a_vals, fv_a);
  Eigen::Matrix<fvar<fvar<var>>, 2, 5> ffv_a;
  ::fill(a_vals, ffv_a);

  Eigen::Matrix<double, 5, 1> b;
  ::fill(b_vals, b);
  Eigen::Matrix<var, 5, 1> v_b;
  ::fill(b_vals, v_b);
  Eigen::Matrix<fvar<var>, 5, 1> fv_b;
  ::fill(b_vals, fv_b);
  Eigen::Matrix<fvar<fvar<var>>, 5, 1> ffv_b;
  ::fill(b_vals, ffv_b);

  std::vector<fvar<fvar<var>>> ffv_a_std_vec(10);
  std::vector<double> a_std_vec(10);
  for (size_t i = 0; i < 10; ++i) {
    a_std_vec[i] = i;
    ffv_a_std_vec[i] = fvar<fvar<var>>(i);
  }
  std::vector<fvar<fvar<var>>> ffv_b_std_vec(5);
  std::vector<double> b_std_vec(5);
  for (size_t i = 0; i < 5; ++i) {
    b_std_vec[i] = 10 + i;
    ffv_b_std_vec[i] = fvar<fvar<var>>(10 + i);
  }
  auto b_tuple_dbl = std::make_tuple(b, b, b, b_std_vec);
  auto a_b_tuple_dbl = std::make_tuple(a, a, a, a_std_vec, b_tuple_dbl);
  std::vector a_b_tuple_vec_dbl{a_b_tuple_dbl, a_b_tuple_dbl, a_b_tuple_dbl};
  auto a_b_tuple_vec_tuple_dbl = std::make_tuple(a, a_b_tuple_vec_dbl, b);
  auto b_tuple_ad = std::make_tuple(v_b, fv_b, ffv_b, ffv_b_std_vec);
  auto a_b_tuple_ad
      = std::make_tuple(v_a, fv_a, ffv_a, ffv_a_std_vec, b_tuple_ad);
  std::vector a_b_tuple_vec_ad{a_b_tuple_ad, a_b_tuple_ad, a_b_tuple_ad};
  // tuple(vector, array[tuple(vec, vec, vec, array[], tuple(mat, mat, mat,
  // array[]))])
  auto a_b_tuple_vec_tuple_ad = std::make_tuple(v_a, a_b_tuple_vec_ad, ffv_b);
  stan::math::test::recursive_for_each(
      [](auto&& x_ad, auto&& x_dbl) {
        EXPECT_FLOAT_EQ(stan::math::test::get_val(x_ad), x_dbl);
      },
      value_of(a_b_tuple_vec_tuple_ad), a_b_tuple_vec_tuple_dbl);
  stan::math::test::recursive_for_each(
      [](auto&& x_value_of, auto&& x_ad) {
        using value_of_t = std::decay_t<decltype(x_value_of)>;
        using ad_t = std::decay_t<decltype(x_ad)>;
        static_assert(std::is_same_v<value_of_t, stan::partials_type_t<ad_t>>,
                      "value_of() type and partials type should be the same!!");
      },
      value_of(a_b_tuple_vec_tuple_ad), a_b_tuple_vec_tuple_ad);
}

TEST(AgradMix, value_of_expr) {
  using stan::math::fvar;
  using stan::math::value_of;
  using stan::math::var;
  Eigen::Matrix<double, -1, -1> x_d = Eigen::MatrixXd::Random(3, 3);
  Eigen::Matrix<var, -1, -1> x_v = x_d;
  Eigen::Matrix<fvar<double>, -1, -1> x_fd = x_d;
  Eigen::Matrix<fvar<var>, -1, -1> x_fv = x_d;

  using stan::math::as_array_or_scalar;
  using stan::math::to_ref;
  using stan::math::value_of;
  auto y_d = value_of(as_array_or_scalar(to_ref(x_d * x_d)));
  auto y_v = value_of(as_array_or_scalar(to_ref(x_v * x_v)));
  auto y_fd = value_of(as_array_or_scalar(to_ref(x_fd * x_fd)));
  auto y_fv = value_of(as_array_or_scalar(to_ref(x_fv * x_fv)));
  stan::math::recover_memory();
}
