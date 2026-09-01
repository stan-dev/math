#include <gtest/gtest.h>
#include <stan/math/rev/core.hpp>
#include <stan/math.hpp>
#include <test/unit/math/rev/util.hpp>
#include <tuple>
#include <type_traits>
#include <vector>

using stan::math::var;
using stan::math::vari;

TEST_F(AgradRev, Rev_deep_copy_vars_int_arg) {
  int arg = 5;

  decltype(stan::math::deep_copy_vars(arg)) out
      = stan::math::deep_copy_vars(arg);

  EXPECT_EQ(&out, &arg);
}

TEST_F(AgradRev, Rev_deep_copy_vars_double_arg) {
  double arg = 5.0;

  decltype(stan::math::deep_copy_vars(arg)) out
      = stan::math::deep_copy_vars(arg);

  EXPECT_EQ(&out, &arg);
}

TEST_F(AgradRev, Rev_deep_copy_vars_std_vector_int_arg) {
  std::vector<int> arg(5, 10);

  decltype(stan::math::deep_copy_vars(arg)) out
      = stan::math::deep_copy_vars(arg);

  EXPECT_EQ(&out, &arg);
}

TEST_F(AgradRev, Rev_deep_copy_vars_std_vector_double_arg) {
  std::vector<double> arg(5, 10.0);

  decltype(stan::math::deep_copy_vars(arg)) out
      = stan::math::deep_copy_vars(arg);

  EXPECT_EQ(&out, &arg);
}

TEST_F(AgradRev, Rev_deep_copy_vars_eigen_vector_arg) {
  Eigen::VectorXd arg = Eigen::VectorXd::Ones(5);

  decltype(stan::math::deep_copy_vars(arg)) out
      = stan::math::deep_copy_vars(arg);

  EXPECT_EQ(&out, &arg);
}

TEST_F(AgradRev, Rev_deep_copy_vars_eigen_row_vector_arg) {
  Eigen::RowVectorXd arg = Eigen::RowVectorXd::Ones(5);

  decltype(stan::math::deep_copy_vars(arg)) out
      = stan::math::deep_copy_vars(arg);

  EXPECT_EQ(&out, &arg);
}

TEST_F(AgradRev, Rev_deep_copy_vars_eigen_matrix_arg) {
  Eigen::MatrixXd arg = Eigen::MatrixXd::Ones(5, 5);

  decltype(stan::math::deep_copy_vars(arg)) out
      = stan::math::deep_copy_vars(arg);

  EXPECT_EQ(&out, &arg);
}

TEST_F(AgradRev, Rev_deep_copy_vars_std_vector_std_vector_double_arg) {
  std::vector<std::vector<double>> arg(5, std::vector<double>(5, 10.0));

  decltype(stan::math::deep_copy_vars(arg)) out
      = stan::math::deep_copy_vars(arg);

  EXPECT_EQ(&out, &arg);
}

TEST_F(AgradRev, Rev_deep_copy_vars_std_vector_eigen_vector_arg) {
  std::vector<Eigen::VectorXd> arg(2, Eigen::VectorXd::Ones(5));

  decltype(stan::math::deep_copy_vars(arg)) out
      = stan::math::deep_copy_vars(arg);

  EXPECT_EQ(&out, &arg);
}

TEST_F(AgradRev, Rev_deep_copy_vars_std_vector_eigen_row_vector_arg) {
  std::vector<Eigen::RowVectorXd> arg(2, Eigen::VectorXd::Ones(5));

  decltype(stan::math::deep_copy_vars(arg)) out
      = stan::math::deep_copy_vars(arg);

  EXPECT_EQ(&out, &arg);
}

TEST_F(AgradRev, Rev_deep_copy_vars_std_vector_eigen_matrix_arg) {
  std::vector<Eigen::MatrixXd> arg(2, Eigen::MatrixXd::Ones(5, 3));

  decltype(stan::math::deep_copy_vars(arg)) out
      = stan::math::deep_copy_vars(arg);

  EXPECT_EQ(&out, &arg);
}

TEST_F(AgradRev, Rev_deep_copy_vars_var_arg) {
  var arg(5.0);
  arg.vi_->adj_ = 2.0;

  decltype(stan::math::deep_copy_vars(arg)) out
      = stan::math::deep_copy_vars(arg);

  out.grad();

  EXPECT_EQ(out.vi_->adj_, 1.0);
  EXPECT_EQ(arg.vi_->adj_, 2.0);
  EXPECT_EQ(out, arg);
  EXPECT_NE(out.vi_, arg.vi_);
}

TEST_F(AgradRev, Rev_deep_copy_vars_std_vector_var_arg) {
  std::vector<var> arg(5);
  for (size_t i = 0; i < arg.size(); ++i)
    arg[i] = i + 1.0;

  decltype(stan::math::deep_copy_vars(arg)) out
      = stan::math::deep_copy_vars(arg);

  for (int i = 0; i < arg.size(); ++i) {
    stan::math::set_zero_all_adjoints();
    arg[i].vi_->adj_ = 2.0;

    out[i].grad();

    EXPECT_EQ(out[i].vi_->adj_, 1.0);
    EXPECT_EQ(arg[i].vi_->adj_, 2.0);
    EXPECT_EQ(out[i], arg[i]);
    EXPECT_NE(out[i].vi_, arg[i].vi_);
  }
}

TEST_F(AgradRev, Rev_deep_copy_vars_eigen_vector_var_arg) {
  Eigen::Matrix<var, Eigen::Dynamic, 1> arg(5);
  for (size_t i = 0; i < arg.size(); ++i) {
    arg(i) = i + 1.0;
    arg(i).vi_->adj_ = i + 1.0;
  }

  decltype(stan::math::deep_copy_vars(arg)) out
      = stan::math::deep_copy_vars(arg);

  for (int i = 0; i < arg.size(); ++i) {
    stan::math::set_zero_all_adjoints();
    arg(i).vi_->adj_ = 2.0;

    out(i).grad();

    EXPECT_EQ(out(i).vi_->adj_, 1.0);
    EXPECT_EQ(arg(i).vi_->adj_, 2.0);
    EXPECT_EQ(out(i), arg(i));
    EXPECT_NE(out(i).vi_, arg(i).vi_);
  }
}

TEST_F(AgradRev, Rev_deep_copy_vars_eigen_row_vector_var_arg) {
  Eigen::Matrix<var, 1, Eigen::Dynamic> arg(5);
  for (size_t i = 0; i < arg.size(); ++i) {
    arg(i) = i + 1.0;
    arg(i).vi_->adj_ = i + 1.0;
  }

  decltype(stan::math::deep_copy_vars(arg)) out
      = stan::math::deep_copy_vars(arg);

  for (int i = 0; i < arg.size(); ++i) {
    stan::math::set_zero_all_adjoints();
    arg(i).vi_->adj_ = 2.0;

    out(i).grad();

    EXPECT_EQ(out(i).vi_->adj_, 1.0);
    EXPECT_EQ(arg(i).vi_->adj_, 2.0);
    EXPECT_EQ(out(i), arg(i));
    EXPECT_NE(out(i).vi_, arg(i).vi_);
  }
}

TEST_F(AgradRev, Rev_deep_copy_vars_eigen_matrix_var_arg) {
  Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> arg(5, 5);
  for (size_t i = 0; i < arg.size(); ++i) {
    arg(i) = i + 1.0;
    arg(i).vi_->adj_ = i + 1.0;
  }

  decltype(stan::math::deep_copy_vars(arg)) out
      = stan::math::deep_copy_vars(arg);

  for (int i = 0; i < arg.size(); ++i) {
    stan::math::set_zero_all_adjoints();
    arg(i).vi_->adj_ = 2.0;

    out(i).grad();

    EXPECT_EQ(out(i).vi_->adj_, 1.0);
    EXPECT_EQ(arg(i).vi_->adj_, 2.0);
    EXPECT_EQ(out(i), arg(i));
    EXPECT_NE(out(i).vi_, arg(i).vi_);
  }
}

TEST_F(AgradRev, Rev_deep_copy_vars_std_vector_std_vector_var_arg) {
  std::vector<var> arg_(5);
  std::vector<std::vector<var>> arg(5, arg_);
  for (size_t i = 0; i < arg.size(); ++i)
    for (size_t j = 0; j < arg[i].size(); ++j)
      arg[i][j] = i * arg[i].size() + j + 5.0;

  decltype(stan::math::deep_copy_vars(arg)) out
      = stan::math::deep_copy_vars(arg);

  for (int i = 0; i < arg.size(); ++i)
    for (int j = 0; j < arg[i].size(); ++j) {
      stan::math::set_zero_all_adjoints();
      arg[i][j].vi_->adj_ = 2.0;

      out[i][j].grad();

      EXPECT_EQ(out[i][j].vi_->adj_, 1.0);
      EXPECT_EQ(arg[i][j].vi_->adj_, 2.0);
      EXPECT_EQ(out[i][j], arg[i][j]);
      EXPECT_NE(out[i][j].vi_, arg[i][j].vi_);
    }
}

TEST_F(AgradRev, Rev_deep_copy_vars_std_vector_eigen_vector_var_arg) {
  Eigen::Matrix<var, Eigen::Dynamic, 1> arg_(5);
  std::vector<Eigen::Matrix<var, Eigen::Dynamic, 1>> arg(2, arg_);
  for (size_t i = 0; i < arg.size(); ++i)
    for (size_t j = 0; j < arg[i].size(); ++j)
      arg[i](j) = i * arg[i].size() + j + 5.0;

  decltype(stan::math::deep_copy_vars(arg)) out
      = stan::math::deep_copy_vars(arg);

  for (int i = 0; i < arg.size(); ++i)
    for (int j = 0; j < arg[i].size(); ++j) {
      stan::math::set_zero_all_adjoints();
      arg[i](j).vi_->adj_ = 2.0;

      out[i](j).grad();

      EXPECT_EQ(out[i](j).vi_->adj_, 1.0);
      EXPECT_EQ(arg[i](j).vi_->adj_, 2.0);
      EXPECT_EQ(out[i](j), arg[i](j));
      EXPECT_NE(out[i](j).vi_, arg[i](j).vi_);
    }
}

TEST_F(AgradRev, Rev_deep_copy_vars_std_vector_eigen_row_vector_var_arg) {
  Eigen::Matrix<var, 1, Eigen::Dynamic> arg_(5);
  std::vector<Eigen::Matrix<var, 1, Eigen::Dynamic>> arg(2, arg_);
  for (size_t i = 0; i < arg.size(); ++i)
    for (size_t j = 0; j < arg[i].size(); ++j)
      arg[i](j) = i * arg[i].size() + j + 5.0;

  decltype(stan::math::deep_copy_vars(arg)) out
      = stan::math::deep_copy_vars(arg);

  for (int i = 0; i < arg.size(); ++i)
    for (int j = 0; j < arg[i].size(); ++j) {
      stan::math::set_zero_all_adjoints();
      arg[i](j).vi_->adj_ = 2.0;

      out[i](j).grad();

      EXPECT_EQ(out[i](j).vi_->adj_, 1.0);
      EXPECT_EQ(arg[i](j).vi_->adj_, 2.0);
      EXPECT_EQ(out[i](j), arg[i](j));
      EXPECT_NE(out[i](j).vi_, arg[i](j).vi_);
    }
}

TEST_F(AgradRev, Rev_deep_copy_vars_std_vector_eigen_matrix_var_arg) {
  Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> arg_(5, 3);
  std::vector<Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>> arg(2, arg_);
  for (size_t i = 0; i < arg.size(); ++i)
    for (size_t j = 0; j < arg[i].size(); ++j)
      arg[i](j) = i * arg[i].size() + j + 5.0;

  decltype(stan::math::deep_copy_vars(arg)) out
      = stan::math::deep_copy_vars(arg);

  for (int i = 0; i < arg.size(); ++i)
    for (int j = 0; j < arg[i].size(); ++j) {
      stan::math::set_zero_all_adjoints();
      arg[i](j).vi_->adj_ = 2.0;

      out[i](j).grad();

      EXPECT_EQ(out[i](j).vi_->adj_, 1.0);
      EXPECT_EQ(arg[i](j).vi_->adj_, 2.0);
      EXPECT_EQ(out[i](j), arg[i](j));
      EXPECT_NE(out[i](j).vi_, arg[i](j).vi_);
    }
}

TEST_F(AgradRev, Rev_deep_copy_vars_tuple_data_arg) {
  const std::tuple<> empty;
  const auto arg = std::make_tuple(5, Eigen::VectorXd::Ones(2).eval());

  auto empty_out = stan::math::deep_copy_vars(empty);
  auto out = stan::math::deep_copy_vars(arg);

  static_assert(std::is_same_v<decltype(empty_out), std::tuple<>>);
  static_assert(std::is_same_v<decltype(out),
                               std::tuple<const int&, const Eigen::VectorXd&>>);
  EXPECT_EQ(&std::get<0>(out), &std::get<0>(arg));
  EXPECT_EQ(&std::get<1>(out), &std::get<1>(arg));
}

TEST_F(AgradRev, Rev_deep_copy_vars_nested_tuple_var_arg) {
  Eigen::Matrix<var, Eigen::Dynamic, 1> vars(2);
  vars << 2.0, 3.0;
  auto arg = std::make_tuple(1.0, var(4.0), vars, std::make_tuple(var(5.0), 6));

  auto out = stan::math::deep_copy_vars(arg);

  static_assert(std::is_reference_v<std::tuple_element_t<0, decltype(out)>>);
  EXPECT_EQ(&std::get<0>(out), &std::get<0>(arg));
  EXPECT_FLOAT_EQ(std::get<1>(out).val(), std::get<1>(arg).val());
  EXPECT_NE(std::get<1>(out).vi_, std::get<1>(arg).vi_);
  for (int i = 0; i < vars.size(); ++i) {
    EXPECT_FLOAT_EQ(std::get<2>(out)(i).val(), std::get<2>(arg)(i).val());
    EXPECT_NE(std::get<2>(out)(i).vi_, std::get<2>(arg)(i).vi_);
  }
  EXPECT_FLOAT_EQ(std::get<0>(std::get<3>(out)).val(),
                  std::get<0>(std::get<3>(arg)).val());
  EXPECT_NE(std::get<0>(std::get<3>(out)).vi_,
            std::get<0>(std::get<3>(arg)).vi_);
  EXPECT_EQ(&std::get<1>(std::get<3>(out)), &std::get<1>(std::get<3>(arg)));

  auto rvalue_out = stan::math::deep_copy_vars(
      std::make_tuple(7.0, var(8.0), std::make_tuple(9)));
  static_assert(std::is_same_v<decltype(rvalue_out),
                               std::tuple<double, var, std::tuple<int>>>);
  EXPECT_FLOAT_EQ(8.0, std::get<1>(rvalue_out).val());
}
