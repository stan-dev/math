#include <gtest/gtest.h>
#include <stan/math/rev/core.hpp>
#include <stan/math.hpp>
#include <vector>

using stan::math::var;
using stan::math::vari;

TEST(AgradRev_deep_copy_vars, int_arg) {
  int arg = 5;

  decltype(stan::math::deep_copy_vars(arg)) out
      = stan::math::deep_copy_vars(arg);

  EXPECT_EQ(&out, &arg);
}

TEST(AgradRev_deep_copy_vars, double_arg) {
  double arg = 5.0;

  decltype(stan::math::deep_copy_vars(arg)) out
      = stan::math::deep_copy_vars(arg);

  EXPECT_EQ(&out, &arg);
}

TEST(AgradRev_deep_copy_vars, std_vector_int_arg) {
  std::vector<int> arg(5, 10);

  decltype(stan::math::deep_copy_vars(arg)) out
      = stan::math::deep_copy_vars(arg);

  EXPECT_EQ(&out, &arg);
}

TEST(AgradRev_deep_copy_vars, std_vector_double_arg) {
  std::vector<double> arg(5, 10.0);

  decltype(stan::math::deep_copy_vars(arg)) out
      = stan::math::deep_copy_vars(arg);

  EXPECT_EQ(&out, &arg);
}

TEST(AgradRev_deep_copy_vars, eigen_vector_arg) {
  Eigen::VectorXd arg = Eigen::VectorXd::Ones(5);

  decltype(stan::math::deep_copy_vars(arg)) out
      = stan::math::deep_copy_vars(arg);

  EXPECT_EQ(&out, &arg);
}

TEST(AgradRev_deep_copy_vars, eigen_row_vector_arg) {
  Eigen::RowVectorXd arg = Eigen::RowVectorXd::Ones(5);

  decltype(stan::math::deep_copy_vars(arg)) out
      = stan::math::deep_copy_vars(arg);

  EXPECT_EQ(&out, &arg);
}

TEST(AgradRev_deep_copy_vars, eigen_matrix_arg) {
  Eigen::MatrixXd arg = Eigen::MatrixXd::Ones(5, 5);

  decltype(stan::math::deep_copy_vars(arg)) out
      = stan::math::deep_copy_vars(arg);

  EXPECT_EQ(&out, &arg);
}

TEST(AgradRev_deep_copy_vars, std_vector_std_vector_double_arg) {
  std::vector<std::vector<double>> arg(5, std::vector<double>(5, 10.0));

  decltype(stan::math::deep_copy_vars(arg)) out
      = stan::math::deep_copy_vars(arg);

  EXPECT_EQ(&out, &arg);
}

TEST(AgradRev_deep_copy_vars, std_vector_eigen_vector_arg) {
  std::vector<Eigen::VectorXd> arg(2, Eigen::VectorXd::Ones(5));

  decltype(stan::math::deep_copy_vars(arg)) out
      = stan::math::deep_copy_vars(arg);

  EXPECT_EQ(&out, &arg);
}

TEST(AgradRev_deep_copy_vars, std_vector_eigen_row_vector_arg) {
  std::vector<Eigen::RowVectorXd> arg(2, Eigen::VectorXd::Ones(5));

  decltype(stan::math::deep_copy_vars(arg)) out
      = stan::math::deep_copy_vars(arg);

  EXPECT_EQ(&out, &arg);
}

TEST(AgradRev_deep_copy_vars, std_vector_eigen_matrix_arg) {
  std::vector<Eigen::MatrixXd> arg(2, Eigen::MatrixXd::Ones(5, 3));

  decltype(stan::math::deep_copy_vars(arg)) out
      = stan::math::deep_copy_vars(arg);

  EXPECT_EQ(&out, &arg);
}

TEST(AgradRev_deep_copy_vars, var_arg) {
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

TEST(AgradRev_deep_copy_vars, std_vector_var_arg) {
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

TEST(AgradRev_deep_copy_vars, eigen_vector_var_arg) {
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

TEST(AgradRev_deep_copy_vars, eigen_row_vector_var_arg) {
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

TEST(AgradRev_deep_copy_vars, eigen_matrix_var_arg) {
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

TEST(AgradRev_deep_copy_vars, std_vector_std_vector_var_arg) {
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

TEST(AgradRev_deep_copy_vars, std_vector_eigen_vector_var_arg) {
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

TEST(AgradRev_deep_copy_vars, std_vector_eigen_row_vector_var_arg) {
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

TEST(AgradRev_deep_copy_vars, std_vector_eigen_matrix_var_arg) {
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
