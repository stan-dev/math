#include <gtest/gtest.h>
#include <stan/math/rev/core.hpp>
#include <stan/math.hpp>
#include <vector>

using stan::math::var;
using stan::math::vari;

Eigen::VectorXd storage(1000);

TEST(AgradRev_count_vars, int_arg) {
  int arg = 5;

  storage.setZero();
  double* ptr = stan::math::accumulate_adjoints(storage.data(), arg);

  size_t num_vars = stan::math::count_vars(arg);

  for (int i = 0; i < num_vars; ++i)
    EXPECT_FLOAT_EQ(storage(i), 1.0);
  for (int i = num_vars; i < storage.size(); ++i)
    EXPECT_FLOAT_EQ(storage(i), 0.0);

  EXPECT_EQ(ptr, storage.data() + num_vars);
}

TEST(AgradRev_count_vars, double_arg) {
  double arg = 5.0;

  storage.setZero();
  double* ptr = stan::math::accumulate_adjoints(storage.data(), arg);

  size_t num_vars = stan::math::count_vars(arg);

  for (int i = 0; i < num_vars; ++i)
    EXPECT_FLOAT_EQ(storage(i), 1.0);
  for (int i = num_vars; i < storage.size(); ++i)
    EXPECT_FLOAT_EQ(storage(i), 0.0);

  EXPECT_EQ(ptr, storage.data() + num_vars);
}

TEST(AgradRev_count_vars, std_vector_int_arg) {
  std::vector<int> arg(5, 10);

  storage.setZero();
  double* ptr = stan::math::accumulate_adjoints(storage.data(), arg);

  size_t num_vars = stan::math::count_vars(arg);

  for (int i = 0; i < num_vars; ++i)
    EXPECT_FLOAT_EQ(storage(i), 1.0);
  for (int i = num_vars; i < storage.size(); ++i)
    EXPECT_FLOAT_EQ(storage(i), 0.0);

  EXPECT_EQ(ptr, storage.data() + num_vars);
}

TEST(AgradRev_count_vars, std_vector_double_arg) {
  std::vector<double> arg(5, 10.0);

  storage.setZero();
  double* ptr = stan::math::accumulate_adjoints(storage.data(), arg);

  size_t num_vars = stan::math::count_vars(arg);

  for (int i = 0; i < num_vars; ++i)
    EXPECT_FLOAT_EQ(storage(i), 1.0);
  for (int i = num_vars; i < storage.size(); ++i)
    EXPECT_FLOAT_EQ(storage(i), 0.0);

  EXPECT_EQ(ptr, storage.data() + num_vars);
}

TEST(AgradRev_count_vars, eigen_vector_arg) {
  Eigen::VectorXd arg = Eigen::VectorXd::Ones(5);

  storage.setZero();
  double* ptr = stan::math::accumulate_adjoints(storage.data(), arg);

  size_t num_vars = stan::math::count_vars(arg);

  for (int i = 0; i < num_vars; ++i)
    EXPECT_FLOAT_EQ(storage(i), 1.0);
  for (int i = num_vars; i < storage.size(); ++i)
    EXPECT_FLOAT_EQ(storage(i), 0.0);

  EXPECT_EQ(ptr, storage.data() + num_vars);
}

TEST(AgradRev_count_vars, eigen_row_vector_arg) {
  Eigen::RowVectorXd arg = Eigen::RowVectorXd::Ones(5);

  storage.setZero();
  double* ptr = stan::math::accumulate_adjoints(storage.data(), arg);

  size_t num_vars = stan::math::count_vars(arg);

  for (int i = 0; i < num_vars; ++i)
    EXPECT_FLOAT_EQ(storage(i), 1.0);
  for (int i = num_vars; i < storage.size(); ++i)
    EXPECT_FLOAT_EQ(storage(i), 0.0);

  EXPECT_EQ(ptr, storage.data() + num_vars);
}

TEST(AgradRev_count_vars, eigen_matrix_arg) {
  Eigen::MatrixXd arg = Eigen::MatrixXd::Ones(5, 5);

  storage.setZero();
  double* ptr = stan::math::accumulate_adjoints(storage.data(), arg);

  size_t num_vars = stan::math::count_vars(arg);

  for (int i = 0; i < num_vars; ++i)
    EXPECT_FLOAT_EQ(storage(i), 1.0);
  for (int i = num_vars; i < storage.size(); ++i)
    EXPECT_FLOAT_EQ(storage(i), 0.0);

  EXPECT_EQ(ptr, storage.data() + num_vars);
}

TEST(AgradRev_count_vars, std_vector_std_vector_double_arg) {
  std::vector<std::vector<double>> arg(5, std::vector<double>(5, 10.0));

  storage.setZero();
  double* ptr = stan::math::accumulate_adjoints(storage.data(), arg);

  size_t num_vars = stan::math::count_vars(arg);

  for (int i = 0; i < num_vars; ++i)
    EXPECT_FLOAT_EQ(storage(i), 1.0);
  for (int i = num_vars; i < storage.size(); ++i)
    EXPECT_FLOAT_EQ(storage(i), 0.0);

  EXPECT_EQ(ptr, storage.data() + num_vars);
}

TEST(AgradRev_count_vars, std_vector_eigen_vector_arg) {
  std::vector<Eigen::VectorXd> arg(2, Eigen::VectorXd::Ones(5));

  storage.setZero();
  double* ptr = stan::math::accumulate_adjoints(storage.data(), arg);

  size_t num_vars = stan::math::count_vars(arg);

  for (int i = 0; i < num_vars; ++i)
    EXPECT_FLOAT_EQ(storage(i), 1.0);
  for (int i = num_vars; i < storage.size(); ++i)
    EXPECT_FLOAT_EQ(storage(i), 0.0);

  EXPECT_EQ(ptr, storage.data() + num_vars);
}

TEST(AgradRev_count_vars, std_vector_eigen_row_vector_arg) {
  std::vector<Eigen::RowVectorXd> arg(2, Eigen::VectorXd::Ones(5));

  storage.setZero();
  double* ptr = stan::math::accumulate_adjoints(storage.data(), arg);

  size_t num_vars = stan::math::count_vars(arg);

  for (int i = 0; i < num_vars; ++i)
    EXPECT_FLOAT_EQ(storage(i), 1.0);
  for (int i = num_vars; i < storage.size(); ++i)
    EXPECT_FLOAT_EQ(storage(i), 0.0);

  EXPECT_EQ(ptr, storage.data() + num_vars);
}

TEST(AgradRev_count_vars, std_vector_eigen_matrix_arg) {
  std::vector<Eigen::MatrixXd> arg(2, Eigen::MatrixXd::Ones(5, 3));

  storage.setZero();
  double* ptr = stan::math::accumulate_adjoints(storage.data(), arg);

  size_t num_vars = stan::math::count_vars(arg);

  for (int i = 0; i < num_vars; ++i)
    EXPECT_FLOAT_EQ(storage(i), 1.0);
  for (int i = num_vars; i < storage.size(); ++i)
    EXPECT_FLOAT_EQ(storage(i), 0.0);

  EXPECT_EQ(ptr, storage.data() + num_vars);
}

TEST(AgradRev_count_vars, var_arg) {
  var arg(5.0);
  arg.vi_->adj_ = 1.0;

  storage.setZero();
  double* ptr = stan::math::accumulate_adjoints(storage.data(), arg);

  size_t num_vars = stan::math::count_vars(arg);

  for (int i = 0; i < num_vars; ++i)
    EXPECT_FLOAT_EQ(storage(i), 1.0);
  for (int i = num_vars; i < storage.size(); ++i)
    EXPECT_FLOAT_EQ(storage(i), 0.0);

  EXPECT_EQ(ptr, storage.data() + num_vars);
}

TEST(AgradRev_count_vars, std_vector_var_arg) {
  std::vector<var> arg(5, var(5.0));
  for (size_t i = 0; i < arg.size(); ++i)
    arg[i].vi_->adj_ = 1.0;

  storage.setZero();
  double* ptr = stan::math::accumulate_adjoints(storage.data(), arg);

  size_t num_vars = stan::math::count_vars(arg);

  for (int i = 0; i < num_vars; ++i)
    EXPECT_FLOAT_EQ(storage(i), 1.0);
  for (int i = num_vars; i < storage.size(); ++i)
    EXPECT_FLOAT_EQ(storage(i), 0.0);

  EXPECT_EQ(ptr, storage.data() + num_vars);
}

TEST(AgradRev_count_vars, eigen_vector_var_arg) {
  Eigen::Matrix<var, Eigen::Dynamic, 1> arg(5);
  for (size_t i = 0; i < arg.size(); ++i) {
    arg(i) = 5.0;
    arg(i).vi_->adj_ = 1.0;
  }

  storage.setZero();
  double* ptr = stan::math::accumulate_adjoints(storage.data(), arg);

  size_t num_vars = stan::math::count_vars(arg);

  for (int i = 0; i < num_vars; ++i)
    EXPECT_FLOAT_EQ(storage(i), 1.0);
  for (int i = num_vars; i < storage.size(); ++i)
    EXPECT_FLOAT_EQ(storage(i), 0.0);

  EXPECT_EQ(ptr, storage.data() + num_vars);
}

TEST(AgradRev_count_vars, eigen_row_vector_var_arg) {
  Eigen::Matrix<var, 1, Eigen::Dynamic> arg(5);
  for (size_t i = 0; i < arg.size(); ++i) {
    arg(i) = 5.0;
    arg(i).vi_->adj_ = 1.0;
  }

  storage.setZero();
  double* ptr = stan::math::accumulate_adjoints(storage.data(), arg);

  size_t num_vars = stan::math::count_vars(arg);

  for (int i = 0; i < num_vars; ++i)
    EXPECT_FLOAT_EQ(storage(i), 1.0);
  for (int i = num_vars; i < storage.size(); ++i)
    EXPECT_FLOAT_EQ(storage(i), 0.0);

  EXPECT_EQ(ptr, storage.data() + num_vars);
}

TEST(AgradRev_count_vars, eigen_matrix_var_arg) {
  Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> arg(5, 5);
  for (size_t i = 0; i < arg.size(); ++i) {
    arg(i) = 5.0;
    arg(i).vi_->adj_ = 1.0;
  }

  storage.setZero();
  double* ptr = stan::math::accumulate_adjoints(storage.data(), arg);

  size_t num_vars = stan::math::count_vars(arg);

  for (int i = 0; i < num_vars; ++i)
    EXPECT_FLOAT_EQ(storage(i), 1.0);
  for (int i = num_vars; i < storage.size(); ++i)
    EXPECT_FLOAT_EQ(storage(i), 0.0);

  EXPECT_EQ(ptr, storage.data() + num_vars);
}

TEST(AgradRev_count_vars, std_vector_std_vector_var_arg) {
  std::vector<var> arg_(5, var(5.0));
  for (size_t i = 0; i < arg_.size(); ++i)
    arg_[i].vi_->adj_ = 1.0;
  std::vector<std::vector<var>> arg(5, arg_);

  storage.setZero();
  double* ptr = stan::math::accumulate_adjoints(storage.data(), arg);

  size_t num_vars = stan::math::count_vars(arg);

  for (int i = 0; i < num_vars; ++i)
    EXPECT_FLOAT_EQ(storage(i), 1.0);
  for (int i = num_vars; i < storage.size(); ++i)
    EXPECT_FLOAT_EQ(storage(i), 0.0);

  EXPECT_EQ(ptr, storage.data() + num_vars);
}

TEST(AgradRev_count_vars, std_vector_eigen_vector_var_arg) {
  Eigen::Matrix<var, Eigen::Dynamic, 1> arg_(5);
  for (size_t i = 0; i < arg_.size(); ++i) {
    arg_(i) = 5.0;
    arg_(i).vi_->adj_ = 1.0;
  }
  std::vector<Eigen::Matrix<var, Eigen::Dynamic, 1>> arg(2, arg_);

  storage.setZero();
  double* ptr = stan::math::accumulate_adjoints(storage.data(), arg);

  size_t num_vars = stan::math::count_vars(arg);

  for (int i = 0; i < num_vars; ++i)
    EXPECT_FLOAT_EQ(storage(i), 1.0);
  for (int i = num_vars; i < storage.size(); ++i)
    EXPECT_FLOAT_EQ(storage(i), 0.0);

  EXPECT_EQ(ptr, storage.data() + num_vars);
}

TEST(AgradRev_count_vars, std_vector_eigen_row_vector_var_arg) {
  Eigen::Matrix<var, 1, Eigen::Dynamic> arg_(5);
  for (size_t i = 0; i < arg_.size(); ++i) {
    arg_(i) = 5.0;
    arg_(i).vi_->adj_ = 1.0;
  }
  std::vector<Eigen::Matrix<var, 1, Eigen::Dynamic>> arg(2, arg_);

  storage.setZero();
  double* ptr = stan::math::accumulate_adjoints(storage.data(), arg);

  size_t num_vars = stan::math::count_vars(arg);

  for (int i = 0; i < num_vars; ++i)
    EXPECT_FLOAT_EQ(storage(i), 1.0);
  for (int i = num_vars; i < storage.size(); ++i)
    EXPECT_FLOAT_EQ(storage(i), 0.0);

  EXPECT_EQ(ptr, storage.data() + num_vars);
}

TEST(AgradRev_count_vars, std_vector_eigen_matrix_var_arg) {
  Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> arg_(5, 3);
  for (size_t i = 0; i < arg_.size(); ++i) {
    arg_(i) = 5.0;
    arg_(i).vi_->adj_ = 1.0;
  }
  std::vector<Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>> arg(2, arg_);

  storage.setZero();
  double* ptr = stan::math::accumulate_adjoints(storage.data(), arg);

  size_t num_vars = stan::math::count_vars(arg);

  for (int i = 0; i < num_vars; ++i)
    EXPECT_FLOAT_EQ(storage(i), 1.0);
  for (int i = num_vars; i < storage.size(); ++i)
    EXPECT_FLOAT_EQ(storage(i), 0.0);

  EXPECT_EQ(ptr, storage.data() + num_vars);
}

TEST(AgradRev_count_vars, sum) {
  int arg1 = 1;
  double arg2 = 1.0;
  std::vector<int> arg3(5, 1);
  std::vector<double> arg4(5, 1.0);
  Eigen::VectorXd arg5 = Eigen::VectorXd::Ones(5);
  Eigen::RowVectorXd arg6 = Eigen::RowVectorXd::Ones(5);
  Eigen::MatrixXd arg7 = Eigen::MatrixXd::Ones(5, 5);
  std::vector<std::vector<double>> arg8(2, arg4);
  std::vector<Eigen::VectorXd> arg9(2, arg5);

  var arg10(new vari(5.0));
  arg10.vi_->adj_ = 1.0;
  std::vector<var> arg11(5, new vari(5.0, 1.0));
  for (size_t i = 0; i < arg11.size(); ++i)
    arg11[i].vi_->adj_ = 1.0;
  Eigen::Matrix<var, Eigen::Dynamic, 1> arg12(3);
  for (size_t i = 0; i < arg12.size(); ++i) {
    arg12(i) = 5.0;
    arg12(i).vi_->adj_ = 1.0;
  }
  Eigen::Matrix<var, 1, Eigen::Dynamic> arg13(4);
  for (size_t i = 0; i < arg13.size(); ++i) {
    arg13(i) = 5.0;
    arg13(i).vi_->adj_ = 1.0;
  }
  Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> arg14(5, 3);
  for (size_t i = 0; i < arg14.size(); ++i) {
    arg14(i) = 5.0;
    arg14(i).vi_->adj_ = 1.0;
  }
  std::vector<std::vector<var>> arg15(2, arg11);
  std::vector<Eigen::Matrix<var, Eigen::Dynamic, 1>> arg16(2, arg12);
  std::vector<Eigen::Matrix<var, 1, Eigen::Dynamic>> arg17(2, arg13);
  std::vector<Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>> arg18(2,
                                                                        arg14);

  storage.setZero();
  double* ptr = stan::math::accumulate_adjoints(
      storage.data(), arg1, arg18, arg17, arg2, arg16, arg3, arg15, arg4, arg14,
      arg5, arg13, arg12, arg6, arg11, arg7, arg10, arg8, arg9);

  size_t num_vars = stan::math::count_vars(
      arg1, arg18, arg17, arg2, arg16, arg3, arg15, arg4, arg14, arg5, arg13,
      arg12, arg6, arg11, arg7, arg10, arg8, arg9);

  for (int i = 0; i < num_vars; ++i)
    EXPECT_FLOAT_EQ(storage(i), 1.0);
  for (int i = num_vars; i < storage.size(); ++i)
    EXPECT_FLOAT_EQ(storage(i), 0.0);

  EXPECT_EQ(ptr, storage.data() + num_vars);
}
