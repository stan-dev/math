#include <stan/math/fwd.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(ProbDistributionsOrdLog, fd_fd) {
  using stan::math::fvar;
  using stan::math::ordered_logistic_lpmf;
  using stan::math::vector_d;
  using stan::math::vector_fd;
  using stan::math::vector_ffd;

  int y = 1;

  fvar<double> lam_fd = -1.32;
  lam_fd.d_ = 1.0;

  vector_fd c_fd(3);
  c_fd << -0.95, -0.10, 0.95;
  for (int i = 0; i < 3; i++)
    c_fd[i].d_ = 1.0;

  fvar<fvar<double>> lam_ffd;
  lam_ffd.val_ = -1.32;
  lam_ffd.d_ = 1.0;
  lam_ffd.val_.d_ = 1.0;

  vector_ffd c_ffd(3);
  c_ffd << -0.95, -0.10, 0.95;
  for (int i = 0; i < 3; i++) {
    c_ffd[i].d_ = 1.0;
    c_ffd[i].val_.d_ = 1.0;
  }

  fvar<double> out_fd = ordered_logistic_lpmf(y, lam_fd, c_fd);

  EXPECT_FLOAT_EQ(out_fd.val_, -0.52516294973063);
  EXPECT_FLOAT_EQ(out_fd.d_ + 1, 0.0 + 1);

  fvar<fvar<double>> out_ffd = ordered_logistic_lpmf(y, lam_ffd, c_ffd);

  EXPECT_FLOAT_EQ(out_ffd.val_.val_, -0.52516294973063);
  EXPECT_FLOAT_EQ(out_ffd.d_.val_ + 1, 0.0 + 1);
}

TEST(ProbDistributionsOrdLog, fd_d) {
  using stan::math::fvar;
  using stan::math::ordered_logistic_lpmf;
  using stan::math::vector_d;
  using stan::math::vector_fd;
  using stan::math::vector_ffd;

  int y = 1;

  fvar<double> lam_fd = -1.32;
  lam_fd.d_ = 1.0;

  double lam_d = -1.32;

  vector_fd c_fd(3);
  c_fd << -0.95, -0.10, 0.95;
  for (int i = 0; i < 3; i++)
    c_fd[i].d_ = 1.0;

  vector_d c_d(3);
  c_d << -0.95, -0.10, 0.95;

  fvar<fvar<double>> lam_ffd;
  lam_ffd.val_ = -1.32;
  lam_ffd.d_ = 1.0;
  lam_ffd.val_.d_ = 1.0;

  vector_ffd c_ffd(3);
  c_ffd << -0.95, -0.10, 0.95;
  for (int i = 0; i < 3; i++) {
    c_ffd[i].d_ = 1.0;
    c_ffd[i].val_.d_ = 1.0;
  }

  fvar<double> out = ordered_logistic_lpmf(y, lam_fd, c_d);

  EXPECT_FLOAT_EQ(out.val_, -0.52516294973063);
  EXPECT_FLOAT_EQ(out.d_, -0.40854102156722);

  out = ordered_logistic_lpmf(y, lam_d, c_fd);

  EXPECT_FLOAT_EQ(out.val_, -0.52516294973063);
  EXPECT_FLOAT_EQ(out.d_, 0.40854102156722);

  fvar<fvar<double>> out_ffd = ordered_logistic_lpmf(y, lam_ffd, c_d);

  EXPECT_FLOAT_EQ(out_ffd.val_.val_, -0.52516294973063);
  EXPECT_FLOAT_EQ(out_ffd.d_.val_, -0.40854102156722);

  out_ffd = ordered_logistic_lpmf(y, lam_d, c_ffd);

  EXPECT_FLOAT_EQ(out_ffd.val_.val_, -0.52516294973063);
  EXPECT_FLOAT_EQ(out_ffd.d_.val_, 0.40854102156722);
}

TEST(ProbDistributionsOrdLog, fd_fd_vec) {
  using stan::math::fvar;
  using stan::math::ordered_logistic_lpmf;
  using stan::math::vector_d;
  using stan::math::vector_fd;
  using stan::math::vector_ffd;

  std::vector<int> y{1, 2, 3, 4};

  vector_fd lam_fd(4);
  lam_fd << -1.32, -0.05, 0.56, 1.55;
  for (int i = 0; i < 4; i++)
    lam_fd[i].d_ = 1.0;

  vector_fd c_fd(3);
  c_fd << -0.95, -0.10, 0.95;
  for (int i = 0; i < 3; i++)
    c_fd[i].d_ = 1.0;

  vector_ffd lam_ffd(4);
  lam_ffd << -1.32, -0.05, 0.56, 1.55;
  for (int i = 0; i < 4; i++) {
    lam_ffd[i].d_ = 1.0;
    lam_ffd[i].val_.d_ = 1.0;
  }

  vector_ffd c_ffd(3);
  c_ffd << -0.95, -0.10, 0.95;
  for (int i = 0; i < 3; i++) {
    c_ffd[i].d_ = 1.0;
    c_ffd[i].val_.d_ = 1.0;
  }

  fvar<double> out_fd = ordered_logistic_lpmf(y, lam_fd, c_fd);

  EXPECT_FLOAT_EQ(out_fd.val_, -3.9442226297351447171906);
  EXPECT_FLOAT_EQ(out_fd.d_ + 1, 0.0 + 1);

  fvar<fvar<double>> out_ffd = ordered_logistic_lpmf(y, lam_ffd, c_ffd);

  EXPECT_FLOAT_EQ(out_ffd.val_.val_, -3.9442226297351447171906);
  EXPECT_FLOAT_EQ(out_ffd.d_.val_ + 1, 0.0 + 1);
}

TEST(ProbDistributionsOrdLog, fd_d_vec) {
  using stan::math::fvar;
  using stan::math::ordered_logistic_lpmf;
  using stan::math::vector_d;
  using stan::math::vector_fd;
  using stan::math::vector_ffd;

  std::vector<int> y{1, 2, 3, 4};

  vector_fd lam_fd(4);
  lam_fd << -1.32, -0.05, 0.56, 1.55;
  for (int i = 0; i < 4; i++)
    lam_fd[i].d_ = 1.0;

  vector_fd c_fd(3);
  c_fd << -0.95, -0.10, 0.95;
  for (int i = 0; i < 3; i++)
    c_fd[i].d_ = 1.0;

  vector_ffd lam_ffd(4);
  lam_ffd << -1.32, -0.05, 0.56, 1.55;
  for (int i = 0; i < 4; i++) {
    lam_ffd[i].d_ = 1.0;
    lam_ffd[i].val_.d_ = 1.0;
  }

  vector_ffd c_ffd(3);
  c_ffd << -0.95, -0.10, 0.95;
  for (int i = 0; i < 3; i++) {
    c_ffd[i].d_ = 1.0;
    c_ffd[i].val_.d_ = 1.0;
  }

  vector_d lam_d(4);
  lam_d << -1.32, -0.05, 0.56, 1.55;

  vector_d c_d(3);
  c_d << -0.95, -0.10, 0.95;

  fvar<double> out = ordered_logistic_lpmf(y, lam_fd, c_d);

  EXPECT_FLOAT_EQ(out.val_, -3.9442226297351447171906);
  EXPECT_FLOAT_EQ(out.d_, -0.340621916056826);

  out = ordered_logistic_lpmf(y, lam_d, c_fd);

  EXPECT_FLOAT_EQ(out.val_, -3.9442226297351447171906);
  EXPECT_FLOAT_EQ(out.d_, 0.340621916056825);

  fvar<fvar<double>> out_ffd = ordered_logistic_lpmf(y, lam_ffd, c_d);

  EXPECT_FLOAT_EQ(out_ffd.val_.val_, -3.9442226297351447171906);
  EXPECT_FLOAT_EQ(out_ffd.d_.val_, -0.340621916056826);

  out_ffd = ordered_logistic_lpmf(y, lam_d, c_ffd);

  EXPECT_FLOAT_EQ(out_ffd.val_.val_, -3.9442226297351447171906);
  EXPECT_FLOAT_EQ(out_ffd.d_.val_, 0.340621916056826);
}

TEST(ProbDistributionsOrdLog, fd_fd_stvec) {
  using stan::math::fvar;
  using stan::math::ordered_logistic_lpmf;
  using stan::math::vector_d;
  using stan::math::vector_fd;
  using stan::math::vector_ffd;

  std::vector<int> y{1, 2, 3, 4};

  vector_fd lam_fd(4);
  lam_fd << -1.32, -0.05, 0.56, 1.55;
  for (int i = 0; i < 4; i++)
    lam_fd[i].d_ = 1.0;

  vector_fd c1_fd(3);
  c1_fd << -0.95, -0.10, 0.95;
  for (int i = 0; i < 3; i++)
    c1_fd[i].d_ = 1.0;

  vector_fd c2_fd(3);
  c2_fd << -2.15, -1.56, -0.23;
  for (int i = 0; i < 3; i++)
    c2_fd[i].d_ = 1.0;

  vector_fd c3_fd(3);
  c3_fd << -1.52, 0.21, 1.32;
  for (int i = 0; i < 3; i++)
    c3_fd[i].d_ = 1.0;

  vector_fd c4_fd(3);
  c4_fd << -2.15, -2.00, -0.51;
  for (int i = 0; i < 3; i++)
    c4_fd[i].d_ = 1.0;

  vector_ffd lam_ffd(4);
  lam_ffd << -1.32, -0.05, 0.56, 1.55;
  for (int i = 0; i < 4; i++) {
    lam_ffd[i].d_ = 1.0;
    lam_ffd[i].val_.d_ = 1.0;
  }

  vector_ffd c1_ffd(3);
  c1_ffd << -0.95, -0.10, 0.95;
  for (int i = 0; i < 3; i++) {
    c1_ffd[i].d_ = 1.0;
    c1_ffd[i].val_.d_ = 1.0;
  }

  vector_ffd c2_ffd(3);
  c2_ffd << -2.15, -1.56, -0.23;
  for (int i = 0; i < 3; i++) {
    c2_ffd[i].d_ = 1.0;
    c2_ffd[i].val_.d_ = 1.0;
  }

  vector_ffd c3_ffd(3);
  c3_ffd << -1.52, 0.21, 1.32;
  for (int i = 0; i < 3; i++) {
    c3_ffd[i].d_ = 1.0;
    c3_ffd[i].val_.d_ = 1.0;
  }

  vector_ffd c4_ffd(3);
  c4_ffd << -2.15, -2.00, -0.51;
  for (int i = 0; i < 3; i++) {
    c4_ffd[i].d_ = 1.0;
    c4_ffd[i].val_.d_ = 1.0;
  }

  std::vector<vector_fd> std_c_fd{c1_fd, c2_fd, c3_fd, c4_fd};
  std::vector<vector_ffd> std_c_ffd{c1_ffd, c2_ffd, c3_ffd, c4_ffd};

  fvar<double> out_fd = ordered_logistic_lpmf(y, lam_fd, std_c_fd);

  EXPECT_FLOAT_EQ(out_fd.val_, -4.59528666878753);
  EXPECT_FLOAT_EQ(out_fd.d_ + 1, 0.0 + 1);

  fvar<fvar<double>> out_ffd = ordered_logistic_lpmf(y, lam_ffd, std_c_ffd);

  EXPECT_FLOAT_EQ(out_ffd.val_.val_, -4.59528666878753);
  EXPECT_FLOAT_EQ(out_ffd.d_.val_ + 1, 0.0 + 1);
}

TEST(ProbDistributionsOrdLog, fd_d_stvec) {
  using stan::math::fvar;
  using stan::math::ordered_logistic_lpmf;
  using stan::math::vector_d;
  using stan::math::vector_fd;
  using stan::math::vector_ffd;

  std::vector<int> y{1, 2, 3, 4};

  vector_fd lam_fd(4);
  lam_fd << -2.55, -1.63, 0.28, 1.09;
  for (int i = 0; i < 4; i++)
    lam_fd[i].d_ = 1.0;

  vector_fd c1_fd(3);
  c1_fd << -1.52, -0.88, -0.11;
  for (int i = 0; i < 3; i++)
    c1_fd[i].d_ = 1.0;

  vector_fd c2_fd(3);
  c2_fd << -5.11, -2.36, -0.55;
  for (int i = 0; i < 3; i++)
    c2_fd[i].d_ = 1.0;

  vector_fd c3_fd(3);
  c3_fd << -0.13, 0.53, 1.56;
  for (int i = 0; i < 3; i++)
    c3_fd[i].d_ = 1.0;

  vector_fd c4_fd(3);
  c4_fd << -3.51, -1.63, -0.51;
  for (int i = 0; i < 3; i++)
    c4_fd[i].d_ = 1.0;

  vector_ffd lam_ffd(4);
  lam_ffd << -2.55, -1.63, 0.28, 1.09;
  for (int i = 0; i < 4; i++) {
    lam_ffd[i].d_ = 1.0;
    lam_ffd[i].val_.d_ = 1.0;
  }

  vector_ffd c1_ffd(3);
  c1_ffd << -1.52, -0.88, -0.11;
  for (int i = 0; i < 3; i++) {
    c1_ffd[i].d_ = 1.0;
    c1_ffd[i].val_.d_ = 1.0;
  }

  vector_ffd c2_ffd(3);
  c2_ffd << -5.11, -2.36, -0.55;
  for (int i = 0; i < 3; i++) {
    c2_ffd[i].d_ = 1.0;
    c2_ffd[i].val_.d_ = 1.0;
  }

  vector_ffd c3_ffd(3);
  c3_ffd << -0.13, 0.53, 1.56;
  for (int i = 0; i < 3; i++) {
    c3_ffd[i].d_ = 1.0;
    c3_ffd[i].val_.d_ = 1.0;
  }

  vector_ffd c4_ffd(3);
  c4_ffd << -3.51, -1.63, -0.51;
  for (int i = 0; i < 3; i++) {
    c4_ffd[i].d_ = 1.0;
    c4_ffd[i].val_.d_ = 1.0;
  }

  vector_d lam_d(4);
  lam_d << -2.55, -1.63, 0.28, 1.09;

  vector_d c1_d(3);
  c1_d << -1.52, -0.88, -0.11;

  vector_d c2_d(3);
  c2_d << -5.11, -2.36, -0.55;

  vector_d c3_d(3);
  c3_d << -0.13, 0.53, 1.56;

  vector_d c4_d(3);
  c4_d << -3.51, -1.63, -0.51;

  std::vector<vector_fd> std_c_fd{c1_fd, c2_fd, c3_fd, c4_fd};
  std::vector<vector_ffd> std_c_ffd{c1_ffd, c2_ffd, c3_ffd, c4_ffd};
  std::vector<vector_d> std_c_d{c1_d, c2_d, c3_d, c4_d};

  fvar<double> out_fd = ordered_logistic_lpmf(y, lam_fd, std_c_d);

  EXPECT_FLOAT_EQ(out_fd.val_, -3.22180483131937);
  EXPECT_FLOAT_EQ(out_fd.d_, -0.395394804770271);

  out_fd = ordered_logistic_lpmf(y, lam_d, std_c_fd);

  EXPECT_FLOAT_EQ(out_fd.val_, -3.22180483131937);
  EXPECT_FLOAT_EQ(out_fd.d_, 0.395394804770271);

  fvar<fvar<double>> out_ffd = ordered_logistic_lpmf(y, lam_ffd, std_c_d);

  EXPECT_FLOAT_EQ(out_ffd.val_.val_, -3.22180483131937);
  EXPECT_FLOAT_EQ(out_ffd.d_.val_, -0.395394804770271);

  out_ffd = ordered_logistic_lpmf(y, lam_d, std_c_ffd);

  EXPECT_FLOAT_EQ(out_ffd.val_.val_, -3.22180483131937);
  EXPECT_FLOAT_EQ(out_ffd.d_.val_, 0.395394804770271);
}
