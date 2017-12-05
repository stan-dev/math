#include <stan/math/mix/mat.hpp>
#include <gtest/gtest.h>
#include <stdexcept>
#include <vector>

using stan::math::fvar;
using stan::math::var;


TEST(AgradMixMatrixAssign, vector_fvar_var) {
  using stan::math::assign;
  using std::vector;

  vector<fvar<var> > y(3);
  y[0] = 1.2;
  y[1] = 100;
  y[2] = -5.1;
  y[0].d_ = 1.0;
  y[1].d_ = 2.0;
  y[2].d_ = 3.0;

  vector<fvar<var> > x(3);
  assign(x, y);
  EXPECT_EQ(3U, x.size());
  EXPECT_EQ(3U, y.size());
  for (size_t i = 0; i < 3; ++i) {
    EXPECT_FLOAT_EQ(y[i].val_.val(), x[i].val_.val());
    EXPECT_FLOAT_EQ(y[i].d_.val(), x[i].d_.val());
  }

  vector<fvar<var> > z(2);
  EXPECT_THROW(assign(x, z), std::invalid_argument);

  std::vector<double> grads;
  std::vector<var> vars;
  vars.push_back(y[0].val_);
  vars.push_back(y[1].val_);
  vars.push_back(y[2].val_);

  x[0].val_.grad(vars, grads);
  EXPECT_FLOAT_EQ(1, grads[0]);
  EXPECT_FLOAT_EQ(0, grads[1]);
  EXPECT_FLOAT_EQ(0, grads[2]);
}

TEST(AgradMixMatrixAssign, eigen_row_vector_fvar_var_to_fvar_var) {
  using stan::math::assign;
  using Eigen::Matrix;
  using Eigen::Dynamic;

  Matrix<fvar<var>, 1, Dynamic> y(3);
  y[0] = 1.2;
  y[1] = 100;
  y[2] = -5.1;
  y[0].d_ = 1.0;
  y[1].d_ = 2.0;
  y[2].d_ = 3.0;

  Matrix<fvar<var>, 1, Dynamic> x(3);
  assign(x, y);
  EXPECT_EQ(3, x.size());
  EXPECT_EQ(3, y.size());
  for (int i = 0; i < 3; ++i) {
    EXPECT_FLOAT_EQ(y[i].val_.val(), x[i].val_.val());
    EXPECT_FLOAT_EQ(y[i].d_.val(), x[i].d_.val());
  }

  std::vector<double> grads;
  std::vector<var> vars;
  vars.push_back(y[0].val_);
  vars.push_back(y[1].val_);
  vars.push_back(y[2].val_);

  x[0].d_.grad(vars, grads);
  EXPECT_FLOAT_EQ(0, grads[0]);
  EXPECT_FLOAT_EQ(0, grads[1]);
  EXPECT_FLOAT_EQ(0, grads[2]);
}
TEST(AgradMixMatrixAssign, eigen_row_vector_fvar_var_shape_mismatch) {
  using stan::math::assign;
  using Eigen::Matrix;
  using Eigen::Dynamic;

  Matrix<fvar<var>, 1, Dynamic> x(3);
  x[0] = 1.2;
  x[1] = 100;
  x[2] = -5.1;

  Matrix<fvar<var>, 1, Dynamic> z(2);
  EXPECT_THROW(assign(x, z), std::invalid_argument);

  Matrix<fvar<var>, Dynamic, 1> zz(3);
  zz << 1, 2, 3;
  EXPECT_THROW(assign(x, zz), std::invalid_argument);

  Matrix<fvar<var>, Dynamic, Dynamic> zzz(3, 1);
  zzz << 1, 2, 3;
  EXPECT_THROW(assign(x, zzz), std::invalid_argument);

  Matrix<fvar<var>, Dynamic, Dynamic> zzzz(1, 3);
  EXPECT_THROW(assign(x, zzzz), std::invalid_argument);
}


TEST(AgradMixMatrixAssign, eigen_matrix_fvar_var_to_fvar_var) {
  using stan::math::assign;
  using Eigen::Matrix;
  using Eigen::Dynamic;

  Matrix<fvar<var>, Dynamic, Dynamic> y(3, 2);
  y << 1.2, 100, -5.1, 12, 1000, -5100;
  y(0, 0).d_ = 1.0;
  y(0, 1).d_ = 2.0;
  y(1, 0).d_ = 3.0;
  y(1, 1).d_ = 4.0;
  y(2, 0).d_ = 5.0;
  y(2, 1).d_ = 6.0;

  Matrix<fvar<var>, Dynamic, Dynamic> x(3, 2);
  assign(x, y);
  EXPECT_EQ(6, x.size());
  EXPECT_EQ(6, y.size());
  EXPECT_EQ(3, x.rows());
  EXPECT_EQ(3, y.rows());
  EXPECT_EQ(2, x.cols());
  EXPECT_EQ(2, y.cols());
  for (size_t i = 0; i < 6; ++i) {
    EXPECT_FLOAT_EQ(y(i).val_.val(), x(i).val_.val());
    EXPECT_FLOAT_EQ(y(i).d_.val(), x(i).d_.val());
  }

  std::vector<double> grads;
  std::vector<var> vars;
  vars.push_back(y(0, 0).val_);
  vars.push_back(y(0, 1).val_);
  vars.push_back(y(1, 0).val_);
  vars.push_back(y(1, 1).val_);
  vars.push_back(y(2, 0).val_);
  vars.push_back(y(2, 1).val_);

  x(0).val_.grad(vars, grads);
  EXPECT_FLOAT_EQ(1, grads[0]);
  EXPECT_FLOAT_EQ(0, grads[1]);
  EXPECT_FLOAT_EQ(0, grads[2]);
  EXPECT_FLOAT_EQ(0, grads[3]);
  EXPECT_FLOAT_EQ(0, grads[4]);
  EXPECT_FLOAT_EQ(0, grads[5]);
}

TEST(AgradMixMatrixAssign, eigen_matrix_fvar_var_shape_mismatch) {
  using stan::math::assign;
  using Eigen::Matrix;
  using Eigen::Dynamic;

  Matrix<fvar<var>, Dynamic, Dynamic> x(2, 3);
  x << 1, 2, 3, 4, 5, 6;

  Matrix<fvar<var>, 1, Dynamic> z(6);
  z << 1, 2, 3, 4, 5, 6;
  EXPECT_THROW(assign(x, z), std::invalid_argument);
  EXPECT_THROW(assign(z, x), std::invalid_argument);

  Matrix<fvar<var>, Dynamic, 1> zz(6);
  zz << 1, 2, 3, 4, 5, 6;
  EXPECT_THROW(assign(x, zz), std::invalid_argument);
  EXPECT_THROW(assign(zz, x), std::invalid_argument);

  Matrix<fvar<var>, Dynamic, Dynamic> zzz(6, 1);
  zzz << 1, 2, 3, 4, 5, 6;
  EXPECT_THROW(assign(x, zzz), std::invalid_argument);
  EXPECT_THROW(assign(zzz, x), std::invalid_argument);
}

TEST(AgradMixMatrixAssign, block_fvar_var) {
  using Eigen::Matrix;
  using Eigen::Dynamic;
  using stan::math::get_base1_lhs;
  using stan::math::assign;

  Matrix<fvar<var>, Dynamic, Dynamic> m(2, 3);
  m << 1, 2, 3, 4, 5, 6;

  Matrix<fvar<var>, 1, Dynamic> rv(3);
  rv << 10, 100, 1000;

  assign(get_base1_lhs(m, 1, "m", 1), rv);
  EXPECT_FLOAT_EQ(10.0, m(0, 0).val_.val());
  EXPECT_FLOAT_EQ(100.0, m(0, 1).val_.val());
  EXPECT_FLOAT_EQ(1000.0, m(0, 2).val_.val());
}


TEST(AgradMixMatrixAssign, vector_vector_fvar_var) {
  using std::vector;
  using stan::math::assign;
  vector<vector<fvar<var> > > x(3, vector<fvar<var> >(2));
  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 2; ++j) {
      x[i][j] = static_cast<double>((i + 1) * (j - 10));
      x[i][j].d_ = static_cast<double>(i+j);
    }

  vector<vector<fvar<var> > > y(3, vector<fvar<var> >(2));

  assign(y, x);
  EXPECT_EQ(3U, y.size());
  for (size_t i = 0; i < 3U; ++i) {
    EXPECT_EQ(2U, y[i].size());
    for (size_t j = 0; j < 2U; ++j) {
      EXPECT_FLOAT_EQ(x[i][j].val_.val(), y[i][j].val_.val());
      EXPECT_FLOAT_EQ(x[i][j].d_.val(), y[i][j].d_.val());
    }
  }

  std::vector<double> grads;
  std::vector<var> vars;
  vars.push_back(x[0][0].val_);
  vars.push_back(x[0][1].val_);
  vars.push_back(x[1][0].val_);
  vars.push_back(x[1][1].val_);
  vars.push_back(x[2][0].val_);
  vars.push_back(x[2][1].val_);

  y[0][0].val_.grad(vars, grads);
  EXPECT_FLOAT_EQ(1, grads[0]);
  EXPECT_FLOAT_EQ(0, grads[1]);
  EXPECT_FLOAT_EQ(0, grads[2]);
  EXPECT_FLOAT_EQ(0, grads[3]);
  EXPECT_FLOAT_EQ(0, grads[4]);
  EXPECT_FLOAT_EQ(0, grads[5]);
}


TEST(AgradMixMatrixAssign, vector_vector_vector_fvar_var) {
  using std::vector;
  using stan::math::assign;
  std::vector<var> vars;
  vector<vector<vector<fvar<var> > > >
    x(4, vector<vector<fvar<var> > >(3, vector<fvar<var> >(2)));
  for (size_t k = 0; k < 4; ++k)
    for (size_t i = 0; i < 3; ++i)
      for (size_t j = 0; j < 2; ++j) {
        x[k][i][j] = static_cast<double>((i + 1) * (j - 10) * (20 * k + 100));
        x[k][i][j].d_ = static_cast<double>(i+j+k);
        vars.push_back(x[k][i][j].val_);
      }

  vector<vector<vector<fvar<var> > > >
    y(4, vector<vector<fvar<var> > >(3, vector<fvar<var> >(2)));

  assign(y, x);
  EXPECT_EQ(4U, y.size());
  for (size_t k = 0; k < 4U; ++k) {
    EXPECT_EQ(3U, y[k].size());
    for (size_t i = 0; i < 3U; ++i) {
      EXPECT_EQ(2U, y[k][i].size());
      for (size_t j = 0; j < 2U; ++j) {
        EXPECT_FLOAT_EQ(x[k][i][j].val_.val(), y[k][i][j].val_.val());
        EXPECT_FLOAT_EQ(x[k][i][j].d_.val(), y[k][i][j].d_.val());
      }
    }
  }

  std::vector<double> grads;

  y[0][0][0].val_.grad(vars, grads);
  EXPECT_FLOAT_EQ(1, grads[0]);
  for (int i = 1; i < 24; i++)
    EXPECT_FLOAT_EQ(0, grads[i]);
}

TEST(AgradMixMatrixAssign, vector_eigen_vector_fvar_var) {
  using std::vector;
  using Eigen::Matrix;
  using Eigen::Dynamic;
  using stan::math::assign;
  std::vector<var> vars;
  vector<Matrix<fvar<var>, Dynamic, 1> > x(2, Matrix<fvar<var>, Dynamic, 1>(3));
  for (size_t i = 0; i < 2; ++i)
    for (int j = 0; j < 3; ++j) {
      x[i](j) = static_cast<double>((i + 1) * (10 * j + 2));
      x[i](j).d_ = static_cast<double>(i+j);
      vars.push_back(x[i](j).val_);
    }
  vector<Matrix<fvar<var>, Dynamic, 1> > y(2, Matrix<fvar<var>, Dynamic, 1>(3));

  assign(y, x);

  EXPECT_EQ(2U, y.size());
  for (size_t i = 0; i < 2U; ++i) {
    EXPECT_EQ(3U, y[i].size());
    for (size_t j = 0; j < 3U; ++j) {
      EXPECT_FLOAT_EQ(x[i](j).val_.val(), y[i](j).val_.val());
      EXPECT_FLOAT_EQ(x[i](j).d_.val(), y[i](j).d_.val());
    }
  }

  std::vector<double> grads;
  y[0](0).val_.grad(vars, grads);
  EXPECT_FLOAT_EQ(1, grads[0]);
  for (int i = 1; i < 6; i++)
    EXPECT_FLOAT_EQ(0, grads[i]);
}

TEST(AgradMixMatrixAssign, get_assign_row_fvar_var) {
  using Eigen::Matrix;
  using Eigen::Dynamic;
  using stan::math::get_base1_lhs;
  using stan::math::assign;

  Matrix<fvar<var>, Dynamic, Dynamic> m(2, 3);
  m << 1, 2, 3, 4, 5, 6;

  Matrix<fvar<var>, 1, Dynamic> rv(3);
  rv << 10, 100, 1000;

  assign(get_base1_lhs(m, 1, "m", 1), rv);
  EXPECT_FLOAT_EQ(10.0, m(0, 0).val_.val());
  EXPECT_FLOAT_EQ(100.0, m(0, 1).val_.val());
  EXPECT_FLOAT_EQ(1000.0, m(0, 2).val_.val());
}

TEST(AgradMixMatrixAssign, vector_fvar_fvar_var) {
  using stan::math::assign;
  using std::vector;

  vector<fvar<fvar<var> > > y(3);
  y[0] = 1.2;
  y[1] = 100;
  y[2] = -5.1;
  y[0].d_.val_ = 1.0;
  y[1].d_.val_ = 2.0;
  y[2].d_.val_ = 3.0;
  y[0].val_.d_ = 1.0;
  y[1].val_.d_ = 2.0;
  y[2].val_.d_ = 3.0;
  y[0].d_.d_ = 1.0;
  y[1].d_.d_ = 2.0;
  y[2].d_.d_ = 3.0;

  vector<fvar<fvar<var> > > x(3);
  assign(x, y);
  EXPECT_EQ(3U, x.size());
  EXPECT_EQ(3U, y.size());
  for (size_t i = 0; i < 3; ++i) {
    EXPECT_FLOAT_EQ(y[i].val_.val_.val(), x[i].val_.val_.val());
    EXPECT_FLOAT_EQ(y[i].d_.val_.val(), x[i].d_.val_.val());
    EXPECT_FLOAT_EQ(y[i].val_.d_.val(), x[i].val_.d_.val());
    EXPECT_FLOAT_EQ(y[i].d_.d_.val(), x[i].d_.d_.val());
  }

  vector<fvar<fvar<var> > > z(2);
  EXPECT_THROW(assign(x, z), std::invalid_argument);

  std::vector<double> grads;
  std::vector<var> vars;
  vars.push_back(y[0].val_.val_);
  vars.push_back(y[1].val_.val_);
  vars.push_back(y[2].val_.val_);

  x[0].val_.val_.grad(vars, grads);
  EXPECT_FLOAT_EQ(1, grads[0]);
  EXPECT_FLOAT_EQ(0, grads[1]);
  EXPECT_FLOAT_EQ(0, grads[2]);
}

TEST(AgradMixMatrixAssign, eigen_row_vector_fvar_fvar_var_to_fvar_fvar_var) {
  using stan::math::assign;
  using Eigen::Matrix;
  using Eigen::Dynamic;

  Matrix<fvar<fvar<var> >, 1, Dynamic> y(3);
  y[0] = 1.2;
  y[1] = 100;
  y[2] = -5.1;
  y[0].d_.val_ = 1.0;
  y[1].d_.val_ = 2.0;
  y[2].d_.val_ = 3.0;
  y[0].val_.d_ = 1.0;
  y[1].val_.d_ = 2.0;
  y[2].val_.d_ = 3.0;
  y[0].d_.d_ = 1.0;
  y[1].d_.d_ = 2.0;
  y[2].d_.d_ = 3.0;

  Matrix<fvar<fvar<var> >, 1, Dynamic> x(3);
  assign(x, y);
  EXPECT_EQ(3, x.size());
  EXPECT_EQ(3, y.size());
  for (int i = 0; i < 3; ++i) {
    EXPECT_FLOAT_EQ(y[i].val_.val_.val(), x[i].val_.val_.val());
    EXPECT_FLOAT_EQ(y[i].d_.val_.val(), x[i].d_.val_.val());
    EXPECT_FLOAT_EQ(y[i].val_.d_.val(), x[i].val_.d_.val());
    EXPECT_FLOAT_EQ(y[i].d_.d_.val(), x[i].d_.d_.val());
  }

  std::vector<double> grads;
  std::vector<var> vars;
  vars.push_back(y[0].val_.val_);
  vars.push_back(y[1].val_.val_);
  vars.push_back(y[2].val_.val_);

  x[0].val_.d_.grad(vars, grads);
  EXPECT_FLOAT_EQ(0, grads[0]);
  EXPECT_FLOAT_EQ(0, grads[1]);
  EXPECT_FLOAT_EQ(0, grads[2]);
}

TEST(AgradMixMatrixAssign, eigen_row_vector_fvar_fvar_var_shape_mismatch) {
  using stan::math::assign;
  using Eigen::Matrix;
  using Eigen::Dynamic;

  Matrix<fvar<fvar<var> >, 1, Dynamic> x(3);
  x[0] = 1.2;
  x[1] = 100;
  x[2] = -5.1;

  Matrix<fvar<fvar<var> >, 1, Dynamic> z(2);
  EXPECT_THROW(assign(x, z), std::invalid_argument);

  Matrix<fvar<fvar<var> >, Dynamic, 1> zz(3);
  zz << 1, 2, 3;
  EXPECT_THROW(assign(x, zz), std::invalid_argument);

  Matrix<fvar<fvar<var> >, Dynamic, Dynamic> zzz(3, 1);
  zzz << 1, 2, 3;
  EXPECT_THROW(assign(x, zzz), std::invalid_argument);

  Matrix<fvar<fvar<var> >, Dynamic, Dynamic> zzzz(1, 3);
  EXPECT_THROW(assign(x, zzzz), std::invalid_argument);
}


TEST(AgradMixMatrixAssign, eigen_matrix_fvar_fvar_var_to_fvar_fvar_var) {
  using stan::math::assign;
  using Eigen::Matrix;
  using Eigen::Dynamic;

  Matrix<fvar<fvar<var> >, Dynamic, Dynamic> y(3, 2);
  y << 1.2, 100, -5.1, 12, 1000, -5100;
  y(0, 0).d_ = 1.0;
  y(0, 1).d_ = 2.0;
  y(1, 0).d_ = 3.0;
  y(1, 1).d_ = 4.0;
  y(2, 0).d_ = 5.0;
  y(2, 1).d_ = 6.0;
  y(0, 0).val_.d_ = 1.0;
  y(0, 1).val_.d_ = 2.0;
  y(1, 0).val_.d_ = 3.0;
  y(1, 1).val_.d_ = 4.0;
  y(2, 0).val_.d_ = 5.0;
  y(2, 1).val_.d_ = 6.0;
  y(0, 0).d_.d_ = 1.0;
  y(0, 1).d_.d_ = 2.0;
  y(1, 0).d_.d_ = 3.0;
  y(1, 1).d_.d_ = 4.0;
  y(2, 0).d_.d_ = 5.0;
  y(2, 1).d_.d_ = 6.0;

  Matrix<fvar<fvar<var> >, Dynamic, Dynamic> x(3, 2);
  assign(x, y);
  EXPECT_EQ(6, x.size());
  EXPECT_EQ(6, y.size());
  EXPECT_EQ(3, x.rows());
  EXPECT_EQ(3, y.rows());
  EXPECT_EQ(2, x.cols());
  EXPECT_EQ(2, y.cols());
  for (size_t i = 0; i < 6; ++i) {
    EXPECT_FLOAT_EQ(y(i).val_.val_.val(), x(i).val_.val_.val());
    EXPECT_FLOAT_EQ(y(i).d_.val_.val(), x(i).d_.val_.val());
    EXPECT_FLOAT_EQ(y(i).val_.d_.val(), x(i).val_.d_.val());
    EXPECT_FLOAT_EQ(y(i).d_.d_.val(), x(i).d_.d_.val());
  }

  std::vector<double> grads;
  std::vector<var> vars;
  vars.push_back(y(0, 0).val_.val_);
  vars.push_back(y(0, 1).val_.val_);
  vars.push_back(y(1, 0).val_.val_);
  vars.push_back(y(1, 1).val_.val_);
  vars.push_back(y(2, 0).val_.val_);
  vars.push_back(y(2, 1).val_.val_);

  x(0).d_.val_.grad(vars, grads);
  EXPECT_FLOAT_EQ(0, grads[0]);
  EXPECT_FLOAT_EQ(0, grads[1]);
  EXPECT_FLOAT_EQ(0, grads[2]);
  EXPECT_FLOAT_EQ(0, grads[3]);
  EXPECT_FLOAT_EQ(0, grads[4]);
  EXPECT_FLOAT_EQ(0, grads[5]);
}

TEST(AgradMixMatrixAssign, eigen_matrix_fvar_fvar_var_shape_mismatch) {
  using stan::math::assign;
  using Eigen::Matrix;
  using Eigen::Dynamic;

  Matrix<fvar<fvar<var> >, Dynamic, Dynamic> x(2, 3);
  x << 1, 2, 3, 4, 5, 6;

  Matrix<fvar<fvar<var> >, 1, Dynamic> z(6);
  z << 1, 2, 3, 4, 5, 6;
  EXPECT_THROW(assign(x, z), std::invalid_argument);
  EXPECT_THROW(assign(z, x), std::invalid_argument);

  Matrix<fvar<fvar<var> >, Dynamic, 1> zz(6);
  zz << 1, 2, 3, 4, 5, 6;
  EXPECT_THROW(assign(x, zz), std::invalid_argument);
  EXPECT_THROW(assign(zz, x), std::invalid_argument);

  Matrix<fvar<fvar<var> >, Dynamic, Dynamic> zzz(6, 1);
  zzz << 1, 2, 3, 4, 5, 6;
  EXPECT_THROW(assign(x, zzz), std::invalid_argument);
  EXPECT_THROW(assign(zzz, x), std::invalid_argument);
}

TEST(AgradMixMatrixAssign, block_fvar_fvar_var) {
  using Eigen::Matrix;
  using Eigen::Dynamic;
  using stan::math::get_base1_lhs;
  using stan::math::assign;

  Matrix<fvar<fvar<var> >, Dynamic, Dynamic> m(2, 3);
  m << 1, 2, 3, 4, 5, 6;

  Matrix<fvar<fvar<var> >, 1, Dynamic> rv(3);
  rv << 10, 100, 1000;

  assign(get_base1_lhs(m, 1, "m", 1), rv);
  EXPECT_FLOAT_EQ(10.0, m(0, 0).val_.val_.val());
  EXPECT_FLOAT_EQ(100.0, m(0, 1).val_.val_.val());
  EXPECT_FLOAT_EQ(1000.0, m(0, 2).val_.val_.val());
}


TEST(AgradMixMatrixAssign, vector_vector_fvar_fvar_var) {
  using std::vector;
  using stan::math::assign;
  vector<vector<fvar<fvar<var> > > > x(3, vector<fvar<fvar<var> > >(2));
  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 2; ++j) {
      x[i][j] = static_cast<double>((i + 1) * (j - 10));
      x[i][j].val_.d_ = static_cast<double>(i+j);
      x[i][j].d_ = static_cast<double>(i+j);
      x[i][j].d_.d_ = static_cast<double>(i+j);
    }

  vector<vector<fvar<fvar<var> > > > y(3, vector<fvar<fvar<var> > >(2));

  assign(y, x);
  EXPECT_EQ(3U, y.size());
  for (size_t i = 0; i < 3U; ++i) {
    EXPECT_EQ(2U, y[i].size());
    for (size_t j = 0; j < 2U; ++j) {
      EXPECT_FLOAT_EQ(x[i][j].val_.val_.val(), y[i][j].val_.val_.val());
      EXPECT_FLOAT_EQ(x[i][j].d_.val_.val(), y[i][j].d_.val_.val());
      EXPECT_FLOAT_EQ(x[i][j].val_.d_.val(), y[i][j].val_.d_.val());
      EXPECT_FLOAT_EQ(x[i][j].d_.d_.val(), y[i][j].d_.d_.val());
    }
  }

  std::vector<double> grads;
  std::vector<var> vars;
  vars.push_back(x[0][0].val_.val_);
  vars.push_back(x[0][1].val_.val_);
  vars.push_back(x[1][0].val_.val_);
  vars.push_back(x[1][1].val_.val_);
  vars.push_back(x[2][0].val_.val_);
  vars.push_back(x[2][1].val_.val_);

  y[0][0].d_.d_.grad(vars, grads);
  EXPECT_FLOAT_EQ(0, grads[0]);
  EXPECT_FLOAT_EQ(0, grads[1]);
  EXPECT_FLOAT_EQ(0, grads[2]);
  EXPECT_FLOAT_EQ(0, grads[3]);
  EXPECT_FLOAT_EQ(0, grads[4]);
  EXPECT_FLOAT_EQ(0, grads[5]);
}


TEST(AgradMixMatrixAssign, vector_vector_vector_fvar_fvar_var) {
  using std::vector;
  using stan::math::assign;
  std::vector<var> vars;
  vector<vector<vector<fvar<fvar<var> > > > >
    x(4, vector<vector<fvar<fvar<var> > > >(3, vector<fvar<fvar<var> > >(2)));
  for (size_t k = 0; k < 4; ++k)
    for (size_t i = 0; i < 3; ++i)
      for (size_t j = 0; j < 2; ++j) {
        x[k][i][j] = static_cast<double>((i + 1) * (j - 10) * (20 * k + 100));
        x[k][i][j].val_.d_ = static_cast<double>(i+j+k);
        x[k][i][j].d_ = static_cast<double>(i+j+k);
        x[k][i][j].d_.d_ = static_cast<double>(i+j+k);
        vars.push_back(x[k][i][j].val_.val_);
      }

  vector<vector<vector<fvar<fvar<var> > > > >
    y(4, vector<vector<fvar<fvar<var> > > >(3, vector<fvar<fvar<var> > >(2)));

  assign(y, x);
  EXPECT_EQ(4U, y.size());
  for (size_t k = 0; k < 4U; ++k) {
    EXPECT_EQ(3U, y[k].size());
    for (size_t i = 0; i < 3U; ++i) {
      EXPECT_EQ(2U, y[k][i].size());
      for (size_t j = 0; j < 2U; ++j) {
        EXPECT_FLOAT_EQ(x[k][i][j].val_.val_.val(), y[k][i][j].val_.val_.val());
        EXPECT_FLOAT_EQ(x[k][i][j].d_.val_.val(), y[k][i][j].d_.val_.val());
        EXPECT_FLOAT_EQ(x[k][i][j].val_.d_.val(), y[k][i][j].val_.d_.val());
        EXPECT_FLOAT_EQ(x[k][i][j].d_.d_.val(), y[k][i][j].d_.d_.val());
      }
    }
  }

  std::vector<double> grads;

  y[0][0][0].val_.val_.grad(vars, grads);
  EXPECT_FLOAT_EQ(1, grads[0]);
  for (int i = 1; i < 24; i++)
    EXPECT_FLOAT_EQ(0, grads[i]);
}

TEST(AgradMixMatrixAssign, vector_eigen_vector_fvar_fvar_var) {
  using std::vector;
  using Eigen::Matrix;
  using Eigen::Dynamic;
  using stan::math::assign;
  std::vector<var> vars;
  vector<Matrix<fvar<fvar<var> >, Dynamic, 1> >
      x(2, Matrix<fvar<fvar<var> >, Dynamic, 1>(3));
  for (size_t i = 0; i < 2; ++i)
    for (int j = 0; j < 3; ++j) {
      x[i](j) = static_cast<double>((i + 1) * (10 * j + 2));
      x[i](j).val_.d_ = static_cast<double>(i+j);
      x[i](j).d_ = static_cast<double>(i+j);
      x[i](j).d_.d_ = static_cast<double>(i+j);
      vars.push_back(x[i](j).val_.val_);
    }
  vector<Matrix<fvar<fvar<var> >, Dynamic, 1> >
      y(2, Matrix<fvar<fvar<var> >, Dynamic, 1>(3));

  assign(y, x);

  EXPECT_EQ(2U, y.size());
  for (size_t i = 0; i < 2U; ++i) {
    EXPECT_EQ(3U, y[i].size());
    for (size_t j = 0; j < 3U; ++j) {
      EXPECT_FLOAT_EQ(x[i](j).val_.val_.val(), y[i](j).val_.val_.val());
      EXPECT_FLOAT_EQ(x[i](j).d_.val_.val(), y[i](j).d_.val_.val());
      EXPECT_FLOAT_EQ(x[i](j).val_.d_.val(), y[i](j).val_.d_.val());
      EXPECT_FLOAT_EQ(x[i](j).d_.d_.val(), y[i](j).d_.d_.val());
    }
  }

  std::vector<double> grads;
  y[0](0).val_.val_.grad(vars, grads);
  EXPECT_FLOAT_EQ(1, grads[0]);
  for (int i = 1; i < 6; i++)
    EXPECT_FLOAT_EQ(0, grads[i]);
}

TEST(AgradMixMatrixAssign, get_assign_row_fvar_fvar_var) {
  using Eigen::Matrix;
  using Eigen::Dynamic;
  using stan::math::get_base1_lhs;
  using stan::math::assign;

  Matrix<fvar<fvar<var> >, Dynamic, Dynamic> m(2, 3);
  m << 1, 2, 3, 4, 5, 6;

  Matrix<fvar<fvar<var> >, 1, Dynamic> rv(3);
  rv << 10, 100, 1000;

  assign(get_base1_lhs(m, 1, "m", 1), rv);
  EXPECT_FLOAT_EQ(10.0, m(0, 0).val_.val_.val());
  EXPECT_FLOAT_EQ(100.0, m(0, 1).val_.val_.val());
  EXPECT_FLOAT_EQ(1000.0, m(0, 2).val_.val_.val());
}
