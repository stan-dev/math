#include <stan/math/fwd/mat.hpp>
#include <gtest/gtest.h>
#include <stdexcept>
#include <vector>

using stan::math::fvar;

TEST(AgradFwdMatrixAssign, vector_fvar_double) {
  using stan::math::assign;
  using std::vector;

  vector<fvar<double> > y(3);
  y[0] = 1.2;
  y[1] = 100;
  y[2] = -5.1;
  y[0].d_ = 1.0;
  y[1].d_ = 2.0;
  y[2].d_ = 3.0;

  vector<fvar<double> > x(3);
  assign(x, y);
  EXPECT_EQ(3U, x.size());
  EXPECT_EQ(3U, y.size());
  for (size_t i = 0; i < 3; ++i) {
    EXPECT_FLOAT_EQ(y[i].val_, x[i].val_);
    EXPECT_FLOAT_EQ(y[i].d_, x[i].d_);
  }

  vector<fvar<double> > z(2);
  EXPECT_THROW(assign(x, z), std::invalid_argument);
}

TEST(AgradFwdMatrixAssign, eigen_row_vector_fvar_double_to_fvar_double) {
  using stan::math::assign;
  using Eigen::Matrix;
  using Eigen::Dynamic;

  Matrix<fvar<double>, 1, Dynamic> y(3);
  y[0] = 1.2;
  y[1] = 100;
  y[2] = -5.1;
  y[0].d_ = 1.0;
  y[1].d_ = 2.0;
  y[2].d_ = 3.0;

  Matrix<fvar<double>, 1, Dynamic> x(3);
  assign(x, y);
  EXPECT_EQ(3, x.size());
  EXPECT_EQ(3, y.size());
  for (int i = 0; i < 3; ++i) {
    EXPECT_FLOAT_EQ(y[i].val_, x[i].val_);
    EXPECT_FLOAT_EQ(y[i].d_, x[i].d_);
  }
}

TEST(AgradFwdMatrixAssign, eigen_row_vector_fvar_double_shape_mismatch) {
  using stan::math::assign;
  using Eigen::Matrix;
  using Eigen::Dynamic;

  Matrix<fvar<double>, 1, Dynamic> x(3);
  x[0] = 1.2;
  x[1] = 100;
  x[2] = -5.1;

  Matrix<fvar<double>, 1, Dynamic> z(2);
  EXPECT_THROW(assign(x, z), std::invalid_argument);

  Matrix<fvar<double>, Dynamic, 1> zz(3);
  zz << 1, 2, 3;
  EXPECT_THROW(assign(x, zz), std::invalid_argument);

  Matrix<fvar<double>, Dynamic, Dynamic> zzz(3, 1);
  zzz << 1,
         2,
         3;
  EXPECT_THROW(assign(x, zzz), std::invalid_argument);

  Matrix<fvar<double>, Dynamic, Dynamic> zzzz(1, 3);
  EXPECT_THROW(assign(x, zzzz), std::invalid_argument);
}

TEST(AgradFwdMatrixAssign, eigen_matrix_fvar_double_to_fvar_double) {
  using stan::math::assign;
  using Eigen::Matrix;
  using Eigen::Dynamic;

  Matrix<fvar<double>, Dynamic, Dynamic> y(3, 2);
  y << 1.2, 100,
      -5.1, 12,
       1000, -5100;
  y(0, 0).d_ = 1.0;
  y(0, 1).d_ = 2.0;
  y(1, 0).d_ = 3.0;
  y(1, 1).d_ = 4.0;
  y(2, 0).d_ = 5.0;
  y(2, 1).d_ = 6.0;

  Matrix<fvar<double>, Dynamic, Dynamic> x(3, 2);
  assign(x, y);
  EXPECT_EQ(6, x.size());
  EXPECT_EQ(6, y.size());
  EXPECT_EQ(3, x.rows());
  EXPECT_EQ(3, y.rows());
  EXPECT_EQ(2, x.cols());
  EXPECT_EQ(2, y.cols());
  for (size_t i = 0; i < 6; ++i) {
    EXPECT_FLOAT_EQ(y(i).val_, x(i).val_);
    EXPECT_FLOAT_EQ(y(i).d_, x(i).d_);
  }
}

TEST(AgradFwdMatrixAssign, eigen_matrix_fvar_double_shape_mismatch) {
  using stan::math::assign;
  using Eigen::Matrix;
  using Eigen::Dynamic;

  Matrix<fvar<double>, Dynamic, Dynamic> x(2, 3);
  x << 1, 2, 3,
       4, 5, 6;

  Matrix<fvar<double>, 1, Dynamic> z(6);
  z << 1, 2, 3, 4, 5, 6;
  EXPECT_THROW(assign(x, z), std::invalid_argument);
  EXPECT_THROW(assign(z, x), std::invalid_argument);

  Matrix<fvar<double>, Dynamic, 1> zz(6);
  zz << 1, 2, 3, 4, 5, 6;
  EXPECT_THROW(assign(x, zz), std::invalid_argument);
  EXPECT_THROW(assign(zz, x), std::invalid_argument);

  Matrix<fvar<double>, Dynamic, Dynamic> zzz(6, 1);
  zzz << 1,
         2,
         3,
         4,
         5,
         6;
  EXPECT_THROW(assign(x, zzz), std::invalid_argument);
  EXPECT_THROW(assign(zzz, x), std::invalid_argument);
}

TEST(AgradFwdMatrixAssign, block) {
  using Eigen::Matrix;
  using Eigen::Dynamic;
  using stan::math::get_base1_lhs;
  using stan::math::assign;

  Matrix<fvar<double>, Dynamic, Dynamic> m(2, 3);
  m << 1, 2, 3,
       4, 5, 6;

  Matrix<fvar<double>, 1, Dynamic> rv(3);
  rv << 10, 100, 1000;

  assign(get_base1_lhs(m, 1, "m", 1), rv);
  EXPECT_FLOAT_EQ(10.0, m(0, 0).val_);
  EXPECT_FLOAT_EQ(100.0, m(0, 1).val_);
  EXPECT_FLOAT_EQ(1000.0, m(0, 2).val_);
}

TEST(AgradFwdMatrixAssign, vector_vector_fvar_double) {
  using std::vector;
  using stan::math::assign;
  vector<vector<fvar<double> > > x(3, vector<fvar<double> >(2));
  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 2; ++j) {
      x[i][j] = (i + 1) * (j - 10);
      x[i][j].d_ = i+j;
    }

  vector<vector<fvar<double> > > y(3, vector<fvar<double> >(2));

  assign(y, x);
  EXPECT_EQ(3U, y.size());
  for (size_t i = 0; i < 3U; ++i) {
    EXPECT_EQ(2U, y[i].size());
    for (size_t j = 0; j < 2U; ++j) {
      EXPECT_FLOAT_EQ(x[i][j].val_, y[i][j].val_);
      EXPECT_FLOAT_EQ(x[i][j].d_, y[i][j].d_);
    }
  }
}

TEST(AgradFwdMatrixAssign, vector_vector_vector_fvar_double) {
  using std::vector;
  using stan::math::assign;
  vector<vector<vector<fvar<double> > > >
    x(4, vector<vector<fvar<double> > >(3, vector<fvar<double> >(2)));
  for (size_t k = 0; k < 4; ++k)
    for (size_t i = 0; i < 3; ++i)
      for (size_t j = 0; j < 2; ++j) {
        x[k][i][j] = (i + 1) * (j - 10) * (20 * k + 100);
        x[k][i][j].d_ = i+j+k;
      }

  vector<vector<vector<fvar<double> > > >
    y(4, vector<vector<fvar<double> > >(3, vector<fvar<double> >(2)));

  assign(y, x);
  EXPECT_EQ(4U, y.size());
  for (size_t k = 0; k < 4U; ++k) {
    EXPECT_EQ(3U, y[k].size());
    for (size_t i = 0; i < 3U; ++i) {
      EXPECT_EQ(2U, y[k][i].size());
      for (size_t j = 0; j < 2U; ++j) {
        EXPECT_FLOAT_EQ(x[k][i][j].val_, y[k][i][j].val_);
        EXPECT_FLOAT_EQ(x[k][i][j].d_, y[k][i][j].d_);
      }
    }
  }
}

TEST(AgradFwdMatrixAssign, vector_eigen_vector_fvar_double) {
  using std::vector;
  using Eigen::Matrix;
  using Eigen::Dynamic;
  using stan::math::assign;

  vector<Matrix<fvar<double>, Dynamic, 1> >
    x(2, Matrix<fvar<double>, Dynamic, 1>(3));
  for (size_t i = 0; i < 2; ++i)
    for (int j = 0; j < 3; ++j) {
      x[i](j) = (i + 1) * (10 * j + 2);
      x[i](j).d_ = i+j;
    }

  vector<Matrix<fvar<double>, Dynamic, 1> >
    y(2, Matrix<fvar<double>, Dynamic, 1>(3));

  assign(y, x);

  EXPECT_EQ(2U, y.size());
  for (size_t i = 0; i < 2U; ++i) {
    EXPECT_EQ(3U, y[i].size());
    for (size_t j = 0; j < 3U; ++j) {
      EXPECT_FLOAT_EQ(x[i](j).val_, y[i](j).val_);
      EXPECT_FLOAT_EQ(x[i](j).d_, y[i](j).d_);
    }
  }
}

TEST(AgradFwdMatrixAssign, get_assign_row_fvar_double) {
  using Eigen::Matrix;
  using Eigen::Dynamic;
  using stan::math::get_base1_lhs;
  using stan::math::assign;

  Matrix<fvar<double>, Dynamic, Dynamic> m(2, 3);
  m << 1, 2, 3,
       4, 5, 6;

  Matrix<fvar<double>, 1, Dynamic> rv(3);
  rv << 10, 100, 1000;

  assign(get_base1_lhs(m, 1, "m", 1), rv);
  EXPECT_FLOAT_EQ(10.0, m(0, 0).val_);
  EXPECT_FLOAT_EQ(100.0, m(0, 1).val_);
  EXPECT_FLOAT_EQ(1000.0, m(0, 2).val_);
}

TEST(AgradFwdMatrixAssign, vector_fvar_fvar_double) {
  using stan::math::assign;
  using std::vector;

  vector<fvar<fvar<double> > > y(3);
  y[0] = 1.2;
  y[1] = 100;
  y[2] = -5.1;
  y[0].val_.d_ = 1.0;
  y[1].val_.d_ = 2.0;
  y[2].val_.d_ = 3.0;
  y[0].d_.d_ = 1.0;
  y[1].d_.d_ = 2.0;
  y[2].d_.d_ = 3.0;
  y[0].d_ = 1.0;
  y[1].d_ = 2.0;
  y[2].d_ = 3.0;

  vector<fvar<fvar<double> > > x(3);
  assign(x, y);
  EXPECT_EQ(3U, x.size());
  EXPECT_EQ(3U, y.size());
  for (size_t i = 0; i < 3; ++i) {
    EXPECT_FLOAT_EQ(y[i].val_.val_, x[i].val_.val_);
    EXPECT_FLOAT_EQ(y[i].d_.val_, x[i].d_.val_);
    EXPECT_FLOAT_EQ(y[i].val_.d_, x[i].val_.d_);
    EXPECT_FLOAT_EQ(y[i].d_.d_, x[i].d_.d_);
  }

  vector<fvar<fvar<double> > > z(2);
  EXPECT_THROW(assign(x, z), std::invalid_argument);
}

TEST(AgradFwdMatrixAssign,
     eigen_row_vector_fvar_fvar_double_to_fvar_fvar_double) {
  using stan::math::assign;
  using Eigen::Matrix;
  using Eigen::Dynamic;

  Matrix<fvar<fvar<double> >, 1, Dynamic> y(3);
  y[0] = 1.2;
  y[1] = 100;
  y[2] = -5.1;
  y[0].val_.d_ = 1.0;
  y[1].val_.d_ = 2.0;
  y[2].val_.d_ = 3.0;
  y[0].d_.d_ = 1.0;
  y[1].d_.d_ = 2.0;
  y[2].d_.d_ = 3.0;
  y[0].d_ = 1.0;
  y[1].d_ = 2.0;
  y[2].d_ = 3.0;

  Matrix<fvar<fvar<double> >, 1, Dynamic> x(3);
  assign(x, y);
  EXPECT_EQ(3, x.size());
  EXPECT_EQ(3, y.size());
  for (int i = 0; i < 3; ++i) {
    EXPECT_FLOAT_EQ(y[i].val_.val_, x[i].val_.val_);
    EXPECT_FLOAT_EQ(y[i].d_.val_, x[i].d_.val_);
    EXPECT_FLOAT_EQ(y[i].val_.d_, x[i].val_.d_);
    EXPECT_FLOAT_EQ(y[i].d_.d_, x[i].d_.d_);
  }
}

TEST(AgradFwdMatrixAssign, eigen_row_vector_fvar_fvar_double_shape_mismatch) {
  using stan::math::assign;
  using Eigen::Matrix;
  using Eigen::Dynamic;

  Matrix<fvar<fvar<double> >, 1, Dynamic> x(3);
  x[0] = 1.2;
  x[1] = 100;
  x[2] = -5.1;

  Matrix<fvar<fvar<double> >, 1, Dynamic> z(2);
  EXPECT_THROW(assign(x, z), std::invalid_argument);

  Matrix<fvar<fvar<double> >, Dynamic, 1> zz(3);
  zz << 1, 2, 3;
  EXPECT_THROW(assign(x, zz), std::invalid_argument);

  Matrix<fvar<fvar<double> >, Dynamic, Dynamic> zzz(3, 1);
  zzz << 1, 2, 3;
  EXPECT_THROW(assign(x, zzz), std::invalid_argument);

  Matrix<fvar<fvar<double> >, Dynamic, Dynamic> zzzz(1, 3);
  EXPECT_THROW(assign(x, zzzz), std::invalid_argument);
}

TEST(AgradFwdMatrixAssign, eigen_matrix_fvar_fvar_double_to_fvar_fvar_double) {
  using stan::math::assign;
  using Eigen::Matrix;
  using Eigen::Dynamic;

  Matrix<fvar<fvar<double> >, Dynamic, Dynamic> y(3, 2);
  y << 1.2, 100,
      -5.1, 12,
      1000, -5100;
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
  y(0, 0).d_.val_ = 1.0;
  y(0, 1).d_.val_ = 2.0;
  y(1, 0).d_.val_ = 3.0;
  y(1, 1).d_.val_ = 4.0;
  y(2, 0).d_.val_ = 5.0;
  y(2, 1).d_.val_ = 6.0;

  Matrix<fvar<fvar<double> >, Dynamic, Dynamic> x(3, 2);
  assign(x, y);
  EXPECT_EQ(6, x.size());
  EXPECT_EQ(6, y.size());
  EXPECT_EQ(3, x.rows());
  EXPECT_EQ(3, y.rows());
  EXPECT_EQ(2, x.cols());
  EXPECT_EQ(2, y.cols());
  for (size_t i = 0; i < 6; ++i) {
    EXPECT_FLOAT_EQ(y(i).val_.val_, x(i).val_.val_);
    EXPECT_FLOAT_EQ(y(i).d_.val_, x(i).d_.val_);
    EXPECT_FLOAT_EQ(y(i).val_.d_, x(i).val_.d_);
    EXPECT_FLOAT_EQ(y(i).d_.d_, x(i).d_.d_);
  }
}

TEST(AgradFwdMatrixAssign, eigen_matrix_fvar_fvar_double_shape_mismatch) {
  using stan::math::assign;
  using Eigen::Matrix;
  using Eigen::Dynamic;

  Matrix<fvar<fvar<double> >, Dynamic, Dynamic> x(2, 3);
  x << 1, 2, 3,
       4, 5, 6;

  Matrix<fvar<fvar<double> >, 1, Dynamic> z(6);
  z << 1, 2, 3, 4, 5, 6;
  EXPECT_THROW(assign(x, z), std::invalid_argument);
  EXPECT_THROW(assign(z, x), std::invalid_argument);

  Matrix<fvar<fvar<double> >, Dynamic, 1> zz(6);
  zz << 1, 2, 3, 4, 5, 6;
  EXPECT_THROW(assign(x, zz), std::invalid_argument);
  EXPECT_THROW(assign(zz, x), std::invalid_argument);

  Matrix<fvar<fvar<double> >, Dynamic, Dynamic> zzz(6, 1);
  zzz << 1,
         2,
         3,
         4,
         5,
         6;
  EXPECT_THROW(assign(x, zzz), std::invalid_argument);
  EXPECT_THROW(assign(zzz, x), std::invalid_argument);
}

TEST(AgradFwdMatrixAssign, block_fvar_fvar_double) {
  using Eigen::Matrix;
  using Eigen::Dynamic;
  using stan::math::get_base1_lhs;
  using stan::math::assign;

  Matrix<fvar<fvar<double> >, Dynamic, Dynamic> m(2, 3);
  m << 1, 2, 3,
       4, 5, 6;

  Matrix<fvar<fvar<double> >, 1, Dynamic> rv(3);
  rv << 10, 100, 1000;

  assign(get_base1_lhs(m, 1, "m", 1), rv);
  EXPECT_FLOAT_EQ(10.0, m(0, 0).val_.val_);
  EXPECT_FLOAT_EQ(100.0, m(0, 1).val_.val_);
  EXPECT_FLOAT_EQ(1000.0, m(0, 2).val_.val_);
}

TEST(AgradFwdMatrixAssign, vector_vector_fvar_fvar_double) {
  using std::vector;
  using stan::math::assign;
  vector<vector<fvar<fvar<double> > > > x(3, vector<fvar<fvar<double> > >(2));
  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 2; ++j) {
      x[i][j] = (i + 1) * (j - 10);
      x[i][j].val_.d_ = i+j;
      x[i][j].d_ = i+j;
      x[i][j].d_.d_ = i+j;
    }

  vector<vector<fvar<fvar<double> > > > y(3, vector<fvar<fvar<double> > >(2));

  assign(y, x);
  EXPECT_EQ(3U, y.size());
  for (size_t i = 0; i < 3U; ++i) {
    EXPECT_EQ(2U, y[i].size());
    for (size_t j = 0; j < 2U; ++j) {
      EXPECT_FLOAT_EQ(x[i][j].val_.val_, y[i][j].val_.val_);
      EXPECT_FLOAT_EQ(x[i][j].d_.val_, y[i][j].d_.val_);
      EXPECT_FLOAT_EQ(x[i][j].val_.d_, y[i][j].val_.d_);
      EXPECT_FLOAT_EQ(x[i][j].d_.d_, y[i][j].d_.d_);
    }
  }
}

TEST(AgradFwdMatrixAssign, vector_vector_vector_fvar_fvar_double) {
  using std::vector;
  using stan::math::assign;
  vector<vector<vector<fvar<fvar<double> > > > >
    x(4,
      vector<vector<fvar<fvar<double> > > >(3,
                                            vector<fvar<fvar<double> > >(2)));
  for (size_t k = 0; k < 4; ++k)
    for (size_t i = 0; i < 3; ++i)
      for (size_t j = 0; j < 2; ++j) {
        x[k][i][j] = (i + 1) * (j - 10) * (20 * k + 100);
        x[k][i][j].val_.d_ = i+j+k;
        x[k][i][j].d_ = i+j+k;
        x[k][i][j].d_.d_ = i+j+k;
      }

  vector<vector<vector<fvar<fvar<double> > > > >
    y(4,
      vector<vector<fvar<fvar<double> > > >(3,
                                            vector<fvar<fvar<double> > >(2)));

  assign(y, x);
  EXPECT_EQ(4U, y.size());
  for (size_t k = 0; k < 4U; ++k) {
    EXPECT_EQ(3U, y[k].size());
    for (size_t i = 0; i < 3U; ++i) {
      EXPECT_EQ(2U, y[k][i].size());
      for (size_t j = 0; j < 2U; ++j) {
        EXPECT_FLOAT_EQ(x[k][i][j].val_.val_, y[k][i][j].val_.val_);
        EXPECT_FLOAT_EQ(x[k][i][j].d_.val_, y[k][i][j].d_.val_);
        EXPECT_FLOAT_EQ(x[k][i][j].val_.d_, y[k][i][j].val_.d_);
        EXPECT_FLOAT_EQ(x[k][i][j].d_.d_, y[k][i][j].d_.d_);
      }
    }
  }
}

TEST(AgradFwdMatrixAssign, vector_eigen_vector_fvar_fvar_double) {
  using std::vector;
  using Eigen::Matrix;
  using Eigen::Dynamic;
  using stan::math::assign;

  vector<Matrix<fvar<fvar<double> >, Dynamic, 1> >
    x(2, Matrix<fvar<fvar<double> >, Dynamic, 1>(3));
  for (size_t i = 0; i < 2; ++i)
    for (int j = 0; j < 3; ++j) {
      x[i](j) = (i + 1) * (10 * j + 2);
      x[i](j).val_.d_ = i+j;
      x[i](j).d_ = i+j;
      x[i](j).d_.d_ = i+j;
    }
  vector<Matrix<fvar<fvar<double> >, Dynamic, 1> >
    y(2, Matrix<fvar<fvar<double> >, Dynamic, 1>(3));

  assign(y, x);

  EXPECT_EQ(2U, y.size());
  for (size_t i = 0; i < 2U; ++i) {
    EXPECT_EQ(3U, y[i].size());
    for (size_t j = 0; j < 3U; ++j) {
      EXPECT_FLOAT_EQ(x[i](j).val_.val_, y[i](j).val_.val_);
      EXPECT_FLOAT_EQ(x[i](j).d_.val_, y[i](j).d_.val_);
      EXPECT_FLOAT_EQ(x[i](j).val_.d_, y[i](j).val_.d_);
      EXPECT_FLOAT_EQ(x[i](j).d_.d_, y[i](j).d_.d_);
    }
  }
}

TEST(AgradFwdMatrixAssign, get_assign_row_fvar_fvar_double) {
  using Eigen::Matrix;
  using Eigen::Dynamic;
  using stan::math::get_base1_lhs;
  using stan::math::assign;

  Matrix<fvar<fvar<double> >, Dynamic, Dynamic> m(2, 3);
  m << 1, 2, 3,
       4, 5, 6;

  Matrix<fvar<fvar<double> >, 1, Dynamic> rv(3);
  rv << 10, 100, 1000;

  assign(get_base1_lhs(m, 1, "m", 1), rv);
  EXPECT_FLOAT_EQ(10.0, m(0, 0).val_.val_);
  EXPECT_FLOAT_EQ(100.0, m(0, 1).val_.val_);
  EXPECT_FLOAT_EQ(1000.0, m(0, 2).val_.val_);
}
