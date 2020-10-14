#include <stan/math.hpp>
#include <gtest/gtest.h>

namespace stan {
    namespace test {
      namespace indexing {
        template <typename T1, typename I, typename T2>
        void test_throw(T1& lhs, const I& idxs, const T2& rhs) {
          EXPECT_THROW(stan::math::assign(lhs, idxs, rhs), std::out_of_range);
        }

        template <typename T1, typename I, typename T2>
        void test_throw_ia(T1& lhs, const I& idxs, const T2& rhs) {
          EXPECT_THROW(stan::math::assign(lhs, idxs, rhs), std::invalid_argument);
        }
      }
  }
}



TEST(ModelIndexing, lvalueNil) {
  using stan::math::assign;
  using stan::math::nil_index_list;
  double x = 3;
  double y = 5;
  assign(x, nil_index_list(), y);
  EXPECT_FLOAT_EQ(5, x);

  std::vector<double> xs;
  xs.push_back(3);
  xs.push_back(5);
  std::vector<double> ys;
  ys.push_back(13);
  ys.push_back(15);
  assign(xs, nil_index_list(), ys);
  EXPECT_FLOAT_EQ(ys[0], xs[0]);
  EXPECT_FLOAT_EQ(ys[1], xs[1]);
}

TEST(ModelIndexing, lvalueUni) {
  using stan::math::assign;
  using stan::math::index_list;
  using stan::math::index_uni;
  std::vector<double> xs;
  xs.push_back(3);
  xs.push_back(5);
  xs.push_back(7);
  double y = 15;
  assign(xs, index_list(index_uni(2)), y);
  EXPECT_FLOAT_EQ(y, xs[1]);

  stan::test::indexing::test_throw(xs, index_list(index_uni(0)), y);
  stan::test::indexing::test_throw(xs, index_list(index_uni(4)), y);
}

TEST(ModelIndexing, lvalueUniEigen) {
  using stan::math::assign;
  using stan::math::index_list;
  using stan::math::index_uni;
  Eigen::VectorXd xs(3);
  xs << 3, 5, 7;
  double y = 15;
  assign(xs, index_list(index_uni(2)), y);
  EXPECT_FLOAT_EQ(y, xs[1]);
  double z = 10;
  assign(xs.segment(0, 3), index_list(index_uni(2)), z);
  EXPECT_FLOAT_EQ(z, xs[1]);

  stan::test::indexing::test_throw(xs, index_list(index_uni(0)), y);
  stan::test::indexing::test_throw(xs, index_list(index_uni(4)), y);
}

TEST(model_indexing, assign_eigvec_scalar_uni_index_segment) {
  using stan::math::assign;
  using stan::math::index_list;
  using stan::math::index_uni;
  Eigen::VectorXd lhs_x(5);
  lhs_x << 0, 1, 2, 3, 4;
  double y = 13;
  assign(lhs_x.segment(0, 5), index_list(index_uni(3)), y);
  EXPECT_FLOAT_EQ(y, lhs_x(2));

  stan::test::indexing::test_throw(lhs_x, index_list(index_uni(0)), y);
  stan::test::indexing::test_throw(lhs_x, index_list(index_uni(6)), y);
}

TEST(model_indexing, assign_eigrowvec_scalar_uni_index_segment) {
  using stan::math::assign;
  using stan::math::index_list;
  using stan::math::index_uni;

  Eigen::RowVectorXd lhs_x(5);
  lhs_x << 0, 1, 2, 3, 4;
  double y = 13;
  assign(lhs_x.segment(0, 5), index_list(index_uni(3)), y);
  EXPECT_FLOAT_EQ(y, lhs_x(2));
  stan::test::indexing::test_throw(lhs_x, index_list(index_uni(0)), y);
  stan::test::indexing::test_throw(lhs_x, index_list(index_uni(6)), y);
}

TEST(ModelIndexing, lvalueUniUni) {
  using stan::math::assign;
  using stan::math::index_list;
  using stan::math::index_max;
  using stan::math::index_min;
  using stan::math::index_min_max;
  using stan::math::index_multi;
  using stan::math::index_omni;
  using stan::math::index_uni;
  using stan::math::nil_index_list;
  std::vector<double> xs0;
  xs0.push_back(0.0);
  xs0.push_back(0.1);
  xs0.push_back(0.2);

  std::vector<double> xs1;
  xs1.push_back(1.0);
  xs1.push_back(1.1);
  xs1.push_back(1.2);

  std::vector<std::vector<double> > xs;
  xs.push_back(xs0);
  xs.push_back(xs1);

  double y = 15;
  assign(xs, index_list(index_uni(2), index_uni(3)), y);
  EXPECT_FLOAT_EQ(y, xs[1][2]);

  stan::test::indexing::test_throw(xs, index_list(index_uni(0), index_uni(3)), y);
  stan::test::indexing::test_throw(xs, index_list(index_uni(2), index_uni(0)), y);
  stan::test::indexing::test_throw(xs, index_list(index_uni(10), index_uni(3)), y);
  stan::test::indexing::test_throw(xs, index_list(index_uni(2), index_uni(10)), y);
}

TEST(ModelIndexing, lvalueUniUniEigen) {
  using stan::math::assign;
  using stan::math::index_list;
  using stan::math::index_uni;
  Eigen::VectorXd xs0(3);
  xs0 << 0.0, 0.1, 0.2;

  Eigen::VectorXd xs1(3);
  xs1 << 1.0, 1.1, 1.2;

  std::vector<Eigen::VectorXd> xs;
  xs.push_back(xs0);
  xs.push_back(xs1);

  double y = 15;
  assign(xs, index_list(index_uni(2), index_uni(3)), y);
  EXPECT_FLOAT_EQ(y, xs[1][2]);

  stan::test::indexing::test_throw(xs, index_list(index_uni(0), index_uni(3)), y);
  stan::test::indexing::test_throw(xs, index_list(index_uni(2), index_uni(0)), y);
  stan::test::indexing::test_throw(xs, index_list(index_uni(10), index_uni(3)), y);
  stan::test::indexing::test_throw(xs, index_list(index_uni(2), index_uni(10)), y);
}

TEST(ModelIndexing, lvalueMulti) {
  using stan::math::assign;
  using stan::math::index_list;
  using stan::math::index_max;
  using stan::math::index_min;
  using stan::math::index_multi;
  std::vector<double> x;
  for (int i = 0; i < 10; ++i)
    x.push_back(i);

  std::vector<double> y;
  y.push_back(8.1);
  y.push_back(9.1);

  assign(x, index_list(index_min(9)), y);
  EXPECT_FLOAT_EQ(y[0], x[8]);
  EXPECT_FLOAT_EQ(y[1], x[9]);
  stan::test::indexing::test_throw_ia(x, index_list(index_min(0)), y);

  assign(x, index_list(index_max(2)), y);
  EXPECT_FLOAT_EQ(y[0], x[0]);
  EXPECT_FLOAT_EQ(y[1], x[1]);
  EXPECT_FLOAT_EQ(2, x[2]);
  stan::test::indexing::test_throw_ia(x, index_list(index_max(10)), y);

  std::vector<int> ns;
  ns.push_back(4);
  ns.push_back(6);
  assign(x, index_list(index_multi(ns)), y);
  EXPECT_FLOAT_EQ(y[0], x[3]);
  EXPECT_FLOAT_EQ(y[1], x[5]);

  ns[0] = 0;
  stan::test::indexing::test_throw(x, index_list(index_multi(ns)), y);

  ns[0] = 11;
  stan::test::indexing::test_throw(x, index_list(index_multi(ns)), y);

  ns.push_back(3);
  stan::test::indexing::test_throw_ia(x, index_list(index_multi(ns)), y);
}

TEST(ModelIndexing, lvalueMultiEigen) {
  using stan::math::assign;
  using stan::math::index_list;
  using stan::math::index_max;
  using stan::math::index_min;
  using stan::math::index_multi;
  Eigen::VectorXd x(10);
  for (int i = 0; i < 10; ++i) {
    x(i) = i;
  }

  Eigen::VectorXd y(2);
  y << 8.1, 9.1;

  assign(x, index_list(index_min(9)), y);
  EXPECT_FLOAT_EQ(y[0], x[8]);
  EXPECT_FLOAT_EQ(y[1], x[9]);
  stan::test::indexing::test_throw_ia(x, index_list(index_min(0)), y);

  assign(x, index_list(index_max(2)), y);
  EXPECT_FLOAT_EQ(y[0], x[0]);
  EXPECT_FLOAT_EQ(y[1], x[1]);
  EXPECT_FLOAT_EQ(2, x[2]);
  stan::test::indexing::test_throw_ia(x, index_list(index_max(10)), y);

  std::vector<int> ns;
  ns.push_back(4);
  ns.push_back(6);
  assign(x, index_list(index_multi(ns)), y);
  EXPECT_FLOAT_EQ(y[0], x[3]);
  EXPECT_FLOAT_EQ(y[1], x[5]);

  ns[0] = 0;
  stan::test::indexing::test_throw(x, index_list(index_multi(ns)), y);

  ns[0] = 11;
  stan::test::indexing::test_throw(x, index_list(index_multi(ns)), y);

  ns.push_back(3);
  stan::test::indexing::test_throw_ia(x, index_list(index_multi(ns)), y);
}

TEST(ModelIndexing, lvalueMultiMulti) {
  using stan::math::assign;
  using stan::math::index_list;
  using stan::math::index_max;
  using stan::math::index_min;
  std::vector<std::vector<double> > xs;
  for (int i = 0; i < 10; ++i) {
    std::vector<double> xsi;
    for (int j = 0; j < 20; ++j)
      xsi.push_back(i + j / 10.0);
    xs.push_back(xsi);
  }

  std::vector<std::vector<double> > ys;
  for (int i = 0; i < 2; ++i) {
    std::vector<double> ysi;
    for (int j = 0; j < 3; ++j)
      ysi.push_back(10 + i + j / 10.0);
    ys.push_back(ysi);
  }

  assign(xs, index_list(index_min(9), index_max(3)), ys);

  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 3; ++j)
      EXPECT_FLOAT_EQ(ys[i][j], xs[8 + i][j]);

  stan::test::indexing::test_throw_ia(xs, index_list(index_min(7), index_max(3)), ys);
  stan::test::indexing::test_throw_ia(xs, index_list(index_min(9), index_max(2)), ys);
}

TEST(ModelIndexing, lvalueMultiMultiEigen) {
  using stan::math::assign;
  using stan::math::index_list;
  using stan::math::index_max;
  using stan::math::index_min;
  std::vector<Eigen::VectorXd> xs;
  for (int i = 0; i < 10; ++i) {
    Eigen::VectorXd xsi(20);
    for (int j = 0; j < 20; ++j) {
      xsi(j) = (i + j / 10.0);
    }
    xs.push_back(xsi);
  }

  std::vector<Eigen::VectorXd> ys;
  for (int i = 0; i < 2; ++i) {
    Eigen::VectorXd ysi(3);
    for (int j = 0; j < 3; ++j) {
      ysi(j) = (10 + i + j / 10.0);
    }
    ys.push_back(ysi);
  }

  assign(xs, index_list(index_min(9), index_max(3)), ys);

  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 3; ++j)
      EXPECT_FLOAT_EQ(ys[i][j], xs[8 + i][j]);

  stan::test::indexing::test_throw_ia(xs, index_list(index_min(7), index_max(3)), ys);
  stan::test::indexing::test_throw_ia(xs, index_list(index_min(9), index_max(2)), ys);
}

TEST(ModelIndexing, lvalueUniMulti) {
  using stan::math::assign;
  using stan::math::index_list;
  using stan::math::index_min_max;
  using stan::math::index_uni;
  std::vector<std::vector<double> > xs;
  for (int i = 0; i < 10; ++i) {
    std::vector<double> xsi;
    for (int j = 0; j < 20; ++j)
      xsi.push_back(i + j / 10.0);
    xs.push_back(xsi);
  }

  std::vector<double> ys;
  for (int i = 0; i < 3; ++i)
    ys.push_back(10 + i);

  assign(xs, index_list(index_uni(4), index_min_max(3, 5)), ys);

  for (int j = 0; j < 3; ++j)
    EXPECT_FLOAT_EQ(ys[j], xs[3][j + 2]);

  stan::test::indexing::test_throw(xs, index_list(index_uni(0), index_min_max(3, 5)), ys);
  stan::test::indexing::test_throw(xs, index_list(index_uni(11), index_min_max(3, 5)), ys);
  stan::test::indexing::test_throw_ia(xs, index_list(index_uni(4), index_min_max(2, 5)), ys);
}

TEST(ModelIndexing, lvalueMultiUni) {
  using stan::math::assign;
  using stan::math::index_list;
  using stan::math::index_min_max;
  using stan::math::index_uni;
  std::vector<std::vector<double> > xs;
  for (int i = 0; i < 10; ++i) {
    std::vector<double> xsi;
    for (int j = 0; j < 20; ++j)
      xsi.push_back(i + j / 10.0);
    xs.push_back(xsi);
  }

  std::vector<double> ys;
  for (int i = 0; i < 3; ++i)
    ys.push_back(10 + i);

  assign(xs, index_list(index_min_max(5, 7), index_uni(8)), ys);

  for (int j = 0; j < 3; ++j)
    EXPECT_FLOAT_EQ(ys[j], xs[j + 4][7]);

  stan::test::indexing::test_throw_ia(xs, index_list(index_min_max(3, 6), index_uni(7)), ys);
  stan::test::indexing::test_throw(xs, index_list(index_min_max(4, 6), index_uni(0)), ys);
  stan::test::indexing::test_throw(xs, index_list(index_min_max(4, 6), index_uni(30)), ys);
}

TEST(ModelIndexing, lvalueVecUni) {
  using stan::math::assign;
  using stan::math::index_list;
  using stan::math::index_uni;
  Eigen::VectorXd xs(5);
  xs << 0, 1, 2, 3, 4;
  double y = 13;
  assign(xs, index_list(index_uni(3)), y);
  EXPECT_FLOAT_EQ(y, xs(2));

  stan::test::indexing::test_throw(xs, index_list(index_uni(0)), y);
  stan::test::indexing::test_throw(xs, index_list(index_uni(6)), y);
}

TEST(ModelIndexing, lvalueRowVecUni) {
  using stan::math::assign;
  using stan::math::index_list;
  using stan::math::index_uni;
  using stan::math::nil_index_list;
  Eigen::RowVectorXd xs(5);
  xs << 0, 1, 2, 3, 4;
  double y = 13;
  assign(xs, index_list(index_uni(3)), y);
  EXPECT_FLOAT_EQ(y, xs(2));
  stan::test::indexing::test_throw(xs, index_list(index_uni(0)), y);
  stan::test::indexing::test_throw(xs, index_list(index_uni(6)), y);
}

TEST(ModelIndexing, lvalueVecMulti) {
  using stan::math::assign;
  using stan::math::index_list;
  using stan::math::index_min;
  using stan::math::index_multi;
  Eigen::VectorXd xs(5);
  xs << 0, 1, 2, 3, 4;
  Eigen::VectorXd ys(3);
  ys << 10, 11, 12;
  assign(xs, index_list(index_min(3)), ys);
  EXPECT_FLOAT_EQ(ys(0), xs(2));
  EXPECT_FLOAT_EQ(ys(1), xs(3));
  EXPECT_FLOAT_EQ(ys(2), xs(4));
  stan::test::indexing::test_throw_ia(xs, index_list(index_min(0)), ys);

  xs << 0, 1, 2, 3, 4;
  std::vector<int> ns;
  ns.push_back(4);
  ns.push_back(1);
  ns.push_back(3);
  assign(xs, index_list(index_multi(ns)), ys);
  EXPECT_FLOAT_EQ(ys(0), xs(3));
  EXPECT_FLOAT_EQ(ys(1), xs(0));
  EXPECT_FLOAT_EQ(ys(2), xs(2));

  ns[ns.size() - 1] = 0;
  stan::test::indexing::test_throw(xs, index_list(index_multi(ns)), ys);

  ns[ns.size() - 1] = 10;
  stan::test::indexing::test_throw(xs, index_list(index_multi(ns)), ys);

  ns[ns.size() - 1] = 3;
  ns.push_back(1);
  stan::test::indexing::test_throw_ia(xs, index_list(index_multi(ns)), ys);
}

TEST(ModelIndexing, lvalueRowVecMulti) {
  using stan::math::assign;
  using stan::math::index_list;
  using stan::math::index_min;
  using stan::math::index_multi;
  Eigen::RowVectorXd xs(5);
  xs << 0, 1, 2, 3, 4;
  Eigen::RowVectorXd ys(3);
  ys << 10, 11, 12;
  assign(xs, index_list(index_min(3)), ys);
  EXPECT_FLOAT_EQ(ys(0), xs(2));
  EXPECT_FLOAT_EQ(ys(1), xs(3));
  EXPECT_FLOAT_EQ(ys(2), xs(4));
  stan::test::indexing::test_throw_ia(xs, index_list(index_min(2)), ys);
  stan::test::indexing::test_throw_ia(xs, index_list(index_min(0)), ys);

  xs << 0, 1, 2, 3, 4;
  std::vector<int> ns;
  ns.push_back(4);
  ns.push_back(1);
  ns.push_back(3);
  assign(xs, index_list(index_multi(ns)), ys);
  EXPECT_FLOAT_EQ(ys(0), xs(3));
  EXPECT_FLOAT_EQ(ys(1), xs(0));
  EXPECT_FLOAT_EQ(ys(2), xs(2));

  ns[ns.size() - 1] = 0;
  stan::test::indexing::test_throw(xs, index_list(index_multi(ns)), ys);

  ns[ns.size() - 1] = 10;
  stan::test::indexing::test_throw(xs, index_list(index_multi(ns)), ys);

  ns[ns.size() - 1] = 3;
  ns.push_back(1);
  stan::test::indexing::test_throw_ia(xs, index_list(index_multi(ns)), ys);
}

TEST(ModelIndexing, lvalueMatrixUni) {
  using stan::math::assign;
  using stan::math::index_list;
  using stan::math::index_uni;
  Eigen::MatrixXd x(3, 4);
  x << 0.0, 0.1, 0.2, 0.3, 1.0, 1.1, 1.2, 1.3, 2.0, 2.1, 2.2, 2.3;

  Eigen::RowVectorXd y(4);
  y << 10.0, 10.1, 10.2, 10.3;

  assign(x, index_list(index_uni(3)), y);
  for (int j = 0; j < 4; ++j)
    EXPECT_FLOAT_EQ(x(2, j), y(j));

  stan::test::indexing::test_throw(x, index_list(index_uni(0)), y);
  stan::test::indexing::test_throw(x, index_list(index_uni(5)), y);
}

TEST(ModelIndexing, lvalueMatrixMulti) {
  using stan::math::assign;
  using stan::math::index_list;
  using stan::math::index_min;
  Eigen::MatrixXd x(3, 4);
  x << 0.0, 0.1, 0.2, 0.3, 1.0, 1.1, 1.2, 1.3, 2.0, 2.1, 2.2, 2.3;

  Eigen::MatrixXd y(2, 4);
  y << 10.0, 10.1, 10.2, 10.3, 11.0, 11.1, 11.2, 11.3;

  assign(x, index_list(index_min(2)), y);
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 4; ++j)
      EXPECT_FLOAT_EQ(y(i, j), x(i + 1, j));
  stan::test::indexing::test_throw_ia(x, index_list(index_min(1)), y);

  Eigen::MatrixXd z(1, 2);
  z << 10, 20;
  stan::test::indexing::test_throw_ia(x, index_list(index_min(1)), z);
  stan::test::indexing::test_throw_ia(x, index_list(index_min(2)), z);
}

TEST(ModelIndexing, lvalueMatrixUniUni) {
  using stan::math::assign;
  using stan::math::index_list;
  using stan::math::index_uni;
  Eigen::MatrixXd x(3, 4);
  x << 0.0, 0.1, 0.2, 0.3, 1.0, 1.1, 1.2, 1.3, 2.0, 2.1, 2.2, 2.3;

  double y = 10.12;
  assign(x, index_list(index_uni(2), index_uni(3)), y);
  EXPECT_FLOAT_EQ(y, x(1, 2));

  stan::test::indexing::test_throw(x, index_list(index_uni(0), index_uni(3)), y);
  stan::test::indexing::test_throw(x, index_list(index_uni(2), index_uni(0)), y);
  stan::test::indexing::test_throw(x, index_list(index_uni(4), index_uni(3)), y);
  stan::test::indexing::test_throw(x, index_list(index_uni(2), index_uni(5)), y);
}

TEST(ModelIndexing, lvalueMatrixUniMulti) {
  using stan::math::assign;
  using stan::math::index_list;
  using stan::math::index_multi;
  using stan::math::index_min_max;
  using stan::math::index_uni;
  Eigen::MatrixXd x(3, 4);
  x << 0.0, 0.1, 0.2, 0.3, 1.0, 1.1, 1.2, 1.3, 2.0, 2.1, 2.2, 2.3;

  Eigen::RowVectorXd y(3);
  y << 10, 11, 12;
  puts("1");
  assign(x, index_list(index_uni(2), index_min_max(2, 4)), y);
  EXPECT_FLOAT_EQ(y(0), x(1, 1));
  EXPECT_FLOAT_EQ(y(1), x(1, 2));
  EXPECT_FLOAT_EQ(y(2), x(1, 3));

  puts("2");
  stan::test::indexing::test_throw(x, index_list(index_uni(0), index_min_max(2, 4)), y);
  puts("3");
  stan::test::indexing::test_throw(x, index_list(index_uni(5), index_min_max(2, 4)), y);
  puts("4");
  stan::test::indexing::test_throw(x, index_list(index_uni(2), index_min_max(0, 2)), y);
  puts("5");
  stan::test::indexing::test_throw_ia(x, index_list(index_uni(2), index_min_max(2, 5)), y);

  std::vector<int> ns;
  ns.push_back(4);
  ns.push_back(1);
  ns.push_back(3);
  puts("6");
  std::cout << "\nx: \n" << x << "\n";
  assign(x, index_list(index_uni(3), index_multi(ns)), y);
  std::cout << "\nnewx: \n" << x << "\n";
  std::cout << "\ny: \n" << y << "\n";
  EXPECT_FLOAT_EQ(y(0), x(2, 3));
  EXPECT_FLOAT_EQ(y(1), x(2, 0));
  EXPECT_FLOAT_EQ(y(2), x(2, 2));

  ns[ns.size() - 1] = 0;
  stan::test::indexing::test_throw(x, index_list(index_uni(3), index_multi(ns)), y);

  ns[ns.size() - 1] = 20;
  stan::test::indexing::test_throw(x, index_list(index_uni(3), index_multi(ns)), y);

  ns.push_back(2);
  stan::test::indexing::test_throw_ia(x, index_list(index_uni(3), index_multi(ns)), y);
}

TEST(ModelIndexing, lvalueMatrixMultiUni) {
  using stan::math::assign;
  using stan::math::index_list;
  using stan::math::index_min_max;
  using stan::math::index_uni;
  using stan::math::index_multi;
  Eigen::MatrixXd x(3, 4);
  x << 0.0, 0.1, 0.2, 0.3, 1.0, 1.1, 1.2, 1.3, 2.0, 2.1, 2.2, 2.3;

  Eigen::VectorXd y(2);
  y << 10, 11;
  puts("get1");
  assign(x, index_list(index_min_max(2, 3), index_uni(4)), y);
  EXPECT_FLOAT_EQ(y(0), x(1, 3));
  EXPECT_FLOAT_EQ(y(1), x(2, 3));

  stan::test::indexing::test_throw(x, index_list(index_min_max(2, 3), index_uni(0)), y);
  stan::test::indexing::test_throw(x, index_list(index_min_max(2, 3), index_uni(5)), y);
  stan::test::indexing::test_throw(x, index_list(index_min_max(0, 1), index_uni(4)), y);
  stan::test::indexing::test_throw_ia(x, index_list(index_min_max(1, 3), index_uni(4)), y);

  std::vector<int> ns;
  ns.push_back(3);
  ns.push_back(1);
  puts("get2");
  assign(x, index_list(index_multi(ns), index_uni(3)), y);
  EXPECT_FLOAT_EQ(y(0), x(2, 2));
  EXPECT_FLOAT_EQ(y(1), x(0, 2));

  ns[ns.size() - 1] = 0;
  stan::test::indexing::test_throw(x, index_list(index_multi(ns), index_uni(3)), y);

  ns[ns.size() - 1] = 20;
  stan::test::indexing::test_throw(x, index_list(index_multi(ns), index_uni(3)), y);

  ns.push_back(2);
  stan::test::indexing::test_throw_ia(x, index_list(index_multi(ns), index_uni(3)), y);
}

TEST(ModelIndexing, lvalueMatrixMultiMulti) {
  using stan::math::assign;
  using stan::math::index_list;
  using stan::math::index_min_max;
  using stan::math::index_min;
  using stan::math::index_multi;
  Eigen::MatrixXd x(3, 4);
  x << 0.0, 0.1, 0.2, 0.3, 1.0, 1.1, 1.2, 1.3, 2.0, 2.1, 2.2, 2.3;

  Eigen::MatrixXd y(2, 3);
  y << 10, 11, 12, 20, 21, 22;
  assign(x, index_list(index_min_max(2, 3), index_min(2)), y);
  EXPECT_FLOAT_EQ(y(0, 0), x(1, 1));
  EXPECT_FLOAT_EQ(y(0, 1), x(1, 2));
  EXPECT_FLOAT_EQ(y(0, 2), x(1, 3));
  EXPECT_FLOAT_EQ(y(1, 0), x(2, 1));
  EXPECT_FLOAT_EQ(y(1, 1), x(2, 2));
  EXPECT_FLOAT_EQ(y(1, 2), x(2, 3));

  stan::test::indexing::test_throw_ia(x, index_list(index_min_max(2, 3), index_min(0)), y);
  stan::test::indexing::test_throw_ia(x, index_list(index_min_max(2, 3), index_min(10)), y);
  stan::test::indexing::test_throw_ia(x, index_list(index_min_max(1, 3), index_min(2)), y);

  x << 0.0, 0.1, 0.2, 0.3, 1.0, 1.1, 1.2, 1.3, 2.0, 2.1, 2.2, 2.3;
  std::vector<int> ms;
  ms.push_back(3);
  ms.push_back(1);

  std::vector<int> ns;
  ns.push_back(2);
  ns.push_back(3);
  ns.push_back(1);
  assign(x, index_list(index_multi(ms), index_multi(ns)), y);
  EXPECT_FLOAT_EQ(y(0, 0), x(2, 1));
  EXPECT_FLOAT_EQ(y(0, 1), x(2, 2));
  EXPECT_FLOAT_EQ(y(0, 2), x(2, 0));
  EXPECT_FLOAT_EQ(y(1, 0), x(0, 1));
  EXPECT_FLOAT_EQ(y(1, 1), x(0, 2));
  EXPECT_FLOAT_EQ(y(1, 2), x(0, 0));

  ms[ms.size() - 1] = 0;
  stan::test::indexing::test_throw(x, index_list(index_multi(ms), index_multi(ns)), y);

  ms[ms.size() - 1] = 10;
  stan::test::indexing::test_throw(x, index_list(index_multi(ms), index_multi(ns)), y);

  ms[ms.size() - 1] = 1;  // back to original valid value
  ns[ns.size() - 1] = 0;
  stan::test::indexing::test_throw(x, index_list(index_multi(ms), index_multi(ns)), y);

  ns[ns.size() - 1] = 10;
  stan::test::indexing::test_throw(x, index_list(index_multi(ms), index_multi(ns)), y);
}

TEST(ModelIndexing, doubleToVar) {
  using stan::math::assign;
  using stan::math::index_list;
  using stan::math::index_uni;
  using stan::math::index_multi;
  using stan::math::index_omni;
  using stan::math::nil_index_list;
  using stan::math::var;
  using Eigen::Dynamic;
  std::vector<double> xs;
  xs.push_back(1);
  xs.push_back(2);
  xs.push_back(3);
  std::vector<std::vector<double> > xss;
  xss.push_back(xs);

  std::vector<var> ys(3);
  std::vector<std::vector<var> > yss;
  yss.push_back(ys);

  assign(yss, index_list(index_omni()), xss, "foo");

  // test both cases where matrix indexed by rows
  // case 1: double matrix with single multi-index on LHS, var matrix on RHS
  Eigen::Matrix<var, Dynamic, Dynamic> a(4, 3);
  for (int i = 0; i < 12; ++i)
    a(i) = -(i + 1);

  Eigen::Matrix<double, Dynamic, Dynamic> b(2, 3);
  b << 1, 2, 3, 4, 5, 6;

  std::vector<int> is;
  is.push_back(2);
  is.push_back(3);
  assign(a, index_list(index_multi(is)), b);
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 3; ++j)
      EXPECT_FLOAT_EQ(a(i + 1, j).val(), b(i, j));

  // case 2: double matrix with single multi-index on LHS, row vector
  // on RHS
  Eigen::Matrix<var, Dynamic, Dynamic> c(4, 3);
  for (int i = 0; i < 12; ++i)
    c(i) = -(i + 1);
  Eigen::Matrix<double, 1, Dynamic> d(3);
  d << 100, 101, 102;
  assign(c, index_list(index_uni(2)), d);
  for (int j = 0; j < 3; ++j)
    EXPECT_FLOAT_EQ(c(1, j).val(), d(j));
}
TEST(ModelIndexing, resultSizeNegIndexing) {
  using stan::math::assign;
  using stan::math::index_list;
  using stan::math::index_min_max;

  std::vector<double> rhs;
  rhs.push_back(2);
  rhs.push_back(5);
  rhs.push_back(-125);

  std::vector<double> lhs;
  assign(rhs, index_list(index_min_max(1, 0)), lhs);
  EXPECT_EQ(0, lhs.size());
}

TEST(ModelIndexing, resultSizeIndexingEigen) {
  using stan::math::assign;
  using stan::math::index_list;
  using stan::math::index_min_max;
  Eigen::VectorXd lhs(5);
  lhs << 1, 2, 3, 4, 5;
  Eigen::VectorXd rhs(4);
  rhs << 4, 3, 2, 1;
  assign(lhs, index_list(index_min_max(1, 4)), rhs);
  EXPECT_FLOAT_EQ(lhs(0), 4);
  EXPECT_FLOAT_EQ(lhs(1), 3);
  EXPECT_FLOAT_EQ(lhs(2), 2);
  EXPECT_FLOAT_EQ(lhs(3), 1);
  EXPECT_FLOAT_EQ(lhs(4), 5);
}

TEST(ModelIndexing, resultSizeNegIndexingEigen) {
  using stan::math::assign;
  using stan::math::index_list;
  using stan::math::index_min_max;
  Eigen::VectorXd lhs(5);
  lhs << 1, 2, 3, 4, 5;
  Eigen::VectorXd rhs(4);
  rhs << 1, 2, 3, 4;
  assign(lhs, index_list(index_min_max(4, 1)), rhs);
  EXPECT_FLOAT_EQ(lhs(0), 4);
  EXPECT_FLOAT_EQ(lhs(1), 3);
  EXPECT_FLOAT_EQ(lhs(2), 2);
  EXPECT_FLOAT_EQ(lhs(3), 1);
  EXPECT_FLOAT_EQ(lhs(4), 5);
}

TEST(ModelIndexing, resultSizePosMinMaxPosMinMaxEigenMatrix) {
  using stan::math::assign;
  using stan::math::cons_list;
  using stan::math::index_list;
  using stan::math::index_min_max;
  using std::vector;
  Eigen::Matrix<double, -1, -1> x(5, 5);
  Eigen::Matrix<double, -1, -1> x_rev(5, 5);
  for (int i = 0; i < x.size(); ++i) {
    x(i) = i;
    x_rev(i) = x.size() - i - 1;
  }

  for (int i = 0; i < x.rows(); ++i) {
    Eigen::MatrixXd x_colwise_rev = x_rev.block(0, 0, i + 1, i + 1);
    assign(x, index_list(index_min_max(1, i + 1), index_min_max(1, i + 1)),
           x_rev.block(0, 0, i + 1, i + 1));
    for (int kk = 0; kk < i; ++kk) {
      for (int jj = 0; jj < i; ++jj) {
        EXPECT_FLOAT_EQ(x(kk, jj), x_rev(kk, jj));
      }
    }
    for (int j = 0; j < x.size(); ++j) {
      x(j) = j;
    }
  }
}

TEST(ModelIndexing, resultSizePosMinMaxNegMinMaxEigenMatrix) {
  using stan::math::assign;
  using stan::math::cons_list;
  using stan::math::index_list;
  using stan::math::index_min_max;
  Eigen::Matrix<double, -1, -1> x(5, 5);
  Eigen::Matrix<double, -1, -1> x_rev(5, 5);
  for (int i = 0; i < x.size(); ++i) {
    x(i) = i;
    x_rev(i) = x.size() - i - 1;
  }

  for (int i = 0; i < x.rows(); ++i) {
    Eigen::MatrixXd x_rowwise_reverse
        = x_rev.block(0, 0, i + 1, i + 1).rowwise().reverse();
    assign(x, index_list(index_min_max(1, i + 1), index_min_max(i + 1, 1)),
           x_rev.block(0, 0, i + 1, i + 1));
    for (int kk = 0; kk < i; ++kk) {
      for (int jj = 0; jj < i; ++jj) {
        EXPECT_FLOAT_EQ(x(kk, jj), x_rowwise_reverse(kk, jj));
      }
    }
    for (int j = 0; j < x.size(); ++j) {
      x(j) = j;
    }
  }
}

TEST(ModelIndexing, resultSizeNigMinMaxPosMinMaxEigenMatrix) {
  using stan::math::assign;
  using stan::math::cons_list;
  using stan::math::index_list;
  using stan::math::index_min_max;
  Eigen::Matrix<double, -1, -1> x(5, 5);
  Eigen::Matrix<double, -1, -1> x_rev(5, 5);
  for (int i = 0; i < x.size(); ++i) {
    x(i) = i;
    x_rev(i) = x.size() - i - 1;
  }

  for (int i = 0; i < x.rows(); ++i) {
    Eigen::MatrixXd x_colwise_reverse
        = x_rev.block(0, 0, i + 1, i + 1).colwise().reverse();
    assign(x, index_list(index_min_max(i + 1, 1), index_min_max(1, i + 1)),
           x_rev.block(0, 0, i + 1, i + 1));
    for (int kk = 0; kk < i; ++kk) {
      for (int jj = 0; jj < i; ++jj) {
        EXPECT_FLOAT_EQ(x(kk, jj), x_colwise_reverse(kk, jj));
      }
    }
    for (int j = 0; j < x.size(); ++j) {
      x(j) = j;
    }
  }
}

TEST(ModelIndexing, resultSizeNegMinMaxNegMinMaxEigenMatrix) {
  using stan::math::assign;
  using stan::math::cons_list;
  using stan::math::index_list;
  using stan::math::index_min_max;
  Eigen::Matrix<double, -1, -1> x(5, 5);
  Eigen::Matrix<double, -1, -1> x_rev(5, 5);
  for (int i = 0; i < x.size(); ++i) {
    x(i) = i;
    x_rev(i) = x.size() - i - 1;
  }

  for (int i = 0; i < x.rows(); ++i) {
    Eigen::MatrixXd x_reverse = x_rev.block(0, 0, i + 1, i + 1).reverse();
    assign(x, index_list(index_min_max(i + 1, 1), index_min_max(i + 1, 1)),
           x_rev.block(0, 0, i + 1, i + 1));
    for (int kk = 0; kk < i; ++kk) {
      for (int jj = 0; jj < i; ++jj) {
        EXPECT_FLOAT_EQ(x(kk, jj), x_reverse(kk, jj));
      }
    }
    for (int j = 0; j < x.size(); ++j) {
      x(j) = j;
    }
  }
}

TEST(modelIndexing, doubleToVarSimple) {
  using stan::math::assign;
  using stan::math::var;
  using stan::math::nil_index_list;

  Eigen::MatrixXd a(2, 2);
  a << 1, 2, 3, 4;
  Eigen::Matrix<var, -1, -1> b;
  assign(b, nil_index_list(), a);
  for (int i = 0; i < a.size(); ++i)
    EXPECT_FLOAT_EQ(a(i), b(i).val());
}

TEST(model_indexing, assign_eigvec_eigvec_index_min) {
  using stan::math::assign;
  using stan::math::index_list;
  using stan::math::index_min;
  Eigen::VectorXd lhs_x(5);
  lhs_x << 0, 1, 2, 3, 4;
  Eigen::VectorXd rhs_y(3);
  rhs_y << 10, 11, 12;
  assign(lhs_x, index_list(index_min(3)), rhs_y);
  EXPECT_FLOAT_EQ(rhs_y(0), lhs_x(2));
  EXPECT_FLOAT_EQ(rhs_y(1), lhs_x(3));
  EXPECT_FLOAT_EQ(rhs_y(2), lhs_x(4));
  stan::test::indexing::test_throw_ia(lhs_x, index_list(index_min(0)), rhs_y);

  assign(lhs_x, index_list(index_min(3)), rhs_y.array() + 1.0);
  EXPECT_FLOAT_EQ(rhs_y(0) + 1.0, lhs_x(2));
  EXPECT_FLOAT_EQ(rhs_y(1) + 1.0, lhs_x(3));
  EXPECT_FLOAT_EQ(rhs_y(2) + 1.0, lhs_x(4));
}

TEST(model_indexing, assign_eigvec_eigvec_index_multi) {
  using stan::math::assign;
  using stan::math::index_list;
  using stan::math::index_multi;
  Eigen::VectorXd lhs_x(5);
  lhs_x << 0, 1, 2, 3, 4;
  Eigen::VectorXd rhs_y(3);
  rhs_y << 10, 11, 12;

  std::vector<int> ns;
  ns.push_back(4);
  ns.push_back(1);
  ns.push_back(3);
  assign(lhs_x, index_list(index_multi(ns)), rhs_y);
  EXPECT_FLOAT_EQ(rhs_y(0), lhs_x(3));
  EXPECT_FLOAT_EQ(rhs_y(1), lhs_x(0));
  EXPECT_FLOAT_EQ(rhs_y(2), lhs_x(2));

  assign(lhs_x, index_list(index_multi(ns)), rhs_y.array() + 4);
  EXPECT_FLOAT_EQ(rhs_y(0) + 4, lhs_x(3));
  EXPECT_FLOAT_EQ(rhs_y(1) + 4, lhs_x(0));
  EXPECT_FLOAT_EQ(rhs_y(2) + 4, lhs_x(2));

  ns[ns.size() - 1] = 0;
  stan::test::indexing::test_throw(lhs_x, index_list(index_multi(ns)), rhs_y);

  ns[ns.size() - 1] = 10;
  stan::test::indexing::test_throw(lhs_x, index_list(index_multi(ns)), rhs_y);

  ns[ns.size() - 1] = 3;
  ns.push_back(1);
  stan::test::indexing::test_throw_ia(lhs_x, index_list(index_multi(ns)), rhs_y);
}

TEST(model_indexing, assign_eigrowvec_eigrowvec_index_min) {
  using stan::math::assign;
  using stan::math::index_list;
  using stan::math::index_min;
  Eigen::RowVectorXd lhs_x(5);
  lhs_x << 0, 1, 2, 3, 4;
  Eigen::RowVectorXd rhs_y(3);
  rhs_y << 10, 11, 12;
  assign(lhs_x, index_list(index_min(3)), rhs_y);
  EXPECT_FLOAT_EQ(rhs_y(0), lhs_x(2));
  EXPECT_FLOAT_EQ(rhs_y(1), lhs_x(3));
  EXPECT_FLOAT_EQ(rhs_y(2), lhs_x(4));
  stan::test::indexing::test_throw_ia(lhs_x, index_list(index_min(0)), rhs_y);

  assign(lhs_x, index_list(index_min(3)), rhs_y.array() + 1.0);
  EXPECT_FLOAT_EQ(rhs_y(0) + 1.0, lhs_x(2));
  EXPECT_FLOAT_EQ(rhs_y(1) + 1.0, lhs_x(3));
  EXPECT_FLOAT_EQ(rhs_y(2) + 1.0, lhs_x(4));
}

TEST(model_indexing, assign_eigrowvec_eigrowvec_index_multi) {
  using stan::math::assign;
  using stan::math::index_list;
  using stan::math::index_multi;
  Eigen::RowVectorXd lhs_x(5);
  lhs_x << 0, 1, 2, 3, 4;
  Eigen::RowVectorXd rhs_y(3);
  rhs_y << 10, 11, 12;

  std::vector<int> ns;
  ns.push_back(4);
  ns.push_back(1);
  ns.push_back(3);
  assign(lhs_x, index_list(index_multi(ns)), rhs_y);
  EXPECT_FLOAT_EQ(rhs_y(0), lhs_x(3));
  EXPECT_FLOAT_EQ(rhs_y(1), lhs_x(0));
  EXPECT_FLOAT_EQ(rhs_y(2), lhs_x(2));

  assign(lhs_x, index_list(index_multi(ns)), rhs_y.array() + 4);
  EXPECT_FLOAT_EQ(rhs_y(0) + 4, lhs_x(3));
  EXPECT_FLOAT_EQ(rhs_y(1) + 4, lhs_x(0));
  EXPECT_FLOAT_EQ(rhs_y(2) + 4, lhs_x(2));

  ns[ns.size() - 1] = 0;
  stan::test::indexing::test_throw(lhs_x, index_list(index_multi(ns)), rhs_y);

  ns[ns.size() - 1] = 10;
  stan::test::indexing::test_throw(lhs_x, index_list(index_multi(ns)), rhs_y);

  ns[ns.size() - 1] = 3;
  ns.push_back(1);
  stan::test::indexing::test_throw_ia(lhs_x, index_list(index_multi(ns)), rhs_y);
}

TEST(model_indexing, assign_densemat_rowvec_uni_index) {
  using stan::math::assign;
  using stan::math::index_list;
  using stan::math::index_uni;
  Eigen::MatrixXd x(3, 4);
  x << 0.0, 0.1, 0.2, 0.3, 1.0, 1.1, 1.2, 1.3, 2.0, 2.1, 2.2, 2.3;

  Eigen::RowVectorXd y(4);
  y << 10.0, 10.1, 10.2, 10.3;

  assign(x, index_list(index_uni(3)), y.array() + 3);
  for (int j = 0; j < 4; ++j)
    EXPECT_FLOAT_EQ(x(2, j), y(j) + 3);

  stan::test::indexing::test_throw(x, index_list(index_uni(0)), y);
  stan::test::indexing::test_throw(x, index_list(index_uni(5)), y);
}

TEST(model_indexing, assign_densemat_densemat_index_min) {
  using stan::math::assign;
  using stan::math::index_list;
  using stan::math::index_min;
  Eigen::MatrixXd x(3, 4);
  x << 0.0, 0.1, 0.2, 0.3, 1.0, 1.1, 1.2, 1.3, 2.0, 2.1, 2.2, 2.3;

  Eigen::MatrixXd y(2, 4);
  y << 10.0, 10.1, 10.2, 10.3, 11.0, 11.1, 11.2, 11.3;

  assign(x, index_list(index_min(2)), y);
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 4; ++j) {
      EXPECT_FLOAT_EQ(y(i, j), x(i + 1, j));
    }
  }
  assign(x, index_list(index_min(2)), y.transpose().transpose());
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 4; ++j) {
      EXPECT_FLOAT_EQ(y(i, j), x(i + 1, j));
    }
  }
  stan::test::indexing::test_throw_ia(x, index_list(index_min(1)), y);

  Eigen::MatrixXd z(1, 2);
  z << 10, 20;
  stan::test::indexing::test_throw_ia(x, index_list(index_min(1)), z);
  stan::test::indexing::test_throw_ia(x, index_list(index_min(2)), z);
}

TEST(model_indexing, assign_densemat_scalar_index_uni) {
  using stan::math::assign;
  using stan::math::index_list;
  using stan::math::index_uni;
  Eigen::MatrixXd x(3, 4);
  x << 0.0, 0.1, 0.2, 0.3, 1.0, 1.1, 1.2, 1.3, 2.0, 2.1, 2.2, 2.3;

  double y = 10.12;
  assign(x, index_list(index_uni(2), index_uni(3)), y);
  EXPECT_FLOAT_EQ(y, x(1, 2));

  stan::test::indexing::test_throw(x, index_list(index_uni(0), index_uni(3)), y);
  stan::test::indexing::test_throw(x, index_list(index_uni(2), index_uni(0)), y);
  stan::test::indexing::test_throw(x, index_list(index_uni(4), index_uni(3)), y);
  stan::test::indexing::test_throw(x, index_list(index_uni(2), index_uni(5)), y);
}

TEST(model_indexing, assign_densemat_eigrowvec_uni_index_min_max_index) {
  using stan::math::assign;
  using stan::math::index_list;
  using stan::math::index_uni;
  using stan::math::index_min_max;
  using stan::math::index_multi;
  Eigen::MatrixXd x(3, 4);
  x << 0.0, 0.1, 0.2, 0.3, 1.0, 1.1, 1.2, 1.3, 2.0, 2.1, 2.2, 2.3;

  Eigen::RowVectorXd y(3);
  y << 10, 11, 12;
  assign(x, index_list(index_uni(2), index_min_max(2, 4)), y);
  EXPECT_FLOAT_EQ(y(0), x(1, 1));
  EXPECT_FLOAT_EQ(y(1), x(1, 2));
  EXPECT_FLOAT_EQ(y(2), x(1, 3));

  assign(x, index_list(index_uni(2), index_min_max(2, 4)), y.array() + 2);
  EXPECT_FLOAT_EQ(y(0) + 2, x(1, 1));
  EXPECT_FLOAT_EQ(y(1) + 2, x(1, 2));
  EXPECT_FLOAT_EQ(y(2) + 2, x(1, 3));

  stan::test::indexing::test_throw(x, index_list(index_uni(0), index_min_max(2, 4)), y);
  stan::test::indexing::test_throw(x, index_list(index_uni(5), index_min_max(2, 4)), y);
  stan::test::indexing::test_throw(x, index_list(index_uni(2), index_min_max(0, 2)), y);
  stan::test::indexing::test_throw_ia(x, index_list(index_uni(2), index_min_max(2, 5)), y);

  std::vector<int> ns;
  ns.push_back(4);
  ns.push_back(1);
  ns.push_back(3);
  assign(x, index_list(index_uni(3), index_multi(ns)), y);
  EXPECT_FLOAT_EQ(y(0), x(2, 3));
  EXPECT_FLOAT_EQ(y(1), x(2, 0));
  EXPECT_FLOAT_EQ(y(2), x(2, 2));

  assign(x, index_list(index_uni(3), index_multi(ns)), y.array() + 2);
  EXPECT_FLOAT_EQ(y(0) + 2, x(2, 3));
  EXPECT_FLOAT_EQ(y(1) + 2, x(2, 0));
  EXPECT_FLOAT_EQ(y(2) + 2, x(2, 2));

  ns[ns.size() - 1] = 0;
  stan::test::indexing::test_throw(x, index_list(index_uni(3), index_multi(ns)), y);

  ns[ns.size() - 1] = 20;
  stan::test::indexing::test_throw(x, index_list(index_uni(3), index_multi(ns)), y);

  ns.push_back(2);
  stan::test::indexing::test_throw_ia(x, index_list(index_uni(3), index_multi(ns)), y);
}

TEST(model_indexing, assign_densemat_eigvec_min_max_index_uni_index) {
  using stan::math::assign;
  using stan::math::index_list;
  using stan::math::index_uni;
  using stan::math::index_min_max;
  using stan::math::index_multi;
  Eigen::MatrixXd x(3, 4);
  x << 0.0, 0.1, 0.2, 0.3, 1.0, 1.1, 1.2, 1.3, 2.0, 2.1, 2.2, 2.3;

  Eigen::VectorXd y(2);
  y << 10, 11;

  assign(x, index_list(index_min_max(2, 3), index_uni(4)), y);
  EXPECT_FLOAT_EQ(y(0), x(1, 3));
  EXPECT_FLOAT_EQ(y(1), x(2, 3));

  assign(x, index_list(index_min_max(2, 3), index_uni(4)), y.array() + 2);
  EXPECT_FLOAT_EQ(y(0) + 2, x(1, 3));
  EXPECT_FLOAT_EQ(y(1) + 2, x(2, 3));

  stan::test::indexing::test_throw(x, index_list(index_min_max(2, 3), index_uni(0)), y);
  stan::test::indexing::test_throw(x, index_list(index_min_max(2, 3), index_uni(5)), y);
  stan::test::indexing::test_throw(x, index_list(index_min_max(0, 1), index_uni(4)), y);
  stan::test::indexing::test_throw_ia(x, index_list(index_min_max(1, 3), index_uni(4)), y);

  std::vector<int> ns;
  ns.push_back(3);
  ns.push_back(1);
  assign(x, index_list(index_multi(ns), index_uni(3)), y);
  EXPECT_FLOAT_EQ(y(0), x(2, 2));
  EXPECT_FLOAT_EQ(y(1), x(0, 2));

  assign(x.block(0, 0, 3, 3), index_list(index_multi(ns), index_uni(3)),
         y.array() + 2);
  EXPECT_FLOAT_EQ(y(0) + 2, x(2, 2));
  EXPECT_FLOAT_EQ(y(1) + 2, x(0, 2));

  ns[ns.size() - 1] = 0;
  stan::test::indexing::test_throw(x, index_list(index_multi(ns), index_uni(3)), y);

  ns[ns.size() - 1] = 20;
  stan::test::indexing::test_throw(x, index_list(index_multi(ns), index_uni(3)), y);

  ns.push_back(2);
  stan::test::indexing::test_throw_ia(x, index_list(index_multi(ns), index_uni(3)), y);
}

TEST(model_indexing, assign_densemat_densemat_min_max_index_min_index) {
  using stan::math::assign;
  using stan::math::index_list;
  using stan::math::index_min;
  using stan::math::index_min_max;
  Eigen::MatrixXd x(3, 4);
  x << 0.0, 0.1, 0.2, 0.3, 1.0, 1.1, 1.2, 1.3, 2.0, 2.1, 2.2, 2.3;

  Eigen::MatrixXd y(2, 3);
  y << 10, 11, 12, 20, 21, 22;

  assign(x, index_list(index_min_max(2, 3), index_min(2)), y);
  EXPECT_FLOAT_EQ(y(0, 0), x(1, 1));
  EXPECT_FLOAT_EQ(y(0, 1), x(1, 2));
  EXPECT_FLOAT_EQ(y(0, 2), x(1, 3));
  EXPECT_FLOAT_EQ(y(1, 0), x(2, 1));
  EXPECT_FLOAT_EQ(y(1, 1), x(2, 2));
  EXPECT_FLOAT_EQ(y(1, 2), x(2, 3));

  assign(x.block(0, 0, 3, 3), index_list(index_min_max(2, 3), index_min(2)),
         y.block(0, 0, 2, 2));
  EXPECT_FLOAT_EQ(y(0, 0), x(1, 1));
  EXPECT_FLOAT_EQ(y(0, 1), x(1, 2));
  EXPECT_FLOAT_EQ(y(0, 2), x(1, 3));
  EXPECT_FLOAT_EQ(y(1, 0), x(2, 1));
  EXPECT_FLOAT_EQ(y(1, 1), x(2, 2));
  EXPECT_FLOAT_EQ(y(1, 2), x(2, 3));

  stan::test::indexing::test_throw_ia(x, index_list(index_min_max(2, 3), index_min(0)), y);
  stan::test::indexing::test_throw_ia(x, index_list(index_min_max(2, 3), index_min(10)), y);
  stan::test::indexing::test_throw_ia(x, index_list(index_min_max(1, 3), index_min(2)), y);
}

TEST(model_indexing, assign_densemat_densemat_multi_index_multi_index) {
  using stan::math::assign;
  using stan::math::index_list;
  using stan::math::index_multi;
  Eigen::MatrixXd x(3, 4);
  x << 0.0, 0.1, 0.2, 0.3, 1.0, 1.1, 1.2, 1.3, 2.0, 2.1, 2.2, 2.3;

  Eigen::MatrixXd y(2, 3);
  y << 10, 11, 12, 20, 21, 22;
  std::vector<int> ms;
  ms.push_back(3);
  ms.push_back(1);

  std::vector<int> ns;
  ns.push_back(2);
  ns.push_back(3);
  ns.push_back(1);
  assign(x, index_list(index_multi(ms), index_multi(ns)), y);
  EXPECT_FLOAT_EQ(y(0, 0), x(2, 1));
  EXPECT_FLOAT_EQ(y(0, 1), x(2, 2));
  EXPECT_FLOAT_EQ(y(0, 2), x(2, 0));
  EXPECT_FLOAT_EQ(y(1, 0), x(0, 1));
  EXPECT_FLOAT_EQ(y(1, 1), x(0, 2));
  EXPECT_FLOAT_EQ(y(1, 2), x(0, 0));

  Eigen::MatrixXd y2 = y.array() + 2;
  assign(x.block(0, 0, 3, 4), index_list(index_multi(ms), index_multi(ns)),
         y.array() + 2);
  EXPECT_FLOAT_EQ(y2(0, 0), x(2, 1));
  EXPECT_FLOAT_EQ(y2(0, 1), x(2, 2));
  EXPECT_FLOAT_EQ(y2(0, 2), x(2, 0));
  EXPECT_FLOAT_EQ(y2(1, 0), x(0, 1));
  EXPECT_FLOAT_EQ(y2(1, 1), x(0, 2));
  EXPECT_FLOAT_EQ(y2(1, 2), x(0, 0));

  ms[ms.size() - 1] = 0;
  stan::test::indexing::test_throw(x, index_list(index_multi(ms), index_multi(ns)), y);

  ms[ms.size() - 1] = 10;
  stan::test::indexing::test_throw(x, index_list(index_multi(ms), index_multi(ns)), y);

  ms[ms.size() - 1] = 1;  // back to original valid value
  ns[ns.size() - 1] = 0;
  stan::test::indexing::test_throw(x, index_list(index_multi(ms), index_multi(ns)), y);

  ns[ns.size() - 1] = 10;
  stan::test::indexing::test_throw(x, index_list(index_multi(ms), index_multi(ns)), y);
}
