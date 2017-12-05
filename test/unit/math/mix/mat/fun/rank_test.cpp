#include <stan/math/mix/mat.hpp>
#include <gtest/gtest.h>
#include <stdexcept>
#include <vector>

template <typename T>
void test_rank() {
  using stan::math::rank;

  T c(1);
  c[0] = 1.7;
  EXPECT_EQ(0U, rank(c, 1));
  EXPECT_THROW(rank(c, 0), std::out_of_range);
  EXPECT_THROW(rank(c, 2), std::out_of_range);

  T e(2);
  e[0] = 5.9; e[1] = -1.2;
  EXPECT_EQ(1U, rank(e, 1));
  EXPECT_EQ(0U, rank(e, 2));
  EXPECT_THROW(rank(e, 0), std::out_of_range);
  EXPECT_THROW(rank(e, 3), std::out_of_range);

  T g(3);
  g[0] = 5.9; g[1] = -1.2; g[2] = 192.13456;
  EXPECT_EQ(1U, rank(g, 1));
  EXPECT_EQ(0U, rank(g, 2));
  EXPECT_EQ(2U, rank(g, 3));
  EXPECT_THROW(rank(g, 0), std::out_of_range);
  EXPECT_THROW(rank(g, 4), std::out_of_range);

  T z;
  EXPECT_THROW(rank(z, 0), std::out_of_range);
  EXPECT_THROW(rank(z, 1), std::out_of_range);
  EXPECT_THROW(rank(z, 2), std::out_of_range);
}


template <typename T>
void test_rank_int() {
  using stan::math::rank;

  T c(1);
  c[0] = 1;
  EXPECT_EQ(0U, rank(c, 1));
  EXPECT_THROW(rank(c, 0), std::out_of_range);
  EXPECT_THROW(rank(c, 2), std::out_of_range);

  T e(2);
  e[0] = 5; e[1] = -1;
  EXPECT_EQ(1U, rank(e, 1));
  EXPECT_EQ(0U, rank(e, 2));
  EXPECT_THROW(rank(e, 0), std::out_of_range);
  EXPECT_THROW(rank(e, 3), std::out_of_range);

  T g(3);
  g[0] = 5; g[1] = -1; g[2] = 192;
  EXPECT_EQ(1U, rank(g, 1));
  EXPECT_EQ(0U, rank(g, 2));
  EXPECT_EQ(2U, rank(g, 3));
  EXPECT_THROW(rank(g, 0), std::out_of_range);
  EXPECT_THROW(rank(g, 4), std::out_of_range);

  T z;
  EXPECT_THROW(rank(z, 0), std::out_of_range);
  EXPECT_THROW(rank(z, 1), std::out_of_range);
  EXPECT_THROW(rank(z, 2), std::out_of_range);
}

TEST(AgradMixMatrixRank, rank) {
  using stan::math::rank;
  using stan::math::fvar;
  using stan::math::var;


  test_rank< std::vector<fvar<var> > >();
  test_rank< Eigen::Matrix<fvar<var> , Eigen::Dynamic, 1> >();
  test_rank< Eigen::Matrix<fvar<var> , 1, Eigen::Dynamic> >();

  test_rank_int< std::vector<fvar<var> > >();
  test_rank_int< Eigen::Matrix<fvar<var> , Eigen::Dynamic, 1> >();
  test_rank_int< Eigen::Matrix<fvar<var> , 1, Eigen::Dynamic> >();

  test_rank< std::vector<fvar<fvar<var> > > >();
  test_rank< Eigen::Matrix<fvar<fvar<var> > , Eigen::Dynamic, 1> >();
  test_rank< Eigen::Matrix<fvar<fvar<var> > , 1, Eigen::Dynamic> >();

  test_rank_int< std::vector<fvar<fvar<var> > > >();
  test_rank_int< Eigen::Matrix<fvar<fvar<var> > , Eigen::Dynamic, 1> >();
  test_rank_int< Eigen::Matrix<fvar<fvar<var> > , 1, Eigen::Dynamic> >();
}
