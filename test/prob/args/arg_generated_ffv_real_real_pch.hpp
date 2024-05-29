#include <gtest/gtest.h>
#include <tuple>
#include <test/prob/test_fixture_distr.hpp>
#include <test/prob/test_fixture_cdf.hpp>
#include <test/prob/test_fixture_cdf_log.hpp>
#include <test/prob/test_fixture_ccdf_log.hpp>

typedef std::tuple<double, fvar<fvar<var> >, empty, empty, empty, empty> type_ffv_real_real_0;
typedef std::tuple<double, std::vector<fvar<fvar<var> >>, empty, empty, empty, empty> type_ffv_real_real_1;
typedef std::tuple<double, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, empty, empty, empty, empty> type_ffv_real_real_2;
typedef std::tuple<std::vector<double>, fvar<fvar<var> >, empty, empty, empty, empty> type_ffv_real_real_3;
typedef std::tuple<std::vector<double>, std::vector<fvar<fvar<var> >>, empty, empty, empty, empty> type_ffv_real_real_4;
typedef std::tuple<std::vector<double>, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, empty, empty, empty, empty> type_ffv_real_real_5;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, fvar<fvar<var> >, empty, empty, empty, empty> type_ffv_real_real_6;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, std::vector<fvar<fvar<var> >>, empty, empty, empty, empty> type_ffv_real_real_7;
typedef std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, empty, empty, empty, empty> type_ffv_real_real_8;
typedef std::tuple<fvar<fvar<var> >, double, empty, empty, empty, empty> type_ffv_real_real_9;
typedef std::tuple<fvar<fvar<var> >, std::vector<double>, empty, empty, empty, empty> type_ffv_real_real_10;
typedef std::tuple<fvar<fvar<var> >, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty, empty, empty> type_ffv_real_real_11;
typedef std::tuple<fvar<fvar<var> >, fvar<fvar<var> >, empty, empty, empty, empty> type_ffv_real_real_12;
typedef std::tuple<fvar<fvar<var> >, std::vector<fvar<fvar<var> >>, empty, empty, empty, empty> type_ffv_real_real_13;
typedef std::tuple<fvar<fvar<var> >, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, empty, empty, empty, empty> type_ffv_real_real_14;
typedef std::tuple<std::vector<fvar<fvar<var> >>, double, empty, empty, empty, empty> type_ffv_real_real_15;
typedef std::tuple<std::vector<fvar<fvar<var> >>, std::vector<double>, empty, empty, empty, empty> type_ffv_real_real_16;
typedef std::tuple<std::vector<fvar<fvar<var> >>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty, empty, empty> type_ffv_real_real_17;
typedef std::tuple<std::vector<fvar<fvar<var> >>, fvar<fvar<var> >, empty, empty, empty, empty> type_ffv_real_real_18;
typedef std::tuple<std::vector<fvar<fvar<var> >>, std::vector<fvar<fvar<var> >>, empty, empty, empty, empty> type_ffv_real_real_19;
typedef std::tuple<std::vector<fvar<fvar<var> >>, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, empty, empty, empty, empty> type_ffv_real_real_20;
typedef std::tuple<Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, double, empty, empty, empty, empty> type_ffv_real_real_21;
typedef std::tuple<Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, std::vector<double>, empty, empty, empty, empty> type_ffv_real_real_22;
typedef std::tuple<Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, Eigen::Matrix<double, Eigen::Dynamic, 1>, empty, empty, empty, empty> type_ffv_real_real_23;
typedef std::tuple<Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, fvar<fvar<var> >, empty, empty, empty, empty> type_ffv_real_real_24;
typedef std::tuple<Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, std::vector<fvar<fvar<var> >>, empty, empty, empty, empty> type_ffv_real_real_25;
typedef std::tuple<Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, Eigen::Matrix<fvar<fvar<var> >, Eigen::Dynamic, 1>, empty, empty, empty, empty> type_ffv_real_real_26;

