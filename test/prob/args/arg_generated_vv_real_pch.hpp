#include <gtest/gtest.h>
#include <tuple>
#include <test/prob/test_fixture_distr.hpp>
#include <test/prob/test_fixture_cdf.hpp>
#include <test/prob/test_fixture_cdf_log.hpp>
#include <test/prob/test_fixture_ccdf_log.hpp>

typedef std::tuple<var, empty, empty, empty, empty, empty> type_vv_real_0;
typedef std::tuple<std::vector<var>, empty, empty, empty, empty, empty> type_vv_real_1;
typedef std::tuple<stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>>, empty, empty, empty, empty, empty> type_vv_real_2;

