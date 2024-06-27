#include <gtest/gtest.h>
#include <tuple>
#include <test/prob/test_fixture_distr.hpp>
#include <test/prob/test_fixture_cdf.hpp>
#include <test/prob/test_fixture_cdf_log.hpp>
#include <test/prob/test_fixture_ccdf_log.hpp>

typedef std::tuple<fvar<fvar<var>>, empty, empty, empty, empty, empty>
    type_ffv_real_0;
typedef std::tuple<std::vector<fvar<fvar<var>>>, empty, empty, empty, empty,
                   empty>
    type_ffv_real_1;
typedef std::tuple<Eigen::Matrix<fvar<fvar<var>>, Eigen::Dynamic, 1>, empty,
                   empty, empty, empty, empty>
    type_ffv_real_2;
