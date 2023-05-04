//
// Copyright (c) 2019 Vinnie Falco (vinnie.falco@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

// Test that header file is self-contained.
#include <boost/url/param.hpp>

#include <boost/core/ignore_unused.hpp>
#include <boost/optional.hpp>

#include "test_suite.hpp"

namespace boost {
namespace urls {

struct param_test
{
    BOOST_STATIC_ASSERT(std::is_copy_constructible<param>::value);
    BOOST_STATIC_ASSERT(std::is_copy_assignable<param>::value);
    BOOST_STATIC_ASSERT(std::is_move_constructible<param>::value);
    BOOST_STATIC_ASSERT(std::is_move_assignable<param>::value);

    BOOST_STATIC_ASSERT(std::is_copy_constructible<param_view>::value);
    BOOST_STATIC_ASSERT(std::is_copy_assignable<param_view>::value);
    BOOST_STATIC_ASSERT(std::is_move_constructible<param_view>::value);
    BOOST_STATIC_ASSERT(std::is_move_assignable<param_view>::value);

    BOOST_STATIC_ASSERT(std::is_copy_constructible<param_pct_view>::value);
    BOOST_STATIC_ASSERT(std::is_copy_assignable<param_pct_view>::value);
    BOOST_STATIC_ASSERT(std::is_move_constructible<param_pct_view>::value);
    BOOST_STATIC_ASSERT(std::is_move_assignable<param_pct_view>::value);

    BOOST_STATIC_ASSERT(std::is_constructible<param_view, param>::value);

    // explicit, expensive
    BOOST_STATIC_ASSERT(std::is_constructible<param, param_pct_view>::value);

    // cheap, loses pct-validation
    BOOST_STATIC_ASSERT(std::is_constructible<param_view, param_pct_view>::value);

    // expensive constructions
    BOOST_STATIC_ASSERT(std::is_constructible<param_pct_view, param>::value);
    BOOST_STATIC_ASSERT(std::is_constructible<param_pct_view, param_view>::value);

    void
    testParam()
    {
        auto const check =
            []( param const& qp,
                string_view key,
                string_view value,
                bool has_value)
        {
            BOOST_TEST_EQ(qp.key, key);
            BOOST_TEST_EQ(qp.value, value);
            BOOST_TEST_EQ(qp.has_value, has_value);
        };

        // param()
        {
            param qp;
            check(qp, "", "", false);
        }

        // param(param&&)
        {
            param qp0("key", "value");
            param qp1(std::move(qp0));
            check(qp0, "", "", false);
            check(qp1, "key", "value", true);

            // gcc 4.8 strings use copy-on-write
            //BOOST_TEST_NE(qp0.key.data(), qp1.key.data());
        }

        // param(param const&)
        {
            param qp0("key", "value");
            param qp1(qp0);
            check(qp1, qp0.key, qp0.value, qp0.has_value);

            // gcc 4.8 strings use copy-on-write
            //BOOST_TEST_NE(qp0.key.data(), qp1.key.data());
        }

        // operator=(param&&)
        {
            param qp0("key", "value");
            param qp1("t", "v");
            qp1 = std::move(qp0);
            check(qp0, "", "", false);
            check(qp1, "key", "value", true);
        }

        // operator=(param const&)
        {
            param qp0("key", "value");
            param qp1("t", "v");
            qp1 = qp0;
            check(qp0, "key", "value", true);
            check(qp1, "key", "value", true);
        }

        //----------------------------------------

        // param(string_view, no_value_t)
        {
            param qp("key", no_value);
            check(qp, "key", "", false);
        }

        // param(string_view, string_view)
        {
            {
                param qp("key", "value");
                check(qp, "key", "value", true);
            }
            {
                param qp("key", "");
                check(qp, "key", "", true);
            }
        }

        // param(string_view, optional<string_view>)
        {
            {
                param qp("key", optional<string_view>("value"));
                check(qp, "key", "value", true);
            }
            {
                param qp("key", optional<string_view>(none));
                check(qp, "key", "", false);
            }
        }

        // operator=(param_view)
        {
            param qp;
            qp.key.reserve(100);
            qp.value.reserve(100);
            qp = param_view("key", "value");
            check(qp, "key", "value", true);
            // capacity preserved on assignment
            BOOST_TEST_GE(qp.key.capacity(), 100);
            BOOST_TEST_GE(qp.value.capacity(), 100);
        }

        // operator=(param_pct_view)
        {
            param qp;
            qp.key.reserve(100);
            qp.value.reserve(100);
            qp = param_pct_view("key", "value");
            check(qp, "key", "value", true);
            // capacity preserved on assignment
            BOOST_TEST_GE(qp.key.capacity(), 100);
            BOOST_TEST_GE(qp.value.capacity(), 100);
        }

        // operator->
        {
            param qp("key", "value");
            BOOST_TEST_EQ(qp->key, "key");
        }
    }

    void
    testParamView()
    {
        auto const check =
            []( param_view const& qp,
                string_view key,
                string_view value,
                bool has_value)
        {
            BOOST_TEST_EQ(qp.key, key);
            BOOST_TEST_EQ(qp.value, value);
            BOOST_TEST_EQ(qp.has_value, has_value);
        };

        // param_view()
        {
            param_view qp;
            check(qp, "", "", false);
        }

        // param_view(string_view, no_value_t)
        {
            param_view qp("key", no_value);
            check(qp, "key", "", false);
        }

        // param_view(string_view, string_view)
        {
            {
                param_view qp("key", "value");
                check(qp, "key", "value", true);
            }
            {
                param_view qp("key", "");
                check(qp, "key", "", true);
            }
        }

        // param_view(string_view, optional<string_view>)
        {
            {
                param_view qp("key", optional<string_view>("value"));
                check(qp, "key", "value", true);
            }
            {
                param_view qp("key", optional<string_view>(none));
                check(qp, "key", "", false);
            }
        }

        // param_view(param)
        {
            param qp("key", "value");
            param_view pv(qp);
            check(pv, "key", "value", true);
            BOOST_TEST_EQ(pv.key.data(), qp.key.data());
            BOOST_TEST_EQ(pv.value.data(), qp.value.data());
        }

        // operator param
        {
            {
                param_view pv;
                param qp(pv);
                check(qp, "", "", false);
                BOOST_TEST_NE(qp.key.data(), pv.key.data());
                BOOST_TEST_NE(qp.value.data(), pv.value.data());
            }
            {
                param_view pv("key", no_value);
                param qp(pv);
                check(qp, "key", "", false);
                BOOST_TEST_NE(qp.key.data(), pv.key.data());
                BOOST_TEST_NE(qp.value.data(), pv.value.data());
            }
            {
                param_view pv("key", "");
                param qp(pv);
                check(qp, "key", "", true);
                BOOST_TEST_NE(qp.key.data(), pv.key.data());
                BOOST_TEST_NE(qp.value.data(), pv.value.data());
            }
            {
                param_view pv("key", "value");
                param qp(pv);
                check(qp, "key", "value", true);
                BOOST_TEST_NE(qp.key.data(), pv.key.data());
                BOOST_TEST_NE(qp.value.data(), pv.value.data());
            }
        }

    }

    void
    testParamPctView()
    {
        auto const check =
            []( param_view const& qp,
                pct_string_view key,
                pct_string_view value,
                bool has_value)
        {
            BOOST_TEST_EQ(qp.key, key);
            BOOST_TEST_EQ(qp.value, value);
            BOOST_TEST_EQ(qp.has_value, has_value);
        };

        // param_pct_view()
        {
            param_view qp;
            check(qp, "", "", false);
        }

        // param_pct_view(pct_string_view)
        {
            param_view qp("key", no_value);
            check(qp, "key", "", false);
        }

        // param_pct_view(pct_string_view, pct_string_view)
        {
            {
                param_view qp("key", "value");
                check(qp, "key", "value", true);
            }
            {
                param_view qp("key", "");
                check(qp, "key", "", true);
            }
        }

        // param_pct_view(pct_string_view, optional<pct_string_view>)
        {
            {
                param_view qp("key", optional<string_view>("value"));
                check(qp, "key", "value", true);
            }
            {
                param_view qp("key", optional<string_view>(none));
                check(qp, "key", "", false);
            }
        }

        // operator param
        {
            param_pct_view pv("key", "value");
            param qp(pv);
            check(qp, "key", "value", true);
            BOOST_TEST_NE(qp.key.data(), pv.key.data());
            BOOST_TEST_NE(qp.value.data(), pv.value.data());
        }
    }

    void
    testNatvis()
    {
        param v0 = {};
        param v1 = { "key", no_value };
        param v2 = { "key", "" };
        param v3 = { "key", "value" };
        boost::ignore_unused(v0,  v1,  v2,  v3);

        param_view pv0 = {};
        param_view pv1 = { "key", no_value };
        param_view pv2 = { "key", "" };
        param_view pv3 = { "key", "value" };
        boost::ignore_unused(pv0, pv1, pv2, pv3);

        param_pct_view d0 = {};
        param_pct_view d1 = { "key", no_value };
        param_pct_view d2 = { "key", "" };
        param_pct_view d3 = { "key", "value" };
        boost::ignore_unused(d0,  d1,  d2,  d3);
    }

    void
    run()
    {
        testParam();
        testParamView();
        testParamPctView();
        testNatvis();
    }
};

TEST_SUITE(
    param_test,
    "boost.url.param");

} // urls
} // boost
