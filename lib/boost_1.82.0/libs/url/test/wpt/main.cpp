//
// Copyright (c) 2021 Vinnie Falco (vinnie dot falco at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

#include <boost/json.hpp>
#include <boost/url.hpp>
#include "test_suite.hpp"
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <fstream>
#include <string>

/*
    Runs tests against the web-platform-tests set
    of input vectors. Paths to zero or more files
    must be provided on the command line. The names
    have to match the original names from this
    repo:

    https://github.com/web-platform-tests/wpt/tree/master/url/resources

    or

    https://github.com/web-platform-tests/wpt/tree/982c7addc45086db44c44e5c442c97703409d675/url/resources
*/

namespace json = boost::json;
using namespace boost::urls;

test_suite::debug_stream Log(std::cout);
int fail_ = 0;
int total_ = 0;

json::value
read_json(char const* path)
{
    std::ostringstream oss;
    std::ifstream f(path);
    oss << f.rdbuf();
    std::string s = oss.str();
    return json::parse(s);
}

string_view
filename(char const* path)
{
    string_view s(path);
#ifdef _MSC_VER
    auto n = s.find_last_of('\\');
#else
    auto n = s.find_last_of('/');
#endif
    if(n == string_view::npos)
        return s;
    return s.substr(n + 1);
}

void
do_setters_scheme(json::array const& ja)
{
    for(auto const& jv : ja)
    {
        ++total_;
        auto href = jv.at("href").as_string();
        url u = parse_uri_reference(href).value();
        auto const& ex = jv.at("expected").as_object();
        try
        {
            u.set_scheme(jv.at("new_value").as_string());
        }
        catch(std::exception const& e)
        {
            if(ex.at("href").as_string() != href)
            {
                Log << "caught exception: " << e.what() << std::endl;
                Log << "set_scheme failed: " << href <<
                    ", " << jv.at("new_value") << std::endl;
                ++fail_;
            }
        }
    }
}

void
do_setters_user(json::array const& ja)
{
    for(auto const& jv : ja)
    {
        ++total_;
        auto href = jv.at("href").as_string();
        url u = parse_uri_reference(href).value();
        // VFALCO TODO
    }
}

void
do_setters_tests(json::value const& jv)
{
    for(auto const& v : jv.as_object())
    {
        if(v.key() == "protocol")
            do_setters_scheme(v.value().as_array());
    }
}

int main(int argc, char** argv)
{
    for(int i = 1; i < argc; ++i)
    {
        try
        {
            auto jv = read_json(argv[i]);
            auto s = filename(argv[i]);
            Log << "file: " << s << std::endl;
            if(s == "setters_tests.json")
                do_setters_tests(jv);
        }
        catch(std::exception const& e)
        {
            Log << "caught exception: " << e.what() << std::endl;
        }
    }
    if(fail_ == 0)
        Log << total_ << " total success" << std::endl;
    else
        Log << fail_ << " of " << total_ << " failures" << std::endl;

    if(fail_ > 0)
        return EXIT_FAILURE;
    return EXIT_SUCCESS;
}
