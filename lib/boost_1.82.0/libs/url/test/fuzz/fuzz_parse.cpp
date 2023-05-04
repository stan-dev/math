//
// Copyright (c) 2023 alandefreitas (alandefreitas@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//

#include <boost/url/parse.hpp>

using namespace boost::urls;

enum class parser {
    absolute_uri,
    origin_form,
    relative_ref,
    uri,
    uri_reference
};

bool
fuzz_parse(parser p, string_view sv)
{
    result<url_view> r;
    switch (p)
    {
        case parser::absolute_uri:
            r = parse_absolute_uri(sv);
            break;
        case parser::origin_form:
            r = parse_origin_form(sv);
            break;
        case parser::relative_ref:
            r = parse_relative_ref(sv);
            break;
        case parser::uri:
            r = parse_uri(sv);
            break;
        case parser::uri_reference:
            r = parse_uri_reference(sv);
            break;
    }
    return r.has_value();
}

extern "C"
int
LLVMFuzzerTestOneInput(
    const uint8_t* data,
    size_t size)
{
    if (size == 0)
        return -1;
    try
    {
        auto p = static_cast<parser>(data[0] % 5);
        string_view s{reinterpret_cast<
            const char*>(data + 1), size - 1};
        fuzz_parse(p, s);
    }
    catch(...)
    {
    }
    return 0;
}

