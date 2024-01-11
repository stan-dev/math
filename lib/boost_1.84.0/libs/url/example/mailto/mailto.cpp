//
// Copyright (c) 2022 alandefreitas (alandefreitas@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//

//[example_mailto

/*
    This example parses a mailto URL into a new
    view type and prints its components to
    standard output.
*/

#include <boost/url/grammar/ci_string.hpp>
#include <boost/url/grammar/parse.hpp>
#include <boost/url/optional.hpp>
#include <boost/url/rfc/absolute_uri_rule.hpp>
#include <boost/url/url.hpp>
#include <boost/url/url_view.hpp>
#include <algorithm>
#include <iostream>
#include "rfc.hpp"

namespace urls = boost::urls;

// fwd-declaration for mailto_view
struct mailto_rule_t;

/// A new url type for mailto URLs
/**
    This class represents a URI with the mailto
    scheme.

    Unlike a urls::url_view, which only represents
    the general syntax of urls, a mailto_view
    represents a reference to fields that are
    relevant to mailto URLs, while ignoring
    elements of the general syntax
    that are not relevant to the scheme.

    This allows us to use the general syntax
    parsers to create a representation that
    is more appropriate for the specified scheme
    syntax.

    @par Specification
    @li <a href="https://www.rfc-editor.org/rfc/rfc6068"
        >The 'mailto' URI Scheme</a>
    @li <a href="https://www.rfc-editor.org/errata/rfc6068"
        >RFC Errata Report</a>

    @par References
    @li <a href="https://en.wikipedia.org/wiki/Mailto"
        >mailto (Wikipedia)</a>

 */
class mailto_view
{
    urls::url_view u_;

public:
    /// Return the specified email address in the URL
    /**
        A mailto URL might contain multiple email
        addresses separated by commas.

        The first addresses are represented in
        the path. Other addresses are in
        any query parameter whose key is "to".

        @param i Address index

        @return The specified address
     */
    std::string
    address(std::size_t i = 0) const;

    /// @copydoc address()
    urls::pct_string_view
    encoded_address(std::size_t i = 0) const noexcept;

    /// Return number of email addresses in the URL
    std::size_t
    size() const noexcept;

    /// Return the specified cc email address in the URL
    /**
        A mailto URL might contain multiple cc
        email addresses separated by commas.

        Addresses can be represented in any query
        parameter whose key is "cc".

        @param i Address index

        @return The specified cc address
     */
    std::string
    cc(std::size_t i) const;

    /// @copydoc cc()
    urls::pct_string_view
    encoded_cc(std::size_t i) const noexcept;

    /// Return number of "cc" email addresses in the URL
    std::size_t
    size_cc() const noexcept;

    /// Return email message subject
    std::string
    subject() const;

    /// @copydoc subject()
    urls::pct_string_view
    encoded_subject() const noexcept;

    /// Return email message body
    std::string
    body() const;

    /// @copydoc body()
    urls::pct_string_view
    encoded_body() const noexcept;

    friend
    std::ostream&
    operator<<(std::ostream& os, mailto_view m)
    {
        return os << m.u_;
    }

private:
    // Count number of addresses in a string
    static
    std::size_t
    addr_in_str(boost::core::string_view s);

    // Get the ith address from a string
    static
    boost::optional<urls::pct_string_view>
    get_nth_address(boost::core::string_view to, std::size_t &i) noexcept;

    // Get param value or empty otherwise
    urls::pct_string_view
    param_or_empty(urls::pct_string_view k) const noexcept;

    friend mailto_rule_t;
};

/** Rule to match a mailto URL
*/
struct mailto_rule_t
{
    /// Value type returned by the rule
    using value_type = mailto_view;

    /// Parse a sequence of characters into a mailto_view
    boost::system::result< value_type >
    parse( char const*& it, char const* end ) const noexcept;
};

constexpr mailto_rule_t mailto_rule{};

/** Return a parsed mailto URL from a string, or error.

    This is a more convenient user-facing function
    to parse mailto URLs.
*/
boost::system::result< mailto_view >
parse_mailto( boost::core::string_view s ) noexcept
{
    return urls::grammar::parse(s, mailto_rule);
}

int main(int argc, char** argv)
{
    // This example shows how to use custom parsing
    // to process alternate URI schemes, in this
    // case "mailto"
    if (argc != 2) {
        std::cout << argv[0] << "\n";
        std::cout << "mailto <URL>\n"
                     "examples:\n"
                     // Single e-mail address
                     "mailto mailto:someone@example.com\n"
                     // Two e-mail addresses
                     "mailto mailto:someone@example.com,someoneelse@example.com\n"
                     // E-mail headers
                     "mailto mailto:someone@example.com?subject=Our%20meeting&cc=someone_else@example.com&body=Hi%21\n"
                     // E-mail headers only
                     "mailto mailto:?to=&subject=mailto%20example&body=https%3A%2F%2Fen.wikipedia.org%2Fwiki%2FMailto\n"
                     // All fields
                     "mailto mailto:someone@example.com,%73omeoneelse@me.com?to=thirdperson@example.com&subject=Our%20meeting&cc=someone_else@example.com,onemore@ex%61mple.com&body=Hi%21\n";
        return EXIT_FAILURE;
    }

    boost::system::result<mailto_view> r =
        parse_mailto(argv[1]);
    if (!r)
        return EXIT_FAILURE;

    mailto_view m = *r;
    std::cout << "link: " << m << "\n";

    for (std::size_t i = 0; i < m.size(); ++i)
        std::cout <<
            "to[" << i << "]: " <<
            m.address(i) << "\n";

    for (std::size_t i = 0; i < m.size_cc(); ++i)
        std::cout <<
            "cc[" << i << "]: " <<
            m.address(i) << "\n";

    std::cout << "subject: " << m.subject() << "\n";
    std::cout << "body: " << m.body() << "\n";

    return EXIT_SUCCESS;
}

std::string
mailto_view::address(std::size_t i) const
{
    return encoded_address(i).decode();
}

urls::pct_string_view
mailto_view::encoded_address(std::size_t i) const noexcept
{
    // Look for ith email address in the path string
    auto s = get_nth_address(u_.encoded_path(), i);
    if (s)
        return *s;

    // Look for ith email address in one of the "to" headers
    auto ps = u_.encoded_params();
    auto it = ps.find("to", urls::ignore_case);
    while (it != ps.end())
    {
        s = get_nth_address((*it++).value, i);
        if (s)
            return *s;
        it = ps.find(it, "to", urls::ignore_case);
    }
    return {};
}

std::size_t
mailto_view::size() const noexcept
{
    // Count addresses in path
    std::size_t n = addr_in_str(u_.encoded_path());

    // Count addresses in "to" headers
    auto ps = u_.encoded_params();
    auto it = ps.find("to", urls::ignore_case);
    while (it != ps.end())
    {
        n += addr_in_str((*it++).value);
        it = ps.find(it, "to", urls::ignore_case);
    }
    return n;
}

std::string
mailto_view::cc(std::size_t i) const
{
    return encoded_cc(i).decode();
}

urls::pct_string_view
mailto_view::encoded_cc(std::size_t i) const noexcept
{
    // Look for ith email address in one of the "to" headers
    auto ps = u_.encoded_params();
    auto it = ps.find("cc", urls::ignore_case);
    while (it != ps.end())
    {
        auto s = get_nth_address((*it++).value, i);
        if (s)
            return *s;
        it = ps.find(it, "cc", urls::ignore_case);
    }
    return {};
}

std::size_t
mailto_view::size_cc() const noexcept
{
    // Count addresses in "to" headers
    std::size_t n = 0;
    auto ps = u_.encoded_params();
    auto it = ps.find("cc", urls::ignore_case);
    while (it != ps.end())
    {
        n += addr_in_str((*it++).value);
        it = ps.find(it, "cc", urls::ignore_case);
    }
    return n;
}

std::string
mailto_view::subject() const
{
    return encoded_subject().decode();
}

urls::pct_string_view
mailto_view::encoded_subject() const noexcept
{
    return param_or_empty("subject");
}

std::string
mailto_view::mailto_view::body() const
{
    return encoded_body().decode();
}

urls::pct_string_view
mailto_view::encoded_body() const noexcept
{
    return param_or_empty("body");
}

std::size_t
mailto_view::addr_in_str(boost::core::string_view s)
{
    std::size_t n = 0;
    bool empty = true;
    for (char c : s)
    {
        if (c == ',')
        {
            n += !empty;
            empty = true;
        }
        else
        {
            empty = false;
        }
    }
    n += !empty;
    return n;
}

boost::optional<urls::pct_string_view>
mailto_view::get_nth_address(boost::core::string_view to, std::size_t &i) noexcept
{
    auto p = to.find(',');
    while (p != boost::core::string_view::npos)
    {
        if (i == 0)
            return urls::pct_string_view(
                to.substr(0, p));
        --i;
        to.remove_prefix(p + 1);
        p = to.find(',');
    }
    if (!to.empty())
    {
        if (i == 0)
            return urls::pct_string_view(
                to.substr(0, p));
        --i;
    }
    return boost::none;
}

urls::pct_string_view
mailto_view::param_or_empty(urls::pct_string_view k) const noexcept
{
    auto ps = u_.encoded_params();
    auto it = ps.find(k, urls::ignore_case);
    if (it != ps.end())
        return (*it).value;
    return {};
}

auto
mailto_rule_t::parse( char const*& it, char const* end ) const noexcept
    -> boost::system::result< value_type >
{
    // Syntax-based rules
    boost::system::result<urls::url_view> r =
        urls::grammar::parse(it, end, urls::absolute_uri_rule);
    if (!r)
        return r.error();

    // Scheme-based rules
    mailto_view m;
    m.u_ = *r;
    auto valid_header = [](urls::param_pct_view p) {
        return
            urls::grammar::parse(p.key, hfname_rule) &&
            urls::grammar::parse(p.value, hfvalue_rule) &&
            p.has_value &&
            (!urls::grammar::ci_is_equal(p.key, "to") ||
             urls::grammar::parse(p.value, addr_spec_rule));
    };
    auto ps = m.u_.encoded_params();
    if (m.u_.scheme() == "mailto" &&
        !m.u_.has_authority() &&
        urls::grammar::parse(m.u_.encoded_path(), to_rule) &&
        std::all_of(ps.begin(), ps.end(), valid_header))
        return m;
    return urls::grammar::error::invalid;
}

//]

