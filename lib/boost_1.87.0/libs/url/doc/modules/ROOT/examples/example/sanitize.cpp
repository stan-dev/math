//
// Copyright (c) 2023 Alan de Freitas (alandefreitas@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

// tag::example_sanitize_url[]

/*
    This example parses a non-strict / invalid URL
    into path components according to its delimiters.
    This pattern can be adapted to the requirements of other
    applications.

    Once the non-strict components are determined, a new URL is
    created and its parts are set with the set_encoded_X
    functions, which will encode any invalid chars accordingly.

    This sort of transformation is useful in applications that are
    extremely loose in what kinds of URLs they accept, such as
    browsers. The sanitized URL can later be used for machine-to-machine
    communication.

    Using non-strict URLs directly is a security concern in
    machine-to-machine communication, is ambiguous, and also
    involve an extra cost for the transformations.

    Different transformations are required by different applications to
    construct a valid URL appropriate for machine-to-machine communication.
    For instance, if an invalid relative reference includes something that
    looks like a host in the first path segment, browsers usually interpret
    that as the host with an implicit "https" scheme. Other applications
    also have other implicit schemes.

    The example also identifies whether the input url is already valid.
    It includes diagnostics that can be used to help the user determine
    if a URL is invalid and why it's invalid.

    Once all transformations are applied, the result is a URL
    appropriate for machine-to-machine communication.
*/

#include <boost/url/url.hpp>
#include <boost/url/parse.hpp>
#include <boost/url/parse_path.hpp>
#include <boost/url/string_view.hpp>
#include <boost/url/grammar/alpha_chars.hpp>
#include <boost/url/grammar/charset.hpp>
#include <boost/url/grammar/digit_chars.hpp>
#include <boost/url/grammar/lut_chars.hpp>
#include <algorithm>
#include <iostream>
#include <array>

namespace urls = boost::urls;
namespace core = boost::core;

struct url_components
{
    core::string_view scheme;
    core::string_view user;
    core::string_view password;
    core::string_view hostname;
    core::string_view port;
    core::string_view path;
    core::string_view query;
    core::string_view fragment;
};

core::string_view
port_of_scheme(core::string_view scheme_str) {
    static std::array<std::pair<core::string_view, core::string_view>, 21> scheme_ports =
    {{
         {"http", "80"},
         {"ftp", "21"},
         {"https", "443"},
         {"gopher", "70"},
         {"ldap", "389"},
         {"nntp", "119"},
         {"snews", "563"},
         {"imap", "143"},
         {"pop", "110"},
         {"sip", "5060"},
         {"rtsp", "554"},
         {"wais", "210"},
         {"z39.50r", "210"},
         {"z39.50s", "210"},
         {"prospero", "191"},
         {"nfs", "2049"},
         {"tip", "3372"},
         {"acap", "674"},
         {"telnet", "23"},
         {"ssh", "22"},
         {"", "65535"}
    }};

    auto iequals = [](core::string_view a, core::string_view b)
    {
        if (b.size() != a.size()) {
            return false;
        }
        for (unsigned int i = 0; i < a.size(); ++i) {
            if (std::tolower(a[i]) != std::tolower(b[i])) {
                return false;
            }
        }
        return true;
    };

    auto const& it = std::find_if(
        scheme_ports.begin(),
        scheme_ports.end(),
        [&](std::pair<core::string_view, core::string_view> const& s) {
        return iequals(s.first, scheme_str);
    });

    if (it != scheme_ports.end()) {
        return it->second;
    } else {
        return {};
    }
}

void
extract_relative_ref(
    core::string_view hostinfo_relative,
    url_components &out)
{
    // split path and query#fragment
    constexpr urls::grammar::lut_chars path_end_chars("?#\0");
    auto it = urls::grammar::find_if(
        hostinfo_relative.begin(),
        hostinfo_relative.end(), path_end_chars);
    core::string_view query_and_frag = hostinfo_relative.substr(it - hostinfo_relative.begin());
    if (query_and_frag != hostinfo_relative)
        out.path = hostinfo_relative.substr(
            0, query_and_frag.data() - hostinfo_relative.data());
    if (query_and_frag.empty())
        return;

    // ?query#fragment
    if (query_and_frag.front() == '?') {
        query_and_frag = query_and_frag.substr(1);
        core::string_view::size_type hash_pos = query_and_frag.find('#');
        if (hash_pos != core::string_view::npos) {
            core::string_view fragment_part = query_and_frag.substr(hash_pos);
            out.fragment = fragment_part.substr(1);
            out.query = query_and_frag.substr(
                0, fragment_part.data() - query_and_frag.data());
        } else {
            out.query = query_and_frag;
        }
        return;
    }

    // fragment
    out.fragment = query_and_frag.substr(1);
}

void
extract_userinfo_relative(
    core::string_view relative_ref,
    core::string_view userinfo_relative,
    core::string_view host_info,
    url_components& out) {
    // We expect userinfo_relative to point to the first character of
    // the hostname.  If there's a port it is the first colon,
    // except with IPv6.
    auto host_end_pos = host_info.find(':');
    if (host_end_pos == core::string_view::npos)
    {
        // definitely no port
        out.hostname = userinfo_relative.substr(
            0, relative_ref.data() - userinfo_relative.data());
        return extract_relative_ref(relative_ref, out);
    }

    // extract hostname and port
    out.hostname = userinfo_relative.substr(0, host_end_pos);
    core::string_view host_relative = userinfo_relative.substr(host_end_pos + 1);
    out.port = host_relative.substr(0, relative_ref.data() - host_relative.data());

    // validate port
    bool const valid_port =
        urls::grammar::find_if_not(
            out.port.begin(),
            out.port.end(),
            urls::grammar::digit_chars)
        == out.port.end();
    if (!valid_port)
    {
        // move port to hostname where it can be encoded
        out.hostname = {out.hostname.begin(), out.port.end()};
        out.port = {};
    }

    extract_relative_ref(relative_ref, out);
    if (out.port.empty() && !out.scheme.empty())
        out.port = port_of_scheme(out.scheme);
}

void
extract_scheme_relative(
    core::string_view scheme_relative,
    url_components &out)
{
    // hostinfo
    constexpr urls::grammar::lut_chars hostinfo_end_chars("/?#\0");
    auto it = urls::grammar::find_if(
        scheme_relative.begin(),
        scheme_relative.end(),
        hostinfo_end_chars);
    auto path_offset = (std::min)(
        scheme_relative.size(),
        static_cast<std::size_t>(it - scheme_relative.begin()));
    core::string_view host_info = scheme_relative.substr(0, path_offset);

    // userinfo
    core::string_view relative_ref = scheme_relative.substr(path_offset);
    auto host_offset = host_info.find_last_of('@');
    if (host_offset == core::string_view::npos)
        return extract_userinfo_relative(
            relative_ref,
            scheme_relative,
            host_info,
            out);

    // password
    core::string_view userinfo_at_relative = scheme_relative.substr(host_offset);
    core::string_view userinfo(host_info.data(), userinfo_at_relative.data() - host_info.data());
    auto password_offset = std::min(userinfo.size(), userinfo.find(':'));
    if (password_offset != userinfo.size()) {
        out.user = scheme_relative.substr(0, password_offset);
        core::string_view password = scheme_relative.substr(password_offset + 1);
        out.password = password.substr(0, userinfo_at_relative.data() - password.data());
    } else {
        out.user = scheme_relative.substr(0, userinfo_at_relative.data() - scheme_relative.data());
    }

    // userinfo-relative
    core::string_view userinfo_relative = userinfo_at_relative.substr(1);
    it = urls::grammar::find_if(
        userinfo_relative.begin(),
        userinfo_relative.end(),
        hostinfo_end_chars);
    path_offset = (std::min)(
        userinfo_relative.size(),
        static_cast<std::size_t>(it - userinfo_relative.begin()));
    host_info = userinfo_relative.substr(0, path_offset);
    extract_userinfo_relative(
        relative_ref,
        userinfo_relative,
        host_info,
        out);
}

void
extract_uri_components(
   core::string_view s,
   url_components &out)
{
    if (s.starts_with("//") && !s.starts_with("///"))
        return extract_scheme_relative(s.substr(2), out);

    if (s.starts_with('/'))
        return extract_relative_ref(s, out);

    // extract scheme
    // first char in a scheme must be letter (we accept uppercase here)
    bool has_scheme = false;
    if (!s.empty() && urls::grammar::alpha_chars(s.front())) {
        constexpr
            urls::grammar::lut_chars scheme_chars(
                "0123456789+-.ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz");
        char const* it = urls::grammar::find_if_not(
            s.begin() + 1, s.end(), scheme_chars);
        size_t scheme_size = (std::min)(
            s.size(), static_cast<std::size_t>(it - s.begin()));
        // scheme must be non-empty and followed by ':'
        if (s.size() > scheme_size && s[scheme_size] == ':') {
            out.scheme = s.substr(0, scheme_size);
            has_scheme = true;
        }
    }

    // The usual route, parse scheme first
    core::string_view scheme_relative = s;
    if (has_scheme)
        scheme_relative = s.substr(out.scheme.size() + 1);

    const bool has_authority = scheme_relative.starts_with("//");
    const bool is_relative_ref = !has_scheme && !has_authority;
    if (is_relative_ref)
    {
        // this is the trick browsers usually apply when 1) there's no
        // authority because the "//" is missing, 2) the scheme is also missing,
        // and 3) the first path segment looks like an authority
        //
        // This behavior is widespread, although it's ambiguous because valid
        // host characters are also valid path characters.
        //
        // It's this rule that allows for things like "www.boost.org" in the
        // browser. This is an invalid URL because it has no "//" to indicate
        // this is the authority and "www.boost.org" is a perfectly valid
        // path segment.
        auto first_seg_offset = (std::min)(s.size(), s.find_first_of('/'));
        core::string_view first_seg = s.substr(0, first_seg_offset);
        auto host_delimiter_pos = first_seg.find_first_of(".:");
        bool const looks_like_authority =
            urls::parse_authority(first_seg) &&
            host_delimiter_pos != core::string_view::npos &&
            host_delimiter_pos != first_seg.size() - 1;
        if (looks_like_authority)
            return extract_scheme_relative(s, out);

        // if the first_seg is really a seg, parse as relative ref
        return extract_relative_ref(s, out);
    }

    if (has_authority)
        scheme_relative = scheme_relative.substr(2);

    // all that's left is a relative path
    return extract_relative_ref(scheme_relative, out);
}

void
sanitize_uri(core::string_view s, urls::url_base& dest) {
    url_components o;
    dest.clear();
    extract_uri_components(s, o);
    if (o.scheme.data())
        dest.set_scheme(o.scheme);
    if (o.user.data())
        dest.set_encoded_user(o.user);
    if (o.password.data())
        dest.set_encoded_password(o.password);
    if (o.hostname.data())
        dest.set_encoded_host(o.hostname);
    if (o.port.data())
        dest.set_port(o.port);
    if (o.path.data())
        dest.set_encoded_path(o.path);
    if (o.query.data())
        dest.set_encoded_query(o.query);
    if (o.fragment.data())
        dest.set_encoded_fragment(o.fragment);
}

urls::url
sanitize_uri(core::string_view s) {
    urls::url u;
    sanitize_uri(s, u);
    return u;
}

void
print_url_components(urls::url_view u)
{
    std::cout << "url: " << u.buffer() << '\n';
    if (u.has_scheme())
        std::cout << "scheme: " << u.scheme() << '\n';
    if (u.has_userinfo())
        std::cout << "user: " << u.encoded_user() << '\n';
    if (u.has_password())
        std::cout << "password: " << u.encoded_password() << '\n';
    if (u.has_authority())
        std::cout << "hostname: " << u.encoded_host() << '\n';
    if (u.has_port())
        std::cout << "port: " << u.port() << '\n';
    std::cout << "path: " << u.encoded_path() << '\n';
    std::cout << "segments:\n";
    for (auto seg: u.encoded_segments())
        std::cout << "- " << seg << '\n';
    if (u.has_query())
        std::cout << "query: " << u.encoded_query() << '\n';
    std::cout << "params:\n";
    for (auto param: u.encoded_params())
    {
        if (param.has_value)
            std::cout << "- " << param.key << ": " << param.value << '\n';
        else
            std::cout << "- " << param.key << '\n';
    }
    if (u.has_fragment())
        std::cout << "fragment: " << u.encoded_fragment() << '\n';
}

int
main(int argc, char **argv)
{
    if (argc != 2)
    {
        core::string_view exec = argv[0];
        auto p = exec.find_last_of("/\\");
        if (p != core::string_view::npos)
            exec = exec.substr(p);
        std::cerr
            << "Usage: " << exec
            << " <url>\n"
               "target: a non-strict url\n";
        return EXIT_FAILURE;
    }

    core::string_view uri_str = argv[1];

    boost::system::result<urls::url_view> ru = urls::parse_uri_reference(uri_str);
    if (ru)
    {
        urls::url_view u = *ru;
        if (u.has_scheme() && u.has_fragment())
            std::cout << "Input is a valid URL\n";
        else if (u.has_scheme())
            std::cout << "Input is a valid absolute URL\n";
        else
            std::cout << "Input is a valid relative URL\n";
        print_url_components(u);
        return EXIT_SUCCESS;
    }

    std::cout << "Sanitizing URL:\n";
    std::cout << "input: " << uri_str << '\n';
    urls::url u = sanitize_uri(uri_str);
    print_url_components(u);
    return EXIT_SUCCESS;
}

// end::example_sanitize_url[]
