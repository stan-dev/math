//
// Copyright (c) 2022 alandefreitas (alandefreitas@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//

#ifndef BOOST_URL_EXAMPLE_MAILTO_MAILTO_GRAMMAR_HPP
#define BOOST_URL_EXAMPLE_MAILTO_MAILTO_GRAMMAR_HPP

#include <boost/url/grammar/alnum_chars.hpp>
#include <boost/url/grammar/ci_string.hpp>
#include <boost/url/grammar/delim_rule.hpp>
#include <boost/url/grammar/digit_chars.hpp>
#include <boost/url/grammar/optional_rule.hpp>
#include <boost/url/grammar/range_rule.hpp>
#include <boost/url/grammar/token_rule.hpp>
#include <boost/url/grammar/tuple_rule.hpp>
#include <boost/url/grammar/variant_rule.hpp>
#include <boost/url/grammar/vchars.hpp>
#include <boost/url/rfc/pct_encoded_rule.hpp>
#include <boost/url/rfc/unreserved_chars.hpp>
#include <algorithm>

namespace urls = boost::urls;
namespace grammar = boost::urls::grammar;

/// The set of dtext_no_obs characters
struct dtext_no_obs_chars_t
{
    constexpr
    bool
    operator()(char c) const noexcept
    {
        // dtext-no-obs = %d33-90 / %d94-126
        return (c >= '!' && c <= 'Z') ||
               (c >= '^' && c <= '~');
    }
};

/// A character set containing dtext_no_obs characters.
constexpr dtext_no_obs_chars_t dtext_no_obs_chars{};

/// A character set containing crlf characters.
constexpr auto crlf_chars =
    grammar::lut_chars("\n\r");

/// A character set containing atext characters.
constexpr auto atext_chars =
    grammar::lut_chars("!#$%&'*+-/=?^_`{|}~") +
    grammar::alnum_chars +
    grammar::digit_chars;

/// A character set containing wsp characters.
constexpr auto wsp_chars =
    grammar::lut_chars(" \t");

/// The set of obs_ctext characters
struct obs_ctext_chars_t
{
    static constexpr unsigned char c1 = char(1);
    static constexpr unsigned char c2 = char(8);
    static constexpr unsigned char c3 = char(11);
    static constexpr unsigned char c4 = char(12);
    static constexpr unsigned char c5 = char(14);
    static constexpr unsigned char c6 = char(31);
    static constexpr unsigned char c7 = char(127);

    constexpr obs_ctext_chars_t() noexcept = default;

    constexpr
    bool
    operator()(char c) const noexcept
    {
        // obs-ctext = obs-NO-WS-CTL
        // obs-NO-WS-CTL = %d1-8 / %d11 / %d12 / %d14-31 / %d127
        return
            (static_cast<unsigned char>(c) >= c1 &&
             static_cast<unsigned char>(c) <= c2) ||
            (static_cast<unsigned char>(c) >= c3 &&
             static_cast<unsigned char>(c) <= c4) ||
            static_cast<unsigned char>(c) == c5 ||
            static_cast<unsigned char>(c) == c6 ||
            static_cast<unsigned char>(c) == c7;
    }
};

/// A character set containing obs_ctext characters.
constexpr obs_ctext_chars_t obs_ctext_chars{};

/// A character set containing obs_qtext characters.
constexpr auto obs_qtext_chars = obs_ctext_chars;

/// The set of ctext characters
struct ctext_chars_t
{
    constexpr ctext_chars_t() noexcept = default;

    constexpr
    bool
    operator()(char c) const noexcept
    {
        // ctext = %d33-39 / %d42-91 / %d93-126 / obs-ctext
        return
            (c >= '!' && c <= '\'') ||
            (c >= '*' && c <= '[') ||
            (c >= ']' && c <= '~') ||
            obs_ctext_chars(c);
    }
};

/// A character set containing ctext characters.
constexpr ctext_chars_t ctext_chars{};

/// The set of qtext characters
struct qtext_chars_t
{
    constexpr qtext_chars_t() noexcept = default;

    constexpr
    bool
    operator()(char c) const noexcept
    {
        // qtext = %d33 / %d35-91 / %d93-126 / obs-qtext
        return c == '!' ||
               (c >= '#' && c <= '[') ||
               (c >= ']' && c <= '~') ||
               obs_qtext_chars(c);
    }
};

/// A character set containing qtext characters.
constexpr qtext_chars_t qtext_chars{};

/// A character set containing qchars
constexpr auto qchars = urls::unreserved_chars + "!$\\()*+,;:@";

constexpr auto atext_token =
    grammar::token_rule(atext_chars);

/// Rule for dot-atom-text = 1*atext *("." 1*atext)
constexpr auto dot_atom_text_rule =
    grammar::range_rule(
        atext_token,
        grammar::tuple_rule(
            grammar::squelch(
                grammar::delim_rule('.')),
            atext_token));

/// Rule for "[" *dtext-no-obs "]"
constexpr auto quoted_dtext_no_obs =
    grammar::tuple_rule(
        grammar::squelch(
            grammar::delim_rule('[')),
        grammar::optional_rule(
            grammar::token_rule(dtext_no_obs_chars)),
        grammar::squelch(
            grammar::delim_rule(']')));

/// Rule for domain = dot-atom-text / "[" *dtext-no-obs "]"
constexpr auto domain_rule =
    grammar::variant_rule(
        dot_atom_text_rule,
        quoted_dtext_no_obs);

/// Rule for obs-qp = "\" (%d0 / obs-NO-WS-CTL / LF / CR)
constexpr auto obs_qp_rule =
    grammar::tuple_rule(
        grammar::delim_rule('\\'),
        grammar::variant_rule(
            grammar::delim_rule('\0'),
            grammar::delim_rule(ctext_chars),
            grammar::delim_rule(crlf_chars)));

/// Rule for quoted-pair = ("\" (VCHAR / WSP)) / obs-qp
constexpr auto quoted_pair_rule =
    grammar::variant_rule(
        grammar::tuple_rule(
            grammar::delim_rule('\\'),
            grammar::variant_rule(
                grammar::delim_rule(grammar::vchars),
                grammar::delim_rule(' '))),
        obs_qp_rule);

/// Rule for obs-FWS = 1*WSP *(CRLF 1*WSP)
constexpr auto obs_fws_rule =
    grammar::tuple_rule(
        grammar::token_rule(wsp_chars),
        grammar::range_rule(
            grammar::delim_rule(crlf_chars),
            grammar::token_rule(wsp_chars)));

/// Rule for FWS = ([*WSP CRLF] 1*WSP) / obs-FWS
constexpr auto fws_rule =
    grammar::variant_rule(
        grammar::tuple_rule(
            grammar::optional_rule(
                grammar::tuple_rule(
                    grammar::optional_rule(
                        grammar::token_rule(wsp_chars)),
                    grammar::delim_rule(crlf_chars))),
            grammar::token_rule(wsp_chars)),
        obs_fws_rule);

namespace detail
{
    // workaround for value-based recursive rules
    struct ccontent_and_comment_rules {
        struct ccontent_rule_t
        {
            using value_type = urls::string_view;

            urls::result< value_type >
            parse(
                char const*& it,
                char const* end
            ) const noexcept
            {
                auto it0 = it;
                bool v = ccontent_and_comment_rules::
                    parse_ccontent(it, end);
                if (v)
                    return urls::string_view(it0, it);
                return grammar::error::invalid;
            }
        };

        static
        bool
        parse_ccontent(char const*& it, char const* end) noexcept
        {
            // ccontent = ctext / quoted-pair / comment
            return
                grammar::parse(
                    it, end,
                    grammar::variant_rule(
                        grammar::delim_rule(ctext_chars),
                        quoted_pair_rule,
                        ccontent_rule_t{})).has_value();
        };

        struct comment_rule_t
        {
            using value_type = urls::string_view;

            urls::result< value_type >
            parse(
                char const*& it,
                char const* end
            ) const noexcept
            {
                auto it0 = it;
                bool v = ccontent_and_comment_rules::
                    parse_comment(it, end);
                if (v)
                    return urls::string_view(it0, it);
                return grammar::error::invalid;
            }
        };

        static
        bool
        parse_comment(char const*& it, char const* end) noexcept
        {
            // comment = "(" *([FWS] ccontent) [FWS] ")"
            return grammar::parse(
                it, end,
                grammar::tuple_rule(
                    grammar::delim_rule('('),
                    grammar::range_rule(
                        grammar::tuple_rule(
                            grammar::optional_rule(fws_rule),
                            ccontent_rule_t{})),
                    grammar::optional_rule(fws_rule),
                    grammar::delim_rule(')'))).has_value();
        }
    };
}

/// Rule for ccontent = ctext / quoted-pair / comment
constexpr auto ccontent_rule =
    detail::ccontent_and_comment_rules::ccontent_rule_t{};

/// Rule for comment = "(" *([FWS] ccontent) [FWS] ")"
constexpr auto comment_rule =
    detail::ccontent_and_comment_rules::comment_rule_t{};

/// Rule for CFWS = (1*([FWS] comment) [FWS]) / FWS
constexpr auto cfws_rule =
    grammar::variant_rule(
        grammar::tuple_rule(
            grammar::range_rule(
                grammar::tuple_rule(
                    grammar::optional_rule(fws_rule),
                    comment_rule), 1),
            grammar::optional_rule(fws_rule)),
        fws_rule);

/// Rule for qcontent = qtext / quoted-pair
constexpr auto qcontent_rule =
    grammar::variant_rule(
        grammar::delim_rule(qtext_chars),
        quoted_pair_rule);

/// Rule for quoted-string = [CFWS] DQUOTE *([FWS] qcontent) [FWS] DQUOTE [CFWS]
constexpr auto quoted_string_rule =
    grammar::tuple_rule(
        grammar::optional_rule(cfws_rule),
        grammar::delim_rule('"'),
        grammar::range_rule(
            grammar::tuple_rule(
                grammar::optional_rule(fws_rule),
                qcontent_rule)),
        grammar::optional_rule(fws_rule),
        grammar::delim_rule('"'),
        grammar::optional_rule(cfws_rule));

/// Rule for local-part = dot-atom-text / quoted-string
constexpr auto local_part_rule =
    grammar::variant_rule(
        dot_atom_text_rule,
        quoted_string_rule);

/// Rule for addr-spec = local-part "@" domain
constexpr auto addr_spec_rule =
    grammar::tuple_rule(
        local_part_rule,
        grammar::squelch(
            grammar::delim_rule('@')),
        domain_rule);

/// Rule for to = addr-spec *("," addr-spec )
constexpr auto to_rule =
    grammar::range_rule(
        addr_spec_rule,
        grammar::tuple_rule(
            grammar::squelch(
                grammar::delim_rule(',')),
            addr_spec_rule));

/// Rule for hfvalue = *qchar
constexpr auto hfvalue_rule = urls::pct_encoded_rule(qchars);

/// Rule for hfname = *qchar
struct hfname_rule_t
{
    using value_type = urls::string_view;

    urls::result<value_type>
    parse(
        char const*& it,
        char const* end
    ) const noexcept
    {
        urls::string_view s(it, end);
        it += s.size();
        auto r = grammar::parse(s, hfvalue_rule);
        if (!r)
            return r;

        // The user agent interpreting a 'mailto' URI SHOULD NOT create a
        // message if any of the header fields are considered dangerous
        static const urls::string_view valid_k[] = {
            "to", "subject", "keywords",
            "cc", "body",    "in-reply-to"
        };
        if (std::any_of(
                std::begin(valid_k), std::end(valid_k),
                [s](urls::string_view valid_k)
            {
                return grammar::ci_is_equal(s, valid_k);
            }))
            return s;
        return grammar::error::invalid;
    }
};

constexpr auto hfname_rule = hfname_rule_t{};

#endif
