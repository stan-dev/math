//
// Copyright (c) 2023 Alan de Freitas (alandefreitas@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

#ifndef BOOST_URL_DETAIL_ROUTER_IPP
#define BOOST_URL_DETAIL_ROUTER_IPP

#include "../router.hpp"
#include <boost/url/decode_view.hpp>
#include <boost/url/grammar/alnum_chars.hpp>
#include <boost/url/grammar/alpha_chars.hpp>
#include <boost/url/grammar/lut_chars.hpp>
#include <boost/url/grammar/token_rule.hpp>
#include <boost/url/grammar/variant_rule.hpp>
#include <boost/url/rfc/detail/path_rules.hpp>
#include <boost/url/detail/replacement_field_rule.hpp>
#include <vector>

namespace boost {
namespace urls {
namespace detail {

// A path segment template
class segment_template
{
    enum class modifier : unsigned char
    {
        none,
        // {id?}
        optional,
        // {id*}
        star,
        // {id+}
        plus
    };

    std::string str_;
    bool is_literal_ = true;
    modifier modifier_ = modifier::none;

    friend struct segment_template_rule_t;
public:
    segment_template() = default;

    bool
    match(pct_string_view seg) const;

    core::string_view
    string() const
    {
        return str_;
    }

    core::string_view
    id() const;

    bool
    empty() const
    {
        return str_.empty();
    }

    bool
    is_literal() const
    {
        return is_literal_;
    }

    bool
    has_modifier() const
    {
        return !is_literal() &&
               modifier_ != modifier::none;
    }

    bool
    is_optional() const
    {
        return modifier_ == modifier::optional;
    }

    bool
    is_star() const
    {
        return modifier_ == modifier::star;
    }

    bool
    is_plus() const
    {
        return modifier_ == modifier::plus;
    }

    friend
    bool operator==(
        segment_template const& a,
        segment_template const& b)
    {
        if (a.is_literal_ != b.is_literal_)
            return false;
        if (a.is_literal_)
            return a.str_ == b.str_;
        return a.modifier_ == b.modifier_;
    }

    // segments have precedence:
    //     - literal
    //     - unique
    //     - optional
    //     - plus
    //     - star
    friend
    bool operator<(
        segment_template const& a,
        segment_template const& b)
    {
        if (b.is_literal())
            return false;
        if (a.is_literal())
            return !b.is_literal();
        return a.modifier_ < b.modifier_;
    }
};

// A segment template is either a literal string
// or a replacement field (as in a format_string).
// Fields cannot contain format specs and might
// have one of the following modifiers:
// - ?: optional segment
// - *: zero or more segments
// - +: one or more segments
struct segment_template_rule_t
{
    using value_type = segment_template;

    system::result<value_type>
    parse(
        char const*& it,
        char const* end
    ) const noexcept;
};

constexpr auto segment_template_rule = segment_template_rule_t{};

constexpr auto path_template_rule =
    grammar::tuple_rule(
        grammar::squelch(
            grammar::optional_rule(
                grammar::delim_rule('/'))),
        grammar::range_rule(
            segment_template_rule,
            grammar::tuple_rule(
                grammar::squelch(grammar::delim_rule('/')),
                segment_template_rule)));

bool
segment_template::
match(pct_string_view seg) const
{
    if (is_literal_)
        return *seg == str_;

    // other nodes match any string
    return true;
}

core::string_view
segment_template::
id() const
{
    // if (is_literal_) return {};
    BOOST_ASSERT(!is_literal());
    core::string_view r = {str_};
    r.remove_prefix(1);
    r.remove_suffix(1);
    if (r.ends_with('?') ||
        r.ends_with('+') ||
        r.ends_with('*'))
        r.remove_suffix(1);
    return r;
}

auto
segment_template_rule_t::
parse(
    char const*& it,
    char const* end) const noexcept
    -> system::result<value_type>
{
    segment_template t;
    if (it != end &&
        *it == '{')
    {
        // replacement field
        auto it0 = it;
        ++it;
        auto send =
            grammar::find_if(
                it, end, grammar::lut_chars('}'));
        if (send != end)
        {
            core::string_view s(it, send);
            static constexpr auto modifiers_cs =
                grammar::lut_chars("?*+");
            static constexpr auto id_rule =
                grammar::tuple_rule(
                    grammar::optional_rule(
                        arg_id_rule),
                    grammar::optional_rule(
                        grammar::delim_rule(modifiers_cs)));
            if (s.empty() ||
                grammar::parse(s, id_rule))
            {
                it = send + 1;

                t.str_ = core::string_view(it0, send + 1);
                t.is_literal_ = false;
                if (s.ends_with('?'))
                    t.modifier_ =
                        segment_template::modifier::optional;
                else if (s.ends_with('*'))
                    t.modifier_ =
                        segment_template::modifier::star;
                else if (s.ends_with('+'))
                    t.modifier_ =
                        segment_template::modifier::plus;
                return t;
            }
        }
        it = it0;
    }
    // literal segment
    auto rv = grammar::parse(
        it, end, urls::detail::segment_rule);
    BOOST_ASSERT(rv);
    rv->decode({}, urls::string_token::assign_to(t.str_));
    t.is_literal_ = true;
    return t;
}

// a small vector for child nodes...
// we shouldn't expect many children per node, and
// we don't want to allocate for that. But we also
// cannot cap the max number of child nodes because
// especially the root nodes can potentially an
// exponentially higher number of child nodes.
class child_idx_vector
{
    static constexpr std::size_t N = 5;
    std::size_t static_child_idx_[N]{};
    std::size_t* child_idx{nullptr};
    std::size_t size_{0};
    std::size_t cap_{0};

public:
    ~child_idx_vector()
    {
        delete[] child_idx;
    }

    child_idx_vector() = default;

    child_idx_vector(child_idx_vector const& other)
        : size_{other.size_}
        , cap_{other.cap_}
    {
        if (other.child_idx)
        {
            child_idx = new std::size_t[cap_];
            std::memcpy(child_idx, other.child_idx, size_ * sizeof(std::size_t));
            return;
        }
        std::memcpy(static_child_idx_, other.static_child_idx_, size_ * sizeof(std::size_t));
    }

    child_idx_vector(child_idx_vector&& other)
        : child_idx{other.child_idx}
        , size_{other.size_}
        , cap_{other.cap_}
    {
        std::memcpy(static_child_idx_, other.static_child_idx_, N);
        other.child_idx = nullptr;
    }

    bool
    empty() const
    {
        return size_ == 0;
    }

    std::size_t
    size() const
    {
        return size_;
    }

    std::size_t*
    begin()
    {
        if (child_idx)
            return child_idx;
        return static_child_idx_;
    }

    std::size_t*
    end()
    {
        return begin() + size_;
    }

    std::size_t const*
    begin() const
    {
        if (child_idx)
            return child_idx;
        return static_child_idx_;
    }

    std::size_t const*
    end() const
    {
        return begin() + size_;
    }

    void
    erase(std::size_t* it)
    {
        BOOST_ASSERT(it - begin() >= 0);
        std::memmove(it - 1, it, end() - it);
        --size_;
    }

    void
    push_back(std::size_t v)
    {
        if (size_ == N && !child_idx)
        {
            child_idx = new std::size_t[N*2];
            cap_ = N*2;
            std::memcpy(child_idx, static_child_idx_, N * sizeof(std::size_t));
        }
        else if (child_idx && size_ == cap_)
        {
            auto* tmp = new std::size_t[cap_*2];
            std::memcpy(tmp, child_idx, cap_ * sizeof(std::size_t));
            delete[] child_idx;
            child_idx = tmp;
            cap_ = cap_*2;
        }
        begin()[size_++] = v;
    }
};

// A node in the resource tree
// Each segment in the resource tree might be
// associated with
struct node
{
    static constexpr std::size_t npos{std::size_t(-1)};

    // literal segment or replacement field
    detail::segment_template seg{};

    // A pointer to the resource
    router_base::any_resource const* resource{nullptr};

    // The complete match for the resource
    std::string path_template;

    // Index of the parent node in the
    // implementation pool of nodes
    std::size_t parent_idx{npos};

    // Index of child nodes in the pool
    detail::child_idx_vector child_idx;
};

class impl
{
    // Pool of nodes in the resource tree
    std::vector<node> nodes_;

public:
    impl()
    {
        // root node with no resource
        nodes_.push_back(node{});
    }

    ~impl()
    {
        for (auto &r: nodes_)
            delete r.resource;
    }

    // include a node for a resource
    void
    insert_impl(
        core::string_view path,
        router_base::any_resource const* v);

    // match a node and return the element
    router_base::any_resource const*
    find_impl(
        segments_encoded_view path,
        core::string_view*& matches,
        core::string_view*& ids) const;

private:
    // try to match from this root node
    node const*
    try_match(
        segments_encoded_view::const_iterator it,
        segments_encoded_view::const_iterator end,
        node const* root,
        int level,
        core::string_view*& matches,
        core::string_view*& ids) const;

    // check if a node has a resource when we
    // also consider optional paths through
    // the child nodes.
    static
    node const*
    find_optional_resource(
        const node* root,
        std::vector<node> const& ns,
        core::string_view*& matches,
        core::string_view*& ids);
};

node const*
impl::
find_optional_resource(
    const node* root,
    std::vector<node> const& ns,
    core::string_view*& matches,
    core::string_view*& ids)
{
    BOOST_ASSERT(root);
    if (root->resource)
        return root;
    BOOST_ASSERT(!root->child_idx.empty());
    for (auto i: root->child_idx)
    {
        auto& c = ns[i];
        if (!c.seg.is_optional() &&
            !c.seg.is_star())
            continue;
        // Child nodes are also
        // potentially optional.
        auto matches0 = matches;
        auto ids0 = ids;
        *matches++ = {};
        *ids++ = c.seg.id();
        auto n = find_optional_resource(
            &c, ns, matches, ids);
        if (n)
            return n;
        matches = matches0;
        ids = ids0;
    }
    return nullptr;
}

void
impl::
insert_impl(
    core::string_view path,
    router_base::any_resource const* v)
{
    // Parse dynamic route segments
    if (path.starts_with("/"))
        path.remove_prefix(1);
    auto segsr =
        grammar::parse(path, detail::path_template_rule);
    if (!segsr)
    {
        delete v;
        segsr.value();
    }
    auto segs = *segsr;
    auto it = segs.begin();
    auto end = segs.end();

    // Iterate existing nodes
    node* cur = &nodes_.front();
    int level = 0;
    while (it != end)
    {
        core::string_view seg = (*it).string();
        if (seg == ".")
        {
            ++it;
            continue;
        }
        if (seg == "..")
        {
            // discount unmatched leaf or
            // keep track of levels behind root
            if (cur == &nodes_.front())
            {
                --level;
                ++it;
                continue;
            }
            // move to parent deleting current
            // if it carries no resource
            std::size_t p_idx = cur->parent_idx;
            if (cur == &nodes_.back() &&
                !cur->resource &&
                cur->child_idx.empty())
            {
                node* p = &nodes_[p_idx];
                std::size_t cur_idx = cur - nodes_.data();
                p->child_idx.erase(
                    std::remove(
                        p->child_idx.begin(),
                        p->child_idx.end(),
                        cur_idx));
                nodes_.pop_back();
            }
            cur = &nodes_[p_idx];
            ++it;
            continue;
        }
        // discount unmatched root parent
        if (level < 0)
        {
            ++level;
            ++it;
            continue;
        }
        // look for child
        auto cit = std::find_if(
            cur->child_idx.begin(),
            cur->child_idx.end(),
            [this, &it](std::size_t ci) -> bool
            {
                return nodes_[ci].seg == *it;
            });
        if (cit != cur->child_idx.end())
        {
            // move to existing child
            cur = &nodes_[*cit];
        }
        else
        {
            // create child if it doesn't exist
            node child;
            child.seg = *it;
            std::size_t cur_id = cur - nodes_.data();
            child.parent_idx = cur_id;
            nodes_.push_back(std::move(child));
            nodes_[cur_id].child_idx.push_back(nodes_.size() - 1);
            if (nodes_[cur_id].child_idx.size() > 1)
            {
                // keep nodes sorted
                auto& cs = nodes_[cur_id].child_idx;
                std::size_t n = cs.size() - 1;
                while (n)
                {
                    if (nodes_[cs.begin()[n]].seg < nodes_[cs.begin()[n - 1]].seg)
                        std::swap(cs.begin()[n], cs.begin()[n - 1]);
                    else
                        break;
                    --n;
                }
            }
            cur = &nodes_.back();
        }
        ++it;
    }
    if (level != 0)
    {
        delete v;
        urls::detail::throw_invalid_argument();
    }
    cur->resource = v;
    cur->path_template = path;
}

node const*
impl::
try_match(
    segments_encoded_view::const_iterator it,
    segments_encoded_view::const_iterator end,
    node const* cur,
    int level,
    core::string_view*& matches,
    core::string_view*& ids) const
{
    while (it != end)
    {
        pct_string_view s = *it;
        if (*s == ".")
        {
            // ignore segment
            ++it;
            continue;
        }
        if (*s == "..")
        {

            // move back to the parent node
            ++it;
            if (level <= 0 &&
                cur != &nodes_.front())
            {
                if (!cur->seg.is_literal())
                {
                    --matches;
                    --ids;
                }
                cur = &nodes_[cur->parent_idx];
            }
            else
                // there's no parent, so we
                // discount that from the implicit
                // tree beyond terminals
                --level;
            continue;
        }

        // we are in the implicit tree above the
        // root, so discount that as a level
        if (level < 0)
        {
            ++level;
            ++it;
            continue;
        }

        // calculate the lower bound on the
        // possible number of branches to
        // determine if we need to branch.
        // We branch when we might have more than
        // one child matching node at this level.
        // If so, we need to potentially branch
        // to find which path leads to a valid
        // resource. Otherwise, we can just
        // consume the node and input without
        // any recursive function calls.
        bool branch = false;
        if (cur->child_idx.size() > 1)
        {
            int branches_lb = 0;
            for (auto i: cur->child_idx)
            {
                auto& c = nodes_[i];
                if (c.seg.is_literal() ||
                    !c.seg.has_modifier())
                {
                    // a literal path counts only
                    // if it matches
                    branches_lb += c.seg.match(s);
                }
                else
                {
                    // everything not matching
                    // a single path counts as
                    // more than one path already
                    branches_lb = 2;
                }
                if (branches_lb > 1)
                {
                    // already know we need to
                    // branch
                    branch = true;
                    break;
                }
            }
        }

        // attempt to match each child node
        node const* r = nullptr;
        bool match_any = false;
        for (auto i: cur->child_idx)
        {
            auto& c = nodes_[i];
            if (c.seg.match(s))
            {
                if (c.seg.is_literal())
                {
                    // just continue from the
                    // next segment
                    if (branch)
                    {
                        r = try_match(
                            std::next(it), end,
                            &c, level,
                            matches, ids);
                        if (r)
                            break;
                    }
                    else
                    {
                        cur = &c;
                        match_any = true;
                        break;
                    }
                }
                else if (!c.seg.has_modifier())
                {
                    // just continue from the
                    // next segment
                    if (branch)
                    {
                        auto matches0 = matches;
                        auto ids0 = ids;
                        *matches++ = *it;
                        *ids++ = c.seg.id();
                        r = try_match(
                            std::next(it), end, &c,
                            level, matches, ids);
                        if (r)
                        {
                            break;
                        }
                        else
                        {
                            // rewind
                            matches = matches0;
                            ids = ids0;
                        }
                    }
                    else
                    {
                        // only path possible
                        *matches++ = *it;
                        *ids++ = c.seg.id();
                        cur = &c;
                        match_any = true;
                        break;
                    }
                }
                else if (c.seg.is_optional())
                {
                    // attempt to match by ignoring
                    // and not ignoring the segment.
                    // we first try the complete
                    // continuation consuming the
                    // input, which is the
                    // longest and most likely
                    // match
                    auto matches0 = matches;
                    auto ids0 = ids;
                    *matches++ = *it;
                    *ids++ = c.seg.id();
                    r = try_match(
                        std::next(it), end,
                        &c, level, matches, ids);
                    if (r)
                        break;
                    // rewind
                    matches = matches0;
                    ids = ids0;
                    // try complete continuation
                    // consuming no segment
                    *matches++ = {};
                    *ids++ = c.seg.id();
                    r = try_match(
                        it, end, &c,
                        level, matches, ids);
                    if (r)
                        break;
                    // rewind
                    matches = matches0;
                    ids = ids0;
                }
                else
                {
                    // check if the next segments
                    // won't send us to a parent
                    // directory
                    auto first = it;
                    std::size_t ndotdot = 0;
                    std::size_t nnondot = 0;
                    auto it1 = it;
                    while (it1 != end)
                    {
                        if (*it1 == "..")
                        {
                            ++ndotdot;
                            if (ndotdot >= (nnondot + c.seg.is_star()))
                                break;
                        }
                        else if (*it1 != ".")
                        {
                            ++nnondot;
                        }
                        ++it1;
                    }
                    if (it1 != end)
                        break;

                    // attempt to match many
                    // segments
                    auto matches0 = matches;
                    auto ids0 = ids;
                    *matches++ = *it;
                    *ids++ = c.seg.id();
                    // if this is a plus seg, we
                    // already consumed the first
                    // segment
                    if (c.seg.is_plus())
                    {
                        ++first;
                    }
                    // {*} is usually the last
                    // match in a path.
                    // try complete continuation
                    // match for every subrange
                    // from {last, last} to
                    // {first, last}.
                    // We also try {last, last}
                    // first because it is the
                    // longest match.
                    auto start = end;
                    while (start != first)
                    {
                        r = try_match(
                            start, end, &c,
                            level, matches, ids);
                        if (r)
                        {
                            core::string_view prev = *std::prev(start);
                            *matches0 = {
                                matches0->data(),
                                prev.data() + prev.size()};
                            break;
                        }
                        matches = matches0 + 1;
                        ids = ids0 + 1;
                        --start;
                    }
                    if (r)
                    {
                        break;
                    }
                    // start == first
                    matches = matches0 + 1;
                    ids = ids0 + 1;
                    r = try_match(
                        start, end, &c,
                        level, matches, ids);
                    if (r)
                    {
                        if (!c.seg.is_plus())
                            *matches0 = {};
                        break;
                    }
                }
            }
        }
        // r represent we already found a terminal
        // node which is a match
        if (r)
            return r;
        // if we couldn't match anything, we go
        // one level up in the implicit tree
        // because the path might still have a
        // "..".
        if (!match_any)
            ++level;
        ++it;
    }
    if (level != 0)
    {
        // the path ended below or above an
        // existing node
        return nullptr;
    }
    if (!cur->resource)
    {
        // we consumed all the input and reached
        // a node with no resource, but it might
        // still have child optional segments
        // with resources we can reach without
        // consuming any input
        return find_optional_resource(
            cur, nodes_, matches, ids);
    }
    return cur;
}

router_base::any_resource const*
impl::
find_impl(
    segments_encoded_view path,
    core::string_view*& matches,
    core::string_view*& ids) const
{
    // parse_path is inconsistent for empty paths
    if (path.empty())
        path = segments_encoded_view("./");

    // Iterate nodes from the root
    node const*p = try_match(
        path.begin(), path.end(),
        &nodes_.front(), 0,
        matches, ids);
    if (p)
        return p->resource;
    return nullptr;
}

router_base::
router_base()
    : impl_(new impl{}) {}

router_base::
~router_base()
{
    delete reinterpret_cast<impl*>(impl_);
}

void
router_base::
insert_impl(
    core::string_view s,
    any_resource const* v)
{
    reinterpret_cast<impl*>(impl_)
        ->insert_impl(s, v);
}

auto
router_base::
find_impl(
    segments_encoded_view path,
    core::string_view*& matches,
    core::string_view*& ids) const noexcept
    -> any_resource const*
{
    return reinterpret_cast<impl*>(impl_)
        ->find_impl(path, matches, ids);
}

} // detail
} // urls
} // boost

#endif
