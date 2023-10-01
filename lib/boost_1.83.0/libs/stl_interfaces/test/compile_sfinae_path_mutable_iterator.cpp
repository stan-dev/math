// Copyright (C) 2022 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
#include <boost/stl_interfaces/iterator_interface.hpp>
#include <type_traits>
#include <utility>

struct node
{
    using value_type = std::pair<int, int>;

    value_type kv_;
    node * next_;
};

template<typename Node>
struct iterator
    : boost::stl_interfaces::iterator_interface<
#if !BOOST_STL_INTERFACES_USE_DEDUCED_THIS
          iterator<Node>,
#endif
          std::forward_iterator_tag, Node>
{
    using value_type = typename Node::value_type;

    constexpr iterator() noexcept = default;
    constexpr explicit iterator(Node * it) : it_{it} {}

    template<
        typename Node2,
        typename Enable = std::enable_if_t<
            std::is_convertible<Node2 *, Node *>{} && !std::is_const<Node2>{} &&
            std::is_const<Node>{}>>
    constexpr iterator(iterator<Node2> other) noexcept : it_{other.it_}
    {}

    constexpr value_type const & operator*() const noexcept { return it_->kv_; }

    template<
        typename T = Node,
        typename std::enable_if_t<!std::is_const<T>{}, bool> = true>
    constexpr value_type & operator*() noexcept
    {
        return it_->kv_;
    }

    constexpr iterator & operator++() noexcept
    {
        it_ = it_->next_;
        return *this;
    }

    using base_type = boost::stl_interfaces::
        iterator_interface<iterator<Node>, std::forward_iterator_tag, Node>;
    using base_type::operator++;

    friend constexpr bool operator==(iterator lhs, iterator rhs) noexcept
    {
        return lhs.it_ == rhs.it_;
    }

private:
    Node * it_;

    template<typename Node2>
    friend struct iterator;
};

void compile_sfinae_path_mutable_iterator()
{
    auto it = iterator<node>{};
    if ((*it).first == 42) {
        (*it).second = 13;
    }
}
