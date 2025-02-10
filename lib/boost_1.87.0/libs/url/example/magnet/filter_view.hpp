//
// Copyright (c) 2022 alandefreitas (alandefreitas@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//

#ifndef BOOST_URL_EXAMPLE_GRAMMAR_MAGNER_FILTERED_VIEW_HPP
#define BOOST_URL_EXAMPLE_GRAMMAR_MAGNER_FILTERED_VIEW_HPP

/// A small filter and transform view
/**
   This is a view that represents filtered
   and transformed members of the underlying
   container.

   All elements in the container are filtered
   with a given predicate. Elements that pass
   the predicate are then transformed into an
   element of a specified value type.

   @tparam Container Type of underlying container
   @tparam Value Type of elements after transform
   @tparam Predicate Type of predicate to filter elements
   @tparam Transform Type of transform function on filtered elements
 */
template <class Container, class Value, class Predicate, class Transform>
class filter_view
{
    Container c_;
    Predicate p_;
    Transform t_;

public:
    using value_type = Value;

    class iterator
    {
        friend filter_view;

        using base_iter = typename Container::iterator;

        base_iter b_;
        filter_view& v_;

        void
        skip_filtered()
        {
            // Increment the base iterator
            // skipping all elements that
            // don't pass the filter.
            auto e = v_.c_.end();
            while (b_ != e && !v_.p_(*b_))
                ++b_;
        }

        iterator(
            filter_view& v,
            base_iter it)
            : b_(it)
            , v_(v)
        {
            skip_filtered();
        }

    public:
        iterator&
        operator++() noexcept
        {
            ++b_;
            skip_filtered();
            return *this;
        }

        iterator
        operator++(int) noexcept
        {
            auto tmp = *this;
            ++*this;
            return tmp;
        }

        value_type
        operator*()
        {
            return v_.t_(*b_);
        }

        friend
        bool
        operator==(
            iterator const& a,
            iterator const& b) noexcept
        {
            return a.b_ == b.b_;
        }

        friend
        bool
        operator!=(
            iterator const& a,
            iterator const& b) noexcept
        {
            return a.b_ != b.b_;
        }
    };

    filter_view(
        Container const& base,
        Predicate pred = {},
        Transform transform = {})
        : c_(base)
        , p_(pred)
        , t_(transform) {};

    iterator begin()
    {
        return {*this, c_.begin()};
    }

    iterator end()
    {
        return {*this, c_.end()};
    }
};

/// The identity function
template <class T>
struct identity_fn
{
    T operator()(T const& t) { return t; }
};

/// Return a container view filtered with Predicate
template <class Container, class Value, class Predicate>
filter_view<
    Container,
    typename Container::value_type,
    Predicate,
    identity_fn<typename Container::value_type>>
filter(Container&& c, Predicate&& p) {
    return {
        std::forward<Container>(c),
        std::forward<Predicate>(p),
        identity_fn<typename Container::value_type>{}
    };
}

/// Return a container view filtered with Predicate and Transform applied
template <class Container, class Value, class Predicate, class Transform>
filter_view<
    Container,
    typename Container::value_type,
    Predicate,
    Transform>
filter(Container&& c, Predicate&& p, Transform&& t) {
    return {
        std::forward<Container>(c),
        std::forward<Predicate>(p),
        std::forward<Transform>(t)
    };
}


#endif
