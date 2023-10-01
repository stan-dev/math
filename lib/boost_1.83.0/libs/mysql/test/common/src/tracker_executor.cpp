//
// Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include "test_common/tracker_executor.hpp"

#include <boost/asio/execution/blocking.hpp>
#include <boost/asio/execution/relationship.hpp>
#include <boost/asio/require.hpp>

using namespace boost::mysql::test;

namespace boost {
namespace mysql {
namespace test {

class tracker_executor
{
public:
    tracker_executor(boost::asio::any_io_executor ex, executor_info* tracked) noexcept
        : ex_(ex), tracked_(tracked)
    {
    }

    tracker_executor(const tracker_executor& rhs) noexcept : ex_(rhs.ex_), tracked_(rhs.tracked_) {}
    tracker_executor(tracker_executor&& rhs) noexcept : ex_(std::move(rhs.ex_)), tracked_(rhs.tracked_) {}
    tracker_executor& operator=(const tracker_executor& rhs) noexcept
    {
        ex_ = rhs.ex_;
        tracked_ = rhs.tracked_;
        return *this;
    }
    tracker_executor& operator=(tracker_executor&& rhs) noexcept
    {
        ex_ = std::move(rhs.ex_);
        tracked_ = std::move(rhs.tracked_);
        return *this;
    }
    ~tracker_executor() = default;

    template <class Property>
    tracker_executor require(
        const Property& p,
        typename std::enable_if<asio::can_require<asio::any_io_executor, Property>::value>::type* = nullptr
    ) const
    {
        return tracker_executor(asio::require(ex_, p), tracked_);
    }

    template <class Property>
    tracker_executor prefer(
        const Property& p,
        typename std::enable_if<asio::can_prefer<asio::any_io_executor, Property>::value>::type* = nullptr
    ) const
    {
        return tracker_executor(asio::prefer(ex_, p), tracked_);
    }

    template <class Property>
    auto query(
        const Property& p,
        typename std::enable_if<asio::can_query<asio::any_io_executor, Property>::value>::type* = nullptr
    ) const -> decltype(asio::query(std::declval<boost::asio::any_io_executor>(), p))
    {
        return boost::asio::query(ex_, p);
    }

    template <typename Function>
    void execute(Function&& f) const
    {
        if (asio::query(ex_, asio::execution::relationship) == asio::execution::relationship.fork &&
            asio::query(ex_, asio::execution::blocking) == asio::execution::blocking.never)
        {
            // This is a post
            ++tracked_->num_posts;
        }
        else
        {
            ++tracked_->num_dispatches;
        }
        ex_.execute(std::forward<Function>(f));
    }

    bool operator==(const tracker_executor& rhs) const noexcept
    {
        return ex_ == rhs.ex_ && tracked_ == rhs.tracked_;
    }
    bool operator!=(const tracker_executor& rhs) const noexcept { return !(*this == rhs); }

    executor_info* get_tracked() const noexcept { return tracked_; }

private:
    boost::asio::any_io_executor ex_;
    executor_info* tracked_;
};

}  // namespace test
}  // namespace mysql
}  // namespace boost

// Asio fails to detect that a type is equaility comparable under MSVC, so we need to do this
namespace boost {
namespace asio {
namespace traits {

template <>
struct equality_comparable<tracker_executor>
{
    static constexpr bool is_valid = true;
    static constexpr bool is_noexcept = true;
};

template <typename Function>
struct execute_member<tracker_executor, Function>
{
    static constexpr bool is_valid = true;
    static constexpr bool is_noexcept = false;
    using result_type = void;
};

template <class Property>
struct require_member<
    tracker_executor,
    Property,
    typename std::enable_if<asio::can_require<asio::any_io_executor, Property>::value>::type>
{
    static constexpr bool is_valid = true;
    static constexpr bool is_noexcept = false;
    using result_type = tracker_executor;
};

template <typename Property>
struct prefer_member<
    tracker_executor,
    Property,
    typename std::enable_if<asio::can_prefer<asio::any_io_executor, Property>::value>::type>
{
    static constexpr bool is_valid = true;
    static constexpr bool is_noexcept = false;
    using result_type = tracker_executor;
};

template <typename Property>
struct query_member<
    tracker_executor,
    Property,
    typename std::enable_if<asio::can_query<asio::any_io_executor, Property>::value>::type>
{
    static constexpr bool is_valid = true;
    static constexpr bool is_noexcept = false;
    using result_type = decltype(boost::asio::query(
        std::declval<const any_io_executor&>(),
        std::declval<const Property&>()
    ));
};

}  // namespace traits
}  // namespace asio
}  // namespace boost

boost::asio::any_io_executor boost::mysql::test::create_tracker_executor(
    asio::any_io_executor inner,
    executor_info* tracked_values
)
{
    return tracker_executor(std::move(inner), tracked_values);
}

executor_info boost::mysql::test::get_executor_info(const asio::any_io_executor& exec)
{
    return *exec.target<tracker_executor>()->get_tracked();
}