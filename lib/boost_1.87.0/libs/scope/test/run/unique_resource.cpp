/*
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * https://www.boost.org/LICENSE_1_0.txt)
 *
 * Copyright (c) 2023-2024 Andrey Semashev
 */
/*!
 * \file   unique_resource.cpp
 * \author Andrey Semashev
 *
 * \brief  This file contains tests for \c unique_resource.
 */

#include <boost/scope/unique_resource.hpp>
#include <boost/core/lightweight_test.hpp>
#include <boost/core/lightweight_test_trait.hpp>
#include <boost/config.hpp>
#include <ostream>
#include <utility>
#include <stdexcept>
#include <type_traits>

template< typename Resource >
struct empty_resource_deleter
{
    void operator() (Resource const&) const noexcept
    {
    }
};

template< typename Resource >
inline void copy_resource(Resource const& from, Resource& to)
{
    to = from;
}

template< typename Resource >
class checking_resource_deleter
{
private:
    Resource* m_deleted;
    int& m_n;

public:
    explicit checking_resource_deleter(int& n) noexcept :
        m_deleted(nullptr),
        m_n(n)
    {
    }

    explicit checking_resource_deleter(Resource& deleted, int& n) noexcept :
        m_deleted(&deleted),
        m_n(n)
    {
    }

    checking_resource_deleter(checking_resource_deleter&& that) noexcept :
        m_deleted(that.m_deleted),
        m_n(that.m_n)
    {
    }

    checking_resource_deleter(checking_resource_deleter const& that) noexcept :
        m_deleted(that.m_deleted),
        m_n(that.m_n)
    {
    }

    // Make sure the deleter is move and copy-assignable
    checking_resource_deleter& operator= (checking_resource_deleter&& that) noexcept
    {
        m_deleted = that.m_deleted;
        return *this;
    }

    checking_resource_deleter& operator= (checking_resource_deleter const& that) noexcept
    {
        m_deleted = that.m_deleted;
        return *this;
    }

    Resource* get_deleted() const noexcept
    {
        return m_deleted;
    }

    void operator() (Resource const& res) const noexcept
    {
        if (m_deleted)
            copy_resource(res, *m_deleted);
        ++m_n;
    }
};

int g_n = 0, g_res1 = 0, g_res2 = 0;

void check_int()
{
    {
        boost::scope::unique_resource< int, empty_resource_deleter< int > > ur;
        BOOST_TEST_EQ(ur.get(), 0);
        BOOST_TEST(!ur.allocated());
        BOOST_TEST(!ur);
    }

    int n = 0;
    {
        boost::scope::unique_resource< int, checking_resource_deleter< int > > ur{ boost::scope::default_resource, checking_resource_deleter< int >(n) };
        BOOST_TEST_EQ(ur.get(), 0);
        BOOST_TEST(!ur.allocated());
        BOOST_TEST(!ur);
    }
    BOOST_TEST_EQ(n, 0);

    {
        boost::scope::unique_resource< int, empty_resource_deleter< int > > ur{ 10 };
        BOOST_TEST_EQ(ur.get(), 10);
        BOOST_TEST(ur.allocated());
        BOOST_TEST(!!ur);
    }

    int deleted_res1 = -1;
    n = 0;
    {
        boost::scope::unique_resource< int, checking_resource_deleter< int > > ur{ 0, checking_resource_deleter< int >(deleted_res1, n) };
        BOOST_TEST_EQ(ur.get(), 0);
        BOOST_TEST(ur.allocated());
    }
    BOOST_TEST_EQ(n, 1);
    BOOST_TEST_EQ(deleted_res1, 0);

    n = 0;
    deleted_res1 = -1;
    {
        boost::scope::unique_resource< int, checking_resource_deleter< int > > ur{ 10, checking_resource_deleter< int >(deleted_res1, n) };
        BOOST_TEST_EQ(ur.get(), 10);
        BOOST_TEST(ur.allocated());
        BOOST_TEST_EQ(ur.get_deleter().get_deleted(), &deleted_res1);
    }
    BOOST_TEST_EQ(n, 1);
    BOOST_TEST_EQ(deleted_res1, 10);

    n = 0;
    deleted_res1 = -1;
    {
        boost::scope::unique_resource< int, checking_resource_deleter< int > > ur{ 10, checking_resource_deleter< int >(deleted_res1, n) };
        BOOST_TEST_EQ(ur.get(), 10);
        BOOST_TEST(ur.allocated());
        ur.release();
        BOOST_TEST(!ur.allocated());
    }
    BOOST_TEST_EQ(n, 0);
    BOOST_TEST_EQ(deleted_res1, -1);

    n = 0;
    deleted_res1 = -1;
    {
        boost::scope::unique_resource< int, checking_resource_deleter< int > > ur{ 10, checking_resource_deleter< int >(deleted_res1, n) };
        BOOST_TEST_EQ(ur.get(), 10);
        BOOST_TEST(ur.allocated());
        BOOST_TEST_EQ(noexcept(ur.reset()), noexcept(std::declval< checking_resource_deleter< int >& >()(10)));
        ur.reset();
        BOOST_TEST(!ur.allocated());
        BOOST_TEST_EQ(n, 1);
        BOOST_TEST_EQ(deleted_res1, 10);
    }
    BOOST_TEST_EQ(n, 1);

    n = 0;
    deleted_res1 = -1;
    {
        boost::scope::unique_resource< int, checking_resource_deleter< int > > ur{ 10, checking_resource_deleter< int >(deleted_res1, n) };
        BOOST_TEST_EQ(ur.get(), 10);
        BOOST_TEST(ur.allocated());
        ur.reset(20);
        BOOST_TEST_EQ(n, 1);
        BOOST_TEST_EQ(deleted_res1, 10);
        deleted_res1 = -1;
        BOOST_TEST_EQ(ur.get(), 20);
        BOOST_TEST(ur.allocated());
    }
    BOOST_TEST_EQ(n, 2);
    BOOST_TEST_EQ(deleted_res1, 20);

    n = 0;
    deleted_res1 = -1;
    {
        boost::scope::unique_resource< int, checking_resource_deleter< int > > ur1{ 10, checking_resource_deleter< int >(deleted_res1, n) };
        BOOST_TEST_EQ(ur1.get(), 10);
        BOOST_TEST(ur1.allocated());
        boost::scope::unique_resource< int, checking_resource_deleter< int > > ur2{ std::move(ur1) };
        BOOST_TEST_EQ(ur2.get(), 10);
        BOOST_TEST(ur2.allocated());
        BOOST_TEST(!ur1.allocated());
    }
    BOOST_TEST_EQ(n, 1);
    BOOST_TEST_EQ(deleted_res1, 10);

    n = 0;
    deleted_res1 = -1;
    int deleted_res2 = -1;
    {
        boost::scope::unique_resource< int, checking_resource_deleter< int > > ur1{ 10, checking_resource_deleter< int >(deleted_res1, n) };
        BOOST_TEST_EQ(ur1.get(), 10);
        BOOST_TEST(ur1.allocated());
        boost::scope::unique_resource< int, checking_resource_deleter< int > > ur2{ 20, checking_resource_deleter< int >(deleted_res2, n) };
        BOOST_TEST_EQ(ur2.get(), 20);
        BOOST_TEST(ur2.allocated());
        ur2 = std::move(ur1);
        BOOST_TEST_EQ(ur2.get(), 10);
        BOOST_TEST(ur2.allocated());
        BOOST_TEST(!ur1.allocated());
        BOOST_TEST_EQ(n, 1);
        BOOST_TEST_EQ(deleted_res1, -1);
        BOOST_TEST_EQ(deleted_res2, 20);
        deleted_res2 = -1;
    }
    BOOST_TEST_EQ(n, 2);
    BOOST_TEST_EQ(deleted_res1, 10);
    BOOST_TEST_EQ(deleted_res2, -1);

    {
        boost::scope::unique_resource< int, empty_resource_deleter< int > > ur1;
        BOOST_TEST_EQ(ur1.get(), 0);
        BOOST_TEST(!ur1.allocated());
        boost::scope::unique_resource< int, empty_resource_deleter< int > > ur2{ 10, empty_resource_deleter< int >() };
        BOOST_TEST_EQ(ur2.get(), 10);
        BOOST_TEST(ur2.allocated());
        using namespace std;
        swap(ur1, ur2);
        BOOST_TEST_EQ(ur1.get(), 10);
        BOOST_TEST(ur1.allocated());
        BOOST_TEST_EQ(ur2.get(), 0);
        BOOST_TEST(!ur2.allocated());
    }

    n = 0;
    deleted_res1 = -1;
    deleted_res2 = -1;
    {
        boost::scope::unique_resource< int, checking_resource_deleter< int > > ur1{ 10, checking_resource_deleter< int >(deleted_res1, n) };
        BOOST_TEST_EQ(ur1.get(), 10);
        BOOST_TEST_EQ(ur1.get_deleter().get_deleted(), &deleted_res1);
        BOOST_TEST(ur1.allocated());
        boost::scope::unique_resource< int, checking_resource_deleter< int > > ur2{ 20, checking_resource_deleter< int >(deleted_res2, n) };
        BOOST_TEST_EQ(ur2.get(), 20);
        BOOST_TEST_EQ(ur2.get_deleter().get_deleted(), &deleted_res2);
        BOOST_TEST(ur2.allocated());
        using namespace std;
        swap(ur1, ur2);
        BOOST_TEST_EQ(n, 0);
        BOOST_TEST_EQ(ur1.get(), 20);
        BOOST_TEST_EQ(ur1.get_deleter().get_deleted(), &deleted_res2);
        BOOST_TEST(ur1.allocated());
        BOOST_TEST_EQ(ur2.get(), 10);
        BOOST_TEST_EQ(ur2.get_deleter().get_deleted(), &deleted_res1);
        BOOST_TEST(ur2.allocated());
    }
    BOOST_TEST_EQ(n, 2);
    BOOST_TEST_EQ(deleted_res1, 10);
    BOOST_TEST_EQ(deleted_res2, 20);

    struct local
    {
        static void raw_func_deleter(int)
        {
            ++g_n;
        }

        static void raw_func_deleter1(int res)
        {
            g_res1 = res;
        }

        static void raw_func_deleter2(int res)
        {
            g_res2 = res;
        }
    };

    g_n = 0;
    {
        boost::scope::unique_resource< int, void (&)(int) > ur(10, local::raw_func_deleter);
        BOOST_TEST_EQ(ur.get(), 10);
        BOOST_TEST(ur.allocated());
        BOOST_TEST(!noexcept(ur.reset()));
    }
    BOOST_TEST_EQ(g_n, 1);

    g_n = 0;
    {
        boost::scope::unique_resource< int, void (&)(int) > ur1(10, local::raw_func_deleter);
        BOOST_TEST_EQ(ur1.get(), 10);
        BOOST_TEST(ur1.allocated());
        boost::scope::unique_resource< int, void (&)(int) > ur2(std::move(ur1));
        BOOST_TEST_EQ(ur2.get(), 10);
        BOOST_TEST(ur2.allocated());
        BOOST_TEST(!ur1.allocated());
    }
    BOOST_TEST_EQ(g_n, 1);

    g_res1 = 0;
    g_res2 = 0;
    {
        boost::scope::unique_resource< int, void (&)(int) > ur1(10, local::raw_func_deleter1);
        BOOST_TEST_EQ(ur1.get(), 10);
        BOOST_TEST(ur1.allocated());
        boost::scope::unique_resource< int, void (&)(int) > ur2(20, local::raw_func_deleter2);
        BOOST_TEST_EQ(ur2.get(), 20);
        BOOST_TEST(ur2.allocated());
        ur2 = std::move(ur1);
        BOOST_TEST_EQ(ur2.get(), 10);
        BOOST_TEST(ur2.allocated());
        BOOST_TEST(!ur1.allocated());
    }
    BOOST_TEST_EQ(g_res1, 10);
    BOOST_TEST_EQ(g_res2, 20);

    g_res1 = 0;
    g_res2 = 0;
    {
        boost::scope::unique_resource< int, void (&)(int) > ur1(10, local::raw_func_deleter1);
        BOOST_TEST_EQ(ur1.get(), 10);
        BOOST_TEST(ur1.allocated());
        boost::scope::unique_resource< int, void (&)(int) > ur2(20, local::raw_func_deleter2);
        BOOST_TEST_EQ(ur2.get(), 20);
        BOOST_TEST(ur2.allocated());
        using namespace std;
        swap(ur1, ur2);
        BOOST_TEST_EQ(ur1.get(), 20);
        BOOST_TEST(ur1.allocated());
        BOOST_TEST_EQ(ur2.get(), 10);
        BOOST_TEST(ur2.allocated());
    }
    BOOST_TEST_EQ(g_res1, 10);
    BOOST_TEST_EQ(g_res2, 20);
}


struct struct_resource
{
    int value;

    struct_resource(int v = 0) noexcept :
        value(v)
    {
    }

    friend bool operator== (struct_resource const& left, struct_resource const& right) noexcept
    {
        return left.value == right.value;
    }

    friend bool operator!= (struct_resource const& left, struct_resource const& right) noexcept
    {
        return !(left == right);
    }

    friend std::ostream& operator<< (std::ostream& strm, struct_resource const& res)
    {
        strm << "{ " << res.value << " }";
        return strm;
    }
};

void check_struct()
{
    {
        boost::scope::unique_resource< struct_resource, empty_resource_deleter< struct_resource > > ur;
        BOOST_TEST_EQ(ur.get(), struct_resource{});
        BOOST_TEST(!ur.allocated());
    }

    int n = 0;
    struct_resource deleted_res1{ -1 };
    {
        boost::scope::unique_resource< struct_resource, checking_resource_deleter< struct_resource > > ur{ struct_resource{ 10 }, checking_resource_deleter< struct_resource >(deleted_res1, n) };
        BOOST_TEST_EQ(ur.get(), struct_resource{ 10 });
        BOOST_TEST(ur.allocated());
    }
    BOOST_TEST_EQ(n, 1);
    BOOST_TEST_EQ(deleted_res1, struct_resource{ 10 });

    n = 0;
    deleted_res1 = struct_resource{ -1 };
    {
        boost::scope::unique_resource< struct_resource, checking_resource_deleter< struct_resource > > ur{ struct_resource{ 10 }, checking_resource_deleter< struct_resource >(deleted_res1, n) };
        BOOST_TEST_EQ(ur.get(), struct_resource{ 10 });
        BOOST_TEST(ur.allocated());
        ur.reset(20);
        BOOST_TEST_EQ(n, 1);
        BOOST_TEST_EQ(deleted_res1, struct_resource{ 10 });
        deleted_res1 = struct_resource{ -1 };
        BOOST_TEST_EQ(ur.get(), struct_resource{ 20 });
        BOOST_TEST(ur.allocated());
    }
    BOOST_TEST_EQ(n, 2);
    BOOST_TEST_EQ(deleted_res1, struct_resource{ 20 });

    n = 0;
    deleted_res1 = struct_resource{ -1 };
    {
        boost::scope::unique_resource< struct_resource, checking_resource_deleter< struct_resource > > ur1{ struct_resource{ 10 }, checking_resource_deleter< struct_resource >(deleted_res1, n) };
        BOOST_TEST_EQ(ur1.get(), struct_resource{ 10 });
        BOOST_TEST(ur1.allocated());
        boost::scope::unique_resource< struct_resource, checking_resource_deleter< struct_resource > > ur2{ std::move(ur1) };
        BOOST_TEST_EQ(ur2.get(), struct_resource{ 10 });
        BOOST_TEST(ur2.allocated());
        BOOST_TEST(!ur1.allocated());
    }
    BOOST_TEST_EQ(n, 1);
    BOOST_TEST_EQ(deleted_res1, struct_resource{ 10 });

    n = 0;
    deleted_res1 = struct_resource{ -1 };
    struct_resource deleted_res2{ -1 };
    {
        boost::scope::unique_resource< struct_resource, checking_resource_deleter< struct_resource > > ur1{ struct_resource{ 10 }, checking_resource_deleter< struct_resource >(deleted_res1, n) };
        BOOST_TEST_EQ(ur1.get(), struct_resource{ 10 });
        BOOST_TEST(ur1.allocated());
        boost::scope::unique_resource< struct_resource, checking_resource_deleter< struct_resource > > ur2{ struct_resource{ 20 }, checking_resource_deleter< struct_resource >(deleted_res2, n) };
        BOOST_TEST_EQ(ur2.get(), struct_resource{ 20 });
        BOOST_TEST(ur2.allocated());
        ur2 = std::move(ur1);
        BOOST_TEST_EQ(n, 1);
        BOOST_TEST_EQ(deleted_res1, struct_resource{ -1 });
        BOOST_TEST_EQ(deleted_res2, struct_resource{ 20 });
        deleted_res2 = struct_resource{ -1 };
        BOOST_TEST_EQ(ur2.get(), struct_resource{ 10 });
        BOOST_TEST(ur2.allocated());
        BOOST_TEST(!ur1.allocated());
    }
    BOOST_TEST_EQ(n, 2);
    BOOST_TEST_EQ(deleted_res1, struct_resource{ 10 });
    BOOST_TEST_EQ(deleted_res2, struct_resource{ -1 });

    n = 0;
    deleted_res1 = struct_resource{ -1 };
    deleted_res2 = struct_resource{ -1 };
    {
        boost::scope::unique_resource< struct_resource, checking_resource_deleter< struct_resource > > ur1{ struct_resource{ 10 }, checking_resource_deleter< struct_resource >(deleted_res1, n) };
        BOOST_TEST_EQ(ur1.get(), struct_resource{ 10 });
        BOOST_TEST(ur1.allocated());
        boost::scope::unique_resource< struct_resource, checking_resource_deleter< struct_resource > > ur2{ struct_resource{ 20 }, checking_resource_deleter< struct_resource >(deleted_res2, n) };
        BOOST_TEST_EQ(ur2.get(), struct_resource{ 20 });
        BOOST_TEST(ur2.allocated());
        using namespace std;
        swap(ur1, ur2);
        BOOST_TEST_EQ(n, 0);
        BOOST_TEST_EQ(ur1.get(), struct_resource{ 20 });
        BOOST_TEST_EQ(ur1.get_deleter().get_deleted(), &deleted_res2);
        BOOST_TEST(ur1.allocated());
        BOOST_TEST_EQ(ur2.get(), struct_resource{ 10 });
        BOOST_TEST_EQ(ur2.get_deleter().get_deleted(), &deleted_res1);
        BOOST_TEST(ur2.allocated());
    }
    BOOST_TEST_EQ(n, 2);
    BOOST_TEST_EQ(deleted_res1, struct_resource{ 10 });
    BOOST_TEST_EQ(deleted_res2, struct_resource{ 20 });
}


void check_ptr()
{
    {
        boost::scope::unique_resource< struct_resource*, empty_resource_deleter< struct_resource* > > ur;
        BOOST_TEST_EQ(ur.get(), static_cast< struct_resource* >(nullptr));
        BOOST_TEST(!ur.allocated());
    }

    int n = 0;
    struct_resource res1{ 10 };
    struct_resource* deleted_res1 = nullptr;
    {
        boost::scope::unique_resource< struct_resource*, checking_resource_deleter< struct_resource* > > ur{ &res1, checking_resource_deleter< struct_resource* >(deleted_res1, n) };
        BOOST_TEST_EQ(ur.get(), &res1);
        BOOST_TEST_EQ(ur.get()->value, 10);
        BOOST_TEST_EQ(ur->value, 10);
        BOOST_TEST_EQ((*ur).value, 10);
        BOOST_TEST(ur.allocated());
    }
    BOOST_TEST_EQ(n, 1);
    BOOST_TEST_EQ(deleted_res1, &res1);
}


void check_ref()
{
    int n = 0;
    struct_resource res1{ 10 };
    struct_resource deleted_res1{ -1 };
    {
        boost::scope::unique_resource< struct_resource&, checking_resource_deleter< struct_resource > > ur{ res1, checking_resource_deleter< struct_resource >(deleted_res1, n) };
        BOOST_TEST_EQ(ur.get(), res1);
        BOOST_TEST_EQ(&ur.get(), &res1);
        BOOST_TEST(ur.allocated());
    }
    BOOST_TEST_EQ(n, 1);
    BOOST_TEST_EQ(deleted_res1, res1);

    n = 0;
    deleted_res1 = struct_resource{ -1 };
    struct_resource res2{ 20 };
    {
        boost::scope::unique_resource< struct_resource&, checking_resource_deleter< struct_resource > > ur{ res1, checking_resource_deleter< struct_resource >(deleted_res1, n) };
        BOOST_TEST_EQ(&ur.get(), &res1);
        BOOST_TEST(ur.allocated());
        ur.reset(res2);
        BOOST_TEST_EQ(n, 1);
        BOOST_TEST_EQ(deleted_res1, res1);
        deleted_res1 = struct_resource{ -1 };
        BOOST_TEST_EQ(&ur.get(), &res2);
        BOOST_TEST(ur.allocated());
    }
    BOOST_TEST_EQ(n, 2);
    BOOST_TEST_EQ(deleted_res1, res2);

    n = 0;
    deleted_res1 = struct_resource{ -1 };
    {
        boost::scope::unique_resource< struct_resource&, checking_resource_deleter< struct_resource > > ur1{ res1, checking_resource_deleter< struct_resource >(deleted_res1, n) };
        BOOST_TEST_EQ(&ur1.get(), &res1);
        BOOST_TEST(ur1.allocated());
        boost::scope::unique_resource< struct_resource&, checking_resource_deleter< struct_resource > > ur2{ std::move(ur1) };
        BOOST_TEST_EQ(&ur2.get(), &res1);
        BOOST_TEST(ur2.allocated());
        BOOST_TEST(!ur1.allocated());
    }
    BOOST_TEST_EQ(n, 1);
    BOOST_TEST_EQ(deleted_res1, res1);

    n = 0;
    deleted_res1 = struct_resource{ -1 };
    struct_resource deleted_res2{ -1 };
    {
        boost::scope::unique_resource< struct_resource&, checking_resource_deleter< struct_resource > > ur1{ res1, checking_resource_deleter< struct_resource >(deleted_res1, n) };
        BOOST_TEST_EQ(&ur1.get(), &res1);
        BOOST_TEST(ur1.allocated());
        boost::scope::unique_resource< struct_resource&, checking_resource_deleter< struct_resource > > ur2{ res2, checking_resource_deleter< struct_resource >(deleted_res2, n) };
        BOOST_TEST_EQ(&ur2.get(), &res2);
        BOOST_TEST(ur2.allocated());
        ur2 = std::move(ur1);
        BOOST_TEST_EQ(n, 1);
        BOOST_TEST_EQ(deleted_res1, struct_resource{ -1 });
        BOOST_TEST_EQ(deleted_res2, struct_resource{ 20 });
        deleted_res2 = struct_resource{ -1 };
        BOOST_TEST_EQ(&ur2.get(), &res1);
        BOOST_TEST(ur2.allocated());
        BOOST_TEST(!ur1.allocated());
    }
    BOOST_TEST_EQ(n, 2);
    BOOST_TEST_EQ(deleted_res1, struct_resource{ 10 });
    BOOST_TEST_EQ(deleted_res2, struct_resource{ -1 });

    n = 0;
    deleted_res1 = struct_resource{ -1 };
    deleted_res2 = struct_resource{ -1 };
    {
        boost::scope::unique_resource< struct_resource&, checking_resource_deleter< struct_resource > > ur1{ res1, checking_resource_deleter< struct_resource >(deleted_res1, n) };
        BOOST_TEST_EQ(&ur1.get(), &res1);
        BOOST_TEST(ur1.allocated());
        struct_resource expected_res2{ 20 };
        boost::scope::unique_resource< struct_resource&, checking_resource_deleter< struct_resource > > ur2{ res2, checking_resource_deleter< struct_resource >(deleted_res2, n) };
        BOOST_TEST_EQ(&ur2.get(), &res2);
        BOOST_TEST(ur2.allocated());
        using namespace std;
        swap(ur1, ur2);
        BOOST_TEST_EQ(n, 0);
        BOOST_TEST_EQ(&ur1.get(), &res2);
        BOOST_TEST_EQ(ur1.get_deleter().get_deleted(), &deleted_res2);
        BOOST_TEST(ur1.allocated());
        BOOST_TEST_EQ(&ur2.get(), &res1);
        BOOST_TEST_EQ(ur2.get_deleter().get_deleted(), &deleted_res1);
        BOOST_TEST(ur2.allocated());
    }
    BOOST_TEST_EQ(n, 2);
    BOOST_TEST_EQ(deleted_res1, struct_resource{ 10 });
    BOOST_TEST_EQ(deleted_res2, struct_resource{ 20 });
}


class throw_on_move_resource
{
private:
    int m_value;

public:
    explicit throw_on_move_resource(int v = 0) noexcept :
        m_value(v)
    {
    }

    throw_on_move_resource(throw_on_move_resource const&) = default;
    throw_on_move_resource& operator= (throw_on_move_resource const&) = default;

    throw_on_move_resource(throw_on_move_resource&&)
    {
        throw std::runtime_error("throw_on_move_resource move ctor");
    }
    throw_on_move_resource& operator= (throw_on_move_resource&&)
    {
        throw std::runtime_error("throw_on_move_resource move assignment");
    }

    int get() const noexcept
    {
        return m_value;
    }

    bool operator== (throw_on_move_resource const& that) const noexcept
    {
        return m_value == that.m_value;
    }

    bool operator!= (throw_on_move_resource const& that) const noexcept
    {
        return !operator==(that);
    }

    friend std::ostream& operator<< (std::ostream& strm, throw_on_move_resource const& res)
    {
        strm << "{ " << res.m_value << " }";
        return strm;
    }
};

class throw_on_copy_resource
{
private:
    int m_value;
    mutable bool m_throw;

public:
    explicit throw_on_copy_resource(int v = 0, bool do_throw = true) noexcept :
        m_value(v),
        m_throw(do_throw)
    {
    }

    throw_on_copy_resource(throw_on_copy_resource const& that) :
        m_value(that.m_value),
        m_throw(that.m_throw)
    {
        if (m_throw)
            throw std::runtime_error("throw_on_copy_resource copy ctor");
    }
    throw_on_copy_resource& operator= (throw_on_copy_resource const& that)
    {
        m_value = that.m_value;
        m_throw = that.m_throw;
        if (m_throw)
            throw std::runtime_error("throw_on_copy_resource copy assignment");
        return *this;
    }

    throw_on_copy_resource(throw_on_copy_resource&&) = delete;
    throw_on_copy_resource& operator= (throw_on_copy_resource&&) = delete;

    int get() const noexcept
    {
        return m_value;
    }

    void set_throw(bool do_throw) const noexcept
    {
        m_throw = do_throw;
    }

    bool operator== (throw_on_copy_resource const& that) const noexcept
    {
        return m_value == that.m_value;
    }

    bool operator!= (throw_on_copy_resource const& that) const noexcept
    {
        return !operator==(that);
    }

    friend std::ostream& operator<< (std::ostream& strm, throw_on_copy_resource const& res)
    {
        strm << "{ " << res.m_value << " }";
        return strm;
    }
};

void check_throw_resource()
{
    int n = 0;
    try
    {
        boost::scope::unique_resource< throw_on_move_resource, checking_resource_deleter< throw_on_move_resource > > ur{ throw_on_move_resource{ 10 }, checking_resource_deleter< throw_on_move_resource >(n) };
        BOOST_TEST_EQ(ur.get().get(), 10);
        BOOST_TEST(ur.allocated());
    }
    catch (...)
    {
        BOOST_ERROR("An exception is not expected to be thrown by throw_on_move_resource (copy ctor should be used)");
    }
    BOOST_TEST_EQ(n, 1);

    n = 0;
    try
    {
        boost::scope::unique_resource< throw_on_copy_resource, checking_resource_deleter< throw_on_copy_resource > > ur{ throw_on_copy_resource{ 10 }, checking_resource_deleter< throw_on_copy_resource >(n) };
        BOOST_ERROR("An exception is expected to be thrown by throw_on_copy_resource");
    }
    catch (...)
    {
    }
    BOOST_TEST_EQ(n, 1);

    n = 0;
    try
    {
        boost::scope::unique_resource< throw_on_copy_resource, checking_resource_deleter< throw_on_copy_resource > > ur1{ throw_on_copy_resource{ 10, false }, checking_resource_deleter< throw_on_copy_resource >(n) };
        ur1.get().set_throw(true);
        try
        {
            boost::scope::unique_resource< throw_on_copy_resource, checking_resource_deleter< throw_on_copy_resource > > ur2 = std::move(ur1);
            BOOST_ERROR("An exception is expected to be thrown by throw_on_copy_resource");
        }
        catch (...)
        {
            BOOST_TEST_EQ(n, 0);
            throw;
        }
    }
    catch (...)
    {
    }
    BOOST_TEST_EQ(n, 1);

    n = 0;
    try
    {
        boost::scope::unique_resource< throw_on_copy_resource, checking_resource_deleter< throw_on_copy_resource > > ur1{ throw_on_copy_resource{ 10, false }, checking_resource_deleter< throw_on_copy_resource >(n) };
        ur1.get().set_throw(true);
        try
        {
            boost::scope::unique_resource< throw_on_copy_resource, checking_resource_deleter< throw_on_copy_resource > > ur2{ throw_on_copy_resource{ 20, false }, checking_resource_deleter< throw_on_copy_resource >(n) };
            ur2 = std::move(ur1);
            BOOST_ERROR("An exception is expected to be thrown by throw_on_copy_resource");
        }
        catch (...)
        {
            BOOST_TEST_EQ(n, 1);
            throw;
        }
    }
    catch (...)
    {
    }
    BOOST_TEST_EQ(n, 2);
}

class moveable_resource
{
private:
    int m_value;

public:
    explicit moveable_resource(int value = 0) noexcept :
        m_value(value)
    {
    }

    moveable_resource(moveable_resource&& that) noexcept :
        m_value(that.m_value)
    {
        that.m_value = -1;
    }

    moveable_resource& operator= (moveable_resource&& that) noexcept
    {
        m_value = that.m_value;
        that.m_value = -1;
        return *this;
    }

    moveable_resource(moveable_resource const&) = delete;
    moveable_resource& operator= (moveable_resource const&) = delete;

    int get() const noexcept
    {
        return m_value;
    }

    bool operator== (moveable_resource const& that) const noexcept
    {
        return m_value == that.m_value;
    }

    bool operator!= (moveable_resource const& that) const noexcept
    {
        return !operator==(that);
    }

    friend std::ostream& operator<< (std::ostream& strm, moveable_resource const& res)
    {
        strm << "{ " << res.m_value << " }";
        return strm;
    }

    friend void copy_resource(moveable_resource const& from, moveable_resource& to)
    {
        to.m_value = from.m_value;
    }
};

class copyable_resource
{
private:
    int m_value;

public:
    explicit copyable_resource(int value = 0) noexcept :
        m_value(value)
    {
    }

    copyable_resource(copyable_resource const& that) noexcept :
        m_value(that.m_value)
    {
    }

    copyable_resource& operator= (copyable_resource const& that) noexcept
    {
        m_value = that.m_value;
        return *this;
    }

    int get() const noexcept
    {
        return m_value;
    }

    bool operator== (copyable_resource const& that) const noexcept
    {
        return m_value == that.m_value;
    }

    bool operator!= (copyable_resource const& that) const noexcept
    {
        return !operator==(that);
    }

    friend std::ostream& operator<< (std::ostream& strm, copyable_resource const& res)
    {
        strm << "{ " << res.m_value << " }";
        return strm;
    }
};

template< typename Resource >
class throwing_resource_deleter :
    public checking_resource_deleter< Resource >
{
private:
    using base_type = checking_resource_deleter< Resource >;

private:
    mutable bool m_throw;

public:
    explicit throwing_resource_deleter(int& n, bool do_throw = false) noexcept :
        base_type(n),
        m_throw(do_throw)
    {
    }

    explicit throwing_resource_deleter(Resource& deleted, int& n, bool do_throw = false) noexcept :
        base_type(deleted, n),
        m_throw(do_throw)
    {
    }

    throwing_resource_deleter(throwing_resource_deleter const& that) :
        base_type(static_cast< base_type const& >(that)),
        m_throw(that.m_throw)
    {
        if (m_throw)
            throw std::runtime_error("throwing_resource_deleter copy ctor");
    }

    throwing_resource_deleter& operator=(throwing_resource_deleter const& that)
    {
        *static_cast< base_type* >(this) = static_cast< base_type const& >(that);
        m_throw = that.m_throw;
        if (m_throw)
            throw std::runtime_error("throwing_resource_deleter copy assignment");
        return *this;
    }

    void set_throw(bool do_throw) const noexcept
    {
        m_throw = do_throw;
    }
};

template< typename Resource >
using default_resource_traits = void;

template< typename Resource >
struct wrapped_int_resource_traits
{
    static Resource make_default() noexcept
    {
        return Resource(-1);
    }

    static bool is_allocated(Resource const& res) noexcept
    {
        return res.get() >= 0;
    }
};

template< template< typename > class Traits >
void check_throw_deleter()
{
    int n = 0;
    {
        try
        {
            using unique_resource_t = boost::scope::unique_resource< copyable_resource, throwing_resource_deleter< copyable_resource >, Traits< copyable_resource > >;
            unique_resource_t ur1{ boost::scope::default_resource, throwing_resource_deleter< copyable_resource >(n) };
            ur1.get_deleter().set_throw(true);
            try
            {
                unique_resource_t ur2 = std::move(ur1);
                BOOST_ERROR("An exception is expected to be thrown by throwing_resource_deleter");
            }
            catch (...)
            {
                BOOST_TEST_EQ(n, 0);
                throw;
            }
        }
        catch (...)
        {
        }
        BOOST_TEST_EQ(n, 0);
    }

    n = 0;
    {
        copyable_resource deleted_res1;
        try
        {
            using unique_resource_t = boost::scope::unique_resource< copyable_resource, throwing_resource_deleter< copyable_resource >, Traits< copyable_resource > >;
            unique_resource_t ur1{ copyable_resource{ 10 }, throwing_resource_deleter< copyable_resource >(deleted_res1, n) };
            ur1.get_deleter().set_throw(true);
            try
            {
                unique_resource_t ur2 = std::move(ur1);
                BOOST_ERROR("An exception is expected to be thrown by throwing_resource_deleter");
            }
            catch (...)
            {
                BOOST_TEST_EQ(n, 0);
                throw;
            }
        }
        catch (...)
        {
        }
        BOOST_TEST_EQ(n, 1);
        BOOST_TEST_EQ(deleted_res1, copyable_resource{ 10 });
    }

    n = 0;
    {
        copyable_resource deleted_res1, deleted_res2;
        try
        {
            using unique_resource_t = boost::scope::unique_resource< copyable_resource, throwing_resource_deleter< copyable_resource >, Traits< copyable_resource > >;
            unique_resource_t ur1{ copyable_resource{ 10 }, throwing_resource_deleter< copyable_resource >(deleted_res1, n) };
            ur1.get_deleter().set_throw(true);
            try
            {
                unique_resource_t ur2{ copyable_resource{ 20 }, throwing_resource_deleter< copyable_resource >(deleted_res2, n) };
                ur2 = std::move(ur1);
                BOOST_ERROR("An exception is expected to be thrown by throwing_resource_deleter");
            }
            catch (...)
            {
                BOOST_TEST_EQ(n, 1);
                BOOST_TEST_EQ(deleted_res2, copyable_resource{ 20 });
                throw;
            }
        }
        catch (...)
        {
        }
        BOOST_TEST_EQ(n, 2);
        BOOST_TEST_EQ(deleted_res1, copyable_resource{ 10 });
    }

    n = 0;
    {
        try
        {
            using unique_resource_t = boost::scope::unique_resource< moveable_resource, throwing_resource_deleter< moveable_resource >, Traits< moveable_resource > >;
            unique_resource_t ur1{ boost::scope::default_resource, throwing_resource_deleter< moveable_resource >(n) };
            ur1.get_deleter().set_throw(true);
            try
            {
                unique_resource_t ur2 = std::move(ur1);
                BOOST_ERROR("An exception is expected to be thrown by throwing_resource_deleter");
            }
            catch (...)
            {
                BOOST_TEST_EQ(n, 0);
                throw;
            }
        }
        catch (...)
        {
        }
        BOOST_TEST_EQ(n, 0);
    }

    n = 0;
    {
        moveable_resource deleted_res1;
        try
        {
            using unique_resource_t = boost::scope::unique_resource< moveable_resource, throwing_resource_deleter< moveable_resource >, Traits< moveable_resource > >;
            unique_resource_t ur1{ moveable_resource{ 10 }, throwing_resource_deleter< moveable_resource >(deleted_res1, n) };
            ur1.get_deleter().set_throw(true);
            try
            {
                unique_resource_t ur2 = std::move(ur1);
                BOOST_ERROR("An exception is expected to be thrown by throwing_resource_deleter");
            }
            catch (...)
            {
                // The resource was moved from ur1, but the deleter could not have been moved. Then the resource was moved back to ur1.
                BOOST_TEST_EQ(n, 0);
                BOOST_TEST(ur1.allocated());
                BOOST_TEST_EQ(ur1.get(), moveable_resource{ 10 });
                throw;
            }
        }
        catch (...)
        {
        }
        BOOST_TEST_EQ(n, 1);
        BOOST_TEST_EQ(deleted_res1, moveable_resource{ 10 });
    }

    n = 0;
    {
        moveable_resource deleted_res1, deleted_res2;
        try
        {
            using unique_resource_t = boost::scope::unique_resource< moveable_resource, throwing_resource_deleter< moveable_resource >, Traits< moveable_resource > >;
            unique_resource_t ur1{ moveable_resource{ 10 }, throwing_resource_deleter< moveable_resource >(deleted_res1, n) };
            ur1.get_deleter().set_throw(true);
            try
            {
                unique_resource_t ur2{ moveable_resource{ 20 }, throwing_resource_deleter< moveable_resource >(deleted_res2, n) };
                ur2 = std::move(ur1);
                BOOST_ERROR("An exception is expected to be thrown by throwing_resource_deleter");
            }
            catch (...)
            {
                // First, ur2 was reset - the deleter was called on the resource 20. Then the deleter could not have been moved from ur1,
                // leaving ur1 intact.
                BOOST_TEST_EQ(n, 1);
                BOOST_TEST_EQ(deleted_res2, moveable_resource{ 20 });
                throw;
            }
        }
        catch (...)
        {
        }
        BOOST_TEST_EQ(n, 2);
        BOOST_TEST_EQ(deleted_res1, moveable_resource{ 10 });
        BOOST_TEST_EQ(deleted_res2, moveable_resource{ 20 });
    }
}

// Specially-crafted resource type that is:
// * nothrow move-constructible
// * NOT-nothrow move-assignable
// * has data layout which favors placing the "allocated" flag in unique_resource in its tail padding
class move_constructible_resource
{
public:
    // A special tag to construct "default" resource value
    struct default_tag { };

private:
    int m_n1;
    signed char m_n2;

public:
    constexpr move_constructible_resource() noexcept : m_n1(0), m_n2(0) { }

    // For compatibility with move_constructible_resource_traits::make_default()
    explicit constexpr move_constructible_resource(default_tag) noexcept : move_constructible_resource() { }

    explicit move_constructible_resource(int n1, signed char n2 = 0) noexcept : m_n1(n1), m_n2(n2) { }

    move_constructible_resource(move_constructible_resource&& that) noexcept : m_n1(that.m_n1), m_n2(that.m_n2)
    {
        that.m_n1 = 0;
        that.m_n2 = 0;
    }

    move_constructible_resource(move_constructible_resource const& that) : m_n1(that.m_n1), m_n2(that.m_n2)
    {
    }

    move_constructible_resource& operator= (move_constructible_resource&& that) // not noexcept
    {
        m_n1 = that.m_n1;
        m_n2 = that.m_n2;
        that.m_n1 = 0;
        that.m_n2 = 0;
        return *this;
    }

    move_constructible_resource& operator= (move_constructible_resource const& that)
    {
        m_n1 = that.m_n1;
        m_n2 = that.m_n2;
        return *this;
    }

    // For compatibility with move_constructible_resource_traits::make_default()
    move_constructible_resource& operator= (default_tag) noexcept
    {
        m_n1 = 0;
        m_n2 = 0;
        return *this;
    }

    bool operator== (move_constructible_resource const& that) const noexcept
    {
        return m_n1 == that.m_n1 && m_n2 == that.m_n2;
    }

    bool operator!= (move_constructible_resource const& that) const noexcept
    {
        return !operator==(that);
    }

    friend std::ostream& operator<< (std::ostream& strm, move_constructible_resource const& res)
    {
        strm << "{ " << res.m_n1 << ", " << static_cast< int >(res.m_n2) << " }";
        return strm;
    }

    friend void copy_resource(move_constructible_resource const& from, move_constructible_resource& to)
    {
        to.m_n1 = from.m_n1;
        to.m_n2 = from.m_n2;
    }
};

//! Resource traits for \c move_constructible_resource
struct move_constructible_resource_traits
{
    static move_constructible_resource::default_tag make_default() noexcept { return move_constructible_resource::default_tag{}; }
    static bool is_allocated(move_constructible_resource const& res) noexcept { return res != move_constructible_resource{}; }
};

void check_throw_deleter_move_constructible_resource()
{
    int n = 0;
    {
        move_constructible_resource deleted_res1;
        try
        {
            using unique_resource_t = boost::scope::unique_resource< move_constructible_resource, throwing_resource_deleter< move_constructible_resource >, move_constructible_resource_traits >;
            unique_resource_t ur1{ move_constructible_resource(10, 5), throwing_resource_deleter< move_constructible_resource >(deleted_res1, n) };
            ur1.get_deleter().set_throw(true);
            try
            {
                unique_resource_t ur2 = std::move(ur1);
                BOOST_ERROR("An exception is expected to be thrown by throwing_resource_deleter");
            }
            catch (...)
            {
                // The resource was moved from ur1, but the deleter could not have been moved. Then the resource was moved back to ur1.
                BOOST_TEST_EQ(n, 0);
                BOOST_TEST(ur1.allocated());
                BOOST_TEST_EQ(ur1.get(), move_constructible_resource(10, 5));
                throw;
            }
        }
        catch (...)
        {
        }
        BOOST_TEST_EQ(n, 1);
        BOOST_TEST_EQ(deleted_res1, move_constructible_resource(10, 5));
    }

    n = 0;
    {
        move_constructible_resource deleted_res1;
        try
        {
            // No resource traits to force using the "allocated" flag
            using unique_resource_t = boost::scope::unique_resource< move_constructible_resource, throwing_resource_deleter< move_constructible_resource > >;
            unique_resource_t ur1{ move_constructible_resource(10, 0), throwing_resource_deleter< move_constructible_resource >(deleted_res1, n) };
            ur1.get_deleter().set_throw(true);
            try
            {
                unique_resource_t ur2 = std::move(ur1);
                BOOST_ERROR("An exception is expected to be thrown by throwing_resource_deleter");
            }
            catch (...)
            {
                // The resource was moved from ur1, but the deleter could not have been moved. Then the resource was moved back to ur1.
                BOOST_TEST_EQ(n, 0);
                BOOST_TEST(ur1.allocated());
                BOOST_TEST_EQ(ur1.get(), move_constructible_resource(10, 0));
                throw;
            }
        }
        catch (...)
        {
        }
        BOOST_TEST_EQ(n, 1);
        BOOST_TEST_EQ(deleted_res1, move_constructible_resource(10, 0));
    }
}


void check_deduction()
{
    struct local
    {
        static void raw_func_deleter(int)
        {
            ++g_n;
        }

#if !defined(BOOST_SCOPE_NO_CXX17_NOEXCEPT_FUNCTION_TYPES)
        static void raw_func_deleter_noexcept(int) noexcept
        {
            ++g_n;
        }
#endif
    };

#if !defined(BOOST_NO_CXX17_DEDUCTION_GUIDES)
    {
        using expected_unique_resource_t = boost::scope::unique_resource< int, empty_resource_deleter< int > >;
        boost::scope::unique_resource ur{ 0, empty_resource_deleter< int >() };
        BOOST_TEST_TRAIT_SAME(decltype(ur), expected_unique_resource_t);
    }
    {
        using expected_unique_resource_t = boost::scope::unique_resource< struct_resource, empty_resource_deleter< struct_resource > >;
        boost::scope::unique_resource ur{ struct_resource(), empty_resource_deleter< struct_resource >() };
        BOOST_TEST_TRAIT_SAME(decltype(ur), expected_unique_resource_t);
    }
    {
        using expected_unique_resource_t = boost::scope::unique_resource< int, void(*)(int) >;
        boost::scope::unique_resource ur{ 0, local::raw_func_deleter };
        BOOST_TEST_TRAIT_SAME(decltype(ur), expected_unique_resource_t);
    }
#if !defined(BOOST_SCOPE_NO_CXX17_NOEXCEPT_FUNCTION_TYPES)
    {
        using expected_unique_resource_t = boost::scope::unique_resource< int, void(*)(int) noexcept >;
        boost::scope::unique_resource ur{ 0, local::raw_func_deleter_noexcept };
        BOOST_TEST_TRAIT_SAME(decltype(ur), expected_unique_resource_t);
    }
#endif
    {
        using expected_unique_resource_t = boost::scope::unique_resource< int, empty_resource_deleter< int > >;
        boost::scope::unique_resource ur1{ 0, empty_resource_deleter< int >() };
        BOOST_TEST_TRAIT_SAME(decltype(ur1), expected_unique_resource_t);
        boost::scope::unique_resource ur2 = std::move(ur1);
        BOOST_TEST_TRAIT_SAME(decltype(ur2), expected_unique_resource_t);
    }
#endif // !defined(BOOST_NO_CXX17_DEDUCTION_GUIDES)

    int n = 0, deleted_res = -1;
    {
        using expected_unique_resource_t = boost::scope::unique_resource< int, checking_resource_deleter< int > >;
        auto ur = boost::scope::make_unique_resource_checked(10, 0, checking_resource_deleter< int >(deleted_res, n));
        BOOST_TEST_TRAIT_SAME(decltype(ur), expected_unique_resource_t);
        BOOST_TEST_EQ(ur.get(), 10);
        BOOST_TEST(ur.allocated());
    }
    BOOST_TEST_EQ(n, 1);
    BOOST_TEST_EQ(deleted_res, 10);

    n = 0;
    deleted_res = -1;
    {
        using expected_unique_resource_t = boost::scope::unique_resource< int, checking_resource_deleter< int > >;
        auto ur = boost::scope::make_unique_resource_checked(0, 0, checking_resource_deleter< int >(deleted_res, n));
        BOOST_TEST_TRAIT_SAME(decltype(ur), expected_unique_resource_t);
        BOOST_TEST_EQ(ur.get(), 0);
        BOOST_TEST(!ur.allocated());
    }
    BOOST_TEST_EQ(n, 0);
    BOOST_TEST_EQ(deleted_res, -1);

    {
        using expected_unique_resource_t = boost::scope::unique_resource< int, void(*)(int) >;
        auto ur = boost::scope::make_unique_resource_checked(0, 0, local::raw_func_deleter);
        BOOST_TEST_TRAIT_SAME(decltype(ur), expected_unique_resource_t);
    }
#if !defined(BOOST_SCOPE_NO_CXX17_NOEXCEPT_FUNCTION_TYPES)
    {
        using expected_unique_resource_t = boost::scope::unique_resource< int, void(*)(int) noexcept >;
        auto ur = boost::scope::make_unique_resource_checked(0, 0, local::raw_func_deleter_noexcept);
        BOOST_TEST_TRAIT_SAME(decltype(ur), expected_unique_resource_t);
    }
#endif

    n = 0;
    try
    {
        // Not required to throw, but either way the deleter should not be called
        auto ur = boost::scope::make_unique_resource_checked(throw_on_copy_resource{ 0 }, throw_on_copy_resource{ 0 }, checking_resource_deleter< throw_on_copy_resource >(n));
    }
    catch (...)
    {
    }
    BOOST_TEST_EQ(n, 0);

    n = 0;
    try
    {
        auto ur = boost::scope::make_unique_resource_checked(0, 0, throwing_resource_deleter< int >(n, true));
        BOOST_ERROR("An exception is expected to be thrown by throwing_resource_deleter");
    }
    catch (...)
    {
    }
    BOOST_TEST_EQ(n, 0);
}

struct int_resource_traits
{
    static int make_default() noexcept
    {
        return -1;
    }

    static bool is_allocated(int res) noexcept
    {
        return res >= 0;
    }
};

void check_resource_traits()
{
    {
        boost::scope::unique_resource< int, empty_resource_deleter< int >, int_resource_traits > ur;
        BOOST_TEST_EQ(ur.get(), int_resource_traits::make_default());
        BOOST_TEST(!ur.allocated());
    }

    int n = 0, deleted_res1 = -1;
    {
        boost::scope::unique_resource< int, checking_resource_deleter< int >, int_resource_traits > ur{ -10, checking_resource_deleter< int >(deleted_res1, n) };
        BOOST_TEST_EQ(ur.get(), -10);
        BOOST_TEST(!ur.allocated());
        BOOST_TEST(!ur);
    }
    BOOST_TEST_EQ(n, 0);
    BOOST_TEST_EQ(deleted_res1, -1);

    n = 0;
    deleted_res1 = -1;
    {
        boost::scope::unique_resource< int, checking_resource_deleter< int >, int_resource_traits > ur{ 0, checking_resource_deleter< int >(deleted_res1, n) };
        BOOST_TEST_EQ(ur.get(), 0);
        BOOST_TEST(ur.allocated());
        BOOST_TEST(!!ur);
    }
    BOOST_TEST_EQ(n, 1);
    BOOST_TEST_EQ(deleted_res1, 0);

    n = 0;
    deleted_res1 = -1;
    {
        boost::scope::unique_resource< int, checking_resource_deleter< int >, int_resource_traits > ur{ 10, checking_resource_deleter< int >(deleted_res1, n) };
        BOOST_TEST_EQ(ur.get(), 10);
        BOOST_TEST(ur.allocated());
        ur.release();
        BOOST_TEST_EQ(ur.get(), int_resource_traits::make_default());
        BOOST_TEST(!ur.allocated());
    }
    BOOST_TEST_EQ(n, 0);
    BOOST_TEST_EQ(deleted_res1, -1);

    n = 0;
    deleted_res1 = -1;
    {
        boost::scope::unique_resource< int, checking_resource_deleter< int >, int_resource_traits > ur{ 10, checking_resource_deleter< int >(deleted_res1, n) };
        BOOST_TEST_EQ(ur.get(), 10);
        BOOST_TEST(ur.allocated());
        ur.reset();
        BOOST_TEST(!ur.allocated());
        BOOST_TEST_EQ(ur.get(), int_resource_traits::make_default());
        BOOST_TEST_EQ(n, 1);
        BOOST_TEST_EQ(deleted_res1, 10);
    }
    BOOST_TEST_EQ(n, 1);

    n = 0;
    deleted_res1 = -1;
    {
        boost::scope::unique_resource< int, checking_resource_deleter< int >, int_resource_traits > ur{ 10, checking_resource_deleter< int >(deleted_res1, n) };
        BOOST_TEST_EQ(ur.get(), 10);
        BOOST_TEST(ur.allocated());
        ur.reset(20);
        BOOST_TEST_EQ(n, 1);
        BOOST_TEST_EQ(deleted_res1, 10);
        deleted_res1 = -1;
        BOOST_TEST_EQ(ur.get(), 20);
        BOOST_TEST(ur.allocated());
    }
    BOOST_TEST_EQ(n, 2);
    BOOST_TEST_EQ(deleted_res1, 20);

    n = 0;
    deleted_res1 = -1;
    {
        boost::scope::unique_resource< int, checking_resource_deleter< int >, int_resource_traits > ur1{ 10, checking_resource_deleter< int >(deleted_res1, n) };
        BOOST_TEST_EQ(ur1.get(), 10);
        BOOST_TEST(ur1.allocated());
        boost::scope::unique_resource< int, checking_resource_deleter< int >, int_resource_traits > ur2{ std::move(ur1) };
        BOOST_TEST_EQ(ur2.get(), 10);
        BOOST_TEST(ur2.allocated());
        BOOST_TEST_EQ(ur1.get(), int_resource_traits::make_default());
        BOOST_TEST(!ur1.allocated());
    }
    BOOST_TEST_EQ(n, 1);
    BOOST_TEST_EQ(deleted_res1, 10);

    n = 0;
    deleted_res1 = -1;
    int deleted_res2 = -1;
    {
        boost::scope::unique_resource< int, checking_resource_deleter< int >, int_resource_traits > ur1{ 10, checking_resource_deleter< int >(deleted_res1, n) };
        BOOST_TEST_EQ(ur1.get(), 10);
        BOOST_TEST(ur1.allocated());
        boost::scope::unique_resource< int, checking_resource_deleter< int >, int_resource_traits > ur2{ 20, checking_resource_deleter< int >(deleted_res2, n) };
        BOOST_TEST_EQ(ur2.get(), 20);
        BOOST_TEST(ur2.allocated());
        ur2 = std::move(ur1);
        BOOST_TEST_EQ(ur2.get(), 10);
        BOOST_TEST(ur2.allocated());
        BOOST_TEST_EQ(ur1.get(), int_resource_traits::make_default());
        BOOST_TEST(!ur1.allocated());
        BOOST_TEST_EQ(n, 1);
        BOOST_TEST_EQ(deleted_res1, -1);
        BOOST_TEST_EQ(deleted_res2, 20);
        deleted_res2 = -1;
    }
    BOOST_TEST_EQ(n, 2);
    BOOST_TEST_EQ(deleted_res1, 10);
    BOOST_TEST_EQ(deleted_res2, -1);

    {
        boost::scope::unique_resource< int, empty_resource_deleter< int >, int_resource_traits > ur1;
        BOOST_TEST_EQ(ur1.get(), int_resource_traits::make_default());
        BOOST_TEST(!ur1.allocated());
        boost::scope::unique_resource< int, empty_resource_deleter< int >, int_resource_traits > ur2{ 10, empty_resource_deleter< int >() };
        BOOST_TEST_EQ(ur2.get(), 10);
        BOOST_TEST(ur2.allocated());
        using namespace std;
        swap(ur1, ur2);
        BOOST_TEST_EQ(ur1.get(), 10);
        BOOST_TEST(ur1.allocated());
        BOOST_TEST_EQ(ur2.get(), int_resource_traits::make_default());
        BOOST_TEST(!ur2.allocated());
    }

    n = 0;
    deleted_res1 = -1;
    deleted_res2 = -1;
    {
        boost::scope::unique_resource< int, checking_resource_deleter< int >, int_resource_traits > ur1{ 10, checking_resource_deleter< int >(deleted_res1, n) };
        BOOST_TEST_EQ(ur1.get(), 10);
        BOOST_TEST_EQ(ur1.get_deleter().get_deleted(), &deleted_res1);
        BOOST_TEST(ur1.allocated());
        boost::scope::unique_resource< int, checking_resource_deleter< int >, int_resource_traits > ur2{ 20, checking_resource_deleter< int >(deleted_res2, n) };
        BOOST_TEST_EQ(ur2.get(), 20);
        BOOST_TEST_EQ(ur2.get_deleter().get_deleted(), &deleted_res2);
        BOOST_TEST(ur2.allocated());
        using namespace std;
        swap(ur1, ur2);
        BOOST_TEST_EQ(n, 0);
        BOOST_TEST_EQ(ur1.get(), 20);
        BOOST_TEST_EQ(ur1.get_deleter().get_deleted(), &deleted_res2);
        BOOST_TEST(ur1.allocated());
        BOOST_TEST_EQ(ur2.get(), 10);
        BOOST_TEST_EQ(ur2.get_deleter().get_deleted(), &deleted_res1);
        BOOST_TEST(ur2.allocated());
    }
    BOOST_TEST_EQ(n, 2);
    BOOST_TEST_EQ(deleted_res1, 10);
    BOOST_TEST_EQ(deleted_res2, 20);
}

#if !defined(BOOST_NO_CXX17_FOLD_EXPRESSIONS) && !defined(BOOST_NO_CXX17_AUTO_NONTYPE_TEMPLATE_PARAMS)

struct global_deleter_int
{
    void operator() (int) noexcept
    {
        ++g_n;
    }
};

struct global_deleter_ptr
{
    void operator() (int*) noexcept
    {
        ++g_n;
    }
};

void check_simple_resource_traits()
{
    g_n = 0;
    {
        boost::scope::unique_resource< int, global_deleter_int, boost::scope::unallocated_resource< -1 > > ur;
        BOOST_TEST_EQ(ur.get(), -1);
        BOOST_TEST(!ur.allocated());
        BOOST_TEST(!ur);
    }
    BOOST_TEST_EQ(g_n, 0);

    g_n = 0;
    {
        boost::scope::unique_resource< int, global_deleter_int, boost::scope::unallocated_resource< -1 > > ur{ -1 };
        BOOST_TEST_EQ(ur.get(), -1);
        BOOST_TEST(!ur.allocated());
        BOOST_TEST(!ur);
    }
    BOOST_TEST_EQ(g_n, 0);

    g_n = 0;
    {
        boost::scope::unique_resource< int, global_deleter_int, boost::scope::unallocated_resource< -1 > > ur{ 10 };
        BOOST_TEST_EQ(ur.get(), 10);
        BOOST_TEST(ur.allocated());
        BOOST_TEST(!!ur);
    }
    BOOST_TEST_EQ(g_n, 1);

    g_n = 0;
    {
        boost::scope::unique_resource< int, global_deleter_int, boost::scope::unallocated_resource< -1, -2, -3 > > ur;
        BOOST_TEST_EQ(ur.get(), -1);
        BOOST_TEST(!ur.allocated());
        BOOST_TEST(!ur);

        ur.reset(-2);
        BOOST_TEST_EQ(g_n, 0);
        BOOST_TEST_EQ(ur.get(), -2);
        BOOST_TEST(!ur.allocated());
        BOOST_TEST(!ur);

        ur.reset(-3);
        BOOST_TEST_EQ(g_n, 0);
        BOOST_TEST_EQ(ur.get(), -3);
        BOOST_TEST(!ur.allocated());
        BOOST_TEST(!ur);

        ur.reset(10);
        BOOST_TEST_EQ(g_n, 0);
        BOOST_TEST_EQ(ur.get(), 10);
        BOOST_TEST(ur.allocated());
        BOOST_TEST(!!ur);
    }
    BOOST_TEST_EQ(g_n, 1);

    g_n = 0;
    {
        boost::scope::unique_resource< int*, global_deleter_ptr, boost::scope::unallocated_resource< nullptr > > ur;
        BOOST_TEST_EQ(ur.get(), static_cast< int* >(nullptr));
        BOOST_TEST(!ur.allocated());
        BOOST_TEST(!ur);
    }
    BOOST_TEST_EQ(g_n, 0);

    g_n = 0;
    {
        boost::scope::unique_resource< int*, global_deleter_ptr, boost::scope::unallocated_resource< nullptr > > ur{ &g_n };
        BOOST_TEST_EQ(ur.get(), &g_n);
        BOOST_TEST(ur.allocated());
        BOOST_TEST(!!ur);
    }
    BOOST_TEST_EQ(g_n, 1);
}

#endif // !defined(BOOST_NO_CXX17_FOLD_EXPRESSIONS) && !defined(BOOST_NO_CXX17_AUTO_NONTYPE_TEMPLATE_PARAMS)

int main()
{
    check_int();
    check_struct();
    check_ptr();
    check_ref();
    check_throw_resource();
    check_throw_deleter< default_resource_traits >();
    check_throw_deleter< wrapped_int_resource_traits >();
    check_throw_deleter_move_constructible_resource();
    check_deduction();
    check_resource_traits();
#if !defined(BOOST_NO_CXX17_FOLD_EXPRESSIONS) && !defined(BOOST_NO_CXX17_AUTO_NONTYPE_TEMPLATE_PARAMS)
    check_simple_resource_traits();
#endif
    return boost::report_errors();
}
