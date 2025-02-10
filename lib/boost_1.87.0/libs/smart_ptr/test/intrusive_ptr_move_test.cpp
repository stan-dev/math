//
//  intrusive_ptr_move_test.cpp
//
//  Copyright (c) 2002-2005 Peter Dimov
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/intrusive_ptr.hpp>
#include <boost/smart_ptr/detail/atomic_count.hpp>
#include <boost/core/lightweight_test.hpp>
#include <boost/config.hpp>
#include <utility>

namespace N
{

class base
{
private:

    mutable boost::detail::atomic_count use_count_;

    base(base const &);
    base & operator=(base const &);

protected:

    base(): use_count_(0)
    {
        ++instances;
    }

    virtual ~base()
    {
        --instances;
    }

public:

    static long instances;

    long use_count() const
    {
        return use_count_;
    }

    inline friend void intrusive_ptr_add_ref(base const * p)
    {
        ++p->use_count_;
    }

    inline friend void intrusive_ptr_release(base const * p)
    {
        if(--p->use_count_ == 0) delete p;
    }
};

long base::instances = 0;

} // namespace N

//

struct X: public virtual N::base
{
};

struct Y: public X
{
};

int main()
{
    BOOST_TEST( N::base::instances == 0 );

    {
        boost::intrusive_ptr<X> p( new X );
        BOOST_TEST( N::base::instances == 1 );

        boost::intrusive_ptr<X> p2( std::move( p ) );
        BOOST_TEST( N::base::instances == 1 );
        BOOST_TEST( p.get() == 0 );

        p2.reset();
        BOOST_TEST( N::base::instances == 0 );
    }

    {
        boost::intrusive_ptr<Y> p( new Y );
        BOOST_TEST( N::base::instances == 1 );

        boost::intrusive_ptr<X> p2( std::move( p ) );
        BOOST_TEST( N::base::instances == 1 );
        BOOST_TEST( p.get() == 0 );

        p2.reset();
        BOOST_TEST( N::base::instances == 0 );
    }

    {
        boost::intrusive_ptr<X> p( new X );
        BOOST_TEST( N::base::instances == 1 );

        boost::intrusive_ptr<X> p2;
        p2 = std::move( p );
        BOOST_TEST( N::base::instances == 1 );
        BOOST_TEST( p.get() == 0 );

        p2.reset();
        BOOST_TEST( N::base::instances == 0 );
    }

    {
        boost::intrusive_ptr<X> p( new X );
        BOOST_TEST( N::base::instances == 1 );

        boost::intrusive_ptr<X> p2( new X );
        BOOST_TEST( N::base::instances == 2 );
        p2 = std::move( p );
        BOOST_TEST( N::base::instances == 1 );
        BOOST_TEST( p.get() == 0 );

        p2.reset();
        BOOST_TEST( N::base::instances == 0 );
    }

    {
        boost::intrusive_ptr<Y> p( new Y );
        BOOST_TEST( N::base::instances == 1 );

        boost::intrusive_ptr<X> p2;
        p2 = std::move( p );
        BOOST_TEST( N::base::instances == 1 );
        BOOST_TEST( p.get() == 0 );

        p2.reset();
        BOOST_TEST( N::base::instances == 0 );
    }

    {
        boost::intrusive_ptr<Y> p( new Y );
        BOOST_TEST( N::base::instances == 1 );

        boost::intrusive_ptr<X> p2( new X );
        BOOST_TEST( N::base::instances == 2 );
        p2 = std::move( p );
        BOOST_TEST( N::base::instances == 1 );
        BOOST_TEST( p.get() == 0 );

        p2.reset();
        BOOST_TEST( N::base::instances == 0 );
    }

    {
        boost::intrusive_ptr<X> px( new Y );

        X * px2 = px.get();

        boost::intrusive_ptr<Y> py = boost::static_pointer_cast<Y>( std::move( px ) );
        BOOST_TEST( py.get() == px2 );
        BOOST_TEST( px.get() == 0 );
        BOOST_TEST( py->use_count() == 1 );
    }

    BOOST_TEST( N::base::instances == 0 );

    {
        boost::intrusive_ptr<X const> px( new X );

        X const * px2 = px.get();

        boost::intrusive_ptr<X> px3 = boost::const_pointer_cast<X>( std::move( px ) );
        BOOST_TEST( px3.get() == px2 );
        BOOST_TEST( px.get() == 0 );
        BOOST_TEST( px3->use_count() == 1 );
    }

    BOOST_TEST( N::base::instances == 0 );

    {
        boost::intrusive_ptr<X> px( new Y );

        X * px2 = px.get();

        boost::intrusive_ptr<Y> py = boost::dynamic_pointer_cast<Y>( std::move( px ) );
        BOOST_TEST( py.get() == px2 );
        BOOST_TEST( px.get() == 0 );
        BOOST_TEST( py->use_count() == 1 );
    }

    BOOST_TEST( N::base::instances == 0 );

    {
        boost::intrusive_ptr<X> px( new X );

        X * px2 = px.get();

        boost::intrusive_ptr<Y> py = boost::dynamic_pointer_cast<Y>( std::move( px ) );
        BOOST_TEST( py.get() == 0 );
        BOOST_TEST( px.get() == px2 );
        BOOST_TEST( px->use_count() == 1 );
    }

    BOOST_TEST( N::base::instances == 0 );

    return boost::report_errors();
}
