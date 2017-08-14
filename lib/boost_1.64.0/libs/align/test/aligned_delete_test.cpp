/*
Copyright 2014 Glen Joseph Fernandes
(glenjofe@gmail.com)

Distributed under the Boost Software License, Version 1.0.
(http://www.boost.org/LICENSE_1_0.txt)
*/
#include <boost/align/aligned_alloc.hpp>
#include <boost/align/aligned_delete.hpp>
#include <boost/align/alignment_of.hpp>
#include <boost/core/lightweight_test.hpp>
#include <new>

template<class T>
class type {
public:
    static unsigned count;

    type()
        : value_() {
        ++count;
    }

    ~type() {
        --count;
    }

private:
    T value_;
};

template<class T>
unsigned type<T>::count = 0;

template<class T>
void test()
{
    typedef type<T> E;
    void* p = boost::alignment::aligned_alloc(boost::
        alignment::alignment_of<E>::value, sizeof(E));
    BOOST_TEST(p != 0);
    E* q = ::new(p) E;
    BOOST_TEST(E::count == 1);
    boost::alignment::aligned_delete()(q);
    BOOST_TEST(E::count == 0);
}

class C { };
union U { };

int main()
{
    test<char>();
    test<bool>();
    test<short>();
    test<int>();
    test<long>();
    test<float>();
    test<double>();
    test<long double>();
    test<void*>();
    test<void(*)()>();
    test<C>();
    test<int C::*>();
    test<int (C::*)()>();
    test<U>();

    return boost::report_errors();
}
