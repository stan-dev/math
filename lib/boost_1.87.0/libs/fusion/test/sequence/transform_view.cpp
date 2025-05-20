/*=============================================================================
    Copyright (c) 2001-2011 Joel de Guzman
    Copyright (c) 2022 Denis Mikhailov

    Distributed under the Boost Software License, Version 1.0. (See accompanying 
    file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#include <boost/detail/lightweight_test.hpp>
#include <boost/fusion/container/vector/vector.hpp>
#include <boost/fusion/container/map/map.hpp>
#include <boost/fusion/container/generation/make_map.hpp>
#include <boost/fusion/adapted/mpl.hpp>
#include <boost/fusion/sequence/io/out.hpp>
#include <boost/fusion/container/generation/make_vector.hpp>
#include <boost/fusion/sequence/comparison/equal_to.hpp>
#include <boost/fusion/view/transform_view/transform_view.hpp>
#include <boost/fusion/sequence/intrinsic/begin.hpp>
#include <boost/fusion/sequence/intrinsic/at.hpp>
#include <boost/fusion/sequence/intrinsic/value_at.hpp>
#include <boost/fusion/sequence/intrinsic/at_key.hpp>
#include <boost/fusion/sequence/intrinsic/value_at_key.hpp>
#include <boost/fusion/sequence/intrinsic/has_key.hpp>
#include <boost/fusion/sequence/intrinsic/begin.hpp>
#include <boost/fusion/iterator/next.hpp>
#include <boost/fusion/iterator/prior.hpp>
#include <boost/fusion/iterator/advance.hpp>
#include <boost/fusion/iterator/deref.hpp>
#include <boost/fusion/iterator/distance.hpp>
#include <boost/fusion/iterator/key_of.hpp>
#include <boost/fusion/iterator/deref_data.hpp>
#include <boost/fusion/iterator/value_of_data.hpp>
#include <boost/fusion/support/pair.hpp>
#include <boost/fusion/support/category_of.hpp>

#include <boost/mpl/range_c.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/mpl/at.hpp>
#include <boost/mpl/map.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/static_assert.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/type_traits/remove_reference.hpp>

struct square
{
    template<typename T>
    struct result;

    template <typename T>
    struct result<square(T)>
    {
        typedef int type;
    };

    template <typename T>
    int operator()(T x) const
    {
        return x * x;
    }
};

struct add
{
    template<typename T>
    struct result;

    template <typename A, typename B>
    struct result<add(A,B)>
    {
        typedef int type;
    };

    template <typename A, typename B>
    int operator()(A a, B b) const
    {
        return a + b;
    }
};

struct abstract
{
    virtual void foo() = 0;
};

struct functor
{
    typedef boost::fusion::pair<short, char> pair_0;
    typedef boost::fusion::pair<double, std::string> pair_1;
    typedef boost::fusion::pair<abstract, long> pair_2;

    typedef boost::mpl::map<
        boost::mpl::pair< boost::fusion::pair<int, char> , pair_0>
      , boost::mpl::pair< pair_1 , pair_1>
      , boost::mpl::pair< boost::fusion::pair<abstract, int> , pair_2>
    > m;

    template<typename Sig>
    struct result;

    template<class F, class T>
    struct result<F(T)> 
        : boost::mpl::at<
              m,
              typename boost::remove_reference<T>::type
          >
    {};

    pair_0 operator() (const boost::fusion::pair<int, char>& arg) const
    {
        return pair_0(arg.second);
    }
    
    pair_1 operator() (const pair_1 & arg) const
    {
        return pair_1(arg.second + "_transformed");
    }

    pair_2 operator() (const boost::fusion::pair<abstract, int>& arg) const
    {
        return pair_2(arg.second);
    }
};

struct simple_identity
{
    template<typename Sig>
    struct result;

    template<typename U>
    struct result<simple_identity(U)>
    {
        typedef U type;
    };

    template<typename T>
    T& operator() (T& arg) const 
    {
        return arg;
    }
};

struct simple_identity_nonref
{
    template<typename Sig>
    struct result;

    template<typename U>
    struct result<simple_identity_nonref(U)>
    : boost::remove_reference<U>
    {};

    template<typename T>
    T operator() (T arg) const 
    {
        return arg;
    }
};

int
main()
{
    using namespace boost::fusion;

    std::cout << tuple_open('[');
    std::cout << tuple_close(']');
    std::cout << tuple_delimiter(", ");

/// Testing the transform_view

    {
        typedef boost::mpl::range_c<int, 5, 9> sequence_type;
        sequence_type sequence;
        square sq;
        typedef transform_view<sequence_type, square> xform_type;
        xform_type xform(sequence, sq);

        std::cout << xform << std::endl;
        BOOST_TEST((xform == make_vector(25, 36, 49, 64)));

        typedef boost::fusion::result_of::begin<xform_type>::type first_type;
        first_type first_it(boost::fusion::begin(xform));

        typedef boost::fusion::result_of::next<first_type>::type next_type;
        next_type next_it(boost::fusion::next(first_it));
        BOOST_TEST((*next_it == 36));
        BOOST_TEST((*boost::fusion::prior(next_it) == 25));
        BOOST_TEST((boost::fusion::distance(first_it, next_it) == 1));

        BOOST_TEST((*boost::fusion::advance_c<3>(boost::fusion::begin(xform)) == 64));
        BOOST_TEST((boost::fusion::at_c<2>(xform) == 49));
        BOOST_MPL_ASSERT((boost::is_same<boost::fusion::result_of::value_at_c<xform_type, 0>::type, int>));
    }
    
    {
        typedef boost::mpl::range_c<int, 5, 9> sequence1_type;
        typedef boost::mpl::range_c<int, 10, 14> sequence2_type;
        sequence1_type sequence1;
        sequence2_type sequence2;
        add f;
        typedef transform_view<sequence1_type, sequence2_type, add> xform_type;
        xform_type xform(sequence1, sequence2, f);

        std::cout << xform << std::endl;
        BOOST_TEST((xform == make_vector(15, 17, 19, 21)));
        BOOST_TEST((boost::fusion::at_c<2>(xform) == 19));
        BOOST_MPL_ASSERT((boost::is_same<boost::fusion::result_of::value_at_c<xform_type, 0>::type, int>));
    }

    /// Associative
    {
        typedef map<
            pair<int, char>
          , pair<double, std::string>
          , pair<abstract, int> >
        map_type;
        typedef transform_view<map_type, functor> transformed_type;

        BOOST_MPL_ASSERT((traits::is_associative<transformed_type>));
        BOOST_MPL_ASSERT((traits::is_random_access<transformed_type>));

        map_type m(
            make_pair<int>('X')
          , make_pair<double>("Men")
          , make_pair<abstract>(2));
        transformed_type t(m, functor());

        std::cout << at_key<short>(t) << std::endl;
        std::cout << at_key<double>(t) << std::endl;
        std::cout << at_key<abstract>(t) << std::endl;

        BOOST_TEST(at_key<short>(t) == 'X');
        BOOST_TEST(at_key<double>(t) == "Men_transformed");
        BOOST_TEST(at_key<abstract>(t) == 2);

        BOOST_STATIC_ASSERT((
            boost::is_same<boost::fusion::result_of::value_at_key<transformed_type, short>::type, char>::value));
        BOOST_STATIC_ASSERT((
            boost::is_same<boost::fusion::result_of::value_at_key<transformed_type, double>::type, std::string>::value));
        BOOST_STATIC_ASSERT((
            boost::is_same<boost::fusion::result_of::value_at_key<transformed_type, abstract>::type, long>::value));

        std::cout << t << std::endl;

        BOOST_STATIC_ASSERT((boost::fusion::result_of::has_key<transformed_type, short>::value));
        BOOST_STATIC_ASSERT((boost::fusion::result_of::has_key<transformed_type, double>::value));
        BOOST_STATIC_ASSERT((boost::fusion::result_of::has_key<transformed_type, abstract>::value));
        BOOST_STATIC_ASSERT((!boost::fusion::result_of::has_key<transformed_type, std::string>::value));

        std::cout << deref_data(begin(t)) << std::endl;
        std::cout << deref_data(boost::fusion::next(begin(t))) << std::endl;

        BOOST_TEST(deref_data(begin(t)) == 'X');
        BOOST_TEST(deref_data(boost::fusion::next(begin(t))) == "Men_transformed");
        BOOST_TEST(deref_data(boost::fusion::next(next(begin(t)))) == 2);

        BOOST_STATIC_ASSERT((boost::is_same<boost::fusion::result_of::key_of<boost::fusion::result_of::begin<transformed_type>::type>::type, short>::value));
        BOOST_STATIC_ASSERT((boost::is_same<boost::fusion::result_of::key_of<boost::fusion::result_of::next<boost::fusion::result_of::begin<transformed_type>::type>::type>::type, double>::value));
        BOOST_STATIC_ASSERT((boost::is_same<boost::fusion::result_of::key_of<boost::fusion::result_of::next<boost::fusion::result_of::next<boost::fusion::result_of::begin<transformed_type>::type>::type>::type>::type, abstract>::value));
        BOOST_STATIC_ASSERT((boost::is_same<boost::fusion::result_of::value_of_data<boost::fusion::result_of::begin<transformed_type>::type>::type, char>::value));
        BOOST_STATIC_ASSERT((boost::is_same<boost::fusion::result_of::value_of_data<boost::fusion::result_of::next<boost::fusion::result_of::begin<transformed_type>::type>::type>::type, std::string>::value));
        BOOST_STATIC_ASSERT((boost::is_same<boost::fusion::result_of::value_of_data<boost::fusion::result_of::next<boost::fusion::result_of::next<boost::fusion::result_of::begin<transformed_type>::type>::type>::type>::type, long>::value));

        // Test random access interface.
        pair<short, char> a = at_c<0>(t); (void) a;
        pair<double, std::string> b = at_c<1>(t);
        pair<abstract, int> c = at_c<2>(t);
        (void)c;

        typedef boost::fusion::result_of::begin<transformed_type>::type first;
        typedef boost::fusion::result_of::next<first>::type second;
        typedef boost::fusion::result_of::next<second>::type third;

        BOOST_MPL_ASSERT((boost::is_same<boost::fusion::result_of::value_of<first>::type, boost::fusion::pair<short, char> >));
        BOOST_MPL_ASSERT((boost::is_same<boost::fusion::result_of::value_of<second>::type, boost::fusion::pair<double, std::string> >));
        BOOST_MPL_ASSERT((boost::is_same<boost::fusion::result_of::value_of<third>::type, boost::fusion::pair<abstract, long> >));
    }

    {
        typedef map<
            pair<int, char>
          , pair<double, std::string>
          , pair<abstract, int> >
        map_type;
        typedef transform_view<map_type, simple_identity> transformed_type;

        BOOST_MPL_ASSERT((traits::is_associative<transformed_type>));
        BOOST_MPL_ASSERT((traits::is_random_access<transformed_type>));

        map_type m(
            make_pair<int>('X')
          , make_pair<double>("Men")
          , make_pair<abstract>(2));
        transformed_type t(m, simple_identity());

        std::cout << at_key<int>(t) << std::endl;
        std::cout << at_key<double>(t) << std::endl;
        std::cout << at_key<abstract>(t) << std::endl;

        BOOST_TEST(at_key<int>(t) == 'X');
        BOOST_TEST(at_key<double>(t) == "Men");
        BOOST_TEST(at_key<abstract>(t) == 2);

        BOOST_STATIC_ASSERT((
            boost::is_same<boost::fusion::result_of::value_at_key<transformed_type, int>::type, char>::value));
        BOOST_STATIC_ASSERT((
            boost::is_same<boost::fusion::result_of::value_at_key<transformed_type, double>::type, std::string>::value));
        BOOST_STATIC_ASSERT((
            boost::is_same<boost::fusion::result_of::value_at_key<transformed_type, abstract>::type, int>::value));

        std::cout << t << std::endl;

        BOOST_STATIC_ASSERT((boost::fusion::result_of::has_key<transformed_type, int>::value));
        BOOST_STATIC_ASSERT((boost::fusion::result_of::has_key<transformed_type, double>::value));
        BOOST_STATIC_ASSERT((boost::fusion::result_of::has_key<transformed_type, abstract>::value));
        BOOST_STATIC_ASSERT((!boost::fusion::result_of::has_key<transformed_type, std::string>::value));

        std::cout << deref_data(begin(t)) << std::endl;
        std::cout << deref_data(boost::fusion::next(begin(t))) << std::endl;

        BOOST_TEST(deref_data(begin(t)) == 'X');
        BOOST_TEST(deref_data(boost::fusion::next(begin(t))) == "Men");
        BOOST_TEST(deref_data(boost::fusion::next(next(begin(t)))) == 2);

        BOOST_STATIC_ASSERT((boost::is_same<boost::fusion::result_of::key_of<boost::fusion::result_of::begin<transformed_type>::type>::type, int>::value));
        BOOST_STATIC_ASSERT((boost::is_same<boost::fusion::result_of::key_of<boost::fusion::result_of::next<boost::fusion::result_of::begin<transformed_type>::type>::type>::type, double>::value));
        BOOST_STATIC_ASSERT((boost::is_same<boost::fusion::result_of::key_of<boost::fusion::result_of::next<boost::fusion::result_of::next<boost::fusion::result_of::begin<transformed_type>::type>::type>::type>::type, abstract>::value));
        BOOST_STATIC_ASSERT((boost::is_same<boost::fusion::result_of::value_of_data<boost::fusion::result_of::begin<transformed_type>::type>::type, char>::value));
        BOOST_STATIC_ASSERT((boost::is_same<boost::fusion::result_of::value_of_data<boost::fusion::result_of::next<boost::fusion::result_of::begin<transformed_type>::type>::type>::type, std::string>::value));
        BOOST_STATIC_ASSERT((boost::is_same<boost::fusion::result_of::value_of_data<boost::fusion::result_of::next<boost::fusion::result_of::next<boost::fusion::result_of::begin<transformed_type>::type>::type>::type>::type, int>::value));

        // Test random access interface.
        pair<int, char> a = at_c<0>(t); (void) a;
        pair<double, std::string> b = at_c<1>(t);
        pair<abstract, int> c = at_c<2>(t);
        (void)c;

        typedef boost::fusion::result_of::begin<transformed_type>::type first;
        typedef boost::fusion::result_of::next<first>::type second;
        typedef boost::fusion::result_of::next<second>::type third;

        // BOOST_MPL_ASSERT((boost::is_same<boost::fusion::result_of::value_of<first>::type, boost::fusion::pair<int, char> >));
        // BOOST_MPL_ASSERT((boost::is_same<boost::fusion::result_of::value_of<second>::type, boost::fusion::pair<double, std::string> >));
        // BOOST_MPL_ASSERT((boost::is_same<boost::fusion::result_of::value_of<third>::type, boost::fusion::pair<abstract, int> >));
    }

    {
        // compile test only
        // make sure result_of::deref_data returns a reference
        typedef map<pair<float, int> > map_type;
        typedef transform_view<map_type, simple_identity> transformed_type;
        typedef boost::fusion::result_of::begin<transformed_type>::type i_type;
        typedef boost::fusion::result_of::deref_data<i_type>::type r_type;
        BOOST_STATIC_ASSERT((boost::is_same<r_type, int&>::value));
    }

    {
        // compile test only
        // make sure result_of::deref_data is const correct
        typedef map<pair<float, int> > const map_type;
        typedef transform_view<map_type, simple_identity> transformed_type;
        typedef boost::fusion::result_of::begin<transformed_type>::type i_type;
        typedef boost::fusion::result_of::deref_data<i_type>::type r_type;
        BOOST_STATIC_ASSERT((boost::is_same<r_type, int const&>::value));
    }

    {
        // compile test only
        // make sure result_of::deref_data will not return a reference to temp object
        typedef map<pair<float, int> > const map_type;
        typedef transform_view<map_type, simple_identity_nonref> transformed_type;
        typedef boost::fusion::result_of::begin<transformed_type>::type i_type;
        typedef boost::fusion::result_of::deref_data<i_type>::type r_type;
        BOOST_STATIC_ASSERT((boost::is_same<r_type, int>::value));
    }

    {
        // compile test only
        // make sure result_of::deref_data will not const for non-constant references
        typedef map<pair<float, int&> > const map_type;
        typedef transform_view<map_type, simple_identity> transformed_type;
        typedef boost::fusion::result_of::begin<transformed_type>::type i_type;
        typedef boost::fusion::result_of::deref_data<i_type>::type r_type;
        BOOST_STATIC_ASSERT((boost::is_same<r_type, int&>::value));
    }

    return boost::report_errors();
}

