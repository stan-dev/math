/*=============================================================================
    Copyright (c) 2022 Denis Mikhailov

    Distributed under the Boost Software License, Version 1.0. (See accompanying 
    file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#include <boost/detail/lightweight_test.hpp>
#include <boost/fusion/container/vector/vector.hpp>
#include <boost/fusion/container/map/map.hpp>
#include <boost/fusion/container/generation/make_map.hpp>
#include <boost/fusion/algorithm/auxiliary/copy.hpp>
#include <boost/fusion/adapted/mpl.hpp>
#include <boost/fusion/sequence/io/out.hpp>
#include <boost/fusion/container/generation/make_vector.hpp>
#include <boost/fusion/sequence/comparison/equal_to.hpp>
#include <boost/fusion/view/identity_view/identity_view.hpp>
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

struct abstract
{
    virtual void foo() = 0;
};

int
main()
{
    using namespace boost::fusion;

    std::cout << tuple_open('[');
    std::cout << tuple_close(']');
    std::cout << tuple_delimiter(", ");

/// Testing the identity_view

    {
        typedef boost::mpl::range_c<int, 5, 9> sequence_type;
        sequence_type sequence;
        typedef identity_view<sequence_type> xform_type;
        xform_type xform(sequence);

        std::cout << xform << std::endl;
        BOOST_TEST((xform == make_vector(5, 6, 7, 8))); 

        typedef boost::fusion::result_of::begin<xform_type>::type first_type;
        first_type first_it(boost::fusion::begin(xform));

        typedef boost::fusion::result_of::next<first_type>::type next_type;
        next_type next_it(boost::fusion::next(first_it));
        BOOST_TEST((*next_it == 6));
        BOOST_TEST((*boost::fusion::prior(next_it) == 5));
        BOOST_TEST((boost::fusion::distance(first_it, next_it) == 1));

        BOOST_TEST((*boost::fusion::advance_c<3>(boost::fusion::begin(xform)) == 8));
        BOOST_TEST((boost::fusion::at_c<2>(xform) == 7));
        BOOST_MPL_ASSERT((boost::is_same<boost::fusion::result_of::value_at_c<xform_type, 0>::type,  boost::mpl::integral_c<int, 5> >));
    }
    
    {
        typedef vector<int, int, int, int, int> sequence_type;
        sequence_type seq;
        identity_view<sequence_type> ident(seq);
        copy(make_vector(1, 2, 3, 4, 5), ident);
        std::cout << seq << std::endl;
        BOOST_TEST((seq == make_vector(1, 2, 3, 4, 5)));
    }
    
    /// Associative
    {
        typedef map<
            pair<int, char>
          , pair<double, std::string>
          , pair<abstract, int> >
        map_type;
        typedef identity_view<map_type> transformed_type;

        BOOST_MPL_ASSERT((traits::is_associative<transformed_type>));
        BOOST_MPL_ASSERT((traits::is_random_access<transformed_type>));

        map_type m(
            make_pair<int>('X')
          , make_pair<double>("Men")
          , make_pair<abstract>(2));
        transformed_type t(m);

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

        BOOST_MPL_ASSERT((boost::is_same<boost::fusion::result_of::value_of<first>::type, boost::fusion::pair<int, char> >));
        BOOST_MPL_ASSERT((boost::is_same<boost::fusion::result_of::value_of<second>::type, boost::fusion::pair<double, std::string> >));
        BOOST_MPL_ASSERT((boost::is_same<boost::fusion::result_of::value_of<third>::type, boost::fusion::pair<abstract, int> >));
    }

    {
        // compile test only
        // make sure result_of::deref_data returns a reference
        typedef map<pair<float, int> > map_type;
        typedef identity_view<map_type> transformed_type;
        typedef boost::fusion::result_of::begin<transformed_type>::type i_type;
        typedef boost::fusion::result_of::deref_data<i_type>::type r_type;
        BOOST_STATIC_ASSERT((boost::is_same<r_type, int&>::value));
    }

    {
        // compile test only
        // make sure result_of::deref_data is const correct
        typedef map<pair<float, int> > const map_type;
        typedef identity_view<map_type> transformed_type;
        typedef boost::fusion::result_of::begin<transformed_type>::type i_type;
        typedef boost::fusion::result_of::deref_data<i_type>::type r_type;
        BOOST_STATIC_ASSERT((boost::is_same<r_type, int const&>::value));
    }

    {
        // compile test only
        // make sure result_of::deref_data will not const for non-constant references
        typedef map<pair<float, int&> > const map_type;
        typedef identity_view<map_type> transformed_type;
        typedef boost::fusion::result_of::begin<transformed_type>::type i_type;
        typedef boost::fusion::result_of::deref_data<i_type>::type r_type;
        BOOST_STATIC_ASSERT((boost::is_same<r_type, int&>::value));
    }

    return boost::report_errors();
}

