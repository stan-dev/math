/*=============================================================================
    Copyright (c) 2022 Denis Mikhailov

    Distributed under the Boost Software License, Version 1.0. (See accompanying 
    file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#include <boost/fusion/support/tag_of.hpp>
#include <boost/static_assert.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/core/enable_if.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/if.hpp>

namespace mpl = boost::mpl;

template<typename T>
struct my_is_implicitly_reflectable;

namespace boost { namespace fusion
{
    struct boost_my_reflectable_tag;
    namespace detail
    {
        template<typename Sequence>
        struct tag_of_fallback<Sequence, typename boost::enable_if<boost::is_same<Sequence, Sequence> >::type>
        {
            typedef typename mpl::if_<my_is_implicitly_reflectable<Sequence>
                         , boost_my_reflectable_tag
                         , non_fusion_tag
            >::type type;
        };
    }
}}

struct reflectable {};
struct non_reflectable {};

template<typename T>
struct my_is_implicitly_reflectable : mpl::false_ {};
template<>
struct my_is_implicitly_reflectable<reflectable> : mpl::true_ {};

typedef boost::fusion::traits::tag_of<reflectable>::type ReflectableTag;
typedef boost::fusion::traits::tag_of<non_reflectable>::type NonReflectableTag;
BOOST_STATIC_ASSERT((boost::is_same<ReflectableTag, boost::fusion::boost_my_reflectable_tag>::value));
BOOST_STATIC_ASSERT((boost::is_same<NonReflectableTag, boost::fusion::non_fusion_tag>::value));

int
main() { }

