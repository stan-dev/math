//  Copyright (c) 2023 Andrey Semashev
//
//  Distributed under the Boost Software License, Version 1.0.
//  See accompanying file LICENSE_1_0.txt or copy at
//  https://www.boost.org/LICENSE_1_0.txt)

#include <boost/generator_iterator.hpp>
#include <boost/indirect_reference.hpp>
#include <boost/next_prior.hpp>
#include <boost/pointee.hpp>
#include <boost/shared_container_iterator.hpp>
#include <boost/iterator/advance.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/iterator/distance.hpp>
#include <boost/iterator/filter_iterator.hpp>
#include <boost/iterator/function_input_iterator.hpp>
#include <boost/iterator/function_output_iterator.hpp>
#include <boost/iterator/indirect_iterator.hpp>
#include <boost/iterator/interoperable.hpp>
#include <boost/iterator/is_iterator.hpp>
#include <boost/iterator/is_lvalue_iterator.hpp>
#include <boost/iterator/is_readable_iterator.hpp>
#include <boost/iterator/iterator_adaptor.hpp>
#include <boost/iterator/iterator_archetypes.hpp>
#include <boost/iterator/iterator_categories.hpp>
#include <boost/iterator/iterator_concepts.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/iterator/iterator_traits.hpp>
#include <boost/iterator/minimum_category.hpp>
#include <boost/iterator/new_iterator_tests.hpp>
#include <boost/iterator/permutation_iterator.hpp>
#include <boost/iterator/reverse_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/iterator/zip_iterator.hpp>

template< typename Iterator >
class adapted_iterator :
    public boost::iterators::iterator_adaptor< adapted_iterator< Iterator >, Iterator >
{
    friend class iterator_core_access;

private:
    typedef boost::iterators::iterator_adaptor< adapted_iterator< Iterator >, Iterator > base_type;

public:
    explicit adapted_iterator(Iterator it) : base_type(it) {}
};

int main()
{
    unsigned char buf[8];
    adapted_iterator< unsigned char* > b(buf), e(buf + sizeof(buf));
    return boost::iterators::distance(b, e) == static_cast< adapted_iterator< unsigned char* >::difference_type >(sizeof(buf));
}
