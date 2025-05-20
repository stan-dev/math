// Copyright (C) 2017 Andrzej Krzemienski.
//
// Use, modification, and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
// See http://www.boost.org/lib/optional for documentation.
//
// You are welcome to contact the author at:
//  akrzemi1@gmail.com

#include "boost/optional/optional.hpp"

#ifdef BOOST_BORLANDC
#pragma hdrstop
#endif

#include "boost/core/lightweight_test.hpp"
#include "boost/core/lightweight_test_trait.hpp"
#include "boost/type_traits/is_base_of.hpp"
#include "boost/optional/detail/experimental_traits.hpp"

#ifndef BOOST_OPTIONAL_DETAIL_NO_DEFAULTED_MOVE_FUNCTIONS

struct PrivDefault
{
  private: PrivDefault() {}
};

struct CustDefault
{
  CustDefault() {}
};

struct CustomizedTrivial
{
  CustomizedTrivial() {}
};

struct DeletedDefault
{
  BOOST_DELETED_FUNCTION(DeletedDefault())
};

namespace boost { namespace optional_config {

template <> struct optional_uses_direct_storage_for<CustomizedTrivial> : boost::true_type {};

}}

struct CustDtor
{
  ~CustDtor() {}
};

struct NoDefault
{
  explicit NoDefault(int) {}
};

struct Empty {};

template <typename T, typename U>
struct Aggregate { T t; U u; };

struct CustAssign
{
  CustAssign& operator=(CustAssign const&) { return *this; }
};

struct CustMove
{
  CustMove(CustMove &&) {}
};

void test_type_traits()
{
  // this only tests if type traits are implemented correctly
  BOOST_TEST_TRAIT_TRUE(( boost::optional_config::optional_uses_direct_storage_for<int> ));
  BOOST_TEST_TRAIT_TRUE(( boost::optional_config::optional_uses_direct_storage_for<double> ));

  BOOST_TEST_TRAIT_TRUE(( boost::optional_config::optional_uses_direct_storage_for<CustomizedTrivial> ));

  BOOST_TEST_TRAIT_FALSE(( boost::optional_config::optional_uses_direct_storage_for<PrivDefault> ));
  BOOST_TEST_TRAIT_FALSE(( boost::optional_config::optional_uses_direct_storage_for<NoDefault> ));
  BOOST_TEST_TRAIT_FALSE(( boost::optional_config::optional_uses_direct_storage_for<CustDefault> ));
  BOOST_TEST_TRAIT_FALSE(( boost::optional_config::optional_uses_direct_storage_for<Aggregate<int, CustDefault> > ));

  BOOST_TEST_TRAIT_FALSE(( boost::optional_config::optional_uses_direct_storage_for<CustDtor> ));
  BOOST_TEST_TRAIT_FALSE(( boost::optional_config::optional_uses_direct_storage_for<CustAssign> ));
  BOOST_TEST_TRAIT_FALSE(( boost::optional_config::optional_uses_direct_storage_for<CustMove> ));
  BOOST_TEST_TRAIT_FALSE(( boost::optional_config::optional_uses_direct_storage_for<Aggregate<int, CustMove> > ));

  BOOST_TEST_TRAIT_TRUE(( boost::optional_detail::is_trivially_semiregular<int> ));
  BOOST_TEST_TRAIT_TRUE(( boost::optional_detail::is_trivially_semiregular<double> ));


#ifndef BOOST_OPTIONAL_DETAIL_NO_SPEC_FOR_TRIVIAL_TYPES
  BOOST_TEST_TRAIT_TRUE(( boost::optional_config::optional_uses_direct_storage_for<Empty> ));
  BOOST_TEST_TRAIT_TRUE(( boost::optional_config::optional_uses_direct_storage_for<Aggregate<int, double> > ));
  BOOST_TEST_TRAIT_TRUE(( boost::optional_config::optional_uses_direct_storage_for<Aggregate<Aggregate<Empty, int>, double> > ));

  BOOST_TEST_TRAIT_TRUE(( boost::optional_detail::is_trivially_semiregular<Empty> ));
  BOOST_TEST_TRAIT_TRUE(( boost::optional_detail::is_trivially_semiregular<Aggregate<int, double> > ));
  BOOST_TEST_TRAIT_TRUE(( boost::optional_detail::is_trivially_semiregular<Aggregate<Aggregate<Empty, int>, double> > ));

  BOOST_TEST_TRAIT_TRUE(( boost::optional_detail::is_trivially_semiregular<PrivDefault> ));
  BOOST_TEST_TRAIT_TRUE(( boost::optional_detail::is_trivially_semiregular<NoDefault> ));
  BOOST_TEST_TRAIT_TRUE(( boost::optional_detail::is_trivially_semiregular<CustDefault> ));
  BOOST_TEST_TRAIT_TRUE(( boost::optional_detail::is_trivially_semiregular<Aggregate<int, CustDefault> > ));
  BOOST_TEST_TRAIT_TRUE(( boost::optional_detail::is_trivially_semiregular<DeletedDefault> ));
  BOOST_TEST_TRAIT_TRUE(( boost::optional_detail::is_trivially_semiregular<Aggregate<int, DeletedDefault> > ));
#endif

  BOOST_TEST_TRAIT_FALSE(( boost::optional_config::optional_uses_direct_storage_for<DeletedDefault> ));
  BOOST_TEST_TRAIT_FALSE(( boost::optional_config::optional_uses_direct_storage_for<Aggregate<int, DeletedDefault> > ));

  BOOST_TEST_TRAIT_FALSE(( boost::optional_detail::is_trivially_semiregular<CustDtor> ));
  BOOST_TEST_TRAIT_FALSE(( boost::optional_detail::is_trivially_semiregular<CustAssign> ));
  BOOST_TEST_TRAIT_FALSE(( boost::optional_detail::is_trivially_semiregular<CustMove> ));
  BOOST_TEST_TRAIT_FALSE(( boost::optional_detail::is_trivially_semiregular<Aggregate<int, CustMove> > ));
}

void test_trivial_copyability()
{
  BOOST_TEST_TRAIT_TRUE((boost::is_base_of<boost::optional_detail::tc_optional_base<int>, boost::optional<int> > ));
  BOOST_TEST_TRAIT_TRUE((boost::is_base_of<boost::optional_detail::tc_optional_base<double>, boost::optional<double> > ));
  BOOST_TEST_TRAIT_TRUE((boost::is_base_of<boost::optional_detail::tc_optional_base<CustomizedTrivial>, boost::optional<CustomizedTrivial> > ));
  BOOST_TEST_TRAIT_FALSE((boost::is_base_of<boost::optional_detail::tc_optional_base<DeletedDefault>, boost::optional<DeletedDefault> > ));

#ifndef BOOST_OPTIONAL_DETAIL_NO_SPEC_FOR_TRIVIAL_TYPES
  BOOST_TEST_TRAIT_TRUE(( boost::optional_detail::is_trivially_semiregular<boost::optional<int> > ));
  BOOST_TEST_TRAIT_TRUE(( boost::optional_detail::is_trivially_semiregular<boost::optional<double> > ));
  BOOST_TEST_TRAIT_TRUE(( boost::optional_detail::is_trivially_semiregular<boost::optional<CustomizedTrivial> > ));
  BOOST_TEST_TRAIT_TRUE(( boost::optional_detail::is_trivially_semiregular<boost::optional<Empty> > ));
  BOOST_TEST_TRAIT_TRUE(( boost::optional_detail::is_trivially_semiregular<boost::optional<Aggregate<int, double> > > ));
  BOOST_TEST_TRAIT_TRUE(( boost::optional_detail::is_trivially_semiregular<boost::optional<Aggregate<Aggregate<Empty, int>, double> > > ));

  BOOST_TEST_TRAIT_FALSE(( boost::optional_detail::is_trivially_semiregular<boost::optional<DeletedDefault> > ));
#endif
}

#endif

int main()
{
#ifndef BOOST_OPTIONAL_DETAIL_NO_DEFAULTED_MOVE_FUNCTIONS
  test_type_traits();
  test_trivial_copyability();
#endif
  return boost::report_errors();
}
