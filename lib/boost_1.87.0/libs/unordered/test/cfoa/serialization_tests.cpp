// Copyright (C) 2023-2024 Joaquin M Lopez Munoz
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "../objects/test.hpp"
#include "../helpers/random_values.hpp"

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/unordered/concurrent_flat_map.hpp>
#include <boost/unordered/concurrent_flat_set.hpp>
#include <boost/unordered/concurrent_node_map.hpp>
#include <boost/unordered/concurrent_node_set.hpp>

namespace {

  template <class Container, typename ArchivePair>
  void serialization_tests(
    Container*, ArchivePair*, test::random_generator generator)
  {
    using output_archive = typename ArchivePair::first_type ;
    using input_archive = typename ArchivePair::second_type;

    BOOST_LIGHTWEIGHT_TEST_OSTREAM << "serialization_tests1\n";
    {
      Container c;
  
      std::ostringstream oss;
      {
        output_archive oa(oss);
        oa << boost::serialization::make_nvp("container", c);
      }

      test::random_values<Container> values(100, generator);
      Container c2(values.begin(), values.end());
      std::istringstream iss(oss.str());
      input_archive ia(iss);
      ia >> boost::serialization::make_nvp("container", c2);
      BOOST_TEST(c2.empty());
    }

    BOOST_LIGHTWEIGHT_TEST_OSTREAM << "serialization_tests2\n";
    {
      test::random_values<Container> values(100, generator);
      Container c(values.begin(), values.end());

      std::ostringstream oss;
      {
        output_archive oa(oss);
        oa << boost::serialization::make_nvp("container", c);
      }

      Container c2;
      std::istringstream iss(oss.str());
      input_archive ia(iss);
      ia >> boost::serialization::make_nvp("container", c2);
      BOOST_TEST(c == c2);
    }
  }

  using test::default_generator;

  std::pair<
    boost::archive::text_oarchive, boost::archive::text_iarchive>*
    text_archive;
  std::pair<
    boost::archive::xml_oarchive, boost::archive::xml_iarchive>*
    xml_archive;

  boost::concurrent_flat_map<
    test::object, test::object, test::hash, test::equal_to>* test_flat_map;
  boost::concurrent_node_map<
    test::object, test::object, test::hash, test::equal_to>* test_node_map;
  boost::concurrent_flat_set<
    test::object, test::hash, test::equal_to>* test_flat_set;
  boost::concurrent_node_set<
    test::object, test::hash, test::equal_to>* test_node_set;

  UNORDERED_TEST(serialization_tests,
    ((test_flat_map)(test_node_map)(test_flat_set)(test_node_set))
    ((text_archive)(xml_archive))
    ((default_generator)))
}

RUN_TESTS()
