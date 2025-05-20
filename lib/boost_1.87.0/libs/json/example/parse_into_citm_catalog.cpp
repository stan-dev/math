//
// Copyright (c) 2021 Peter Dimov
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/json
//

//
// An example that compares the performance of json::parse and
// json::parse_into on citm_catalog.json
//

#include <boost/json.hpp>
#if !defined(BOOST_DESCRIBE_CXX14)

#include <boost/config/pragma_message.hpp>

BOOST_PRAGMA_MESSAGE( "This example requires C++14" )

int main() {}

#else

#include <boost/describe.hpp>

#include <chrono>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <utility>
#include <vector>

struct event
{
    std::nullptr_t description;
    unsigned long long id;
    boost::variant2::variant<std::nullptr_t, std::string> logo;
    std::string name;
    std::vector<unsigned long long> subTopicIds;
    std::nullptr_t subjectCode;
    std::nullptr_t subtitle;
    std::vector<unsigned long long> topicIds;
};

BOOST_DESCRIBE_STRUCT(event, (), (description, id, logo, name, subTopicIds, subjectCode, subtitle, topicIds))

struct price
{
    unsigned amount;
    unsigned long long audienceSubCategoryId;
    unsigned long long seatCategoryId;
};

BOOST_DESCRIBE_STRUCT(price, (), (amount, audienceSubCategoryId, seatCategoryId))

struct area
{
    unsigned long long areaId;
    std::vector<unsigned long long> blockIds;
};

BOOST_DESCRIBE_STRUCT(area, (), (areaId, blockIds))

struct seatCategory
{
    std::vector<area> areas;
    unsigned long long seatCategoryId;
};

BOOST_DESCRIBE_STRUCT(seatCategory, (), (areas, seatCategoryId))

struct performance
{
    unsigned long long eventId;
    unsigned long long id;
    boost::variant2::variant<std::nullptr_t, std::string> logo;
    std::nullptr_t name;
    std::vector<price> prices;
    std::vector<seatCategory> seatCategories;
    std::nullptr_t seatMapImage;
    unsigned long long start;
    std::string venueCode;
};

BOOST_DESCRIBE_STRUCT(performance, (), (eventId, id, logo, name, prices, seatCategories, seatMapImage, start, venueCode))

struct citm_catalog
{
    std::map<std::string, std::string> areaNames;
    std::map<std::string, std::string> audienceSubCategoryNames;
    std::map<std::string, std::string> blockNames;
    std::map<std::string, event> events;
    std::vector<performance> performances;
    std::map<std::string, std::string> seatCategoryNames;
    std::map<std::string, std::string> subTopicNames;
    std::map<std::string, std::string> subjectNames;
    std::map<std::string, std::string> topicNames;
    std::map<std::string, std::vector<unsigned long long>> topicSubTopics;
    std::map<std::string, std::string> venueNames;
};

BOOST_DESCRIBE_STRUCT(citm_catalog, (),
    (areaNames, audienceSubCategoryNames, blockNames, events, performances,
    seatCategoryNames, subTopicNames, subjectNames, topicNames, topicSubTopics, venueNames))

using namespace std::chrono_literals;

int main()
{
    std::ifstream is( "citm_catalog.json" );
    std::string json( std::istreambuf_iterator<char>( is ), std::istreambuf_iterator<char>{} );

    std::cout << "citm_catalog.json: " << json.size() << " bytes\n";

    {
        auto tp1 = std::chrono::steady_clock::now();

        boost::json::value jv = boost::json::parse( json );

        auto tp2 = std::chrono::steady_clock::now();
        std::cout << "boost::json::parse: " << (tp2 - tp1) / 1ms << " ms\n";
    }

    {
        auto tp1 = std::chrono::steady_clock::now();

        citm_catalog w;

        boost::system::error_code ec;
        boost::json::parse_into( w, json, ec );

        if( ec.failed() )
        {
            std::cout << "Error: " << ec.what() << std::endl;
        }

        auto tp2 = std::chrono::steady_clock::now();
        std::cout << "parse_into<citm_catalog>: " << (tp2 - tp1) / 1ms << " ms\n";
    }
}

#endif
