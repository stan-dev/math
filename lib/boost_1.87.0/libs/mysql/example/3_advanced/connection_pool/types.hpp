//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_EXAMPLE_3_ADVANCED_CONNECTION_POOL_TYPES_HPP
#define BOOST_MYSQL_EXAMPLE_3_ADVANCED_CONNECTION_POOL_TYPES_HPP

//[example_connection_pool_types_hpp
//
// File: types.hpp
//

#include <boost/core/span.hpp>
#include <boost/describe/class.hpp>
#include <boost/optional/optional.hpp>

#include <string>
#include <vector>

// Contains type definitions used in the REST API and database code.
// We use Boost.Describe (BOOST_DESCRIBE_STRUCT) to add reflection
// capabilities to our types. This allows using Boost.MySQL
// static interface (i.e. static_results<T>) to parse query results,
// and Boost.JSON automatic serialization/deserialization.

namespace notes {

struct note_t
{
    // The unique database ID of the object.
    std::int64_t id;

    // The note's title.
    std::string title;

    // The note's actual content.
    std::string content;
};
BOOST_DESCRIBE_STRUCT(note_t, (), (id, title, content))

//
// REST API requests.
//

// Used for creating and replacing notes
struct note_request_body
{
    // The title that the new note should have.
    std::string title;

    // The content that the new note should have.
    std::string content;
};
BOOST_DESCRIBE_STRUCT(note_request_body, (), (title, content))

//
// REST API responses.
//

// Used by endpoints returning several notes (like GET /notes).
struct multi_notes_response
{
    // The retrieved notes.
    std::vector<note_t> notes;
};
BOOST_DESCRIBE_STRUCT(multi_notes_response, (), (notes))

// Used by endpoints returning a single note (like GET /notes/<id>)
struct single_note_response
{
    // The retrieved note.
    note_t note;
};
BOOST_DESCRIBE_STRUCT(single_note_response, (), (note))

// Used by DELETE /notes/<id>
struct delete_note_response
{
    // true if the note was found and deleted, false if the note didn't exist.
    bool deleted;
};
BOOST_DESCRIBE_STRUCT(delete_note_response, (), (deleted))

}  // namespace notes

//]

#endif
