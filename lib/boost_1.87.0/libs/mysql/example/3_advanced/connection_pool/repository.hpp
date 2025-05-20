//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_MYSQL_EXAMPLE_3_ADVANCED_CONNECTION_POOL_REPOSITORY_HPP
#define BOOST_MYSQL_EXAMPLE_3_ADVANCED_CONNECTION_POOL_REPOSITORY_HPP

//[example_connection_pool_repository_hpp
//
// File: repository.hpp
//

#include <boost/mysql/connection_pool.hpp>
#include <boost/mysql/string_view.hpp>

#include <boost/asio/error.hpp>
#include <boost/asio/spawn.hpp>
#include <boost/optional/optional.hpp>

#include <cstdint>

#include "types.hpp"

namespace notes {

using boost::optional;
using boost::mysql::string_view;

// A lightweight wrapper around a connection_pool that allows
// creating, updating, retrieving and deleting notes in MySQL.
// This class encapsulates the database logic.
// All operations are async, and use stackful coroutines (boost::asio::yield_context).
// If the database can't be contacted, or unexpected database errors are found,
// an exception of type boost::mysql::error_with_diagnostics is thrown.
class note_repository
{
    boost::mysql::connection_pool& pool_;

public:
    // Constructor (this is a cheap-to-construct object)
    note_repository(boost::mysql::connection_pool& pool) : pool_(pool) {}

    // Retrieves all notes present in the database
    std::vector<note_t> get_notes(boost::asio::yield_context yield);

    // Retrieves a single note by ID. Returns an empty optional
    // if no note with the given ID is present in the database.
    optional<note_t> get_note(std::int64_t note_id, boost::asio::yield_context yield);

    // Creates a new note in the database with the given components.
    // Returns the newly created note, including the newly allocated ID.
    note_t create_note(string_view title, string_view content, boost::asio::yield_context yield);

    // Replaces the note identified by note_id, setting its components to the
    // ones passed. Returns the updated note. If no note with ID matching
    // note_id can be found, an empty optional is returned.
    optional<note_t> replace_note(
        std::int64_t note_id,
        string_view title,
        string_view content,
        boost::asio::yield_context yield
    );

    // Deletes the note identified by note_id. Returns true if
    // a matching note was deleted, false otherwise.
    bool delete_note(std::int64_t note_id, boost::asio::yield_context yield);
};

}  // namespace notes

//]

#endif
