//
// Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/mysql/static_results.hpp>

#include "log_error.hpp"

#ifdef BOOST_MYSQL_CXX14

//[example_connection_pool_handle_request_cpp
//
// File: handle_request.cpp
//

#include <boost/mysql/error_code.hpp>
#include <boost/mysql/error_with_diagnostics.hpp>

#include <boost/asio/cancel_after.hpp>
#include <boost/asio/spawn.hpp>
#include <boost/json/parse.hpp>
#include <boost/json/serialize.hpp>
#include <boost/json/value_from.hpp>
#include <boost/json/value_to.hpp>
#include <boost/optional/optional.hpp>
#include <boost/url/parse.hpp>
#include <boost/variant2/variant.hpp>

#include <chrono>
#include <cstdint>
#include <exception>
#include <string>

#include "handle_request.hpp"
#include "repository.hpp"
#include "types.hpp"

// This file contains all the boilerplate code to dispatch HTTP
// requests to API endpoints. Functions here end up calling
// note_repository fuctions.

namespace asio = boost::asio;
namespace http = boost::beast::http;
using boost::mysql::error_code;
using boost::mysql::string_view;
using namespace notes;

namespace {

// Attempts to parse a numeric ID from a string.
// If you're using C++17, you can use std::from_chars, instead
static boost::optional<std::int64_t> parse_id(const std::string& from)
{
    try
    {
        std::size_t consumed = 0;
        int res = std::stoi(from, &consumed);
        if (consumed != from.size())
            return {};
        else if (res < 0)
            return {};
        return res;
    }
    catch (const std::exception&)
    {
        return {};
    }
}

// Encapsulates the logic required to match a HTTP request
// to an API endpoint, call the relevant note_repository function,
// and return an HTTP response.
class request_handler
{
    // The HTTP request we're handling. Requests are small in size,
    // so we use http::request<http::string_body>
    const http::request<http::string_body>& request_;

    // The repository to access MySQL
    note_repository repo_;

    // Creates an error response
    http::response<http::string_body> error_response(http::status code, string_view msg) const
    {
        http::response<http::string_body> res;

        // Set the status code
        res.result(code);

        // Set the keep alive option
        res.keep_alive(request_.keep_alive());

        // Set the body
        res.body() = msg;

        // Adjust the content-length field
        res.prepare_payload();

        // Done
        return res;
    }

    // Used when the request's Content-Type header doesn't match what we expect
    http::response<http::string_body> invalid_content_type() const
    {
        return error_response(http::status::bad_request, "Invalid content-type");
    }

    // Used when the request body didn't match the format we expect
    http::response<http::string_body> invalid_body() const
    {
        return error_response(http::status::bad_request, "Invalid body");
    }

    // Used when the request's method didn't match the ones allowed by the endpoint
    http::response<http::string_body> method_not_allowed() const
    {
        return error_response(http::status::method_not_allowed, "Method not allowed");
    }

    // Used when the request target couldn't be matched to any API endpoint
    http::response<http::string_body> endpoint_not_found() const
    {
        return error_response(http::status::not_found, "The requested resource was not found");
    }

    // Used when the user requested a note (e.g. using GET /note/<id> or PUT /note/<id>)
    // but the note doesn't exist
    http::response<http::string_body> note_not_found() const
    {
        return error_response(http::status::not_found, "The requested note was not found");
    }

    // Creates a response with a serialized JSON body.
    // T should be a type with Boost.Describe metadata containing the
    // body data to be serialized
    template <class T>
    http::response<http::string_body> json_response(const T& body) const
    {
        http::response<http::string_body> res;

        // A JSON response is always a 200
        res.result(http::status::ok);

        // Set the content-type header
        res.set("Content-Type", "application/json");

        // Set the keep-alive option
        res.keep_alive(request_.keep_alive());

        // Serialize the body data into a string and use it as the response body.
        // We use Boost.JSON's automatic serialization feature, which uses Boost.Describe
        // reflection data to generate a serialization function for us.
        res.body() = boost::json::serialize(boost::json::value_from(body));

        // Adjust the content-length header
        res.prepare_payload();

        // Done
        return res;
    }

    // Returns true if the request's Content-Type is set to JSON
    bool has_json_content_type() const
    {
        auto it = request_.find("Content-Type");
        return it != request_.end() && it->value() == "application/json";
    }

    // Attempts to parse the request body as a JSON into an object of type T.
    // T should be a type with Boost.Describe metadata.
    // We use boost::system::result, which may contain a result or an error.
    template <class T>
    boost::system::result<T> parse_json_request() const
    {
        error_code ec;

        // Attempt to parse the request into a json::value.
        // This will fail if the provided body isn't valid JSON.
        auto val = boost::json::parse(request_.body(), ec);
        if (ec)
            return ec;

        // Attempt to parse the json::value into a T. This will
        // fail if the provided JSON doesn't match T's shape.
        return boost::json::try_value_to<T>(val);
    }

    http::response<http::string_body> handle_request_impl(boost::asio::yield_context yield)
    {
        // Parse the request target. We use Boost.Url to do this.
        auto url = boost::urls::parse_origin_form(request_.target());
        if (url.has_error())
            return error_response(http::status::bad_request, "Invalid request target");

        // We will be iterating over the target's segments to determine
        // which endpoint we are being requested
        auto segs = url->segments();
        auto segit = segs.begin();
        auto seg = *segit++;

        // All endpoints start with /notes
        if (seg != "notes")
            return endpoint_not_found();

        if (segit == segs.end())
        {
            if (request_.method() == http::verb::get)
            {
                // GET /notes: retrieves all the notes.
                // The request doesn't have a body.
                // The response has a JSON body with multi_notes_response format
                auto res = repo_.get_notes(yield);
                return json_response(multi_notes_response{std::move(res)});
            }
            else if (request_.method() == http::verb::post)
            {
                // POST /notes: creates a note.
                // The request has a JSON body with note_request_body format.
                // The response has a JSON body with single_note_response format.

                // Parse the request body
                if (!has_json_content_type())
                    return invalid_content_type();
                auto args = parse_json_request<note_request_body>();
                if (args.has_error())
                    return invalid_body();

                // Actually create the note
                auto res = repo_.create_note(args->title, args->content, yield);

                // Return the newly crated note as response
                return json_response(single_note_response{std::move(res)});
            }
            else
            {
                return method_not_allowed();
            }
        }
        else
        {
            // The URL has the form /notes/<note-id>. Parse the note ID.
            auto note_id = parse_id(*segit++);
            if (!note_id.has_value())
            {
                return error_response(
                    http::status::bad_request,
                    "Invalid note_id specified in request target"
                );
            }

            // /notes/<note-id>/<something-else> is not a valid endpoint
            if (segit != segs.end())
                return endpoint_not_found();

            if (request_.method() == http::verb::get)
            {
                // GET /notes/<note-id>: retrieves a single note.
                // The request doesn't have a body.
                // The response has a JSON body with single_note_response format

                // Get the note
                auto res = repo_.get_note(*note_id, yield);

                // If we didn't find it, return a 404 error
                if (!res.has_value())
                    return note_not_found();

                // Return it as response
                return json_response(single_note_response{std::move(*res)});
            }
            else if (request_.method() == http::verb::put)
            {
                // PUT /notes/<note-id>: replaces a note.
                // The request has a JSON body with note_request_body format.
                // The response has a JSON body with single_note_response format.

                // Parse the JSON body
                if (!has_json_content_type())
                    return invalid_content_type();
                auto args = parse_json_request<note_request_body>();
                if (args.has_error())
                    return invalid_body();

                // Perform the update
                auto res = repo_.replace_note(*note_id, args->title, args->content, yield);

                // Check that it took effect. Otherwise, it's because the note wasn't there
                if (!res.has_value())
                    return note_not_found();

                // Return the updated note as response
                return json_response(single_note_response{std::move(*res)});
            }
            else if (request_.method() == http::verb::delete_)
            {
                // DELETE /notes/<note-id>: deletes a note.
                // The request doesn't have a body.
                // The response has a JSON body with delete_note_response format.

                // Attempt to delete the note
                bool deleted = repo_.delete_note(*note_id, yield);

                // Return whether the delete was successful in the response.
                // We don't fail DELETEs for notes that don't exist.
                return json_response(delete_note_response{deleted});
            }
            else
            {
                return method_not_allowed();
            }
        }
    }

public:
    // Constructor
    request_handler(const http::request<http::string_body>& req, note_repository repo)
        : request_(req), repo_(repo)
    {
    }

    // Generates a response for the request passed to the constructor
    http::response<http::string_body> handle_request(boost::asio::yield_context yield)
    {
        try
        {
            // Attempt to handle the request. We use cancel_after to set
            // a timeout to the overall operation
            return asio::spawn(
                yield.get_executor(),
                [this](asio::yield_context yield2) { return handle_request_impl(yield2); },
                asio::cancel_after(std::chrono::seconds(30), yield)
            );
        }
        catch (const boost::mysql::error_with_diagnostics& err)
        {
            // A Boost.MySQL error. This will happen if you don't have connectivity
            // to your database, your schema is incorrect or your credentials are invalid.
            // Log the error, including diagnostics, and return a generic 500
            log_error(
                "Uncaught exception: ",
                err.what(),
                "\nServer diagnostics: ",
                err.get_diagnostics().server_message()
            );
            return error_response(http::status::internal_server_error, "Internal error");
        }
        catch (const std::exception& err)
        {
            // Another kind of error. This indicates a programming error or a severe
            // server condition (e.g. out of memory). Same procedure as above.
            log_error("Uncaught exception: ", err.what());
            return error_response(http::status::internal_server_error, "Internal error");
        }
    }
};

}  // namespace

// External interface
boost::beast::http::response<boost::beast::http::string_body> notes::handle_request(
    const boost::beast::http::request<boost::beast::http::string_body>& request,
    note_repository repo,
    boost::asio::yield_context yield
)
{
    return request_handler(request, repo).handle_request(yield);
}

//]

#endif
