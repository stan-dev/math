// Copyright (C) 2020 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
#include <boost/parser/parser.hpp>

#include <iostream>
#include <fstream>
#include <string>


namespace bp = boost::parser;

//[ logging_error_handler
struct logging_error_handler
{
    logging_error_handler() {}
    logging_error_handler(std::string_view filename) :
        filename_(filename), ofs_(filename_)
    {
        if (!ofs_)
            throw std::runtime_error("Could not open file.");
    }

    // This is the function called by Boost.Parser after a parser fails the
    // parse at an expectation point and throws a parse_error.  It is expected
    // to create a diagnostic message, and put it where it needs to go.  In
    // this case, we're writing it to a log file.  This function returns a
    // bp::error_handler_result, which is an enum with two enumerators -- fail
    // and rethrow.  Returning fail fails the top-level parse; returning
    // rethrow just re-throws the parse_error exception that got us here in
    // the first place.
    template<typename Iter, typename Sentinel>
    bp::error_handler_result
    operator()(Iter first, Sentinel last, bp::parse_error<Iter> const & e) const
    {
        bp::write_formatted_expectation_failure_error_message(
            ofs_, filename_, first, last, e);
        return bp::error_handler_result::fail;
    }

    // This function is for users to call within a semantic action to produce
    // a diagnostic.
    template<typename Context, typename Iter>
    void diagnose(
        bp::diagnostic_kind kind,
        std::string_view message,
        Context const & context,
        Iter it) const
    {
        bp::write_formatted_message(
            ofs_,
            filename_,
            bp::_begin(context),
            it,
            bp::_end(context),
            message);
    }

    // This is just like the other overload of diagnose(), except that it
    // determines the Iter parameter for the other overload by calling
    // _where(ctx).
    template<typename Context>
    void diagnose(
        bp::diagnostic_kind kind,
        std::string_view message,
        Context const & context) const
    {
        diagnose(kind, message, context, bp::_where(context).begin());
    }

    std::string filename_;
    mutable std::ofstream ofs_;
};
//]

//[ using_logging_error_handler
int main()
{
    std::cout << "Enter a list of integers, separated by commas. ";
    std::string input;
    std::getline(std::cin, input);

    constexpr auto parser = bp::int_ >> *(',' > bp::int_);
    logging_error_handler error_handler("parse.log");
    auto const result = bp::parse(input, bp::with_error_handler(parser, error_handler));

    if (result) {
        std::cout << "It looks like you entered:\n";
        for (int x : *result) {
            std::cout << x << "\n";
        }
    }
}
//]
