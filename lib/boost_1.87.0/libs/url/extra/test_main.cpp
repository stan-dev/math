//
// Copyright (c) 2019 Vinnie Falco (vinnie.falco@gmail.com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/url
//

#include "test_suite.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <vector>

#ifdef _MSC_VER
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <debugapi.h>
#include <crtdbg.h>
#include <sstream>
#endif

namespace test_suite {

//------------------------------------------------
//
// debug_stream
//
//------------------------------------------------

#ifdef _MSC_VER

namespace detail {

class debug_streambuf
    : public std::stringbuf
{
    std::ostream& os_;
    bool dbg_;

    void
    write(char const* s)
    {
        if(dbg_)
            ::OutputDebugStringA(s);
        os_ << s;
    }

public:
    explicit
    debug_streambuf(
        std::ostream& os)
        : os_(os)
        , dbg_(::IsDebuggerPresent() != 0)
    {
    }

    ~debug_streambuf()
    {
        sync();
    }

    int sync() override
    {
        write(this->str().c_str());
        this->str("");
        return 0;
    }
};

} // detail

//------------------------------------------------

/** std::ostream with Visual Studio IDE redirection.

    Instances of this stream wrap a specified
    `ostream` (such as `cout` or `cerr`). If a
    debugger is attached when the stream is
    created, output will additionally copied
    to the debugger's output window.
*/
class debug_stream : public std::ostream
{
    detail::debug_streambuf buf_;

public:
    /** Construct a stream.

        @param os The output stream to wrap.
    */
    explicit
    debug_stream(
        std::ostream& os)
        : std::ostream(&buf_)
        , buf_(os)
    {
        if(os.flags() & std::ios::unitbuf)
            std::unitbuf(*this);
    }

};

#else

using debug_stream  = std::ostream&;

#endif

//------------------------------------------------
//
// suite
//
//------------------------------------------------

any_suite::
~any_suite() = default;

//------------------------------------------------

suites&
suites::
instance() noexcept
{
    class suites_impl : public suites
    {
        std::vector<any_suite const*> v_;

    public:
        virtual ~suites_impl() = default;

        void
        insert(any_suite const& t) override
        {
            v_.push_back(&t);
        }

        iterator
        begin() const noexcept override
        {
            if(! v_.empty())
                return &v_[0];
            return nullptr;
        }

        iterator
        end() const noexcept override
        {
            if(! v_.empty())
                return &v_[0] + v_.size();
            return nullptr;
        }

        void
        sort() override
        {
            std::sort(
                v_.begin(),
                v_.end(),
                []( any_suite const* p0,
                    any_suite const* p1)
                {
                    auto s0 = p0->name();
                    auto n0 = std::strlen(s0);
                    auto s1 = p1->name();
                    auto n1 = std::strlen(s1);
                    return std::lexicographical_compare(
                        s0, s0 + n0, s1, s1 + n1);
                });
        }
    };
    static suites_impl impl;
    return impl;
}

//------------------------------------------------
//
// runner
//
//------------------------------------------------

any_runner*&
any_runner::
instance_impl() noexcept
{
    static any_runner* p = nullptr;
    return p;
}

any_runner&
any_runner::
instance() noexcept
{
    return *instance_impl();
}

any_runner::
any_runner() noexcept
    : prev_(instance_impl())
{
    instance_impl() = this;
}

any_runner::
~any_runner()
{
    instance_impl() = prev_;
}

//------------------------------------------------
//
// implementation
//
//------------------------------------------------

namespace detail {

bool
test_impl(
    bool cond,
    char const* expr,
    char const* func,
    char const* file,
    int line)
{
    return any_runner::instance().test(
        cond, expr, func, file, line);
}

void
throw_failed_impl(
    const char* expr,
    char const* excep,
    char const* func,
    char const* file,
    int line)
{
    std::stringstream ss;
    ss <<
        "expression '" << expr <<
        "' didn't throw '" << excep <<
        "' in function '" << func <<
        "'";
   any_runner::instance().test(false, ss.str().c_str(),
       func, file, line);
}

void
no_throw_failed_impl(
    const char* expr,
    char const* excep,
    char const* func,
    char const* file,
    int line)
{
    std::stringstream ss;
    ss <<
        "expression '" << expr <<
        "' threw '" << excep <<
        "' in function '" << func <<
        "'";
   any_runner::instance().test(false, ss.str().c_str(),
       func, file, line);
}

void
no_throw_failed_impl(
    const char* expr,
    char const* func,
    char const* file,
    int line)
{
    std::stringstream ss;
    ss <<
        "expression '" << expr <<
        "' threw in function '" << func << "'";
   any_runner::instance().test(false, ss.str().c_str(),
       func, file, line);
}

//------------------------------------------------
//
// simple_runner
//
//------------------------------------------------

using clock_type =
    std::chrono::steady_clock;

struct elapsed
{
    clock_type::duration d;
};

std::ostream&
operator<<(
    std::ostream& os,
    elapsed const& e)
{
    using namespace std::chrono;
    auto const ms = duration_cast<
        milliseconds>(e.d);
    if(ms < seconds{1})
    {
        os << ms.count() << "ms";
    }
    else
    {
        std::streamsize width{
            os.width()};
        std::streamsize precision{
            os.precision()};
        os << std::fixed <<
            std::setprecision(1) <<
            (ms.count() / 1000.0) << "s";
        os.precision(precision);
        os.width(width);
    }
    return os;
}

} // detail
} // test_suite

//------------------------------------------------

namespace test_suite {
namespace detail {

class simple_runner : public any_runner
{
    struct summary
    {
        char const* name;
        clock_type::time_point start;
        clock_type::duration elapsed;
        std::atomic<std::size_t> failed;
        std::atomic<std::size_t> total;

        summary(summary const& other) noexcept
            : name(other.name)
            , start(other.start)
            , elapsed(other.elapsed)
            , failed(other.failed.load())
            , total(other.total.load())
        {
        }

        explicit
        summary(char const* name_) noexcept
            : name(name_)
            , start(clock_type::now())
            , failed(0)
            , total(0)
        {
        }
    };

    summary all_;
    std::ostream& log_;
    std::vector<summary> v_;

public:
    explicit
    simple_runner(
        std::ostream& os)
        : all_("(all)")
        , log_(os)
    {
        std::unitbuf(log_);
        v_.reserve(256);
    }

    virtual ~simple_runner()
    {
        log_ <<
            elapsed{clock_type::now() -
                all_.start } << ", " <<
            v_.size() << " suites, " <<
            all_.failed.load() << " failures, " <<
            all_.total.load() << " total." <<
                std::endl;
    }

    // true if no failures
    bool
    success() const noexcept
    {
        return all_.failed.load() == 0;
    }

    void
    run(any_suite const& test) override
    {
        v_.emplace_back(test.name());
        log_ << test.name() << "\n";
        test.run();
        v_.back().elapsed =
            clock_type::now() -
            v_.back().start;
    }

    void
    note(char const* msg) override
    {
        log_ << msg << "\n";
    }

    std::ostream&
    log() noexcept override
    {
        return log_;
    }

    char const*
    filename(
        char const* file)
    {
        auto const p0 = file;
        auto p = p0 + std::strlen(file);
        while(p-- != p0)
        {
        #ifdef _MSC_VER
            if(*p == '\\')
        #else
            if(*p == '/')
        #endif
                break;
        }
        return p + 1;
    }

    bool test(bool, char const*, char const*, char const*, int) override;
};

//------------------------------------------------

bool
simple_runner::
test(
    bool cond,
    char const* expr,
    char const* func,
    char const* file,
    int line)
{
    auto const id = ++all_.total;
    ++v_.back().total;
    if(cond)
        return true;
    ++all_.failed;
    ++v_.back().failed;
    (void)func;
    log_ <<
        "#" << id << " " <<
        filename(file) << "(" << line << ") "
        "failed: " << expr <<
        //" in " << func <<
        "\n";
    log_.flush();
    return false;
}

//------------------------------------------------

int
run(std::ostream& out,
    int argc, char const* const* argv)
{
    if(argc == 2)
    {
        std::string arg(argv[1]);
        if(arg == "-h" || arg == "--help")
        {
            log <<
                "Usage:\n"
                "  " << argv[0] << ": { <suite-name>... }" <<
                std::endl;
            return EXIT_SUCCESS;
        }
    }

    simple_runner any_runner(out);
    suites::instance().sort();
    if(argc == 1)
    {
        for(any_suite const* sp :
                suites::instance())
            any_runner.run(*sp);
    }
    else
    {
        std::vector<std::string> args;
        args.reserve(argc - 1);
        for(int i = 0; i < argc - 1; ++i)
            args.push_back(argv[i + 1]);
        for(auto const& e : suites::instance())
        {
            std::string s(e->name());
            if(std::find_if(
                args.begin(), args.end(),
                [&](std::string const& arg)
                {
                    if(arg.size() > s.size())
                        return false;
                    return s.compare(
                        s.size() - arg.size(),
                        arg.size(),
                        arg.data(),
                        arg.size()) == 0;
                }) != args.end())
            {
                any_runner.run(*e);
            }
        }
    }
    return any_runner.success() ?
        EXIT_SUCCESS : EXIT_FAILURE;
}

} // detail
} // test_suite

//------------------------------------------------

// Simple main used to produce stand
// alone executables that run unit tests.
int main(int argc, char const* const* argv)
{
#if defined(_MSC_VER) && !defined(__clang__)
    int flags = _CrtSetDbgFlag(_CRTDBG_REPORT_FLAG);
    flags |= _CRTDBG_LEAK_CHECK_DF;
    _CrtSetDbgFlag(flags);
#endif

    ::test_suite::debug_stream log(std::cerr);
    return ::test_suite::detail::run(log, argc, argv);
}
