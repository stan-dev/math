// Copyright Antony Polukhin, 2016-2024.
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include <boost/array.hpp>
BOOST_NOINLINE void foo(int i);
BOOST_NOINLINE void bar(int i);
 
BOOST_NOINLINE void bar(int i) {
    boost::array<int, 5> a = {{-1, -231, -123, -23, -32}};
    if (i >= 0) {
        foo(a[i]);
    } else {
        std::terminate();
    }
}

BOOST_NOINLINE void foo(int i) {
    bar(--i);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//[getting_started_signal_handlers

#include <signal.h>     // ::signal, ::raise
#include <boost/stacktrace.hpp>

void my_signal_handler(int signum) {
    ::signal(signum, SIG_DFL);

    // Outputs nothing or trash on majority of platforms
    boost::stacktrace::safe_dump_to("./backtrace.dump");

    ::raise(SIGABRT);
}
//]

void setup_handlers() {
//[getting_started_setup_signel_handlers
    ::signal(SIGSEGV, &my_signal_handler);
    ::signal(SIGABRT, &my_signal_handler);
//]
}


//[getting_started_terminate_handlers
#include <cstdlib>       // std::abort
#include <exception>     // std::set_terminate
#include <iostream>      // std::cerr

#include <boost/stacktrace.hpp>

void my_terminate_handler() {
    try {
        std::cerr << boost::stacktrace::stacktrace();
    } catch (...) {}
    std::abort();
}
//]

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>     // std::cerr
#include <fstream>     // std::ifstream
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>

#ifndef BOOST_WINDOWS
inline void copy_and_run(const char* exec_name, char param, bool not_null) {
    std::cout << "Running with param " << param << std::endl;
    boost::filesystem::path command = exec_name;
    command = command.parent_path() / (command.stem().string() + param + command.extension().string());
    boost::filesystem::copy_file(exec_name, command, boost::filesystem::copy_options::overwrite_existing);

    boost::filesystem::path command_args = command;
    command_args += ' ';
    command_args += param;
    const int ret = std::system(command_args.string().c_str());

    std::cout << "End Running with param " << param << "; ret code is " << ret << std::endl;
    boost::system::error_code ignore;
    boost::filesystem::remove(command, ignore);
    if (not_null && !ret) {
        std::exit(97);
    } else if (!not_null && ret) {
        std::exit(ret);
    }
}
#endif

int run_0(const char* /*argv*/[]) {
//[getting_started_setup_terminate_handlers
    std::set_terminate(&my_terminate_handler);
//]
    foo(5);
    return 1;
}


int run_1(const char* /*argv*/[]) {
    setup_handlers();
    foo(5);
    return 11;
}

int run_2(const char* argv[]) {
    if (!boost::filesystem::exists("./backtrace.dump")) {
        if (std::string(argv[0]).find("noop") == std::string::npos) {
            return 21;
        }

        boost::stacktrace::stacktrace st = boost::stacktrace::stacktrace::from_dump(std::cin);
        if (st) {
            return 22;
        }
        return 0;
    }

//[getting_started_on_program_restart
    if (boost::filesystem::exists("./backtrace.dump")) {
        // there is a backtrace
        std::ifstream ifs("./backtrace.dump");

        boost::stacktrace::stacktrace st = boost::stacktrace::stacktrace::from_dump(ifs);
        std::cout << "Previous run crashed:\n" << st << std::endl; /*<-*/

        if (!st) {
            return 23;
        } /*->*/

        // cleaning up
        ifs.close();
        boost::filesystem::remove("./backtrace.dump");
    }
//]

    return 0;
}

#include <sstream>

int test_inplace() {
    const bool is_noop = !boost::stacktrace::stacktrace();

    {
        // This is very dependent on compiler and link flags. No sane way to make it work, because:
        // * BOOST_NOINLINE could be ignored by MSVC compiler if link-time optimization is enabled.
        // * BOOST_FORCEINLINE could be ignored by GCC depending on the std::vector default constructor length.
        const std::size_t frames_ss1 = boost::stacktrace::safe_dump_to("./backtrace2.dump");
        boost::stacktrace::stacktrace ss2;
        std::ifstream ifs("./backtrace2.dump");
        boost::stacktrace::stacktrace ss1 = boost::stacktrace::stacktrace::from_dump(ifs);
        ifs.close();
        boost::filesystem::remove("./backtrace2.dump");

        if (ss1.size() + 1 != frames_ss1 || ss2.size() != ss1.size()) {
            std::cerr << "51: Stacktraces differ. Dumped size == " << frames_ss1 << ".\n" << ss1 << "\n vs \n" << ss2 << '\n';
        } else if (ss1.size() > 1 && ss1[1].name() != ss2[1].name()) {
            std::cerr << "52: Stacktraces differ:\n" << ss1 << "\n vs \n" << ss2 << '\n';
        }
    }

    {
        // This is very dependent on compiler and link flags. No sane way to make it work, because:
        // * BOOST_NOINLINE could be ignored by MSVC compiler if link-time optimization is enabled.
        // * BOOST_FORCEINLINE could be ignored by GCC depending on the std::vector default constructor length.
        void* data[1024];
        const std::size_t frames_ss1 = boost::stacktrace::safe_dump_to(data, sizeof(data));
        boost::stacktrace::stacktrace ss2;
        boost::stacktrace::stacktrace ss1 = boost::stacktrace::stacktrace::from_dump(data, sizeof(data));

        if (ss1.size() + 1 != frames_ss1 || ss1.size() != ss2.size()) {
            std::cerr << "53: Stacktraces differ. Dumped size == " << frames_ss1 << ".\n" << ss1 << "\n vs \n" << ss2 << '\n';
        } else if (ss1.size() > 1 && ss1[1].name() != ss2[1].name()) {
            std::cerr << "54: Stacktraces differ:\n" << ss1 << "\n vs \n" << ss2 << '\n';
        }
    }

    {
        void* data[1024];
        boost::stacktrace::safe_dump_to(1024, data, sizeof(data));
        if (boost::stacktrace::stacktrace::from_dump(data, sizeof(data))) {
            std::cerr << "Stacktrace not empty!\n";
            return 55;
        }
    }

    {
        void* data[1024];
        boost::stacktrace::safe_dump_to(1, data, sizeof(data));
        if (!is_noop && !boost::stacktrace::stacktrace::from_dump(data, sizeof(data))) {
            std::cerr << "Stacktrace empty!\n";
            return 56;
        }
        const std::size_t size_1_skipped = boost::stacktrace::stacktrace::from_dump(data, sizeof(data)).size();
        boost::stacktrace::safe_dump_to(0, data, sizeof(data));
        const std::size_t size_0_skipped = boost::stacktrace::stacktrace::from_dump(data, sizeof(data)).size();

        if (!is_noop && (size_1_skipped + 1 != size_0_skipped)) {
            std::cerr << "failed to skip 1 frame!\n";
            return 57;
        }
    }

    {
        boost::stacktrace::safe_dump_to(0, 1, "./backtrace3.dump");
        std::ifstream ifs("./backtrace3.dump");
        boost::stacktrace::stacktrace ss1 = boost::stacktrace::stacktrace::from_dump(ifs);
        ifs.close();

        boost::stacktrace::safe_dump_to(1, 1, "./backtrace3.dump");
        ifs.open("./backtrace3.dump");
        boost::stacktrace::stacktrace ss2 = boost::stacktrace::stacktrace::from_dump(ifs);
        ifs.close();

        boost::filesystem::remove("./backtrace3.dump");

#ifdef BOOST_WINDOWS
        // `ss2` could be empty on some combinations of Windows+MSVC.
        if (!ss2) {
            return 0;
        }
#endif

        if (ss1.size() != ss2.size()) {
            std::cerr << "Stacktraces differ:\n" << ss1 << "\n vs \n" << ss2 << '\n';
            return 58;
        }

        if (!is_noop && ss1.size() != 1) {
            std::cerr << "Stacktraces does not have size 1:\n" << ss1 << '\n';
            return 59;
        }

        if (ss1 && ss1[0].address() == ss2[0].address()) {
            std::cerr << "Stacktraces must differ:\n" << ss1 << "\n vs \n" << ss2 << '\n';
            return 60;
        }
    }

    return 0;
}


int main(int argc, const char* argv[]) {
    if (argc < 2) {
        // On Windows the debugger could be active. In that case tests hang and the CI run fails.
#ifndef BOOST_WINDOWS
        copy_and_run(argv[0], '0', true);

        // We are copying files to make sure that stacktrace printing works independently from executable name
        copy_and_run(argv[0], '1', true);
        copy_and_run(argv[0], '2', false);
#endif

        return test_inplace();
    }

    switch (argv[1][0]) {
    case '0': return run_0(argv);
    case '1': return run_1(argv);
    }

    return 404;
}


