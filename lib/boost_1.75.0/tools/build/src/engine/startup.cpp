/*
Copyright 2020 René Ferdinand Rivera Morell
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE_1_0.txt or http://www.boost.org/LICENSE_1_0.txt)
*/

#include "startup.h"
#include "rules.h"
#include "frames.h"
#include "object.h"
#include "pathsys.h"
#include "cwd.h"
#include "filesys.h"
#include "output.h"
#include "variable.h"

#include <string>
#include <algorithm>

namespace
{
    void bind_builtin(
        char const *name_, LIST *(*f)(FRAME *, int flags),
        int flags, char const **args)
    {
        FUNCTION *func;
        OBJECT *name = object_new(name_);
        func = function_builtin(f, flags, args);
        new_rule_body(root_module(), name, func, 1);
        function_free(func);
        object_free(name);
    }
} // namespace

void b2::startup::load_builtins()
{
    {
        char const *args[] = {"dir", "?", 0};
        bind_builtin("boost-build", builtin_boost_build, 0, args);
    }
}

LIST *b2::startup::builtin_boost_build(FRAME *frame, int flags)
{
    b2::jam::list dir_arg{lol_get(frame->args, 0)};
    std::string dir;
    if (!dir_arg.empty()) dir = b2::jam::object(*dir_arg.begin());

    b2::jam::variable dot_bootstrap_file{".bootstrap-file"};
    if (dot_bootstrap_file)
    {
        err_printf(
            "Error: Illegal attempt to re-bootstrap the build system by invoking\n"
            "\n"
            "   'boost-build '%s' ;\n"
            "\n"
            "Please consult the documentation at "
            "'https://boostorg.github.io/build/'.\n\n",
           dir.c_str());
        return L0;
    }

    // # Add the given directory to the path so we can find the build system. If
    // # dir is empty, has no effect.
    b2::jam::variable dot_boost_build_file{".boost-build-file"};
    b2::jam::list dot_boost_build_file_val{static_cast<b2::jam::list>(dot_boost_build_file)};
    std::string boost_build_jam = b2::jam::object{*dot_boost_build_file_val.begin()};
    std::string boost_build_dir;
    if (b2::paths::is_rooted(dir))
        boost_build_dir = dir;
    else
        boost_build_dir = b2::paths::normalize(
            std::string{boost_build_jam}+"/../"+dir);
    b2::jam::list search_path{b2::jam::object{boost_build_dir}};
    b2::jam::variable BOOST_BUILD_PATH{"BOOST_BUILD_PATH"};
    search_path.append(BOOST_BUILD_PATH);

    // We set the global, and env, BOOST_BUILD_PATH so that the loading of the
    // build system finds the initial set of modules needed for starting it up.
    BOOST_BUILD_PATH = search_path;

    // The code that loads the rest of B2, in particular the site-config.jam
    // and user-config.jam configuration files uses os.environ, so we need to
    // update the value there.
    b2::jam::variable dot_ENVIRON__BOOST_BUILD_PATH{".ENVIRON", "BOOST_BUILD_PATH"};
    dot_ENVIRON__BOOST_BUILD_PATH = search_path;

    // # Try to find the build system bootstrap file 'bootstrap.jam'.
    std::string bootstrap_file;
    for (auto path: search_path)
    {
        std::string file = b2::jam::object{path};
        file = b2::paths::normalize(file+"/bootstrap.jam");
        if (b2::filesys::is_file(file))
        {
            bootstrap_file = file;
            break;
        }
    }

    // # There is no bootstrap.jam we can find, exit with an error.
    if (bootstrap_file.empty())
    {
        err_printf(
            "Unable to load B2: could not find build system.\n"
            "-----------------------------------------------\n"
            "%s attempted to load the build system by invoking\n"
            "\n"
            "   'boost-build %s ;'\n"
            "\n"
            "but we were unable to find 'bootstrap.jam' in the specified directory "
            "or in BOOST_BUILD_PATH:\n",
            boost_build_jam.c_str(), dir.c_str());
        for (auto path: search_path)
        {
            std::string file = b2::jam::object{path};
            err_printf("    %s\n", file.c_str());
        }
        err_puts(
            "Please consult the documentation at "
            "'https://boostorg.github.io/build/'.\n\n");
        return L0;
    }

    // Set the bootstrap=file var as it's used by the build system to refer to
    // the rest of the build system files.
    dot_bootstrap_file = b2::jam::list{b2::jam::object{bootstrap_file}};

    // Show where we found it, if asked.
    b2::jam::variable dot_OPTION__debug_configuration{".OPTION", "debug-configration"};
    if (dot_OPTION__debug_configuration)
    {
        out_printf("notice: loading B2 from %s\n", bootstrap_file.c_str());
    }

    // # Load the build system, now that we know where to start from.
    parse_file(b2::jam::object{bootstrap_file}, frame);

    return L0;
}

extern char const *saved_argv0;

bool b2::startup::bootstrap(FRAME *frame)
{
    b2::jam::list ARGV = b2::jam::variable{"ARGV"};
    b2::jam::object opt_debug_configuration{"--debug-configuration"};
    b2::jam::variable dot_OPTION__debug_configuration{".OPTION", "debug-configration"};
    for (auto arg: ARGV)
    {
        if (opt_debug_configuration == arg)
        {
            dot_OPTION__debug_configuration = b2::jam::list{b2::jam::object{"true"}};
            break;
        }
    }

    const std::string b2_exe_path{executable_path(saved_argv0)};
    const std::string boost_build_jam{"boost-build.jam"};
    std::string b2_file_path;

    // Attempt to find the `boost-build.jam` boot file in work directory tree.
    if (b2_file_path.empty())
    {
        std::string work_dir{b2::paths::normalize(b2::cwd_str()) + "/"};
        while (b2_file_path.empty() && !work_dir.empty())
        {
            if (b2::filesys::is_file(work_dir + boost_build_jam))
                b2_file_path = work_dir + boost_build_jam;
            else if (work_dir.length() == 1 && work_dir[0] == '/')
                work_dir.clear();
            else
            {
                auto parent_pos = work_dir.rfind('/', work_dir.length() - 2);
                if (parent_pos != std::string::npos)
                    work_dir.erase(parent_pos + 1);
                else
                    work_dir.clear();
            }
        }
    }

    // Check relative to the executable for portable install location.
    if (b2_file_path.empty())
    {
        const std::string path{
            b2::paths::normalize(
                b2_exe_path + "/../.b2/kernel/" + boost_build_jam)};
        if (b2::filesys::is_file(path))
            b2_file_path = path;
    }

    // Check relative to the executable for portable install location.
    if (b2_file_path.empty())
    {
        const std::string path{
            b2::paths::normalize(
                b2_exe_path + "/../../share/boost-build/" + boost_build_jam)};
        if (b2::filesys::is_file(path))
            b2_file_path = path;
    }

    // Check the BOOST_BUILD_PATH paths.
    if (b2_file_path.empty())
    {
        b2::jam::list BOOST_BUILD_PATH = b2::jam::variable{"BOOST_BUILD_PATH"};
        for (auto search_path: BOOST_BUILD_PATH)
        {
            std::string path = b2::jam::object{search_path};
            path = b2::paths::normalize(path+"/"+boost_build_jam);
            if (b2::filesys::is_file(path))
            {
                b2_file_path = path;
                break;
            }
        }
    }

    // Indicate a load failure when we can't find the build file.
    if (b2_file_path.empty())
    {
        const char * not_found_error =
            "Unable to load B2: could not find 'boost-build.jam'\n"
            "---------------------------------------------------\n"
            "Attempted search from '%s' up to the root "
            "at '%s'\n"
            "Please consult the documentation at "
            "'https://boostorg.github.io/build/'.\n\n";
        err_printf(not_found_error, b2::cwd_str().c_str(), b2_exe_path.c_str());
        return false;
    }

    // Show where we found it, if asked.
    if (dot_OPTION__debug_configuration)
    {
        out_printf("notice: found boost-build.jam at %s\n", b2_file_path.c_str());
    }

    // Load the build system bootstrap file we found. But check we did that.
    b2::jam::variable dot_boost_build_file{".boost-build-file"};
    dot_boost_build_file = b2_file_path;
    b2::jam::object b2_file_path_sym{b2_file_path};
    parse_file(b2_file_path_sym, frame);
    b2::jam::list dot_dot_bootstrap_file_val = b2::jam::variable{".bootstrap-file"};
    if (dot_dot_bootstrap_file_val.empty())
    {
        err_printf(
            "Unable to load B2\n"
            "-----------------\n"
            "'%s' was found by searching from %s up to the root.\n"
            "\n"
            "However, it failed to call the 'boost-build' rule to indicate "
            "the location of the build system.\n"
            "\n"
            "Please consult the documentation at "
            "'https://boostorg.github.io/build/'.\n\n",
            b2_file_path.c_str(), b2::cwd_str().c_str());
        return false;
    }

    return true;
}
