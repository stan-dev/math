/*
Copyright 2020-2022 René Ferdinand Rivera Morell
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE.txt or https://www.bfgroup.xyz/b2/LICENSE.txt)
*/

#include "startup.h"
#include "bindjam.h"
#include "compile.h"
#include "cwd.h"
#include "filesys.h"
#include "frames.h"
#include "mod_jam_modules.h"
#include "modules.h"
#include "object.h"
#include "output.h"
#include "pathsys.h"
#include "rules.h"
#include "value.h"
#include "variable.h"

#include <algorithm>
#include <cstdlib>
#include <string>
#include <vector>

namespace {
void bind_builtin(char const * name_,
	LIST * (*f)(FRAME *, int flags),
	int flags,
	char const ** args)
{
	FUNCTION * func;
	OBJECT * name = object_new(name_);
	func = function_builtin(f, flags, args);
	new_rule_body(root_module(), name, func, 1);
	function_free(func);
	object_free(name);
}
} // namespace

void b2::startup::load_builtins()
{
	{
		char const * args[] = { "dir", "?", 0 };
		bind_builtin("boost-build", builtin_boost_build, 0, args);
	}
}

LIST * b2::startup::builtin_boost_build(FRAME * frame, int flags)
{
	// Do nothing, but keep the rule, for backwards compatability.
	// But do record the path passed in as a fallback to the loading.
	b2::jam::variable(".boost-build-dir")
		= b2::list_ref(lol_get(frame->args, 0));
	return L0;
}

extern char const * saved_argv0;

void bootstrap_dirscan(void * dirs,
	OBJECT * path,
	int found,
	timestamp const * const)
{
	if (file_is_file(path) == 1) return;
	_pathname p(object_str(path));
	std::string basename = p.base() + p.suffix();
	if (basename == "." || basename == "..") return;
	static_cast<std::vector<std::string> *>(dirs)->push_back(object_str(path));
}

bool b2::startup::bootstrap(FRAME * frame)
{
	b2::value_ref opt_debug_configuration { "--debug-configuration" };
	b2::jam::variable dot_OPTION__debug_configuration { ".OPTION",
		"debug-configration" };
	for (auto arg : *b2::jam::variable("ARGV"))
	{
		if (opt_debug_configuration == arg)
		{
			dot_OPTION__debug_configuration = "true";
			break;
		}
	}

	// We use the executable path as a root for searches.
	char * b2_exe_path_pchar = executable_path(saved_argv0);
	const std::string b2_exe_path { b2_exe_path_pchar };
	if (b2_exe_path_pchar)
	{
		std::free(b2_exe_path_pchar);
	}

	/**
	 * Search for a boost-build.jam file to load in various locations. This is
	 * solely for backwards compatiblity. The boost-build.jam file found is
	 * loaded, but the `boost-build` invocation in it is ignored.
	 */

	const std::string boost_build_jam { "boost-build.jam" };
	std::string b2_file_path;

	// Attempt to find the `boost-build.jam` boot file in work directory tree.
	// I.e. in current directory and ancestor directories.
	if (b2_file_path.empty())
	{
		std::string work_dir { b2::paths::normalize(b2::cwd_str()) + "/" };
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
	// ~b2(.exe)/.b2/kernel/boost-build.jam
	if (b2_file_path.empty())
	{
		const std::string path { b2::paths::normalize(
			b2_exe_path + "/../.b2/kernel/" + boost_build_jam) };
		if (b2::filesys::is_file(path)) b2_file_path = path;
	}

	// Check relative to the executable for portable install location.
	// /usr/share/b2/src/kernel/boost-build.jam
	if (b2_file_path.empty())
	{
		const std::string path { b2::paths::normalize(
			b2_exe_path + "/../../share/b2/src/kernel/" + boost_build_jam) };
		if (b2::filesys::is_file(path)) b2_file_path = path;
	}

	// Check the BOOST_BUILD_PATH (and BOOST_ROOT) paths.
	// $(BOOST_BUILD_PATH)/boost-build.jam
	// $(BOOST_ROOT)/boost-build.jam
	if (b2_file_path.empty())
	{
		b2::list_ref BOOST_BUILD_PATH = *b2::jam::variable("BOOST_BUILD_PATH");
		// For back-compat with Boost we also search in the BOOST_ROOT location.
		BOOST_BUILD_PATH.push_back("BOOST_ROOT");
		for (auto search_path : BOOST_BUILD_PATH)
		{
			std::string path = b2::value_ref { search_path };
			path = b2::paths::normalize(path + "/" + boost_build_jam);
			if (b2::filesys::is_file(path))
			{
				b2_file_path = path;
				break;
			}
		}
	}

	// Show where we found it, if asked.
	if (!b2_file_path.empty() && dot_OPTION__debug_configuration)
	{
		out_printf(
			"notice: found boost-build.jam at %s\n", b2_file_path.c_str());
	}

	// Load the boost-build file if we find it for backwards compatability. We
	// ignore the `boost-build ..` invocation it does. Preferring to find our
	// own bootstrap engine file.
	if (!b2_file_path.empty())
	{
		b2::jam::variable dot_boost_build_file { ".boost-build-file" };
		dot_boost_build_file = b2_file_path;
		b2::value_ref b2_file_path_sym { b2_file_path };
		parse_file(b2_file_path_sym, frame);
	}

	/**
	 * Search for a build-system.jam file to load in various locations. The
	 * `build-system.jam` is the starting point of loading the build system.
	 */

	const std::string buildsystem_jam { "build-system.jam" };
	std::string buildsystem_file;
	std::string buildsystem_files_searched;

	// Check various locations relative to executable.
	if (buildsystem_file.empty())
	{
		const char * dirs[] = {
			// Check relative to the executable for portable install location.
			".b2/",
			// Check relative to the exec for system install location.
			"../share/b2/",
			// Check relative to the exec for legacy install location.
			"../share/b2/src/",
			// Check development location relative to executable in src/engine.
			"../",
			// Check development location relative to executable at root.
			"src/",
			// Check for special Boost location.
			"tools/build/src/"
		};
		for (auto dir : dirs)
		{
			const std::string path { b2::paths::normalize(
				b2_exe_path + "/../" + dir + buildsystem_jam) };
			if (b2::filesys::is_file(path))
			{
				buildsystem_file = path;
				break;
			}
			buildsystem_files_searched += "  " + path + "\n";
		}
	}

	// Check the development tree for the file to support not-installed
	// b2 executable locations. I.e. when building b2 with b2 for engine
	// development.
	if (buildsystem_file.empty())
	{
		std::string work_dir(b2::paths::normalize(b2_exe_path + "/..") + "/");
		while (buildsystem_file.empty() && !work_dir.empty())
		{
			buildsystem_files_searched
				+= "  " + work_dir + "src/" + buildsystem_jam + "\n";
			if (b2::filesys::is_file(work_dir + "src/" + buildsystem_jam))
				buildsystem_file = work_dir + "src/" + buildsystem_jam;
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

	// Last resort, search in the directory referenced by the boost-build rule.
	if (buildsystem_file.empty())
	{
		b2::list_cref dot_boost_build_dir
			= b2::jam::variable(".boost-build-dir");
		if (!dot_boost_build_dir.empty())
		{
			std::string dir = b2::value_ref(*dot_boost_build_dir.begin());
			if (!b2_file_path.empty() && !b2::paths::is_rooted(dir))
				dir = b2_file_path + "/../" + dir;
			const std::string path
				= b2::paths::normalize(dir + "/" + buildsystem_jam);
			buildsystem_files_searched += "  " + path + "\n";
			if (b2::filesys::is_file(path)) buildsystem_file = path;
		}
	}

	// Failed to find the build files to load.
	if (buildsystem_file.empty())
	{
		err_printf(
			"Unable to load B2\n"
			"-----------------\n"
			"No '%s' was found by searching for:\n"
			"%s\n"
			"Please consult the documentation at "
			"'https://www.bfgroup.xyz/b2/'.\n\n",
			buildsystem_file.c_str(), buildsystem_files_searched.c_str());
		return false;
	}

	// Show where we found the bootstrap, if asked.
	if (!buildsystem_file.empty() && dot_OPTION__debug_configuration)
	{
		out_printf("notice: loading B2 from %s\n", buildsystem_file.c_str());
	}

	// // Load the build system bootstrap file we found.
	// parse_file(b2::value_ref { buildsystem_file }, frame);

	// Add any subdirs to the build-system.jam to search for the jam files.
	std::string buildsystem_dir
		= b2::paths::normalize(buildsystem_file + "/..");
	std::vector<std::string> buildsystem_subdirs;
	file_dirscan(value_ref(buildsystem_dir), &bootstrap_dirscan,
		&buildsystem_subdirs);
	b2::jam::variable boost_build_path_v("BOOST_BUILD_PATH");
	for (auto subdir : buildsystem_subdirs) boost_build_path_v += subdir;
	boost_build_path_v += buildsystem_dir;
	b2::jam::variable(".ENVIRON", "BOOST_BUILD_PATH") = boost_build_path_v;

	jam::jam_context context(frame);

	// The modules module is required globally for everyone for the `ìmport`
	// rule.
	list_ref modules_imports("modules");
	jam::modules::import(
		list_cref(*modules_imports), list_cref(), list_cref(), context.ref());

	// Process option plugins first to allow them to prevent loading the rest of
	// the build system.
	list_ref option_imports("option");
	jam::modules::import(
		list_cref(*option_imports), list_cref(), list_cref(), context.ref());
	list_ref dont_build(jam::run_rule(context.frame, "option.process"));

	// Should we skip building, i.e. loading the build system, according to the
	// options processed?
	if (!dont_build.empty()) return true;

	// Check `--build-system` option.
	list_ref build_system_imports;
	for (auto arg : *jam::variable("ARGV"))
	{
		std::string opt = arg->str();
		if (opt.substr(0, 15) == "--build-system=")
			build_system_imports.push_back(opt.substr(15));
	}

	// Import the build-system module to start off the regular build process.
	if (build_system_imports.empty())
		build_system_imports.push_back("build-system");
	jam::modules::import(list_cref(*build_system_imports), list_cref(),
		list_cref(), context.ref());

	return true;
}
