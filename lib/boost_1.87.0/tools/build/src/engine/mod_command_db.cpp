/*
Copyright 2024 René Ferdinand Rivera Morell
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE.txt or https://www.bfgroup.xyz/b2/LICENSE.txt)
*/

#include "mod_command_db.h"

#include "command.h"
#include "cwd.h"
#include "events.h"
#include "lists.h"
#include "mod_db.h"
#include "output.h"
#include "pathsys.h"
#include "regexp.h"

#include "ext_bfgroup_lyra.h"

#include <unordered_map>

namespace b2 { namespace command_db {

struct database
{
	bool output_flag = false;
	std::string output_format = "json";
	std::string output_filename = "compile_commands.json";
	std::string output_directory;
	std::string db_directory;
	std::unordered_map<std::string, std::unique_ptr<regex::program>>
		regex_cache;
	std::unique_ptr<property_db> prop_db;
	uint64_t command_index = 0;

	database()
	{
		// B2 doesn't change directories. And runs everything relative to CWD.
		// So we can cache the value to fill into the database.
		db_directory = b2::cwd_str();
		output_directory = db_directory;
	}

	static database & get()
	{
		static database db;
		return db;
	}

	void set_output_format(const std::string & f)
	{
		// Set up for the format and create the output database.
		output_flag = true;
		output_format = f;
		prop_db.reset(new property_db);

		// Default to all targets and no regular action output.
		++globs.noexec;
		for (int i = 0; i < DEBUG_MAX; ++i) globs.debug[i] = 0;
		++anyhow;

		// Events to track the commands and exit to generate the output.
		add_event_callback(event_tag::pre_exec_cmd,
			std::function<void(TARGET *)>(
				[](TARGET * t) { database::get().pre_exec_cmd(t); }));
		add_event_callback(event_tag::exit_main,
			std::function<void(int)>(
				[](int s) { database::get().exit_main(s); }));
	}

	void set_output_filename(const std::string & f) { output_filename = f; }

	void pre_exec_cmd(TARGET * t)
	{
		CMD * cmd = (CMD *)t->cmds;
		auto args_out = list_cref(lol_get((LOL *)&cmd->args, 0));
		auto args_in = list_cref(lol_get((LOL *)&cmd->args, 1));
		auto settings = t->settings;
		for (; settings; settings = settings->next)
		{
			value_ref symbol(settings->symbol);
			list_cref value(settings->value);
			if (symbol == "COMMAND_DATABASE" && !value.empty())
			{
				auto db_file = args_in[0]->str();
				auto db_output = args_out[0]->str();
				auto db_command
					= extract_command(value[0]->str(), cmd->buf->value);
				list_ref db_node;
				db_node.push_back("[]").push_back(double(command_index++));
				prop_db->emplace(
					list_ref(db_node).push_back("output").cref(), db_output);
				prop_db->emplace(
					list_ref(db_node).push_back("file").cref(), db_file);
				prop_db->emplace(
					list_ref(db_node).push_back("directory").cref(),
					db_directory);
				prop_db->emplace(
					list_ref(db_node).push_back("command").cref(), db_command);
			}
		}
	}

	std::string extract_command(const char * cmd_regex, const char * cmd_text)
	{
		regex::program * regex_prog = nullptr;
		auto regex_prog_i = regex_cache.find(cmd_regex);
		if (regex_prog_i == regex_cache.end())
		{
			auto regex_prog_e = std::unique_ptr<regex::program>(
				new regex::program(cmd_regex));
			regex_prog = regex_prog_e.get();
			regex_cache.emplace(cmd_regex, std::move(regex_prog_e));
		}
		else
		{
			regex_prog = regex_prog_i->second.get();
		}
		auto regex_sub = regex_prog->search(cmd_text)[0];
		return std::string(
			regex_sub.cbegin(), regex_sub.cend() - regex_sub.begin());
	}

	void exit_main(int status)
	{
		if (status == EXIT_FAIL) return;
		std::string filename = output_filename;
		if (!b2::paths::is_rooted(output_filename))
		{
			if (!b2::paths::is_rooted(output_directory))
				filename = b2::cwd_str() + "/";
			filename += output_directory + "/" + output_filename;
			filename = b2::paths::normalize(filename);
		}
		if (prop_db->write_file(filename, output_format))
			out_printf("...wrote command database '%s'...\n", filename.c_str());
		else
			err_printf("...writing to command database '%s' FAILED...\n",
				filename.c_str());
	}
};

void declare_args(lyra::cli & cli)
{
	cli |= lyra::opt(
		[](const std::string & f) { database::get().set_output_format(f); },
		"format")
			   .name("--command-database")
			   .help("Output a compile commands database as format.");
	cli |= lyra::opt(
		[](const std::string & f) { database::get().set_output_filename(f); },
		"filename")
			   .name("--command-database-out")
			   .help(
				   "Filename to output the command database to. "
				   "A relative path for the filename is rooted to the project "
				   "build-dir.");
}

void set_output_dir(value_ref dirname)
{
	database::get().output_directory = dirname;
}

}} // namespace b2::command_db
