/*
Copyright 2019-2023 René Ferdinand Rivera Morell
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE.txt or https://www.bfgroup.xyz/b2/LICENSE.txt)
*/

#include "mod_regex.h"

#include "filesys.h"
#include "mem.h"
#include "regexp.h"
#include "tasks.h"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <numeric>
#include <string>
#include <vector>

using namespace b2;
using namespace b2::regex;

list_ref b2::regex_split(
	const std::tuple<value_ref, value_ref> & string_separator)
{
	list_ref result;
	string_t string = std::get<0>(string_separator);
	auto re = program(std::get<1>(string_separator)->str());
	auto pos = string.c_str();
	auto prev = pos;
	for (auto re_i = re.search(pos); re_i; ++re_i)
	{
		result.push_back(prev, re_i->begin());
		prev = re_i->end();
		/* Handle empty matches */
		if (*pos == '\0')
			break;
		else if (pos == re_i->end())
			pos++;
		else
			pos = re_i->end();
	}
	result.push_back(pos);
	return result;
}

list_ref b2::regex_split_each(list_cref to_split, value_ref separator)
{
	list_ref result;
	auto re = program(separator->str());
	for (auto string : to_split)
	{
		auto pos = string->str();
		auto prev = pos;
		for (auto re_i = re.search(pos); re_i; ++re_i)
		{
			result.push_back(prev, re_i->begin());
			prev = re_i->end();
			/* Handle empty matches */
			if (*pos == '\0')
				break;
			else if (pos == re_i->end())
				pos++;
			else
				pos = re_i->end();
		}
		result.push_back(pos);
	}
	return result;
}

list_ref b2::regex_match(
	value_ref pattern, value_ref string, const std::vector<int_t> & indices)
{
	list_ref result;
	auto re = program(pattern->str());
	if (auto re_i = re.search(string->str()))
	{
		if (!indices.empty())
		{
			for (int_t index : indices)
			{
				if (index < NSUBEXP)
					result.push_back(re_i[index].begin(), re_i[index].end());
			}
		}
		else
		{
			for (int_t index = 1;
				 index < NSUBEXP && re_i[index].end() != nullptr; ++index)
				result.push_back(re_i[index].begin(), re_i[index].end());
		}
	}
	return result;
}

list_ref b2::regex_transform(
	list_cref list, value_ref pattern, const std::vector<int_t> & indices)
{
	list_ref result;
	auto re = program(pattern->str());
	if (!indices.empty())
	{
		for (auto string : list)
		{
			if (auto re_i = re.search(string->str()))
				for (int_t index : indices)
				{
					if (index < NSUBEXP && re_i[index].end() != nullptr)
						result.push_back(
							re_i[index].begin(), re_i[index].end());
				}
		}
	}
	else
	{
		for (auto string : list)
		{
			if (auto re_i = re.search(string->str()))
				if (re_i[1].end() != nullptr)
					result.push_back(re_i[1].begin(), re_i[1].end());
		}
	}
	return result;
}

value_ref b2::regex_escape(
	value_ref string, value_ref symbols, value_ref escape_symbol)
{
	string_t result = string;
	string_t::size_type i = 0;
	for (i = result.find_first_of(symbols, i); i != string_t::npos;
		 i = result.find_first_of(symbols, i))
	{
		result.insert(i, escape_symbol);
		i += escape_symbol->as_string().size + 1;
	}
	return value_ref(result);
}

namespace b2 { namespace {
string_t regex_replace(
	const char * string, program & re, const char * replacement)
{
	std::string result;
	auto pos = string;
	for (auto re_i = re.search(pos); re_i; ++re_i)
	{
		result.append(pos, re_i->begin());
		result.append(replacement);
		/* Handle empty matches */
		if (*pos == '\0')
			break;
		else if (pos == re_i->end())
			pos++;
		else
			pos = re_i->end();
	}
	result.append(pos);
	return result;
}
}} // namespace b2

value_ref b2::regex_replace(const std::tuple<value_ref, value_ref, value_ref> &
		string_match_replacement)
{
	value_ref string = std::get<0>(string_match_replacement);
	auto re = program(std::get<1>(string_match_replacement)->str());
	value_ref replacement = std::get<2>(string_match_replacement);
	return value_ref(regex_replace(string->str(), re, replacement->str()));
}

list_ref b2::regex_replace_each(
	list_cref list, value_ref match, value_ref replacement)
{
	list_ref result;
	auto re = program(match->str());
	for (auto string : list)
	{
		result.push_back(regex_replace(string->str(), re, replacement->str()));
	}
	return result;
}

int glob(char const * s, char const * c);

namespace b2 { namespace {

template <typename Int>
struct flags
{
	using value_type = Int;
	static const int size = sizeof(value_type) * 8;

	value_type value = 0;
	inline bool get(int i) const { return (value & (1 << i)) != 0; };
	inline void set(int i, bool v = true)
	{
		value = v ? (value | (1 << i)) : (value & ~(1 << 1));
	}
};

struct regex_grep_task
{
	list_cref file_glob_patterns;
	std::vector<b2::regex::program> text_grep_prog;
	std::shared_ptr<b2::task::group> grep_tasks;
	flags<std::uint16_t> text_grep_result_expressions;
	std::vector<std::unique_ptr<std::vector<std::string>>> intermediate;
	mutex_t mx;
	bool recursive_glob = false;

	regex_grep_task(list_cref files,
		list_cref patterns,
		flags<std::uint16_t> expressions,
		list_cref options)
		: file_glob_patterns(files)
		, text_grep_result_expressions(expressions)
	{
		text_grep_prog.reserve(patterns.length());
		for (auto p : patterns)
		{
			text_grep_prog.emplace_back(p->str());
		}
		for (auto option : options)
		{
			std::string opt = option->str();
			if (opt == "recursive") recursive_glob = true;
		}
		grep_tasks = b2::task::executor::get().make();
	}

	void dir_scan(value_ref dir)
	{
		grep_tasks->queue([this, dir] {
			file_dirscan(dir,
				reinterpret_cast<scanback>(&regex_grep_task::dirscan_callback),
				this);
		});
	}

	static void dirscan_callback(regex_grep_task * self,
		OBJECT * path,
		int found,
		timestamp const * const)
	{
		self->dirscan_file(path);
	}

	void dirscan_file(const value_ref & file)
	{
		// We only care about the filename of the path for globing. So split
		// the path to only get the file.
		std::string filepath = file->str();
		PATHNAME filepathparts(file->str());
		std::string filename(filepathparts.base() + filepathparts.suffix());
		// Ignore meta-dir paths.
		if (filename == "." || filename == "..") return;
		// Regular file.
		if (file_is_file(file) != 0)
		{
			// Match the full `path` to the set of glob patterns we are looking
			// for.
			for (auto glob_pattern : file_glob_patterns)
			{
				if (glob(glob_pattern->str(), filename.c_str()) == 0)
				{
					// We have a glob match for this file.
					auto do_grep = [this, filepath]() {
						auto filedata = file_preload(filepath);
						file_grep(filepath, *filedata);
					};
					do_grep();
					// grep_tasks->queue(do_grep);
				}
			}
			return;
		}
		// Subdir, or equivalent. If indicated, recurse scan subdirectories.
		if (recursive_glob)
		{
			dir_scan(file);
			return;
		}
	}

	std::unique_ptr<b2::filesys::file_buffer> file_preload(std::string filepath)
	{
		return std::unique_ptr<b2::filesys::file_buffer>(
			new b2::filesys::file_buffer(filepath));
	}

	void file_grep(std::string filepath, b2::filesys::file_buffer & filedata)
	{
		// WARNING: We need to avoid Jam operations in this. As we are getting
		// called from different threads. And the Jam memory is not thread-safe.

		// The match results are tuples of filepath+expressions. Collect all
		// those tuples for this file here.
		std::unique_ptr<std::vector<std::string>> result(
			new std::vector<std::string>);

		if (filedata.size() > 0)
		{
			// We have some data to regex search in.
			for (auto & prog : text_grep_prog)
			{
				auto grep_i = prog.search(filedata.begin(), filedata.end());
				for (; grep_i; ++grep_i)
				{
					// We need to add the file to the result (which is
					// followed by the match expressions).
					result->push_back(filepath);
					for (int i = 0; i < text_grep_result_expressions.size; ++i)
					{
						if (text_grep_result_expressions.get(i))
						{
							result->push_back(grep_i[i]);
						}
					}
				}
			}
		}

		// Append this file's results to the global results.
		if (!result->empty())
		{
			scope_lock_t lock(mx);
			intermediate.emplace_back(result.release());
		}
	}

	void wait() { grep_tasks->wait(); }
};

}} // namespace b2

list_ref b2::regex_grep(list_cref directories,
	list_cref files,
	list_cref patterns,
	list_cref result_expressions,
	list_cref options)
{
	// For the glob we always do a case insensitive compare. So we need
	// the globs as lower-case.
	list_ref globs;
	for (auto glob : files)
	{
		std::string g(glob->str());
		std::transform(g.begin(), g.end(), g.begin(),
			[](unsigned char c) { return std::tolower(c); });
		globs.push_back(value_ref(g));
	}

	// Extract the result expressions. Defaulting to the full match.
	flags<std::uint16_t> sub_expr;
	if (result_expressions.empty())
		sub_expr.set(0);
	else
	{
		for (auto result_expression : result_expressions)
		{
			if (result_expression->get_type() == b2::value::type::number)
				sub_expr.set(int(result_expression->as_number()));
			else if (result_expression->get_type() == b2::value::type::string)
				sub_expr.set(std::atoi(result_expression->str()));
		}
	}

	regex_grep_task task(list_cref(*globs), patterns, sub_expr, options);
	for (auto dir : directories)
	{
		task.dir_scan(dir);
	}
	task.wait();

	// Done with the search. Collect the intermediate results into our single
	// linear result.
	list_ref result;
	for (auto & intermediate : task.intermediate)
	{
		for (auto & r : *intermediate)
		{
			result.push_back(r);
		}
	}
	return result;
}

const char * b2::regex_module::init_code = R"jam(
rule __test__ ( )
{
    import assert ;

    assert.result a b c : split "a/b/c" / ;
    assert.result "" a b c : split "/a/b/c" / ;
    assert.result "" "" a b c : split "//a/b/c" / ;
    assert.result "" a "" b c : split "/a//b/c" / ;
    assert.result "" a "" b c "" : split "/a//b/c/" / ;
    assert.result "" a "" b c "" "" : split "/a//b/c//" / ;
    # assert.result "" a b c "" : split "abc" "" ;
    # assert.result "" "" : split "" "" ;
    assert.result "<x>y" "z" "<a>b" "c" "<d>e" "f" : split "<x>y\\z\\<a>b\\c\\<d>e\\f" "[\\/]" ;

    assert.result a b c "" a "" b c "" "" : split-list "a/b/c" "/a//b/c//" : / ;

    assert.result a c b d
        : match (.)(.)(.)(.) : abcd : 1 3 2 4 ;
    assert.result a b c d
        : match (.)(.)(.)(.) : abcd ;
    assert.result ababab cddc
        : match "((ab)*)([cd]+)" : abababcddc : 1 3 ;
    assert.result "a:"
        : match "(^.:)" : "a:/b" ;

    assert.result a.h c.h
        : transform <a.h> \"b.h\" <c.h> : <(.*)> ;

    assert.result a.h b.h c.h
        : transform <a.h> \"b.h\" <c.h> : "<([^>]*)>|\"([^\"]*)\"" : 1 2 ;

    assert.result "^<?xml version=\"1.0\"^>"
        : escape "<?xml version=\"1.0\">" : "&|()<>^" : "^" ;

    assert.result "<?xml version=\\\"1.0\\\">"
        : escape "<?xml version=\"1.0\">" : "\\\"" : "\\" ;

    assert.result "string&nbsp;string&nbsp;" : replace "string string " " " "&nbsp;" ;
    assert.result "&nbsp;string&nbsp;string" : replace " string string" " " "&nbsp;" ;
    assert.result "string&nbsp;&nbsp;string" : replace "string  string" " " "&nbsp;" ;
    assert.result "-" : replace "&" "&" "-" ;
    # assert.result "x" : replace "" "" "x" ;
    # assert.result "xax" : replace "a" "" "x" ;
    # assert.result "xaxbx" : replace "ab" "" "x" ;
    assert.result "a-b-c-d" : replace "a_b.c/d" "[_./]" "-" ;

    assert.result "-" "a-b" : replace-list "&" "a&b" : "&" : "-" ;
}
)jam";
