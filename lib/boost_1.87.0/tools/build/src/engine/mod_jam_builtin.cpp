/*
Copyright 2022 René Ferdinand Rivera Morell
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE.txt or https://www.bfgroup.xyz/b2/LICENSE.txt)
*/

#include "mod_jam_builtin.h"

#include "bindjam.h"
#include "mod_jam_errors.h"
#include "patchlevel.h"
#include "startup.h"

#include <string>
#include <tuple>

namespace b2 { namespace jam { namespace builtin {

namespace {
inline int get_version(const char *& s)
{
	try
	{
		std::size_t len = 0;
		auto val = std::stoi(s, &len);
		s += len;
		if (*s == '.') s += 1;
		return val < 0 ? 0 : val;
	}
	catch (...)
	{}
	return 0;
}
inline std::tuple<int, int, int> get_version_tuple(const char * s)
{
	std::tuple<int, int, int> v;
	std::get<0>(v) = get_version(s);
	std::get<1>(v) = get_version(s);
	std::get<2>(v) = get_version(s);
	return v;
}
inline std::string get_version_str(const char * v)
{
	auto version = get_version_tuple(v);
	return std::to_string(std::get<0>(version)) + "."
		+ std::to_string(std::get<1>(version)) + "."
		+ std::to_string(std::get<2>(version));
}
} // namespace

void require_b2(
	value_ref minimum, value_ref maximum, bind::context_ref_ context_ref)
{
	const std::tuple<int, int, int> b2_version { VERSION_MAJOR, VERSION_MINOR,
		VERSION_PATCH };
	bool success = true;
	if (maximum.has_value())
	{
		auto min_v = get_version_tuple(minimum->str());
		auto max_v = get_version_tuple(maximum->str());
		success = min_v <= b2_version && b2_version < max_v;
	}
	else
	{
		auto min_v = get_version_tuple(minimum->str());
		success = min_v <= b2_version;
	}
	if (!success)
	{
		std::string error = "B2 version ";
		if (maximum.has_value())
			error += std::string("[") + get_version_str(minimum->str()) + ","
				+ get_version_str(maximum->str()) + ")";
		else
			error += get_version_str(minimum->str()) + " or greater";
		error += " required. But this is B2 ";
		error += std::to_string(VERSION_MAJOR) + ".";
		error += std::to_string(VERSION_MINOR) + ".";
		error += std::to_string(VERSION_PATCH) + ".";
		errors::error(lists() | list_ref(error), context_ref);
	}
}

const char * root_module::init_code = R"jam(
)jam";

}}} // namespace b2::jam::builtin
