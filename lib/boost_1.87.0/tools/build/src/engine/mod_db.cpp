/*
Copyright 2024 René Ferdinand Rivera Morell
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE.txt or https://www.bfgroup.xyz/b2/LICENSE.txt)
*/

#include "mod_db.h"

#include "ext_nlohmann_json.hpp"

#include <cstdio>
#include <cstdlib>

void b2::property_db::emplace(list_cref k, value_ref v)
{
	list_ref nk;
	for (list_cref::size_type i = 0; i < k.size(); ++i)
	{
		value_ref a = k[i];
		if (a == array_tag)
		{
			nk.push_back(array_tag);
			i += 1;
			if (i < k.size())
			{
				// Normalize array indices to number type for proper ordering.
				value_ref n = k[i];
				if (n->get_type() != value::type::number)
				{
					double index = list_cref::size_type(
						std::atof(n->to_string()->str()));
					n = value::make(index);
				}
				nk.push_back(n->as_number());
			}
		}
		else
		{
			nk.push_back(a);
		}
	}
	db[nk] = v;
}

bool b2::property_db::write_file(value_ref filename, value_ref format)
{
	if (!format->has_value())
	{
		format = b2::value::make("json");
	}
	if (format == "json") return write_file_json(filename);
	return false;
}

std::string b2::property_db::dump(value_ref format)
{
	if (!format->has_value())
	{
		format = b2::value::make("json");
	}
	if (format == "json") return dump_json();
	return "";
}

namespace {
template <class T>
void build_json_from_db(const T & db, nlohmann::json & out)
{
	b2::value_ref array_tag { "[]" };
	for (auto & item : db)
	{
		auto & k = item.first;
		auto & v = item.second;
		if (v->get_type() != b2::value::type::number
			&& v->get_type() != b2::value::type::string)
			continue;
		nlohmann::json * current = &out;
		for (b2::list_cref::size_type i = 0; i < k.size(); ++i)
		{
			b2::value_ref a = k[i];
			if (a == array_tag)
			{
				if (++i >= k.size()) break;
				auto index = b2::list_cref::size_type(k[i]->as_number());
				if (current->is_null()) *current = nlohmann::json::array();
				if (!current->is_array()) break;
				while (current->size() <= nlohmann::json::size_type(index))
					current->push_back(nlohmann::json());
				current = &(*current)[index];
			}
			else
			{
				if (current->is_null()) *current = nlohmann::json::object();
				current = &(*current)[a->str()];
			}
		}
		if (v->get_type() == b2::value::type::number)
			*current = nlohmann::json(v->as_number());
		else if (v->get_type() == b2::value::type::string)
			*current = nlohmann::json(v->str());
	}
}
} // namespace

bool b2::property_db::write_file_json(value_ref filename)
{
	nlohmann::json out;
	build_json_from_db(db, out);
	FILE * file = std::fopen(filename->str(), "w");
	if (!file) return false;
	bool result = false;
	try
	{
		auto data = out.dump(0);
		result = std::fwrite(data.c_str(), data.size(), 1, file) == 1;
	}
	catch (const std::exception &)
	{}
	std::fclose(file);
	return result;
}

std::string b2::property_db::dump_json()
{
	nlohmann::json out;
	build_json_from_db(db, out);
	try
	{
		return out.dump(-1);
	}
	catch (const std::exception &)
	{}
	return "";
}

const char * b2::db_module::init_code = R"jam(
rule __test__ ( )
{
	import assert ;
	import "class" : new ;

	local pdb = [ new property-db ] ;
	assert.result "null" : $(pdb).dump "json" ;
	$(pdb).emplace "[]" 1 "one" : 1 ;
	assert.result "[null,{\"one\":\"1\"}]" : $(pdb).dump "json" ;
	$(pdb).emplace "[]" 0 "zero" "sub" : "qwerty" ;
	assert.result "[{\"zero\":{\"sub\":\"qwerty\"}},{\"one\":\"1\"}]" : $(pdb).dump "json" ;
}
)jam";
