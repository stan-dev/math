/*
Copyright 2022 René Ferdinand Rivera Morell
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE.txt or https://www.bfgroup.xyz/b2/LICENSE.txt)
*/

#include "mod_path.h"

#include "constants.h"
#include "filesys.h"
#include "pathsys.h"

namespace b2 { namespace paths {

bool exists(value_ref location) { return file_query(location) != nullptr; }

list_ref normalize_all(list_cref paths)
{
	list_ref result;
	for (auto path : paths) result.push_back(normalize(path->str()));
	return result;
}

}} // namespace b2::paths
