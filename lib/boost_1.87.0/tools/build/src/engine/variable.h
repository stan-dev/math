/*
 * Copyright 1993, 2000 Christopher Seiwald.
 *
 * This file is part of Jam - see jam.c for Copyright information.
 */
/*
Copyright 2022 René Ferdinand Rivera Morell
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE.txt or https://www.bfgroup.xyz/b2/LICENSE.txt)
*/

/*
 * variable.h - handle jam multi-element variables
 */

#ifndef VARIABLE_SW20111119_H
#define VARIABLE_SW20111119_H

#include "config.h"
#include "lists.h"
#include "modules.h"
#include "object.h"
#include "value.h"

#include <string>

struct module_t;

void var_defines(struct module_t *, const char * const * e, int preprocess);
LIST * var_get(struct module_t *, OBJECT * symbol);
void var_set(struct module_t *, OBJECT * symbol, LIST * value, int flag);
LIST * var_swap(struct module_t *, OBJECT * symbol, LIST * value);
void var_done(struct module_t *);

/*
 * Defines for var_set().
 */

#define VAR_SET 0 /* override previous value */
#define VAR_APPEND 1 /* append to previous value */
#define VAR_DEFAULT 2 /* set only if no previous value */

namespace b2 { namespace jam {
struct variable
{
	inline variable(const variable & v)
		: var_module(v.var_module)
		, var_symbol(object_copy(v.var_symbol))
	{}
	inline ~variable()
	{
		if (var_symbol) object_free(var_symbol);
	}

	// Refer to variable in given module `m`.
	inline variable(struct module_t * m, const char * v)
		: var_module(m)
		, var_symbol(object_new(v))
	{}
	inline variable(const char * m, const char * v)
	{
		var_module = root_module();
		if (m && m[0] != '\0')
		{
			OBJECT * module_sym = object_new(m);
			var_module = bindmodule(module_sym);
			object_free(module_sym);
		}
		var_symbol = object_new(v);
	}
	inline variable(const std::string & m, const std::string & v)
		: variable(m.c_str(), v.c_str())
	{}

	// Refer to variable in global/root module.
	inline explicit variable(const char * v)
		: variable(root_module(), v)
	{}
	inline explicit variable(const std::string & v)
		: variable(v.c_str())
	{}

	inline list_cref get() const
	{
		return list_cref { var_get(var_module, var_symbol) };
	}
	inline list_cref operator*() const { return this->get(); }
	inline operator list_cref() const { return this->get(); }

	// Get the value at index `i`.
	inline value_ref operator[](int i) const { return this->get()[i]; }

	// Assign, i.e. replace, variable value.
	inline variable & operator=(list_ref && v)
	{
		list_ref vv(std::move(v));
		var_set(var_module, var_symbol, vv.release(), VAR_SET);
		return *this;
	}
	inline variable & operator=(const list_ref & v)
	{
		return *this = list_ref { v };
	}
	inline variable & operator=(const char * v)
	{
		return *this = list_ref { value_ref { v } };
	}
	inline variable & operator=(const std::string & v)
	{
		return *this = list_ref { value_ref { v } };
	}
	inline variable & operator=(const variable & other)
	{
		return *this = *other;
	}

	// Append to current variable value.
	inline variable & operator+=(list_ref & v)
	{
		var_set(var_module, var_symbol, v.release(), VAR_APPEND);
		return *this;
	}
	inline variable & operator+=(list_ref && v)
	{
		list_ref l = std::move(v);
		return *this += l;
	}
	template <class T>
	inline variable & operator+=(T v)
	{
		return *this += list_ref { value_ref { v } };
	}

	// Set variable value if empty.
	inline variable & operator|=(list_ref & v)
	{
		var_set(var_module, var_symbol, v.release(), VAR_DEFAULT);
		return *this;
	}
	inline variable & operator|=(list_ref && v)
	{
		var_set(var_module, var_symbol, v.release(), VAR_DEFAULT);
		return *this;
	}
	template <class T>
	inline variable & operator|=(T v)
	{
		return *this |= list_ref { value_ref { v } };
	}

	inline operator bool() const
	{
		LIST * l = var_get(var_module, var_symbol);
		return (!list_empty(l)) && (list_length(l) > 0);
	}

	private:
	struct module_t * var_module = nullptr;
	OBJECT * var_symbol = nullptr;
};
}} // namespace b2::jam

#endif
