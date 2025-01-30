/*
Copyright 2001-2004 David Abrahams.
Copyright 2022 René Ferdinand Rivera Morell
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE.txt or https://www.bfgroup.xyz/b2/LICENSE.txt)
*/

#ifndef B2_FRAMES_H
#define B2_FRAMES_H

#include "config.h"
#include "lists.h"
#include "modules.h"
#include "object.h"
#include "value.h"

typedef struct frame FRAME;

void frame_init(FRAME *);
void frame_free(FRAME *);

struct frame
{
	FRAME * prev;
	FRAME * prev_user; /* The nearest enclosing frame for which
						  module->user_module is true. */
	LOL args[1];
	module_t * module;
	OBJECT * file;
	int line;
	char const * rulename;
#ifdef JAM_DEBUGGER
	void * function;
#endif

	inline frame() { frame_init(this); }
	inline ~frame() { frame_free(this); }
	inline frame(const frame & other)
		: prev(other.prev)
		, prev_user(other.prev_user)
		, module(other.module)
		, file(other.file)
		, line(other.line)
		, rulename(other.rulename)
#ifdef JAM_DEBUGGER
		, function(other.function)
#endif
	{
		lol_init(args);
		for (int32_t a = 0; a < other.args->count; ++a)
		{
			lol_add(args, list_copy(other.args->list[a]));
		}
	}
};

/* When a call into Python is in progress, this variable points to the bjam
 * frame that was current at the moment of the call. When the call completes,
 * the variable is not defined. Furthermore, if Jam calls Python which calls Jam
 * and so on, this variable only keeps the most recent Jam frame.
 */
extern FRAME * frame_before_python_call;

namespace b2 { namespace jam {

struct module_scope
{
	inline module_scope(FRAME * frame_, const char * module_name_)
		: module_frame(frame_)
	{
		if (module_frame)
		{
			module_saved = module_frame->module;
			module_frame->module
				= b2::ensure_valid(bindmodule(value_ref(module_name_)));
		}
	}

	inline ~module_scope()
	{
		if (module_frame && module_saved)
		{
			module_frame->module = b2::ensure_valid(module_saved);
		}
	}

	inline value_ref name() const { return module_frame->module->name; }
	inline FRAME * frame() const { return module_frame; }
	inline operator module_ptr() const { return module_frame->module; }

	private:
	FRAME * module_frame = nullptr;
	module_t * module_saved = nullptr;
};

struct module_scope_in_function
{
	module_scope_in_function(FRAME * frame_, const char * module_name_);
	~module_scope_in_function();
	inline value_ref name() const { return module_frame->module->name; }
	inline FRAME * frame() const { return module_frame; }
	inline operator module_ptr() const { return module_frame->module; }

	private:
	FRAME * module_frame = nullptr;
    module_t * saved_module = nullptr;
};

}} // namespace b2::jam

#endif
