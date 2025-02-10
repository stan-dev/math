/*
 *  Copyright 2022-2023 René Ferdinand Rivera Morell
 *  Copyright 2011 Steven Watanabe
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE.txt or
 * https://www.bfgroup.xyz/b2/LICENSE.txt)
 */

#ifndef B2_FUNCTION_H
#define B2_FUNCTION_H

#include "config.h"
#include "frames.h"
#include "jam_strings.h"
#include "lists.h"
#include "object.h"
#include "parse.h"

#include <functional>
#include <string>

typedef struct _function FUNCTION;

typedef std::function<LIST *(FRAME *, int32_t flags)> function_builtin_t;

typedef FUNCTION * function_ptr;

FUNCTION * function_compile(PARSE * parse);
FUNCTION * function_builtin(
	function_builtin_t func, int32_t flags, const char ** args);
void function_refer(FUNCTION *);
void function_free(FUNCTION *);
OBJECT * function_rulename(FUNCTION *);
void function_set_rulename(FUNCTION *, OBJECT *);
void function_location(FUNCTION *, OBJECT **, int32_t *);
LIST * function_run(FUNCTION * function, FRAME * frame);

FUNCTION * function_compile_actions(
	const char * actions, OBJECT * file, int32_t line);
void function_run_actions(FUNCTION * function, FRAME * frame, string * out);

FUNCTION * function_bind_variables(
	FUNCTION * f, module_t * module, int32_t * counter);
FUNCTION * function_unbind_variables(FUNCTION * f);

LIST * function_get_variables(FUNCTION * f);

void function_done(void);

namespace b2 { namespace jam {

struct backtrace
{
	static std::string to_string(frame * f);
	static list_ref to_list(frame * f);
};

}} // namespace b2::jam

#endif
