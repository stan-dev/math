/*
 *  Copyright 2011 Steven Watanabe
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE.txt or https://www.bfgroup.xyz/b2/LICENSE.txt)
 */

#ifndef FUNCTION_SW20111123_H
#define FUNCTION_SW20111123_H

#include "config.h"
#include "object.h"
#include "frames.h"
#include "lists.h"
#include "parse.h"
#include "jam_strings.h"

typedef struct _function FUNCTION;
typedef struct _stack STACK;

STACK * stack_global( void );
void stack_push( STACK * s, LIST * l );
LIST * stack_pop( STACK * s );

FUNCTION * function_compile( PARSE * parse );
FUNCTION * function_builtin( LIST * ( * func )( FRAME * frame, int32_t flags ), int32_t flags, const char * * args );
void function_refer( FUNCTION * );
void function_free( FUNCTION * );
OBJECT * function_rulename( FUNCTION * );
void function_set_rulename( FUNCTION *, OBJECT * );
void function_location( FUNCTION *, OBJECT * *, int32_t * );
LIST * function_run( FUNCTION * function, FRAME * frame, STACK * s );

FUNCTION * function_compile_actions( const char * actions, OBJECT * file, int32_t line );
void function_run_actions( FUNCTION * function, FRAME * frame, STACK * s, string * out );

FUNCTION * function_bind_variables( FUNCTION * f, module_t * module, int32_t * counter );
FUNCTION * function_unbind_variables( FUNCTION * f );

LIST * function_get_variables( FUNCTION * f );

void function_done( void );

#ifdef HAVE_PYTHON

FUNCTION * function_python( PyObject * function, PyObject * bjam_signature );

#endif

#endif
