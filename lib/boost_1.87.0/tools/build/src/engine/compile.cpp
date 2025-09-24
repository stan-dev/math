/*
 * Copyright 1993, 2000 Christopher Seiwald.
 *
 * This file is part of Jam - see jam.c for Copyright information.
 */

/*  This file is ALSO:
 *  Copyright 2022 René Ferdinand Rivera Morell
 *  Copyright 2001-2004 David Abrahams.
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE.txt or https://www.bfgroup.xyz/b2/LICENSE.txt)
 */

/*
 * compile.c - compile parsed jam statements
 *
 * External routines:
 *  evaluate_rule() - execute a rule invocation
 *
 * Internal routines:
 *  debug_compile() - printf with indent to show rule expansion
 */

#include "jam.h"
#include "compile.h"

#include "builtins.h"
#include "class.h"
#include "constants.h"
#include "hash.h"
#include "hdrmacro.h"
#include "make.h"
#include "modules.h"
#include "parse.h"
#include "rules.h"
#include "search.h"
#include "jam_strings.h"
#include "variable.h"
#include "output.h"
#include "startup.h"
#include "value.h"

#include <assert.h>
#include <stdarg.h>
#include <string.h>

#include <string>


static void debug_compile( int which, char const * s, FRAME * );

/* Internal functions from builtins.c */
void backtrace( FRAME * );
void backtrace_line( FRAME * );
void print_source_line( FRAME * );
void unknown_rule( FRAME *, char const * key, module_t *, OBJECT * rule_name );


/*
 * evaluate_rule() - execute a rule invocation
 */

LIST * evaluate_rule( RULE * rule, OBJECT * rulename, FRAME * frame )
{
    b2::list_ref result;
    profile_frame   prof[ 1 ];
    module_t      * prev_module = frame->module;

    if ( DEBUG_COMPILE )
    {
        /* Try hard to indicate in which module the rule is going to execute. */
        char buf[ 256 ] = "";
        if ( rule->module->name )
        {
            strncat( buf, object_str( rule->module->name ), sizeof( buf ) -
                1 );
            strncat( buf, ".", sizeof( buf ) - 1 );
            if ( strncmp( buf, object_str( rule->name ), strlen( buf ) ) == 0 )
            {
                buf[ 0 ] = 0;
            }
        }
        strncat( buf, object_str( rule->name ), sizeof( buf ) - 1 );
        debug_compile( 1, buf, frame );

        lol_print( frame->args );
        out_printf( "\n" );
    }

    if ( rule->procedure && rule->module != prev_module )
    {
        /* Propagate current module to nested rule invocations. */
        frame->module = b2::ensure_valid(rule->module);
    }

    /* Record current rule name in frame. */
    if ( rule->procedure )
    {
        frame->rulename = object_str( rulename );
        /* And enter record profile info. */
        if ( DEBUG_PROFILE )
            profile_enter( function_rulename( rule->procedure ), prof );
    }

    /* Check traditional targets $(<) and sources $(>). */
    if ( !rule->actions && !rule->procedure )
        unknown_rule( frame, NULL, frame->module, rulename );

    /* If this rule will be executed for updating the targets then construct the
     * action for make().
     */
    if ( rule->actions )
    {
        targets_ptr t;

        /* The action is associated with this instance of this rule. */
        ACTION * const action = b2::jam::make_ptr<ACTION>();

        action->rule = rule;
        action->targets.reset(); targetlist( action->targets, lol_get( frame->args, 0 ) );
        action->sources.reset(); targetlist( action->sources, lol_get( frame->args, 1 ) );
        action->refs = 1;

        /* If we have a group of targets all being built using the same action
         * and any of these targets is updated, then we have to consider them
         * all to be out-dated.  We do this by adding a REBUILDS in both directions
         * between the first target and all the other targets.
         */
        if ( action->targets )
        {
            TARGET * const t0 = action->targets->target;
            for ( t = action->targets->next.get(); t; t = t->next.get() )
            {
                targetentry( t->target->rebuilds, t0 );
                targetentry( t0->rebuilds, t->target );
            }
        }

        /* Append this action to the actions of each target. */
        for ( t = action->targets.get(); t; t = t->next.get() )
            t->target->actions = actionlist( t->target->actions, action );

        action_free( action );
    }

    /* Now recursively compile any parse tree associated with this rule.
     * function_refer()/function_free() call pair added to ensure the rule does
     * not get freed while in use.
     */
    if ( rule->procedure )
    {
        auto function = b2::jam::make_unique_bare_jptr( rule->procedure, function_refer, function_free );
        result.reset( function_run( function.get(), frame ) );
    }

    if ( DEBUG_PROFILE && rule->procedure )
        profile_exit( prof );

    if ( DEBUG_COMPILE )
        debug_compile( -1, 0, frame );

    return result.release();
}


/*
 * Call the given rule with the specified parameters. The parameters should be
 * of type LIST* and end with a NULL pointer. This differs from 'evaluate_rule'
 * in that frame for the called rule is prepared inside 'call_rule'.
 *
 * This function is useful when a builtin rule (in C) wants to call another rule
 * which might be implemented in Jam.
 */

LIST * call_rule( OBJECT * rulename, FRAME * caller_frame, LOL * args )
{
    FRAME inner[ 1 ];
    frame_init( inner );
    inner->prev = caller_frame;
    inner->prev_user = caller_frame->module->user_module
        ? caller_frame
        : caller_frame->prev_user;
    inner->module = b2::ensure_valid(caller_frame->module);

    for ( int32_t a = 0; a < args->count; ++a)
    {
        lol_add( inner->args, list_copy(lol_get(args, a)) );
    }

    b2::list_ref result( evaluate_rule( bindrule( rulename, inner->module ), rulename, inner ), true );

    frame_free( inner );

    return result.release();
}

LIST * call_rule( OBJECT * rulename, FRAME * caller_frame, LIST * arg, ... )
{
    b2::lists args;
    if ( arg )
    {
        args.push_back( b2::list_cref( arg ) );
        va_list va;
        va_start( va, arg );
        for ( ; ; )
        {
            LIST * l = va_arg( va, LIST * );
            if ( !l )
                break;
            args.push_back( b2::list_cref( l ) );
        }
        va_end( va );
    }

    return call_rule( rulename, caller_frame, args );
}


LIST * call_member_rule(
	OBJECT * rulename, FRAME * caller_frame, b2::list_ref && self_, b2::lists && args_)
{
    b2::list_ref self(std::move(self_));
    b2::lists args(std::move(args_));
    if (self.empty())
    {
        backtrace_line(caller_frame);
        out_printf("warning: object is empty\n");
        backtrace(caller_frame);
        return L0;
    }

    /* FIXME: handle generic case */
    assert( self.length() == 1 );

    module_ptr module = bindmodule(self[0]);
    rule_ptr rule = nullptr;
    b2::value_ref real_rulename;

    if (module->class_module)
    {
        rule = bindrule(rulename, module);
        if (rule->procedure)
        {
            real_rulename = b2::value_ref(function_rulename(rule->procedure));
        }
        else
        {
            real_rulename = std::string(module->name->str())+"."+rulename->str();
        }
    }
    else
    {
        real_rulename = std::string(self[0]->str())+"."+rulename->str();
        rule = bindrule(real_rulename, caller_frame->module);
    }

    FRAME inner[ 1 ];
    frame_init( inner );
    inner->prev = caller_frame;
    inner->prev_user = caller_frame->module->user_module
        ? caller_frame : caller_frame->prev_user;
    inner->module = b2::ensure_valid(caller_frame->module);

    args.swap( inner->args[0] );

    if (self.length() > 1)
    {
        b2::list_ref trailing;
        for (int i = 1; i < self.length(); ++i)
            trailing.push_back(std::string(self[i]->str())+"."+rulename->str());
        if (inner->args->count == 0)
            lol_add(inner->args, trailing.release());
        else
        {
            trailing.append(b2::list_cref(inner->args->list[0]));
            list_free(inner->args->list[0]);
            inner->args->list[0] = trailing.release();
        }
    }

    return evaluate_rule(rule, real_rulename, inner);
}


/*
 * debug_compile() - printf with indent to show rule expansion
 */

static void debug_compile( int which, char const * s, FRAME * frame )
{
    static int level = 0;
    static char indent[ 36 ] = ">>>>|>>>>|>>>>|>>>>|>>>>|>>>>|>>>>|";

    if ( which >= 0 )
    {
        int i;

        print_source_line( frame );

        i = ( level + 1 ) * 2;
        while ( i > 35 )
        {
            out_puts( indent );
            i -= 35;
        }

        out_printf( "%*.*s ", i, i, indent );
    }

    if ( s )
        out_printf( "%s ", s );

    level += which;
}
