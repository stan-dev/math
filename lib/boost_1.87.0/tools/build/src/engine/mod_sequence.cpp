/*
 * Copyright Vladimir Prus 2003.
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE.txt or copy at
 * https://www.bfgroup.xyz/b2/LICENSE.txt)
 */

#include "native.h"
#include "object.h"
#include "lists.h"
#include "compile.h"

#include <stdlib.h>
#include <algorithm>
#include <iterator>
#include <unordered_set>
#include <set>


#ifndef max
# define max(a,b) ((a)>(b)?(a):(b))
#endif


LIST * sequence_select_highest_ranked( FRAME * frame, int flags )
{
   /* Returns all of 'elements' for which corresponding element in parallel */
   /* list 'rank' is equal to the maximum value in 'rank'.                  */

    LIST * const elements = lol_get( frame->args, 0 );
    LIST * const rank = lol_get( frame->args, 1 );

    LIST * result = L0;
    int highest_rank = -1;

    {
        LISTITER iter = list_begin( rank );
        LISTITER const end = list_end( rank );
        for ( ; iter != end; iter = list_next( iter ) )
        {
            int const current = atoi( object_str( list_item( iter ) ) );
            highest_rank = max( highest_rank, current );
        }
    }

    {
        LISTITER iter = list_begin( rank );
        LISTITER const end = list_end( rank );
        LISTITER elements_iter = list_begin( elements );
        for ( ; iter != end; iter = list_next( iter ), elements_iter =
            list_next( elements_iter ) )
            if ( atoi( object_str( list_item( iter ) ) ) == highest_rank )
                result = list_push_back( result, object_copy( list_item(
                    elements_iter ) ) );
    }

    return result;
}

LIST * sequence_transform( FRAME * frame, int flags )
{
    LIST * function = lol_get( frame->args, 0 );
    LIST * sequence = lol_get( frame->args, 1 );
    LIST * result = L0;
    OBJECT * function_name = list_front( function );
    LISTITER args_begin = list_next( list_begin( function ) ), args_end = list_end( function );
    LISTITER iter = list_begin( sequence ), end = list_end( sequence );
    RULE * rule = bindrule( function_name, frame->prev->module );

    for ( ; iter != end; iter = list_next( iter ) )
    {
        FRAME inner[ 1 ];

        frame_init( inner );
        inner->prev = frame;
        inner->prev_user = frame->prev_user;
        inner->module = b2::ensure_valid(frame->prev->module);

        lol_add( inner->args, list_push_back( list_copy_range( function, args_begin, args_end ), object_copy( list_item( iter ) ) ) );
        result = list_append( result, evaluate_rule( rule, function_name, inner ) );

        frame_free( inner );
    }

    return result;
}

namespace {

using list_value_ref = b2::list_cref::value_type;
using list_value_cref = const b2::list_cref::value_type;

struct hash_functor
{
	inline std::size_t operator()(list_value_cref a) const
	{
		return std::size_t(-1) & a->hash64;
	}
};

struct equal_functor
{
	inline bool operator()(list_value_cref a, list_value_cref b) const
	{
		return (a->hash64 == b->hash64) && a->equal_to(*b);
	}
};

struct less_functor
{
	inline bool operator()(list_value_cref a, list_value_cref b) const
	{
		return std::strcmp(a->str(), b->str()) < 0;
	}
};

} // namespace

static b2::list_ref sorted_unique(b2::list_cref values)
{
	std::set<list_value_ref, less_functor> unique_values { values.begin(),
		values.end() };
	b2::list_ref result;
	// result.reserve(unique_values.size()); // TODO: no reserve in list_ref
	for (auto & v : unique_values) result.push_back(v);
	return result;
}

static b2::list_ref stable_unique(b2::list_cref values)
{
	b2::list_ref result;
	// result.reserve(values.size()); // TODO: no reserve in list_ref
	std::unordered_set<list_value_ref, hash_functor, equal_functor>
		unique_values(values.size());
	std::copy_if(values.begin(), values.end(), std::back_inserter(result),
		[&unique_values](
			list_value_cref ptr) { return unique_values.insert(ptr).second; });
	return result;
}

static LIST * sequence_unique(FRAME * frame, int /*flags*/)
{
	b2::list_cref values(lol_get(frame->args, 0));
	LIST * stable = lol_get(frame->args, 1);

	if (!list_empty(stable))
		return stable_unique(values).release();
	else
		// Sorting is done because Boost unintentionally relies on it :-(
		// See https://github.com/boostorg/fiber/issues/305
		return sorted_unique(values).release();
}

void init_sequence()
{
    {
        char const * args[] = { "elements", "*", ":", "rank", "*", 0 };
        declare_native_rule( "sequence", "select-highest-ranked", args,
                            sequence_select_highest_ranked, 1 );
    }
    {
        char const * args[] = { "function", "+", ":", "sequence", "*", 0 };
        declare_native_rule( "sequence", "transform", args,
                            sequence_transform, 1 );
    }
    {
        char const * args[] = { "list", "*", ":", "stable", "?", 0 };
        declare_native_rule( "sequence", "unique", args,
                            sequence_unique, 1 );
    }
}
