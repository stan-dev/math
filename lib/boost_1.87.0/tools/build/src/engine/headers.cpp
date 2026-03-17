/*
 * Copyright 1993, 2000 Christopher Seiwald.
 *
 * This file is part of Jam - see jam.c for Copyright information.
 */
/*  This file is ALSO:
 *  Copyright 2001-2004 David Abrahams.
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE.txt or https://www.bfgroup.xyz/b2/LICENSE.txt)
 */

/*
 * headers.c - handle #includes in source files
 *
 * Using regular expressions provided as the variable $(HDRSCAN), headers()
 * searches a file for #include files and phonies up a rule invocation:
 *    $(HDRRULE) <target> : <include files> ;
 *
 * External routines:
 *    headers() - scan a target for include files and call HDRRULE
 *
 * Internal routines:
 *    headers1() - using regexp, scan a file and build include LIST
 */

#include "jam.h"
#include "headers.h"

#include "compile.h"
#include "hdrmacro.h"
#include "lists.h"
#include "modules.h"
#include "object.h"
#include "parse.h"
#include "regexp.h"
#include "rules.h"
#include "variable.h"
#include "output.h"

#ifdef OPT_HEADER_CACHE_EXT
# include "hcache.h"
#endif

#include <errno.h>
#include <string.h>


/*
 * headers() - scan a target for include files and call HDRRULE
 */

#define MAXINC 10

void headers( TARGET * t )
{
    LIST   * hdrscan;
    LIST   * hdrrule;
    #ifndef OPT_HEADER_CACHE_EXT
    LIST   * headlist = L0;
    #endif
    b2::regex::program re_prog[MAXINC];
    int rec = 0;
    LISTITER iter;
    LISTITER end;

    hdrscan = var_get( root_module(), constant_HDRSCAN );
    if ( list_empty( hdrscan ) )
        return;

    hdrrule = var_get( root_module(), constant_HDRRULE );
    if ( list_empty( hdrrule ) )
        return;

    if ( DEBUG_HEADER )
        out_printf( "header scan %s\n", object_str( t->name ) );

    /* Compile all regular expressions in HDRSCAN */
    iter = list_begin( hdrscan );
    end = list_end( hdrscan );
    for ( ; ( rec < MAXINC ) && iter != end; iter = list_next( iter ) )
    {
        re_prog[ rec ].reset( list_item( iter )->str() );
        rec += 1;
    }

    /* Doctor up call to HDRRULE rule */
    /* Call headers1() to get LIST of included files. */
    {
        FRAME frame[ 1 ];
        frame_init( frame );
        lol_add( frame->args, list_new( object_copy( t->name ) ) );
#ifdef OPT_HEADER_CACHE_EXT
        lol_add( frame->args, hcache( t, rec, re_prog, hdrscan ) );
#else
        lol_add( frame->args, headers1( headlist, t->boundname, rec, re_prog ) );
#endif

        if ( lol_get( frame->args, 1 ) )
        {
            OBJECT * rulename = list_front( hdrrule );
            /* The third argument to HDRRULE is the bound name of $(<). */
            lol_add( frame->args, list_new( object_copy( t->boundname ) ) );
            list_free( evaluate_rule( bindrule( rulename, frame->module ), rulename, frame ) );
        }

        /* Clean up. */
        frame_free( frame );
    }
}


/*
 * headers1() - using regexp, scan a file and build include LIST.
 */

LIST * headers1( LIST * l, OBJECT * file, int rec, b2::regex::program re[] )
{
    FILE * f;
    char buf[ 1024 ];
    int i;

    /* The following regexp is used to detect cases where a file is included
     * through a line like "#include MACRO".
     */
    b2::regex::program re_macros(
        "#[ \t]*include[ \t]*([A-Za-z][A-Za-z0-9_]*).*$" );

#ifdef OPT_IMPROVED_PATIENCE_EXT
    static int count = 0;
    ++count;
    if ( ( ( count == 100 ) || !( count % 1000 ) ) && DEBUG_MAKE )
    {
        out_printf( "...patience...\n" );
        out_flush();
    }
#endif

    if ( file->as_string().size == 0 )
    {
        /* If the scanning was fed empty file names we just ignore them. */
        return l;
    }

    if ( !( f = fopen( object_str( file ), "r" ) ) )
    {
        /* No source files will be generated when -n flag is passed */
        if ( !globs.noexec || errno != ENOENT )
            err_printf( "[errno %d] failed to scan file '%s': %s",
                errno, object_str( file ), strerror(errno) );
        return l;
    }

    while ( fgets( buf, sizeof( buf ), f ) )
    {
        for ( i = 0; i < rec; ++i )
        {
            auto re_i = re [ i ].search( buf );
            if ( re_i && re_i[ 1 ].begin() )
            {
                std::string header(re_i[ 1 ].begin(), re_i[ 1 ].end());
                if ( DEBUG_HEADER )
                    out_printf( "header found: %s\n", header.c_str() );
                l = list_push_back( l, object_new( header.c_str() ) );
            }
        }

        /* Special treatment for #include MACRO. */
        auto re_macros_i = re_macros.search( buf );
        if ( re_macros_i && re_macros_i[ 1 ].end() != nullptr )
        {
            std::string macro_name(re_macros_i[ 1 ].begin(), re_macros_i[ 1 ].end());

            if ( DEBUG_HEADER )
                out_printf( "macro header found: %s", macro_name.c_str() );

            b2::value_ref macro_name_v(macro_name);
            b2::value_ref header_filename_v(macro_header_get( macro_name_v ));
            if ( header_filename_v.has_value() )
            {
                if ( DEBUG_HEADER )
                    out_printf( " resolved to '%s'\n", header_filename_v->str()
                        );
                l = list_push_back( l, header_filename_v );
            }
            else
            {
                if ( DEBUG_HEADER )
                    out_printf( " ignored !!\n" );
            }
        }
    }

    fclose( f );
    return l;
}
