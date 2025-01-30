/*
 * Copyright 2022 René Ferdinand Rivera Morell
 * Copyright 2011 Steven Watanabe
 *
 * This file is part of Jam - see jam.c for Copyright information.
 */

#ifndef BOOST_JAM_OBJECT_H
#define BOOST_JAM_OBJECT_H

#include "config.h"
#include "value.h"

using OBJECT = b2::value;
inline b2::value_ptr object_new_range(const char * string, int32_t size)
{
	return b2::value::make(string, size);
}
inline b2::value_ptr object_new(const char * string)
{
	return b2::value::make(string);
}
inline b2::value_ptr object_copy(b2::value_ptr obj)
{
	return b2::value::copy(obj);
}
inline void object_free(b2::value_ptr & obj) { b2::value::free(obj); }
inline int object_equal(b2::value_ptr a, b2::value_ptr b)
{
	return a->equal_to(*b);
}
inline unsigned int object_hash(b2::value_ptr obj)
{
	return (unsigned int)(obj->hash64);
}
inline const char * object_str(b2::value_ptr obj) { return obj->str(); }
inline void object_done(void) { b2::value::done(); }

#endif
