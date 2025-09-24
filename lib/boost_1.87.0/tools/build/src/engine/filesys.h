/*
 * Copyright 1993-2002 Christopher Seiwald and Perforce Software, Inc.
 *
 * This file is part of Jam - see jam.c for Copyright information.
 */

/* This file is ALSO:
 * Copyright 2023 René Ferdinand Rivera Morell.
 * Copyright 2001-2004 David Abrahams.
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE.txt or https://www.bfgroup.xyz/b2/LICENSE.txt)
 */

/*
 * filesys.h - OS specific file routines
 */

#ifndef FILESYS_DWA20011025_H
#define FILESYS_DWA20011025_H

#include "config.h"
#include "hash.h"
#include "lists.h"
#include "object.h"
#include "pathsys.h"
#include "timestamp.h"

#include <memory>
#include <string>

struct file_info_t
{
	OBJECT * name = nullptr;
	char is_file = 0;
	char is_dir = 0;
	char exists = 0;
	timestamp time;
	LIST * files = L0;
	file_info_t() = default;
	file_info_t(const file_info_t &) = default;
	file_info_t(OBJECT * n)
		: name(n)
	{}
	inline ~file_info_t()
	{
		if (name != nullptr) object_free(name);
		list_free(files);
	}
};

typedef struct file_item FILEITEM;
struct file_item
{
	file_info_t * value; /* expected to be equivalent with &FILEITEM */
	FILEITEM * next;
};

typedef struct file_list
{
	FILEITEM * head;
	FILEITEM * tail;
	int size;
} FILELIST;

typedef file_info_t ** FILELISTITER; /*  also &FILEITEM equivalent */

struct file_archive_info_t
{
	OBJECT * name = nullptr;
	file_info_t * file = nullptr;
	FILELIST * members = nullptr;
	file_archive_info_t() = default;
	file_archive_info_t(const file_archive_info_t &) = default;
	file_archive_info_t(OBJECT * n)
		: name(n)
	{}
	inline ~file_archive_info_t();
};

typedef void (*archive_scanback)(void * closure,
	OBJECT * path,
	LIST * symbols,
	int found,
	timestamp const * const);
typedef void (*scanback)(
	void * closure, OBJECT * path, int found, timestamp const * const);

void file_archscan(char const * arch, scanback func, void * closure);
void file_archivescan(OBJECT * path, archive_scanback func, void * closure);
void file_build1(PATHNAME * const f, string * file);
void file_dirscan(OBJECT * dir, scanback func, void * closure);
file_info_t * file_info(OBJECT * const path, int * found);
int file_is_file(OBJECT * const path);
int file_mkdir(char const * const path);
file_info_t * file_query(OBJECT * const path);
void file_remove_atexit(OBJECT * const path);
void file_supported_fmt_resolution(timestamp * const);
int file_time(OBJECT * const path, timestamp * const);

namespace b2 { namespace filesys {

inline bool is_file(const std::string & path)
{
	OBJECT * path_o = object_new(path.c_str());
	bool result = file_is_file(path_o) == 1;
	object_free(path_o);
	return result;
}

}} // namespace b2::filesys

/*  Archive/library file support */
file_archive_info_t * file_archive_info(OBJECT * const path, int * found);
file_archive_info_t * file_archive_query(OBJECT * path);

/* FILELIST linked-list */
FILELIST * filelist_new(OBJECT * path);
FILELIST * filelist_push_back(FILELIST * list, OBJECT * path);
FILELIST * filelist_push_front(FILELIST * list, OBJECT * path);
FILELIST * filelist_pop_front(FILELIST * list);
int filelist_length(FILELIST * list);
void filelist_free(FILELIST * list);

FILELISTITER filelist_begin(FILELIST * list);
FILELISTITER filelist_end(FILELIST * list);
FILELISTITER filelist_next(FILELISTITER it);
file_info_t * filelist_item(FILELISTITER it);
file_info_t * filelist_front(FILELIST * list);
file_info_t * filelist_back(FILELIST * list);

int filelist_empty(FILELIST * list);

inline file_archive_info_t::~file_archive_info_t()
{
	if (members != nullptr) filelist_free(members);
}

/* Internal utility worker functions. */
void file_query_posix_(file_info_t * const);

void file_done();

namespace b2 { namespace filesys {
class file_buffer
{
	public:
	bool is_memory_mapped = false;

	file_buffer(const std::string & filepath);
	~file_buffer();

	inline std::size_t size() const { return data_size; }
	inline const char * begin() const { return data_c.get(); }
	inline const char * end() const { return data_c.get() + size(); }

	private:
	std::unique_ptr<char[]> data_c;
	std::size_t data_size = 0;
	FILE * file = nullptr;
};
}} // namespace b2::filesys

#endif
