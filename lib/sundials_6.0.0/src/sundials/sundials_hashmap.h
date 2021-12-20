/* -----------------------------------------------------------------
 * Programmer: Cody J. Balos @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * A simple header-only hashmap implementation for char* keys and
 * void* values. Uses linear probing to resolve collisions.
 * The values can be anything, but will be freed by
 * the hash map upon its destruction.
 * -----------------------------------------------------------------*/

#ifndef _SUNDIALS_HASHMAP_H
#define _SUNDIALS_HASHMAP_H

#include <string.h>
#include <stdlib.h>

/*
  For a nice discussion on popular hashing algorithms see:
  https://softwareengineering.stackexchange.com/questions/49550/which-hashing-algorithm-is-best-for-uniqueness-and-speed/145633#145633

  This is a 64-bit implementation of the 'a' modification of the
  Fowler–Noll–Vo hash (i.e., FNV1-a).
 */
static unsigned long fnv1a_hash(const char* str)
{
  const unsigned long prime = 14695981039346656037U;   /* prime */
  unsigned long hash        = 1099511628211U;          /* offset basis */
  char c;
  while ((c = *str++))
    hash = (hash^c) * prime;
  return hash;
}

typedef struct _SUNHashMapKeyValue* SUNHashMapKeyValue;

struct _SUNHashMapKeyValue
{
  const char* key;
  void*       value;
};

typedef struct _SUNHashMap *SUNHashMap;

struct _SUNHashMap
{
  int                 size;     /* current number of entries */
  int                 max_size; /* max number of entries */
  SUNHashMapKeyValue* buckets;
};

/*
  This function creates a new SUNHashMap object allocated to hold
  up to 'max_size' entries.

  **Arguments:**
    * ``max_size`` -- the max number of entries in the hashmap
    * ``map`` -- on input, a SUNHasMap pointer, on output the SUNHashMap will be
                 allocated

  **Returns:**
    * ``0`` -- success
    * ``-1`` -- an error occurred
 */
static int SUNHashMap_New(int max_size, SUNHashMap* map)
{
  int i;

  if (max_size <= 0)
    return(-1);

  *map = NULL;
  *map = (SUNHashMap) malloc(sizeof(struct _SUNHashMap));

  if (map == NULL)
    return(-1);

  (*map)->size = 0;
  (*map)->max_size = max_size;

  (*map)->buckets = NULL;
  (*map)->buckets = (SUNHashMapKeyValue*) malloc(max_size*sizeof(SUNHashMapKeyValue));

  if ((*map)->buckets == NULL)
  {
    free(*map);
    return(-1);
  }

  /* Initialize all buckets to NULL */
  for (i = 0; i < max_size; i++)
    (*map)->buckets[i] = NULL;

  return(0);
}

/*
  This function frees the SUNHashMap object.

  **Arguments:**
    * ``map`` -- on input, a SUNHasMap pointer, on output the SUNHashMap will be
                 deallocated and set to ``NULL``
    * ``freevalue`` -- callback function that should free the value object

  **Returns:**
    * ``0`` -- success
    * ``-1`` -- an error occurred
 */
static int SUNHashMap_Free(SUNHashMap* map, void (*freevalue)(void* ptr))
{
  int i;

  if (map == NULL || freevalue == NULL)
    return(-1);

  for (i = 0; i < (*map)->max_size; i++)
  {
    if ((*map)->buckets[i] && (*map)->buckets[i]->value)
      freevalue((*map)->buckets[i]->value);

    if ((*map)->buckets[i]) free((*map)->buckets[i]);
  }
  if ((*map)->buckets) free((*map)->buckets);
  if (*map) free(*map);
  *map = NULL;

  return(0);
}

/*
  This function iterates the map over the range [start, N]. N is either the
  index at which ``yieldfn`` indicates the iteration should stop, or the max
  entries in the map.

  **Arguments:**
    * ``map`` -- the ``SUNHashMap`` object to operate on
    * ``start`` -- the start of the iteration range
    * ``yieldfn`` -- the callback function to call every iteration
                     this should return -1 to continue the iteration, or >= 0 to
                     stop; the first argument is the current index, the second
                     argument is the current key-value pair, and the final
                     argument is the same pointer ``ctx`` as the final argument
                     to SUNHashMapIterate.
    * ``ctx`` -- a pointer to pass on to ``yieldfn``

  **Returns:**
    * ``max_size`` -- iterated the whole map
    * ``>=0`` -- the index at which the iteration stopped
    * ``<-1`` -- an error occurred
 */
static int SUNHashMap_Iterate(SUNHashMap map, int start,
                              int (*yieldfn)(int, SUNHashMapKeyValue, void*),
                              void* ctx)
{
  int i;

  if (map == NULL || yieldfn == NULL)
    return(-2);

  for (i = start; i < map->max_size; i++)
  {
    int retval = yieldfn(i, map->buckets[i], ctx);
    if (retval >= 0) return(retval); /* yieldfn indicates the loop should break */
    if (retval < -1) return(retval); /* error occurred */
  }

  return(map->max_size);
}

static int sunHashMapLinearProbeInsert(int idx, SUNHashMapKeyValue kv, void* ctx)
{
  /* find the next open spot */
  if (kv == NULL)
    return(idx); /* open spot found at idx */
  return(-1); /* keep looking */
}

/*
  This function creates a key-value pair and attempts to insert it into the map.
  Will use linear probing if there is a collision.

  **Arguments:**
    * ``map`` -- the ``SUNHashMap`` object to operate on
    * ``key`` -- the key to store
    * ``value`` -- the value associated with the key

  **Returns:**
    * ``0`` -- success
    * ``-1`` -- an error occurred
    * ``-2`` -- the map is full
 */
static int SUNHashMap_Insert(SUNHashMap map, const char* key, void* value)
{
  int                idx, retval;
  SUNHashMapKeyValue kvp;

  if (map == NULL || key == NULL || value == NULL)
    return(-1);

  /* We want the index to be in (0, map->max_size) */
  idx = (int) (fnv1a_hash(key) % map->max_size);

  /* Check if the bucket is already filled */
  if (map->buckets[idx] != NULL)
  {
    /* Find the next open spot */
    retval = SUNHashMap_Iterate(map, idx, sunHashMapLinearProbeInsert, NULL);
    if (retval < 0)
      return(-1); /* error occurred */
    else if (retval == map->max_size)
      return(-2); /* no open entry */
    else
      idx = retval;
  }

  /* Create the key-value pair */
  kvp = (SUNHashMapKeyValue) malloc(sizeof(struct _SUNHashMapKeyValue));
  if (kvp == NULL)
    return(-1);

  kvp->key = key;
  kvp->value = value;

  /* Insert the key-value pair */
  map->buckets[idx] = kvp;
  map->size++;

  return(0);
}

static int sunHashMapLinearProbeGet(int idx, SUNHashMapKeyValue kv, void* key)
{
  /* target key cannot be NULL */
  if (key == NULL)
    return(-2);

  /* find the matching entry */
  if (kv == NULL)
    return(-1); /* keep looking since this bucket is empty */
  else if (!strcmp((const char*) kv->key, (const char*) key))
    return(idx); /* found it at idx */
  else
    return(-1); /* keep looking */
}

/*
  This function gets the value for the given key.

  **Arguments:**
    * ``map`` -- the ``SUNHashMap`` object to operate on
    * ``key`` -- the key to look up
    * ``value`` -- the value associated with the key

  **Returns:**
    * ``0`` -- success
    * ``-1`` -- an error occurred
    * ``-2`` -- key not found
 */
static int SUNHashMap_GetValue(SUNHashMap map, const char* key, void** value)
{
  int idx, retval;

  if (map == NULL || key == NULL || value == NULL)
    return(-1);

  /* We want the index to be in (0, map->max_size) */
  idx = (int) (fnv1a_hash(key) % map->max_size);

  /* Check if the key exists */
  if (map->buckets[idx] == NULL)
    return(-2);

  /* Check to see if this is a collision */
  if (strcmp((const char*) map->buckets[idx]->key, (const char*) key))
  {
    /* Keys did not match, so we have a collision and need to probe */
    retval = SUNHashMap_Iterate(map, idx+1, sunHashMapLinearProbeGet, (void*) key);
    if (retval < 0)
      return(-1); /* error occurred */
    else if (retval == map->max_size)
      return(-2); /* not found */
  }

  /* Return a reference to the value only */
  *value = map->buckets[idx]->value;

  return(0);
}

/*
  This function allocates a new array the same max_size as the map,
  then it sorts map into a new array of key-value pairs leaving
  the map unchanged.

  **Arguments:**
    * ``map`` -- the ``SUNHashMap`` object to operate on
    * ``sorted`` -- pointer to the sorted array of key-value pairs, this
                    function will allocate the array
    * ``compar`` -- comparator function that is passed to the C standard qsort
                    function

  **Returns:**
    * ``0`` -- success
    * ``-1`` -- an error occurred
 */
static int SUNHashMap_Sort(SUNHashMap map, SUNHashMapKeyValue** sorted,
                           int (*compar)(const void *, const void*))
{
  int i;

  if (map == NULL || compar == NULL)
    return(-1);

  *sorted = (SUNHashMapKeyValue*) malloc(map->max_size*sizeof(SUNHashMapKeyValue));
  if (*sorted == NULL)
    return(-1);

  /* Copy the buckets into a new array */
  for (i = 0; i < map->max_size; i++)
    (*sorted)[i] = map->buckets[i];

  qsort(*sorted, map->max_size, sizeof(SUNHashMapKeyValue), compar);

  return(0);
}

/*
  This function allocates a new array with just they values of the map.

  **Arguments:**
    * ``map`` -- the ``SUNHashMap`` object to operate on
    * ``values`` -- pointer to the array of keys
    * ``value_size`` -- the size of the values in bytes

  **Returns:**
    * ``0`` -- success
    * ``-1`` -- an error occurred
 */
#if SUNDIALS_MPI_ENABLED
static int SUNHashMap_Values(SUNHashMap map, void*** values, size_t value_size)
{
  int i, count = 0;

  if (map == NULL)
    return(-1);

  *values = (void**) malloc(map->size*sizeof(value_size));
  if (values == NULL)
    return(-1);

  /* Copy the values into a new array */
  for (i = 0; i < map->max_size; i++)
  {
    if (map->buckets[i])
      (*values)[count++] = map->buckets[i]->value;
  }

  return(0);
}
#endif

#endif
