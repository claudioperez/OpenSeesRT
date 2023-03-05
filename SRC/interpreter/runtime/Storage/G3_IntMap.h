/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

// Simple hash table implemented in C.

#ifndef G3_INTMAP_H
#define G3_INTMAP_H

#include <stdbool.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

// Hash table structure: create with G3_NewIntMap, free with G3_DeleteIntMap.
typedef unsigned long int G3_MapTag;
#define G3_FMT_TAG "%lu"
typedef struct G3_IntMap G3_IntMap;

// Create hash table and return pointer to it, or NULL if out of memory.
G3_IntMap *G3_NewIntMap(void);

// Free memory allocated for hash table, including allocated keys.
void G3_DeleteIntMap(G3_IntMap *);

// Get item with given key (NUL-terminated) from hash table. Return
// value (which was set with G3_SetIntMapEntry), or NULL if key not found.
void *G3_GetIntMapEntry(G3_IntMap *, G3_MapTag);

// Set item with given key (NUL-terminated) to value (which must not
// be NULL). If not already present in table, key is copied to newly
// allocated memory (keys are freed automatically when G3_DeleteIntMap is
// called). Return address of copied key, or NULL if out of memory.
const char *G3_SetIntMapEntry(G3_IntMap *table, G3_MapTag tag, void *value);

// Return number of items in hash table.
size_t G3_GetIntMapLength(G3_IntMap *table);


#ifdef __cplusplus
}
#endif
#endif // G3_INTMAP_H
