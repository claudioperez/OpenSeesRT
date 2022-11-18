// Simple hash table implemented in C.

#ifndef G3_STRINGMAP_H
#define G3_STRINGMAP_H

#include <stdbool.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

// Hash table structure: create with ht_create, free with G3_DeleteStringMap.
typedef struct G3_StringMap G3_StringMap;

// Create hash table and return pointer to it, or NULL if out of memory.
G3_StringMap* G3_NewStringMap(void);
G3_StringMap* G3_InitStringMap(G3_StringMap*, size_t);

// Free memory allocated for hash table, including allocated keys.
void G3_DeleteStringMap(G3_StringMap*);

// Get item with given key (NUL-terminated) from hash table. Return
// value (which was set with G3_SetStringMapEntry), or NULL if key not found.
void* G3_GetStringMapEntry(G3_StringMap* table, const char* key);

// Set item with given key (NUL-terminated) to value (which must not
// be NULL). If not already present in table, key is copied to newly
// allocated memory (keys are freed automatically when G3_DeleteStringMap is
// called). Return address of copied key, or NULL if out of memory.
const char* G3_SetStringMapEntry(G3_StringMap* table, const char* key, void* value);

// Return number of items in hash table.
size_t G3_GetStringMapLength(G3_StringMap* table);

// Hash table iterator: create with ht_iterator, iterate with ht_next.
typedef struct {
    const char* key;  // current key
    void* value;      // current value

    // Don't use these fields directly.
    G3_StringMap* _table;       // reference to hash table being iterated
    size_t _index;    // current index into ht._entries
} hti;

// Return new hash table iterator (for use with ht_next).
hti ht_iterator(G3_StringMap* table);

// Move iterator to next item in hash table, update iterator's key
// and value to current item, and return true. If there are no more
// items, return false. Don't call G3_SetStringMapEntry during iteration.
bool ht_next(hti* it);


#ifdef __cplusplus
}
#endif
#endif // G3_STRINGMAP_H
