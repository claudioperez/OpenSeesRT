/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

// Simple hash table iterator implemented in C.

#ifndef G3_TABLEITERATOR_H
#define G3_TABLEITERATOR_H

#include <stdbool.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif


//
// Iterator
//
// Hash table iterator: create with ht_iterator, iterate with ht_next.
typedef struct {
    union {
      int tag;
      const char* key;
    }; // current key
    int _key_type;

    void* value;      // current value

    // Don't use these fields directly.
    void* _table;       // reference to hash table being iterated
    size_t _index;      // current index into ht._entries
} G3_TableIterator;

bool G3_NextTableEntry(G3_TableIterator* it);


#ifdef __cplusplus
}
#endif
#endif // G3_TABLEITERATOR_H
