/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

// G3_Table is a nested data structure which allows integer-tagged void
// pointers to be organized into string-labeled bins, or "namespaces".

/* Claudio Perez */

#ifndef G3_TABLE_H
#define G3_TABLE_H

#include <stddef.h>
#include <stdbool.h>
#include "G3_TableIterator.h"

#ifdef __cplusplus
extern "C" {
#endif

// typedef unsigned long int G3_MapTag;
typedef struct G3_Table G3_Table;

G3_Table *G3_NewTable(void);

#if 0
const char *G3_SetTableEntry(G3_Table *, const char *partition, int tag,
                             void *value);
#endif
int G3_AddTableEntry(G3_Table *, const char *partition, int tag, void *value);

void *G3_GetTableEntry(G3_Table *, const char *partition, unsigned long int tag);

//
// Iterator
//

// Return new hash table iterator (for use with ht_next).
G3_TableIterator G3_IteratePartition(G3_Table* table, const char* partition);

#ifdef __cplusplus
}
#endif

#endif // G3_TABLE_H
