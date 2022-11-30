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

#ifdef __cplusplus
extern "C" {
#endif


typedef unsigned long int G3_Tag;
typedef struct G3_Table G3_Table;

G3_Table* G3_NewTable(void);

const char *G3_SetTableEntry(G3_Table*, const char* partition, int tag, void* value);

int   G3_AddTableEntry(G3_Table*, const char* partition, int tag, void* value);

void *G3_GetTableEntry(G3_Table*, const char* partition, G3_Tag);

#ifdef __cplusplus
}
#endif

#endif // G3_TABLE_H
