/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

#include <stdlib.h>
#include "G3_Table.h"
#include "G3_IntMap.h"
#include "G3_StringMap.h"

#define INITIAL_CAPACITY 4

typedef unsigned long int G3_MapTag;

struct G3_Table {
  G3_StringMap *partitions;
};

G3_Table *
G3_NewTable(void)
{
  G3_Table *table = malloc(sizeof(G3_Table));
  if (table == NULL)
    return NULL;
  if ((table->partitions = G3_NewStringMap()) == NULL) {
    free(table);
    return NULL;
  }
  return table;
}

void
G3_DeleteTable(G3_Table *table)
{
  // TODO
}

#if 0
const char *
G3_SetTableEntry(G3_Table *table, const char *partition, int tag, void *value)
{
  // TODO
}
#endif


int
G3_AddTableEntry(G3_Table *table, const char *partition, int tag, void *value)
{
  G3_IntMap *imap =
      (G3_IntMap *)G3_GetStringMapEntry(table->partitions, partition);
  if ((imap == NULL) && !(imap = G3_NewIntMap())) {
    return 0;
  }
  if (!G3_SetStringMapEntry(table->partitions, partition, (void *)imap)) {
    G3_DeleteIntMap(imap);
    return 0;
  }

  G3_SetIntMapEntry(imap, tag, (void *)value);

  // TODO; send a meaningful return
  return 1;
}

void *
G3_GetTableEntry(G3_Table *table, const char *partition, G3_MapTag tag)
{
  G3_IntMap *imap;
  if ((imap = G3_GetStringMapEntry(table->partitions, partition)) == NULL) {
    return NULL;
  }

  void *value;
  if ((value = G3_GetIntMapEntry(imap, tag)) == NULL) {
    return NULL;
  }
  return value;
}

G3_TableIterator G3_IteratePartition(G3_Table* table, const char* partition)
{
  G3_IntMap *imap;
  G3_TableIterator iter;

  if ((imap = G3_GetStringMapEntry(table->partitions, partition)) == NULL)
    iter._table = NULL;
  else
    iter._table = imap;

//iter._key_type = 'i';
  iter._index = 0;
  
  return iter;
}

