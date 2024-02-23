#ifndef opserr
#include <OPS_Stream.h>
extern OPS_Stream *opserrPtr;
#define opserr (*opserrPtr)
#define endln "\n"
#endif
