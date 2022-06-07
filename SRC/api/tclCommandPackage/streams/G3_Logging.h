#include "color.h"
#define opsdbg (*opsdbgPtr)
#define opswrn ((*opswrnPtr) << G3_WARN_PROMPT)
#define opsmrd opserr

extern OPS_Stream *opswrnPtr;
extern OPS_Stream *opsdbgPtr;
extern const char *G3_WARN_PROMPT;

enum G3_Stream {
  G3_StdOut, G3_StdIn, G3_StdErr, G3_Null
};
enum G3_StreamLevel {
  G3_Error, G3_Debug, G3_Log, G3_Warn
};

int G3_setStreamColor(G3_Runtime* rt, int strm, int flag);
