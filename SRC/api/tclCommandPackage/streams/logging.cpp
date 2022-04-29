#include <G3_Runtime.h>
#include <OPS_Globals.h>

// the following is a little kludgy but it works!
#ifdef _USING_STL_STREAMS
#  include <iomanip>
   using std::ios;
#  include <iostream>
   using std::ofstream;
#else
#  include <StandardStream.h>
#  include <FileStream.h>
#  include <DummyStream.h>
   StandardStream sserr;
   DummyStream    ssnul;
   OPS_Stream *opserrPtr = &sserr;
   OPS_Stream *opsdbgPtr = &ssnul;
   OPS_Stream *opswrnPtr = &sserr;
   OPS_Stream *opsmrdPtr = &sserr;

#endif

#include <G3_Logging.h>

const char * G3_WarnPromptColor   = RED "WARNING " COLOR_RESET;
const char * G3_WarnPromptNoColor = "WARNING ";

const char * G3_WARN_PROMPT = G3_WarnPromptNoColor;

int
G3_setStreamLevel(G3_Runtime* rt, int stream, int level)
{
  OPS_Stream **theStream;
  switch (stream) {
    case G3_Error: theStream = &opserrPtr; break;
    case G3_Debug: theStream = &opsdbgPtr; break;
    case G3_Warn : theStream = &opswrnPtr; break;
  }

  switch (level) {
    case G3_Null: *theStream = &ssnul;
    case G3_Log : *theStream = &sserr;
  }
  return 0;
}

int G3_setStreamColor(G3_Runtime* rt, int strm, int flag)
{
  if (flag == 1)
    G3_WARN_PROMPT = G3_WarnPromptColor;
  if (flag == 0)
    G3_WARN_PROMPT = G3_WarnPromptNoColor;

  return 0;
}

