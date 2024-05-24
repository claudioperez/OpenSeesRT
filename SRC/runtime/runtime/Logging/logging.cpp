#include <cstdarg>
#include <OPS_Globals.h>
#include <tcl.h>
#include <elementAPI.h>
class G3_Runtime;


#include <StandardStream.h>
#include <FileStream.h>
#include <DummyStream.h>
StandardStream sserr;
DummyStream    ssnul;
OPS_Stream *opserrPtr = &sserr;
OPS_Stream *opsdbgPtr = &ssnul;
OPS_Stream *opswrnPtr = &sserr;
OPS_Stream *opsmrdPtr = &sserr;


#include <G3_Logging.h>

const char * G3_WarnPromptColor   = RED "WARNING " COLOR_RESET;
const char * G3_WarnPromptNoColor = "WARNING ";

const char * G3_ErrorPromptColor   = BRED "ERROR " COLOR_RESET;
const char * G3_ErrorPromptNoColor = "ERROR ";

const char * G3_DebugPromptColor   = GRN "DEBUG " COLOR_RESET;
const char * G3_DebugPromptNoColor = "DEBUG ";

const char * G3_WARN_PROMPT  = G3_WarnPromptNoColor;
const char * G3_ERROR_PROMPT = G3_ErrorPromptNoColor;
const char * G3_DEBUG_PROMPT = G3_DebugPromptNoColor;

int
G3_SetStreamLevel(int stream, bool on)
{
  OPS_Stream **theStream;
  switch (stream) {
    case G3_LevelError: theStream = &opserrPtr; break;
    case G3_LevelDebug: theStream = &opsdbgPtr; break;
    case G3_LevelWarn : theStream = &opswrnPtr; break;
    default:
      return -1;
  }

  if (on) {
    *theStream = &sserr;
  } else {
    *theStream = &ssnul;
  }
  return 0;
}

int G3_SetStreamColor(G3_Runtime* rt, int strm, int flag)
{
  if (flag == 1) {
    G3_WARN_PROMPT = G3_WarnPromptColor;
    G3_ERROR_PROMPT = G3_ErrorPromptColor;
    G3_DEBUG_PROMPT = G3_DebugPromptColor;

  } else if (flag == 0) {
    G3_WARN_PROMPT = G3_WarnPromptNoColor;
    G3_ERROR_PROMPT = G3_ErrorPromptNoColor;
    G3_DEBUG_PROMPT = G3_DebugPromptNoColor;
  }

  return 0;
}



int
G3_Raise(G3_Runtime *rt, const char *msg, ...)
{
  va_list ap;

  va_start(ap, msg);
  int n = vsnprintf(NULL, 0, msg, ap);
  va_end(ap);

  if (n < 0)
    return -1;

  size_t size = (size_t)n + 1 + 8;
  char *new_str = (char*)malloc(size);
  if (new_str == NULL)
    return -1;

  strcpy(new_str, "error {");
  va_start(ap, msg);
  n = vsnprintf(new_str+7, size, msg, ap);
  va_end(ap);
  strcpy(new_str+7+n, "}\n");


  Tcl_Interp *tcl_interp = G3_getInterpreter(rt);
  Tcl_Eval(tcl_interp, new_str);
  Tcl_Obj *infoObj = Tcl_GetVar2Ex(tcl_interp, "errorInfo", NULL, TCL_GLOBAL_ONLY);
  const char * error_str = Tcl_GetString(infoObj);
  opserr << error_str;

  /*
  Tcl_Obj *top_interpInfoName ;
  Tcl_Obj *top_interpInfo ;
    top_interpInfoName = Tcl_NewStringObj("errorInfo", -1) ;
    Tcl_IncrRefCount(top_interpInfoName) ;
    top_interpInfo =  Tcl_ObjGetVar2(tcl_interp,
                                     top_interpInfoName,
                                     NULL,
                                     TCL_LEAVE_ERR_MSG) ;
    Tcl_IncrRefCount(top_interpInfo) ;
    const char *error_str = Tcl_GetString(top_interpInfo);
    opserr << "ERROR -- " << msg << "\n\n" << error_str;
    Tcl_DecrRefCount(top_interpInfoName) ;
    Tcl_DecrRefCount(top_interpInfo);
    */

    return TCL_ERROR;
}
