#include <tcl.h>

// formats.cpp
Tcl_CmdProc convertBinaryToText;
Tcl_CmdProc convertTextToBinary;
Tcl_CmdProc stripOpenSeesXML;
Tcl_CmdProc Tcl_Peri;

struct char_cmd {
  const char* name; Tcl_CmdProc*  func;
}  const InterpreterCommands[] =  {

  {"peri",         Tcl_Peri},

  {"stripXML",             stripOpenSeesXML    },
  {"convertBinaryToText",  convertBinaryToText },
  {"convertTextToBinary",  convertTextToBinary },
};
