#include <tcl.h>

// formats.cpp
Tcl_CmdProc convertBinaryToText;
Tcl_CmdProc convertTextToBinary;
Tcl_CmdProc stripOpenSeesXML;
Tcl_CmdProc Tcl_PeridynamicsCommands;

struct char_cmd {
  const char* name; Tcl_CmdProc*  func;
}  const InterpreterCommands[] =  {

  {"peridynamics",         Tcl_PeridynamicsCommands},

  {"stripXML",             stripOpenSeesXML    },
  {"convertBinaryToText",  convertBinaryToText },
  {"convertTextToBinary",  convertTextToBinary },
};
