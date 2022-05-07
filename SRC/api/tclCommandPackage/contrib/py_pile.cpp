// Added by Scott J. Brandenberg (sjbrandenberg@ucdavis.edu)
#include <PySimple1Gen.h>
#include <TzSimple1Gen.h>
// End added by SJB
////////////////////////////////////////////////////////////////////////////////////////////
// Added by Scott J. Brandenberg, UC Davis, sjbrandenberg@ucdavis.edu
int
TclCommand_doPySimple1Gen(ClientData clientData, Tcl_Interp *interp, int argc,
                               TCL_Char **argv)
{
   if(argc < 6 || argc > 7){
           opserr << "WARNING PySimple1Gen file1? file2? file3? file4? file5? <file6?>";
           opserr << "Must have either 5 or 6 arguments." << endln;
   }

   PySimple1Gen *thePySimple1Gen;
   thePySimple1Gen = new PySimple1Gen;

   if(argc==6)
     thePySimple1Gen->WritePySimple1(argv[1], argv[2], argv[3], argv[4], argv[5]);
   if(argc==7) 
     thePySimple1Gen->WritePySimple1(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6]);

   delete thePySimple1Gen;

   return TCL_OK;
}
int
TclCommand_doTzSimple1Gen(ClientData clientData, Tcl_Interp *interp, int argc,
                               TCL_Char **argv)
{
        if(argc < 6 || argc > 7){
                opserr << "WARNING TzSimple1Gen file1? file2? file3? file4? file5? <file6?>"; opserr << "Must have either 5 or 6 arguments." << endln;
        }

        TzSimple1Gen *theTzSimple1Gen;
        theTzSimple1Gen = new TzSimple1Gen;

        if(argc==6)
                theTzSimple1Gen->WriteTzSimple1(argv[1], argv[2], argv[3], argv[4], argv[5]); if(argc==7) theTzSimple1Gen->WriteTzSimple1(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6]);

        delete theTzSimple1Gen;

        return TCL_OK;
}
// End Added by Scott J. Brandenberg
///////////////////////////////////////////////////////////////////////////////////////////////////
