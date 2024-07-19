#include <tcl.h>
#include <assert.h>
#include <Logging.h>
#include <Parsing.h>
#include <BasicModelBuilder.h>
#include <element/Shell/ASDShellQ4.h>
#include <element/Shell/ShellANDeS.h>
#include <element/Shell/ShellDKGQ.h>
#include <element/Shell/ShellDKGT.h>
#include <element/Shell/ShellMITC4.h>
#include <element/Shell/ShellMITC9.h>
#include <element/Shell/ShellNLDKGQ.h>
#include <element/Shell/ShellMITC4Thermal.h>
#include <element/Shell/ShellNLDKGQThermal.h>
#include <element/Shell/ShellNLDKGT.h>

Element*
TclDispatch_newASDShellQ4(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** const argv)
{
  assert(clientData != nullptr);

  bool corotational = false;


  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;

  if (argc < 6) {
    opserr << "Want: element ASDShellQ4 $tag $iNode $jNode $kNode $lNode "
              "$secTag <-corotational>";
    return nullptr;
  }

  int iData[6];
  int numData = 6;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer tag: element ASDShellQ4 \n";
    return nullptr;
  }

  if (argc == 7) {
    const char *type = OPS_GetString();
    if ((strcmp(type, "-corotational") == 0) ||
        (strcmp(type, "-Corotational") == 0))
      corotational = true;
  }

  SectionForceDeformation *section = builder->getTypedObject<SectionForceDeformation>(iData[5]);
  if (section == nullptr)
    return nullptr;

  return new ASDShellQ4(iData[0], iData[1], iData[2], iData[3], iData[4],
                        section, corotational);
}


Element*
TclDispatch_newShellANDeS(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** const argv)
{

  if (argc < 6) {
    opserr << "Want: element ShellANDeS $tag $iNode $jNode $kNode $thick $E "
              "$nu $rho";
    return nullptr;
  }

  int numArgs = OPS_GetNumRemainingInputArgs();

  int iData[4];
  int numData = 4;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer tag: element ShellANDeS \n";
    return nullptr;
  }

  double dData[11];
  numArgs = OPS_GetNumRemainingInputArgs();
  if (OPS_GetDoubleInput(&numArgs, dData) != 0) {
    opserr << "WARNING invalid double thickness: element ShellANDeS \n";
    return nullptr;
  }

  Element *theElement = nullptr;

  if (numArgs == 4) {
    theElement = new ShellANDeS(iData[0], iData[1], iData[2], iData[3],
                                dData[0], dData[1], dData[2], dData[3]);
  } else if (numArgs == 11) {
    theElement =
        new ShellANDeS(iData[0], iData[1], iData[2], iData[3], dData[0],
                       dData[1], dData[2], dData[3], dData[4], dData[5],
                       dData[6], dData[7], dData[8], dData[9], dData[10]);
  }

  return theElement;
}

Element*
TclDispatch_newShellDKGQ(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;

  if (argc < 6) {
    opserr << "Want: element ShellDKGQ $tag $iNode $jNoe $kNode $lNode $secTag";
    return nullptr;
  }

  int iData[6];
  int numData = 6;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer tag: element ShellDKGQ \n";
    return nullptr;
  }

  SectionForceDeformation *theSection = builder->getTypedObject<SectionForceDeformation>(iData[5]);
  if (theSection == nullptr)
    return nullptr;

  return new ShellDKGQ(iData[0], iData[1], iData[2], iData[3], iData[4],
                             *theSection);
}

Element*
TclDispatch_newShellDKGT(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;

  int numArgs = OPS_GetNumRemainingInputArgs();

  if (numArgs < 5) {
    opserr << "Want: element ShellDKGT $tag $iNode $jNoe $kNode $secTag";
    return nullptr;
  }

  int iData[5];
  int numData = 5;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer tag: element ShellDKGT \n";
    return nullptr;
  }

  SectionForceDeformation *theSection = builder->getTypedObject<SectionForceDeformation>(iData[4]);
  if (theSection == nullptr)
    return nullptr;


  double b_data[3] = {0, 0, 0};

  int num_remaining_args = OPS_GetNumRemainingInputArgs();

  if (num_remaining_args > 3) {
    num_remaining_args = 3;
  }
  if (num_remaining_args > 0) {
    if (OPS_GetDoubleInput(&num_remaining_args, b_data) < 0) {
      opserr << "WARNING: invalid double b_data\n";
      return nullptr;
    }
  }

  return new ShellDKGT(iData[0], iData[1], iData[2], iData[3],
                       *theSection, b_data[0], b_data[1], b_data[2]);
}



Element*
TclDispatch_newShellMITC4(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;

  bool updateBasis = false;
  Element *theElement = nullptr;

  if (argc < 6) {
    opserr << "Want: element ShellMITC4 $tag $iNode $jNode $kNode $lNode "
              "$secTag<-updateBasis>";
    return nullptr;
  }

  int iData[6];
  int numData = 6;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer tag: element ShellMITC4 \n";
    return nullptr;
  }

  if (argc == 7) {
    const char *type = OPS_GetString();
    if (strcmp(type, "-updateBasis") == 0)
      updateBasis = true;
  }

  SectionForceDeformation *theSection = builder->getTypedObject<SectionForceDeformation>(iData[5]);
  if (theSection == nullptr)
    return nullptr;

  theElement = new ShellMITC4(iData[0], iData[1], iData[2], iData[3], iData[4],
                              *theSection, updateBasis);

  return theElement;
}



Element*
TclDispatch_newShellMITC4Thermal(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;

  int numArgs = OPS_GetNumRemainingInputArgs();

  if (numArgs < 6) {
    opserr << "Want: element ShellMITC4Thermal $tag $iNode $jNoe $kNode $lNode "
              "$secTag";
    return nullptr;
  }

  int iData[6];
  int numData = 6;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer tag: element ShellMITC4Thermal \n";
    return nullptr;
  }

  SectionForceDeformation *theSection = builder->getTypedObject<SectionForceDeformation>(iData[5]);
  if (theSection == nullptr)
    return nullptr;

  return new ShellMITC4Thermal(iData[0], iData[1], iData[2], iData[3], iData[4], *theSection);
}


Element*
TclDispatch_newShellMITC9(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;

  int numArgs = OPS_GetNumRemainingInputArgs();

  if (numArgs < 11) {
    opserr << "Want: element ShellMITC9 $tag $node1 $node2 .... $node9 $secTag";
    return nullptr;
  }

  int iData[11];
  int numData = 11;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer tag: element ShellMITC9\n";
    return nullptr;
  }

  SectionForceDeformation *theSection = builder->getTypedObject<SectionForceDeformation>(iData[10]);
  if (theSection == nullptr)
    return nullptr;

  return
      new ShellMITC9(iData[0], iData[1], iData[2], iData[3], iData[4], iData[5],
                     iData[6], iData[7], iData[8], iData[9], *theSection);
}



Element*
TclDispatch_newShellNLDKGQ(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;

  int numArgs = OPS_GetNumRemainingInputArgs();

  if (numArgs < 6) {
    opserr
        << "Want: element ShellNLDKGQ $tag $iNode $jNoe $kNode $lNode $secTag";
    return nullptr;
  }

  int iData[6];
  int numData = 6;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer tag: element ShellNLDKGQ \n";
    return nullptr;
  }

  SectionForceDeformation *theSection = builder->getTypedObject<SectionForceDeformation>(iData[5]);
  if (theSection == nullptr)
    return nullptr;

  return new ShellNLDKGQ(iData[0], iData[1], iData[2], iData[3], iData[4],
                               *theSection);

}


Element*
TclDispatch_newShellNLDKGQThermal(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;

  int numArgs = OPS_GetNumRemainingInputArgs();

  if (numArgs < 6) {
    opserr << "Want: element ShellNLDKGQThermal $tag $iNode $jNoe $kNode "
              "$lNode $secTag";
    return nullptr;
  }

  int iData[6];
  int numData = 6;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer tag: element ShellNLDKGQThermal \n";
    return nullptr;
  }

  SectionForceDeformation *theSection = builder->getTypedObject<SectionForceDeformation>(iData[5]);
  if (theSection == nullptr)
    return nullptr;

  return new ShellNLDKGQThermal(iData[0], iData[1], iData[2], iData[3],
                                      iData[4], *theSection);

}


Element*
TclDispatch_newShellNLDKGT(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;

  int numArgs = OPS_GetNumRemainingInputArgs();

  if (numArgs < 5) {
    opserr << "Want: element ShellNLDKGT $tag $iNode $jNoe $kNode $secTag";
    return nullptr;
  }

  int iData[5];
  int numData = 5;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer tag: element ShellNLDKGT \n";
    return nullptr;
  }

  SectionForceDeformation *theSection = builder->getTypedObject<SectionForceDeformation>(iData[4]);
  if (theSection == nullptr)
    return nullptr;

  return
      new ShellNLDKGT(iData[0], iData[1], iData[2], iData[3], *theSection);

}
