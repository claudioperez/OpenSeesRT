

#include <g3_api.h>


#include <SRC/element/pyMacro/PY_Macro2D.h>
OPS_Export void *OPS_PY_Macro2D(void)
{
  if (numPY_Macro2D == 0) {
    opserr << "PY_Macro2D element - Written by V.Varun and A.Shafiee, Georgia "
              "Tech Copyright 2009\n";
    numPY_Macro2D++;
  }
  // get the id and end nodes
  int iData[4]; //!!!!!!!!!!
  double dData[14];
  int numData;

  numData = 3;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid element data for PY_Macro2D\n";
    return 0;
  }

  int eleTag = iData[0];

  numData = 14;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING error reading element data for PY_Macro2D element with "
              "tag: "
           << eleTag << endln;
    return 0;
  }

  // now create the truss and add it to the Domain
  Element *theEle =
      new PY_Macro2D(eleTag, iData[1], iData[2], dData[0], dData[1], dData[2],
                     dData[3], dData[4], dData[5], dData[6], dData[7], dData[8],
                     dData[9], dData[10], dData[11], dData[12], iData[13]);

  if (theEle == 0) {
    opserr << "WARNING ran out of memory creating PY_Macro2D element with tag "
           << eleTag << endln;
    return 0;
  }

  return theEle;
}
