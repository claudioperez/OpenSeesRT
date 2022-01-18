

#include <g3_api.h>


#include <SRC/element/dispBeamColumnInt/DispBeamColumn2dInt.h>
void *OPS_DispBeamColumn2dInt()
{
  int ndm = OPS_GetNDM();
  int ndf = OPS_GetNDF();

  int ok = 0;
  if (ndm == 2 && ndf == 3)
    ok = 1;

  if (ok == 0) {
    opserr << "WARNING -- NDM = " << ndm << " and NDF = " << ndf
           << " not compatible with dispBeamColumn element" << endln;
    return 0;
  }

  if (OPS_GetNumRemainingInputArgs() < 7) { // 8
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: element dispBeamColumn eleTag? iNode? jNode? nIP? secTag? "
              "transfTag? C1? t1? NStrip1? t2? NStrip2? t3? NStrip3?\n";
    return 0;
  }

  // get the id and end nodes
  int eleTag, iNode, jNode, nIP, transfTag;
  double C1;
  int secTag[10]; // Max size of integration rule ... can change if needed

  int idata[4];
  int numdata = 4;
  if (OPS_GetIntInput(&numdata, idata) < 0) {
    opserr << "WARNING invalid dispBeamColumn int inputs" << endln;
    return 0;
  }
  eleTag = idata[0];
  iNode = idata[1];
  jNode = idata[2];
  nIP = idata[3];

  const char *type = OPS_GetString();
  if (strcmp(type, "-sections") == 0) {
    if (nIP > OPS_GetNumRemainingInputArgs()) {
      opserr
          << "WARNING insufficient number of section tags - element "
             "dispBeamColumn eleTag? iNode? jNode? nIP? secTag? transfTag?\n";
      return 0;
    }
    int section;
    numdata = 1;
    for (int i = 0; i < nIP; i++) {
      if (OPS_GetIntInput(&numdata, &section) < 0) {
        opserr << "WARNING invalid secTag - element dispBeamColumn eleTag? "
                  "iNode? jNode? nIP? secTag? transfTag?\n";
        return 0;
      }
      secTag[i] = section;
    }
  }

  else {
    OPS_ResetCurrentInputArg(-1);
    int section;
    numdata = 1;
    if (OPS_GetIntInput(&numdata, &section) < 0) {
      opserr << "WARNING invalid secTag - element dispBeamColumn eleTag? "
                "iNode? jNode? nIP? secTag? transfTag?\n";
      return 0;
    }
    for (int i = 0; i < nIP; i++)
      secTag[i] = section;
  }

  if (OPS_GetNumRemainingInputArgs() > 0) {
    numdata = 1;
    if (OPS_GetIntInput(&numdata, &transfTag) < 0) {
      opserr << "WARNING invalid transfTag? - element dispBeamColumn eleTag? "
                "iNode? jNode? nIP? secTag? transfTag?\n";
      return 0;
    }
  }

  numdata = 1;
  if (OPS_GetDoubleInput(&numdata, &C1) < 0) {
    opserr << "WARNING invalid dispBeamColumn C1" << endln;
    return 0;
  }

  double massDens = 0.0;

  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char *massarg = OPS_GetString();
    if (strcmp(massarg, "-mass") == 0 && OPS_GetNumRemainingInputArgs() > 0) {
      numdata = 1;
      if (OPS_GetDoubleInput(&numdata, &massDens) < 0) {
        opserr << "WARNING invalid massDens - element dispBeamColumn eleTag? "
                  "iNode? jNode? nIP? secTag? transfTag? C1? t? NStrip?\n";
        return 0;
      }
    }
  }

  SectionForceDeformation **sections = new SectionForceDeformation *[nIP];

  if (!sections) {
    opserr << "WARNING TclElmtBuilder - addFrameElement - Insufficient memory "
              "to create sections\n";
    return 0;
  }

  for (int j = 0; j < nIP; j++) {
    SectionForceDeformation *theSection =
        OPS_getSectionForceDeformation(secTag[j]);

    if (theSection == 0) {
      opserr << "WARNING TclElmtBuilder - frameElement - no Section found with "
                "tag ";
      opserr << secTag[j] << endln;
      delete[] sections;
      return 0;
    }

    sections[j] = theSection;
  }

  Element *theElement = 0;

  if (ndm == 2) {

    CrdTransf *theTransf = OPS_getCrdTransf(transfTag);

    if (theTransf == 0) {
      opserr << "WARNING transformation not found\n";
      opserr << "transformation: " << transfTag;
      opserr << "\ndispBeamColumn element: " << eleTag << endln;
      return 0;
    }

    // now create the DispBeamColumn and add it to the Domain
    theElement = new DispBeamColumn2dInt(eleTag, iNode, jNode, nIP, sections,
                                         *theTransf, C1, massDens);

    delete[] sections;
  }

  // if get here we have successfully created the element and added it to the
  // domain
  return theElement;
}
