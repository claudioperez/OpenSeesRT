

#include <g3_api.h>


#include <SRC/element/dispBeamColumn/AxEqDispBeamColumn2d.h>
void *OPS_AxEqDispBeamColumn2d(void)
{
  // print out a message about who wrote this element & any copyright info
  // wanted
  if (numMyDBEle == 0) {
    opserr
        << "AxEqDispBeamColumn2d element - Written by Danilo Tarquini 2017 \n";
    numMyDBEle++;
  }

  // initializing a pointer to the element class
  Element *theEle = 0;

  // added by DANILO: initialization of the beam integration
  BeamIntegration *beamIntegr = 0;

  // creation of an empty element, required for parallel processing only
  int numRemainingArgs = OPS_GetNumRemainingInputArgs();
  if (numRemainingArgs == 0) { // parallel processing
    theEle = new AxEqDispBeamColumn2d();
    return theEle;
  }

  // check the required input parameters are given
  // for this element it is not possible to define different integration
  // sections along the element, to specify element masses, and to form
  // consistent mass matrix Differently from a classical DB element, 1) a
  // tolerance on the axial force unbalance $tol in the sections MUST be
  // specified; It represents the unbalance in axial force at each IP that is
  // deemed acceptable. It depends on the analysis performed and on the employed
  // units 2) a max num of internal element iterations can be specified. Default
  // value is : $maxIters=20
  if (numRemainingArgs < 7) {
    opserr << "insufficient arguments: 1)eleTag? 2)iNode? 3)jNode? "
              "4)numIntgrPts? 5)-$secTag? 6)$transfTag? 7)$tol optionals: "
              "<-integration $intType> <-iter $maxIters>";
    numMyDBEle++;
  }

  // get the input parameters
  int iData[6];
  int numData = 6;

  // check the that 6 integer are given as input data
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid element data\n";
    return 0;
  }

  // variables storing the main quantities required to create the element
  int eleTag = iData[0];    // element tag
  int iNode = iData[1];     // initial node
  int jNode = iData[2];     // final node
  int nIP = iData[3];       // number of integration points
  int secTag = iData[4];    // section tag
  int transfTag = iData[5]; // transformation tag

  // added by DANILO: reading of the specified axial force tolerance
  double tolerance;
  numData = 1;
  if (OPS_GetDoubleInput(&numData, &tolerance) != 0) {
    opserr << "WARNING error reading tolerance" << eleTag << endln;
    return 0;
  }

  // optional quantities that can be defined for the element
  double mass = 0.0;
  int cMass = 0;
  int numdata = 1;
  int maxNumIters = 20;

  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char *type = OPS_GetString();
    if (strcmp(type, "-cMass") == 0) {
      opserr << "WARNING: Consistent mass matrix not available for this "
                "element, Lumped mass matrix is used \n";
      cMass = 0; // changed by Danilo with respect to DB element
    } else if (strcmp(type, "-mass") == 0) {
      if (OPS_GetNumRemainingInputArgs() > 0) {
        if (OPS_GetDoubleInput(&numData, &mass) < 0) {
          opserr << "WARNING: invalid mass\n";
          return 0;
        }
        opserr << "WARNING: Element mass cannot be defined for this element\n";
        return 0;
      }
    }
    // added by DANILO: optional parameter that can be specified is the
    // integration type (LOBATTO or LEGENDRE)
    else if (strcmp(type, "-integration") == 0) {
      if (OPS_GetNumRemainingInputArgs() > 0) {
        const char *typeIntegration = OPS_GetString();
        if (strcmp(typeIntegration, "Lobatto") == 0)
          beamIntegr = new LobattoBeamIntegration();
        else if (strcmp(typeIntegration, "Legendre") == 0)
          beamIntegr = new LegendreBeamIntegration();
        else {
          opserr << "WARNING: invalid integration type\n";
          return 0;
        }
      }
    }
    // added by DANILO: optional parameter that can be specified is the number
    // of internal element iterations
    else if (strcmp(type, "-iter") == 0) {
      if (OPS_GetNumRemainingInputArgs() > 0) {
        if (OPS_GetIntInput(&numData, &maxNumIters) < 0 && maxNumIters <= 0) {
          opserr << "WARNING: invalid max number of iterations\n";
          return 0;
        }
      }
    }
  }

  // added by DANILO: if Beam Integration is not defined, use gauss legendere
  if (beamIntegr == 0) {
    beamIntegr = new LegendreBeamIntegration();
  }

  // get the beam transformation
  CrdTransf *theTransf = OPS_GetCrdTransf(transfTag); // get the transformation
  if (theTransf == 0) { // return error if the transformation is not found
    opserr << "coord transfomration not found\n";
    return 0;
  }

  // pointer to pointers of sections of type secTag
  SectionForceDeformation *theSection = OPS_getSectionForceDeformation(
      secTag); // get the section of type sectionTag
  SectionForceDeformation **sections = new SectionForceDeformation
      *[nIP]; // pointer to vector of pointers to allocate nIP sections
  // check that the section tag exists
  if (theSection == 0) {
    opserr << "WARNING section not found\n";
    opserr << "Section: " << secTag;
    opserr << " element: " << eleTag << endln;
    return 0;
  }
  // definition of the vector of section pointers
  for (int i = 0; i < nIP; i++)
    sections[i] = theSection;

  // now creation of the element
  theEle =
      new AxEqDispBeamColumn2d(eleTag, iNode, jNode, nIP, sections, *beamIntegr,
                               *theTransf, tolerance, mass, cMass, maxNumIters);

  delete[] sections;
  delete beamIntegr;

  return theEle;
}
