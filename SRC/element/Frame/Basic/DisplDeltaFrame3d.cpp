//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
//
#include <array>
#include <math.h>
#include <string.h>

#include <DisplDeltaFrame3d.h>

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

#include <Node.h>
#include <FrameSection.h>
#include <FrameTransform.h>
#include <Domain.h>
#include <Information.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ElementResponse.h>
#include <CompositeResponse.h>
#include <ElementalLoad.h>
#include <BeamIntegration.h>
#include <Parameter.h>

#define ELE_TAG_DisplDeltaFrame3d 0 // TODO

DisplDeltaFrame3d::DisplDeltaFrame3d(int tag, std::array<int,2>& nodes, int numSec,
                                     FrameSection **s,
                                     BeamIntegration &bi, FrameTransform3d &coordTransf,
                                     double r, int cm)

    : FiniteElement(tag, ELE_TAG_DisplDeltaFrame3d, nodes),
      numSections(numSec), 
      theSections(nullptr),
      theCoordTransf(nullptr), beamInt(nullptr), 
      rho(r), cMass(cm), parameterID(0)

{
  // Allocate arrays of pointers to FrameSections
  theSections = new FrameSection *[numSections];

  for (int i = 0; i < numSections; i++) {
    // Get copies of the material model for each integration point
    theSections[i] = s[i]->getFrameCopy();
  }

  beamInt = bi.getCopy();

  theCoordTransf = coordTransf.getCopy();

  q0.zero();
}

DisplDeltaFrame3d::DisplDeltaFrame3d()
    : FiniteElement(0, ELE_TAG_DisplDeltaFrame3d),
      numSections(0), theSections(0),
      beamInt(0),
      rho(0.0), cMass(0), parameterID(0)
{
  q0.zero();
}

DisplDeltaFrame3d::~DisplDeltaFrame3d()
{
  for (int i = 0; i < numSections; i++) {
    if (theSections[i])
      delete theSections[i];
  }

  // Delete the array of pointers to FrameSection pointer arrays
  if (theSections)
    delete[] theSections;

  if (theCoordTransf)
    delete theCoordTransf;

  if (beamInt != nullptr)
    delete beamInt;
}


int
DisplDeltaFrame3d::commitState()
{
  int retVal = this->Element::commitState();

  // Loop over the integration points and commit the material states
  for (int i = 0; i < numSections; i++)
    retVal += theSections[i]->commitState();

  return retVal;
}

int
DisplDeltaFrame3d::revertToLastCommit()
{
  int retVal = 0;

  // Loop over the integration points and revert to last committed state
  for (int i = 0; i < numSections; i++)
    retVal += theSections[i]->revertToLastCommit();

  return retVal;
}

int
DisplDeltaFrame3d::revertToStart()
{
  int retVal = 0;

  // Loop over the integration points and revert states to start
  for (int i = 0; i < numSections; i++)
    retVal += theSections[i]->revertToStart();

  return retVal;
}

int
DisplDeltaFrame3d::update()
{
  int err = 0;

  // Update the transformation
  theCoordTransf->update();

  // Get basic deformations
  const Vector &v = theCoordTransf->getBasicTrialDisp();

  double L        = theCoordTransf->getInitialLength();
  double oneOverL = 1.0 / L;

  //const Matrix &pts = quadRule.getIntegrPointCoords(numSections);
  double xi[maxNumSections];
  beamInt->getSectionLocations(numSections, L, xi);

  // Loop over the integration points
  for (int i = 0; i < numSections; i++) {

    double xi1   = xi[i];
    double dNv1  = 1.0 + 3.0*xi1*xi1 - 4.0*xi1;
    double dNv2  = 3.0 * xi1 * xi1 - 2.0 * xi1;
    double ddNv1 = 6.0 * xi1 * oneOverL - 4.0 * oneOverL;
    double ddNv2 = 6.0 * xi1 * oneOverL - 2.0 * oneOverL;
    double dNw1  = -dNv1;
    double dNw2  = -dNv2;
    double ddNw1 = -ddNv1;
    double ddNw2 = -ddNv2;
    double Nf1   =  xi1;

    double dx   = oneOverL*v(0);            // u'
    double dy   =  dNw1*v(3) +  dNw2*v(4);  // y'
    double dz   =  dNv1*v(1) +  dNv2*v(2);  // z'
    double phi  = Nf1*v(5);                 // phi
                                            //
    double ddy  = ddNw1*v(3) + ddNw2*v(4);  // y"
    double ddz  = ddNv1*v(1) + ddNv2*v(2);  // z"
    double dphi = oneOverL*v(5);            // phi'
    double s7   = v(1);                     // theta_Iz
    double s8   = v(3);                     // theta_Iy
    double s9   = v(2);                     // theta_Jz
    double ddz0 = v(4);                     // theta_Jy


    VectorND<4> e {
            dx,
            ddz + ddy*phi,
           -ddy + ddz*phi,
            dphi
    };

    // Update the section state
    err += theSections[i]->setTrialState<nsr,scheme>(e);
  }

  if (err != 0) {
    opserr << "DisplDeltaFrame3d::update() - failed setTrialSectionDeformations()\n";
    return err;
  }
  return 0;
}

int
DisplDeltaFrame3d::setNodes()
{
  double L  = theCoordTransf->getInitialLength();
  double xi[maxNumSections];
  beamInt->getSectionLocations(numSections, L, xi);
  double wt[maxNumSections];
  beamInt->getSectionWeights(numSections, L, wt);
}

const Matrix &
DisplDeltaFrame3d::getTangentStiff()
{
  static Matrix kb(ndf*nen,  ndf*nen);
  static MatrixND<  8,nsr> A;
  static MatrixND<  8, 12> B;
  static MatrixND<  8,  8> ks;
  static MatrixND< 12, 12> Gm;

  const Vector &v = theCoordTransf->getBasicTrialDisp();

  double L        = theCoordTransf->getInitialLength();
  double oneOverL = 1.0 / L;

  // Zero for integral
  kb.Zero();
  p.zero();

  A.zero();
  B.zero();
  Gm.zero();

  // Loop over the integration points
  for (int i = 0; i < numSections; i++) {

    double xi1   = xi[i];
    double dNv1  = 1.0 + 3.0 * xi1 * xi1 - 4.0 * xi1;
    double ddNv1 = 6.0 * xi1 * oneOverL - 4.0 * oneOverL;
    double dNv2  = 3.0 * xi1 * xi1 - 2.0 * xi1;
    double ddNv2 = 6.0 * xi1 * oneOverL - 2.0 * oneOverL;
    double dNw1  = -dNv1;
    double ddNw1 = -ddNv1;
    double dNw2  = -dNv2;
    double ddNw2 = -ddNv2;
    double Nf1   = xi1;

    double dv  =  dNv1*v(1) +  dNv2*v(2); // v'
    double ddv = ddNv1*v(1) + ddNv2*v(2); // v"
    double dw  =  dNw1*v(3) +  dNw2*v(4); // w'
    double ddw = ddNw1*v(3) + ddNw2*v(4); // w"
    double f   = Nf1*v(5);                // phi
    double df  = oneOverL*v(5);           // phi'

    A( 0, 0)  = 1.0;
    A( 1, 0)  = (4.0 * v(1) - v(2)) / 30.0;
    A( 2, 0)  = (4.0 * v(3) - v(4)) / 30.0;
    A( 3, 0)  = (4.0 * v(2) - v(1)) / 30.0;
    A( 4, 0)  = (4.0 * v(4) - v(3)) / 30.0;
    A( 7, 1)  = 1.0;

    B( 0, 0)  = oneOverL;
    B( 1, 1)  = 1.0;
    B( 2, 3)  = 1.0;
    B( 3, 2)  = 1.0;
    B( 4, 4)  = 1.0;
    B( 5, 1)  = dNv1;
    B( 5, 2)  = dNv2;
    B( 6, 3)  = dNw1;
    B( 6, 4)  = dNw2;
    B( 7, 1)  = ddNv1;
    B( 7, 2)  = ddNv2;
    B( 8, 3)  = ddNw1;
    B( 8, 4)  = ddNw2;
    B( 9, 5)  = Nf1;
    B(10, 5) = oneOverL;

    // Get the section tangent stiffness and stress resultant
    MatrixND<4,4> Ks = theSections[i]->getTangent<nsr,scheme>(State::Pres);
    VectorND<4>   s  = theSections[i]->getResultant<nsr,scheme>();

    ks.addMatrixTripleProduct(0.0, A, Ks, 1.0);

    // Add material stiffness matrix
    kb.addMatrixTripleProduct(1.0, B, ks, L*wt[i]);

    // Beam geometric stiffness matrix
    Gm( 1,  1) = Gm(2,  2) = Gm(3, 3) = Gm(4, 4) =  s[0] * 4.0 / 30.0; // 4/30*N
    Gm( 1,  3) = Gm(2,  4) = Gm(3, 1) = Gm(4, 2) = -s[0] / 30.0;      // -1/30*N
    Gm( 9,  8) = Gm(8,  9) =  s[1];                                        //Mz
    Gm( 9,  7) = Gm(7,  9) =  s[2];                                        //My
    Gm(10, 10)             =  s[3];                                      //W

    // Add geometric stiffness
    kb.addMatrixTripleProduct(1.0, B, Gm, L*wt[i]);


    // assemble internal force vector p
    p.addMatrixTransposeVector(1.0, B, A*s, L*wt[i]);
  }

//q[0] += q0[0];
//q[1] += q0[1];
//q[2] += q0[2];
//q[3] += q0[3];
//q[4] += q0[4];

  return kb;
}


void
DisplDeltaFrame3d::zeroLoad()
{
  q0.zero();
  return;
}

int
DisplDeltaFrame3d::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  opserr
      << "DisplDeltaFrame3d::addLoad() -- load type unknown for element with tag: "
      << this->getTag() << "\n";
  return -1;
}



void
DisplDeltaFrame3d::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_CURRENTSTATE) {
    s << "\nDisplDeltaFrame3d, element id:  " << this->getTag() << "\n";
    s << "\tConnected external nodes:  " << connectedExternalNodes;
    s << "\tCoordTransf: " << theCoordTransf->getTag() << "\n";
    s << "\tmass density:  " << rho << ", cMass: " << cMass << "\n";

    beamInt->Print(s, flag);

    for (int i = 0; i < numSections; i++) {
      theSections[i]->Print(s, flag);
    }
  }

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{";
    s << "\"name\": " << this->getTag() << ", ";
    s << "\"type\": \"DisplDeltaFrame3d\", ";
    s << "\"nodes\": [" << connectedExternalNodes(0) << ", " << connectedExternalNodes(1)
      << "], ";
    s << "\"sections\": [";
    for (int i = 0; i < numSections - 1; i++)
      s << "\"" << theSections[i]->getTag() << "\", ";
    s << "\"" << theSections[numSections - 1]->getTag() << "\"], ";
    s << "\"integration\": ";
    beamInt->Print(s, flag);
    s << ", \"massperlength\": " << rho << ", ";
    s << "\"crdTransformation\": \"" << theCoordTransf->getTag() << "\"}";
  }
}


Response *
DisplDeltaFrame3d::setResponse(const char **argv, int argc, OPS_Stream &output)
{

  Response *theResponse = 0;

  output.tag("ElementOutput");
  output.attr("eleType", "DisplDeltaFrame3d");
  output.attr("eleTag", this->getTag());
  output.attr("node1", connectedExternalNodes[0]);
  output.attr("node2", connectedExternalNodes[1]);

  //
  // we compare argv[0] for known response types
  //

  // global force -
  if (strcmp(argv[0], "forces") == 0 || strcmp(argv[0], "force") == 0 ||
      strcmp(argv[0], "globalForce") == 0 || strcmp(argv[0], "globalForces") == 0) {

    output.tag("ResponseType", "Px_1");
    output.tag("ResponseType", "Py_1");
    output.tag("ResponseType", "Pz_1");
    output.tag("ResponseType", "Mx_1");
    output.tag("ResponseType", "My_1");
    output.tag("ResponseType", "Mz_1");
    output.tag("ResponseType", "Px_2");
    output.tag("ResponseType", "Py_2");
    output.tag("ResponseType", "Pz_2");
    output.tag("ResponseType", "Mx_2");
    output.tag("ResponseType", "My_2");
    output.tag("ResponseType", "Mz_2");

// TODO(cmp)
//  theResponse = new ElementResponse(this, 1, P);

    // local force -
  } else if (strcmp(argv[0], "localForce") == 0 || strcmp(argv[0], "localForces") == 0) {

    output.tag("ResponseType", "N_1");
    output.tag("ResponseType", "Vy_1");
    output.tag("ResponseType", "Vz_1");
    output.tag("ResponseType", "T_1");
    output.tag("ResponseType", "My_1");
    output.tag("ResponseType", "Mz_1");
    output.tag("ResponseType", "N_2");
    output.tag("ResponseType", "Vy_2");
    output.tag("ResponseType", "Vz_2");
    output.tag("ResponseType", "T_2");
    output.tag("ResponseType", "My_2");
    output.tag("ResponseType", "Mz_2");

// TODO(cmp)
//  theResponse = new ElementResponse(this, 2, P);

    // chord rotation -
  } else if (strcmp(argv[0], "chordRotation") == 0 ||
             strcmp(argv[0], "chordDeformation") == 0 ||
             strcmp(argv[0], "basicDeformation") == 0) {

    output.tag("ResponseType", "eps");
    output.tag("ResponseType", "thetaZ_1");
    output.tag("ResponseType", "thetaZ_2");
    output.tag("ResponseType", "thetaY_1");
    output.tag("ResponseType", "thetaY_2");
    output.tag("ResponseType", "thetaX");

    theResponse = new ElementResponse(this, 3, Vector(6));

    // plastic rotation -
  } else if (strcmp(argv[0], "plasticRotation") == 0 ||
             strcmp(argv[0], "plasticDeformation") == 0) {

    output.tag("ResponseType", "epsP");
    output.tag("ResponseType", "thetaZP_1");
    output.tag("ResponseType", "thetaZP_2");
    output.tag("ResponseType", "thetaYP_1");
    output.tag("ResponseType", "thetaYP_2");
    output.tag("ResponseType", "thetaXP");

    theResponse = new ElementResponse(this, 4, Vector(6));

  } else if (strcmp(argv[0], "RayleighForces") == 0 ||
             strcmp(argv[0], "rayleighForces") == 0) {

    theResponse = new ElementResponse(this, 12, P);

  } else if (strcmp(argv[0], "integrationPoints") == 0)
    theResponse = new ElementResponse(this, 10, Vector(numSections));

  else if (strcmp(argv[0], "integrationWeights") == 0)
    theResponse = new ElementResponse(this, 11, Vector(numSections));

  else if (strcmp(argv[0], "sectionTags") == 0)
    theResponse = new ElementResponse(this, 110, ID(numSections));

  // section response -
  else if (strstr(argv[0], "sectionX") != 0) {
    if (argc > 2) {
      float sectionLoc = atof(argv[1]);

      double xi[maxNumSections];
      double L = theCoordTransf->getInitialLength();
      beamInt->getSectionLocations(numSections, L, xi);

      sectionLoc /= L;

      float minDistance = fabs(xi[0] - sectionLoc);
      int sectionNum    = 0;
      for (int i = 1; i < numSections; i++) {
        if (fabs(xi[i] - sectionLoc) < minDistance) {
          minDistance = fabs(xi[i] - sectionLoc);
          sectionNum  = i;
        }
      }

      output.tag("GaussPointOutput");
      output.attr("number", sectionNum + 1);
      output.attr("eta", xi[sectionNum] * L);

      theResponse = theSections[sectionNum]->setResponse(&argv[2], argc - 2, output);
    }
  }

  else if (strcmp(argv[0], "section") == 0) {
    if (argc > 1) {

      int sectionNum = atoi(argv[1]);

      if (sectionNum > 0 && sectionNum <= numSections && argc > 2) {

        double xi[maxNumSections];
        double L = theCoordTransf->getInitialLength();
        beamInt->getSectionLocations(numSections, L, xi);

        output.tag("GaussPointOutput");
        output.attr("number", sectionNum);
        output.attr("eta", xi[sectionNum - 1] * L);

        theResponse =
            theSections[sectionNum - 1]->setResponse(&argv[2], argc - 2, output);

        output.endTag();
      } else if (sectionNum == 0) { // argv[1] was not an int, we want all sections,

        CompositeResponse *theCResponse = new CompositeResponse();
        int numResponse                 = 0;
        double xi[maxNumSections];
        double L = theCoordTransf->getInitialLength();
        beamInt->getSectionLocations(numSections, L, xi);

        for (int i = 0; i < numSections; i++) {

          output.tag("GaussPointOutput");
          output.attr("number", i + 1);
          output.attr("eta", xi[i] * L);

          Response *theSectionResponse =
              theSections[i]->setResponse(&argv[1], argc - 1, output);

          output.endTag();

          if (theSectionResponse != 0) {
            numResponse = theCResponse->addResponse(theSectionResponse);
          }
        }

        if (numResponse == 0) // no valid responses found
          delete theCResponse;
        else
          theResponse = theCResponse;
      }
    }
  }

  output.endTag();
  return theResponse;
}

int
DisplDeltaFrame3d::getResponse(int responseID, Information &eleInfo)
{
  double N, V, M1, M2, T;
  double L        = theCoordTransf->getInitialLength();
  double oneOverL = 1.0 / L;

  if (responseID == 1)
    return eleInfo.setVector(this->getResistingForce());

  else if (responseID == 12)
    return eleInfo.setVector(this->getRayleighDampingForces());

  else if (responseID == 2) {
    // TODO(cmp)
    P.Zero();
    return eleInfo.setVector(P);
  }

  // Chord rotation
  else if (responseID == 3)
    return eleInfo.setVector(theCoordTransf->getBasicTrialDisp());

  // Plastic rotation
  else if (responseID == 4) {
    static Vector vp(6);
    static Vector ve(6);
    // TODO(cmp)
//  const Matrix &kb = this->getInitialBasicStiff();
//  kb.Solve(q, ve);
//  vp = theCoordTransf->getBasicTrialDisp();
//  vp -= ve;
    return eleInfo.setVector(vp);
  }

  else if (responseID == 10) {
    double L = theCoordTransf->getInitialLength();
    double pts[maxNumSections];
    beamInt->getSectionLocations(numSections, L, pts);
    Vector locs(numSections);
    for (int i = 0; i < numSections; i++)
      locs(i) = pts[i] * L;
    return eleInfo.setVector(locs);
  }

  else if (responseID == 11) {
    double L = theCoordTransf->getInitialLength();
    double wts[maxNumSections];
    beamInt->getSectionWeights(numSections, L, wts);
    Vector weights(numSections);
    for (int i = 0; i < numSections; i++)
      weights(i) = wts[i] * L;
    return eleInfo.setVector(weights);
  }

  else if (responseID == 110) {
    ID tags(numSections);
    for (int i = 0; i < numSections; i++)
      tags(i) = theSections[i]->getTag();
    return eleInfo.setID(tags);
  }

  else
    return -1;
}

int
DisplDeltaFrame3d::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;

  if (strstr(argv[0], "sectionX") != 0) {
    if (argc < 3)
      return -1;

    float sectionLoc = atof(argv[1]);

    double xi[maxNumSections];
    double L = theCoordTransf->getInitialLength();
    beamInt->getSectionLocations(numSections, L, xi);

    sectionLoc /= L;

    float minDistance = fabs(xi[0] - sectionLoc);
    int sectionNum    = 0;
    for (int i = 1; i < numSections; i++) {
      if (fabs(xi[i] - sectionLoc) < minDistance) {
        minDistance = fabs(xi[i] - sectionLoc);
        sectionNum  = i;
      }
    }
    return theSections[sectionNum]->setParameter(&argv[2], argc - 2, param);
  }

  // If the parameter belongs to a section or lower
  if (strstr(argv[0], "section") != 0) {

    if (argc < 3)
      return -1;

    // Get section number
    int sectionNum = atoi(argv[1]);

    if (sectionNum > 0 && sectionNum <= numSections)
      return theSections[sectionNum - 1]->setParameter(&argv[2], argc - 2, param);
    else
      return -1;
  }

  else if (strstr(argv[0], "integration") != 0) {

    if (argc < 2)
      return -1;

    return beamInt->setParameter(&argv[1], argc - 1, param);
  }

  // Default, send to every object
  int ok     = 0;
  int result = 0;

  for (int i = 0; i < numSections; i++) {
    ok = theSections[i]->setParameter(argv, argc, param);
    if (ok != -1)
      result = ok;
  }

  ok = beamInt->setParameter(argv, argc, param);
  if (ok != -1)
    result = ok;

  return result;
}

int
DisplDeltaFrame3d::updateParameter(int parameterID, Information &info)
{
  if (parameterID == 1) {
    rho = info.theDouble;
    return 0;
  } else
    return -1;
}

int
DisplDeltaFrame3d::activateParameter(int passedParameterID)
{
  parameterID = passedParameterID;
  return 0;
}


int
DisplDeltaFrame3d::sendSelf(int commitTag, Channel &theChannel)
{
  // place the integer data into an ID

  int dbTag = this->getDbTag();
  int i, j;
  int loc = 0;

  static Vector data(16);
  data(0)            = this->getTag();
  data(1)            = connectedExternalNodes(0);
  data(2)            = connectedExternalNodes(1);
  data(3)            = numSections;
  data(4)            = theCoordTransf->getClassTag();
  int crdTransfDbTag = theCoordTransf->getDbTag();
  if (crdTransfDbTag == 0) {
    crdTransfDbTag = theChannel.getDbTag();
    if (crdTransfDbTag != 0)
      theCoordTransf->setDbTag(crdTransfDbTag);
  }
  data(5)          = crdTransfDbTag;
  data(6)          = beamInt->getClassTag();
  int beamIntDbTag = beamInt->getDbTag();
  if (beamIntDbTag == 0) {
    beamIntDbTag = theChannel.getDbTag();
    if (beamIntDbTag != 0)
      beamInt->setDbTag(beamIntDbTag);
  }
  data(7)  = beamIntDbTag;
  data(8)  = rho;
  data(9)  = cMass;
  data(10) = alphaM;
  data(11) = betaK;
  data(12) = betaK0;
  data(13) = betaKc;

  if (theChannel.sendVector(dbTag, commitTag, data) < 0) {
    opserr << "DisplDeltaFrame3d::sendSelf() - failed to send data Vector\n";
    return -1;
  }

  // send the coordinate transformation
  if (theCoordTransf->sendSelf(commitTag, theChannel) < 0) {
    opserr << "DisplDeltaFrame3d::sendSelf() - failed to send crdTranf\n";
    return -1;
  }

  // send the beam integration
  if (beamInt->sendSelf(commitTag, theChannel) < 0) {
    opserr << "DisplDeltaFrame3d::sendSelf() - failed to send beamInt\n";
    return -1;
  }

  //
  // send an ID for the sections containing each sections dbTag and classTag
  // if section ha no dbTag get one and assign it
  //

  ID idSections(2 * numSections);
  loc = 0;
  for (i = 0; i < numSections; i++) {
    int sectClassTag = theSections[i]->getClassTag();
    int sectDbTag    = theSections[i]->getDbTag();
    if (sectDbTag == 0) {
      sectDbTag = theChannel.getDbTag();
      theSections[i]->setDbTag(sectDbTag);
    }

    idSections(loc)     = sectClassTag;
    idSections(loc + 1) = sectDbTag;
    loc += 2;
  }

  if (theChannel.sendID(dbTag, commitTag, idSections) < 0) {
    opserr << "DisplDeltaFrame3d::sendSelf() - failed to send ID data\n";
    return -1;
  }

  //
  // send the sections
  //

  for (int j = 0; j < numSections; j++) {
    if (theSections[j]->sendSelf(commitTag, theChannel) < 0) {
      opserr << "DisplDeltaFrame3d::sendSelf() - section " << j
             << "failed to send itself\n";
      return -1;
    }
  }

  return 0;
}

int
DisplDeltaFrame3d::recvSelf(int commitTag, Channel &theChannel,
                               FEM_ObjectBroker &theBroker)
{
  //
  // receive the integer data containing tag, numSections and coord transformation info
  //
  int dbTag = this->getDbTag();
  int i;

  static Vector data(16);

  if (theChannel.recvVector(dbTag, commitTag, data) < 0) {
    opserr << "DisplDeltaFrame3d::recvSelf() - failed to recv data Vector\n";
    return -1;
  }

  this->setTag((int)data(0));
  connectedExternalNodes(0) = (int)data(1);
  connectedExternalNodes(1) = (int)data(2);
  int nSect                 = (int)data(3);
  int crdTransfClassTag     = (int)data(4);
  int crdTransfDbTag        = (int)data(5);

  int beamIntClassTag = (int)data(6);
  int beamIntDbTag    = (int)data(7);

  rho   = data(8);
  cMass = (int)data(9);

  alphaM = data(10);
  betaK  = data(11);
  betaK0 = data(12);
  betaKc = data(13);

  // create a new crdTransf object if one needed
  if (theCoordTransf == 0 || theCoordTransf->getClassTag() != crdTransfClassTag) {
    if (theCoordTransf != 0)
      delete theCoordTransf;

    // TODO(cmp) - add FrameTransform to ObjBroker
    theCoordTransf = nullptr; // theBroker.getNewFrameTransform3d(crdTransfClassTag);

    if (theCoordTransf == nullptr) {
      opserr << "DisplDeltaFrame3d::recvSelf() - "
             << "failed to obtain a CrdTrans object with classTag" << crdTransfClassTag
             << "\n";
      return -2;
    }
  }

  theCoordTransf->setDbTag(crdTransfDbTag);

  // invoke recvSelf on the crdTransf object
  if (theCoordTransf->recvSelf(commitTag, theChannel, theBroker) < 0) {
    opserr << "DisplDeltaFrame3d::sendSelf() - failed to recv crdTranf\n";
    return -3;
  }

  // create a new beamInt object if one needed
  if (beamInt == 0 || beamInt->getClassTag() != beamIntClassTag) {
    if (beamInt != 0)
      delete beamInt;

    beamInt = theBroker.getNewBeamIntegration(beamIntClassTag);

    if (beamInt == 0) {
      opserr << "DisplDeltaFrame3d::recvSelf() - failed to obtain the beam "
                "integration object with classTag"
             << beamIntClassTag << "\n";
      return -3;
    }
  }

  beamInt->setDbTag(beamIntDbTag);

  // invoke recvSelf on the beamInt object
  if (beamInt->recvSelf(commitTag, theChannel, theBroker) < 0) {
    opserr << "DisplDeltaFrame3d::sendSelf() - failed to recv beam integration\n";
    return -3;
  }

  //
  // recv an ID for the sections containing each sections dbTag and classTag
  //

  ID idSections(2 * nSect);
  int loc = 0;

  if (theChannel.recvID(dbTag, commitTag, idSections) < 0) {
    opserr << "DisplDeltaFrame3d::recvSelf() - failed to recv ID data\n";
    return -1;
  }

  //
  // now receive the sections
  //

  if (numSections != nSect) {

    //
    // we do not have correct number of sections, must delete the old and create
    // new ones before can recvSelf on the sections
    //

    // delete the old
    if (numSections != 0) {
      for (int i = 0; i < numSections; i++)
        delete theSections[i];
      delete[] theSections;
    }

    // create a new array to hold pointers
    theSections = new FrameSection *[nSect];

    // create a section and recvSelf on it
    numSections = nSect;
    loc         = 0;

    for (i = 0; i < numSections; i++) {
      int sectClassTag = idSections(loc);
      int sectDbTag    = idSections(loc + 1);
      loc += 2;
      // TODO(cmp) add FrameSection to broker
//    theSections[i] = theBroker.getNewSection(sectClassTag);
      if (theSections[i] == nullptr) {
        opserr << "DisplDeltaFrame3d::recvSelf() - Broker could not create Section of "
                  "class type"
               << sectClassTag << "\n";
        return -1;
      }
      theSections[i]->setDbTag(sectDbTag);
      if (theSections[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
        opserr << "DisplDeltaFrame3d::recvSelf() - section " << i
               << "failed to recv itself\n";
        return -1;
      }
    }

  } else {

    //
    // for each existing section, check it is of correct type
    // (if not delete old & create a new one) then recvSelf on it
    //

    loc = 0;
    for (i = 0; i < numSections; i++) {
      int sectClassTag = idSections(loc);
      int sectDbTag    = idSections(loc + 1);
      loc += 2;

      // check of correct type
      if (theSections[i]->getClassTag() != sectClassTag) {
        // delete the old section[i] and create a new one
        delete theSections[i];
      // TODO(cmp) add FrameSection to broker
//      theSections[i] = theBroker.getNewSection(sectClassTag);
        if (theSections[i] == nullptr) {
          opserr << "DisplDeltaFrame3d::recvSelf() - Broker could not create Section "
                    "of class type"
                 << sectClassTag << "\n";
          return -1;
        }
      }

      // recvSelf on it
      theSections[i]->setDbTag(sectDbTag);
      if (theSections[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
        opserr << "DisplDeltaFrame3d::recvSelf() - section " << i
               << "failed to recv itself\n";
        return -1;
      }
    }
  }

  return 0;
}

