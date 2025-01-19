//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
// Description: This file contains the class implementation of FrameSolidSection3d.
//
// Adapted from NDFiberSection3d
//
// Written: CMP,MHS
// Created: 2012
//
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <Channel.h>
#include <Vector.h>
#include <VectorND.h>
#include <Matrix.h>
#include <MatrixND.h>
#include <classTags.h>
#include <FrameSolidSection3d.h>
#include <ID.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <SensitiveResponse.h>
typedef SensitiveResponse<FrameSection> SectionResponse;
#include <NDMaterial.h>
#include <Parameter.h>
// #include <elementAPI.h>

#define SEC_TAG_FrameSolidSection3d 0

using OpenSees::VectorND;
using OpenSees::MatrixND;

ID FrameSolidSection3d::code(6);

FrameSolidSection3d::FrameSolidSection3d(int tag, int num, double a, bool compCentroid): 
    FrameSection(tag, SEC_TAG_FrameSolidSection3d),
    Abar(0.0), QyBar(0.0), QzBar(0.0), yBar(0.0), zBar(0.0), 
    computeCentroid(compCentroid),
    alpha(a), e(nsr), s(0), ks(0), 
    parameterID(0), dedh(nsr)
{
    s = new Vector(sData, nsr);
    ks = new Matrix(kData, nsr, nsr);

    code(0) = SECTION_RESPONSE_P;
    code(1) = SECTION_RESPONSE_MZ;
    code(2) = SECTION_RESPONSE_MY;
    code(3) = SECTION_RESPONSE_VY;
    code(4) = SECTION_RESPONSE_VZ;
    code(5) = SECTION_RESPONSE_T;
}

// constructor for blank object that recvSelf needs to be invoked upon
FrameSolidSection3d::FrameSolidSection3d():
  FrameSection(0, SEC_TAG_FrameSolidSection3d),
  Abar(0.0), QyBar(0.0), QzBar(0.0), 
  yBar(0.0), zBar(0.0), 
  computeCentroid(true),
  alpha(1.0), e(6), s(0), ks(0),
  parameterID(0), dedh(6)
{
  s = new Vector(sData, 6);
  ks = new Matrix(kData, 6, 6);

  code(0) = SECTION_RESPONSE_P;
  code(1) = SECTION_RESPONSE_MZ;
  code(2) = SECTION_RESPONSE_MY;
  code(3) = SECTION_RESPONSE_VY;
  code(4) = SECTION_RESPONSE_VZ;
  code(5) = SECTION_RESPONSE_T;
}

int
FrameSolidSection3d::addFiber(NDMaterial& theMat, 
                              double Area, 
                              double yLoc, 
                              double zLoc)
{
  FiberData fiber {theMat.getCopy("BeamFiber"), yLoc, zLoc, Area};
  fibers.emplace_back(fiber);

  if (fibers[fibers.size()-1].material == nullptr) {
    opserr <<"FrameSolidSection3d::addFiber -- failed to get copy of a Material\n";
    return -1;
  }

  // Recompute centroid
  if (computeCentroid) {
    Abar  += Area;
    QzBar += yLoc*Area;
    QyBar += zLoc*Area;
    
    yBar = QzBar/Abar;
    zBar = QyBar/Abar;
  }
  
  return 0;
}


// destructor:
FrameSolidSection3d::~FrameSolidSection3d()
{
  for (auto fiber : fibers) {
    if (fiber.material != nullptr)
      delete fiber.material;
  }

  if (s != 0)
    delete s;

  if (ks != 0)
    delete ks;
}

//      0  5 4       1       2  3
// a = [1 -y z       0       0  0
//      0  0 0 sqrt(a)       0 -z
//      0  0 0       0 sqrt(a)  y]
int
FrameSolidSection3d::setTrialSectionDeformation(const Vector &deforms)
{
  e = deforms;

  s->Zero();
  ks->Zero();

  double e0 = deforms(0),
         e1 = deforms(1),
         e2 = deforms(2),
         e3 = deforms(3),
         e4 = deforms(4),
         e5 = deforms(5);


  double rootAlpha = 1.0;
  if (alpha != 1.0)
    rootAlpha = sqrt(alpha);

  int res = 0;
  for (int i = 0; i < fibers.size(); i++) {
    auto & fiber = fibers[i];
    NDMaterial *theMat = fibers[i].material;
    const double y  = fiber.y - yBar;
    const double z  = fiber.z - zBar;
    const double A  = fiber.wgt;

    // determine material strain and set it
    VectorND<3> eps;
    eps[0] = e0 - y*e1 + z*e2;
    eps[1] = rootAlpha*e3 - z*e5;
    eps[2] = rootAlpha*e4 + y*e5;

    res += theMat->setTrialStrain(eps);
    const Vector &stress  = theMat->getStress();
    const Matrix &tangent = theMat->getTangent();

    double d00 = tangent(0,0)*A;
    double d01 = tangent(0,1)*A;
    double d02 = tangent(0,2)*A;
    double d10 = tangent(1,0)*A;
    double d11 = tangent(1,1)*A;
    double d12 = tangent(1,2)*A;
    double d20 = tangent(2,0)*A;
    double d21 = tangent(2,1)*A;
    double d22 = tangent(2,2)*A;

    Matrix &ksi = *ks;
    Vector &si  = *s;

    double y2 = y*y;
    double z2 = z*z;
    double yz = y*z;
    double tmp;
    // Bending terms
    ksi(0,0) +=    d00;
    ksi(5,5) += y2*d00;
    ksi(4,4) += z2*d00;
    tmp = -y*d00;
    ksi(0,5) += tmp;
    ksi(5,0) += tmp;
    tmp = z*d00;
    ksi(0,4) += tmp;
    ksi(4,0) += tmp;
    tmp = -yz*d00;
    ksi(5,4) += tmp;
    ksi(4,5) += tmp;
    
    // Shear terms
    ksi(1,1) += alpha*d11;
    ksi(1,2) += alpha*d12;
    ksi(2,1) += alpha*d21;
    ksi(2,2) += alpha*d22;
    
    // Torsion term
    ksi(3,3) += z2*d11 - yz*(d12+d21) + y2*d22;
    
    // Bending-torsion coupling terms
    tmp = -z*d01 + y*d02;
    ksi(0,3) += tmp;
    ksi(5,3) -= y*tmp;
    ksi(4,3) += z*tmp;
    tmp = -z*d10 + y*d20;
    ksi(3,0) += tmp;
    ksi(3,5) -= y*tmp;
    ksi(3,4) += z*tmp;
    
    // Hit tangent terms with rootAlpha
    d01 *= rootAlpha; d02 *= rootAlpha;
    d10 *= rootAlpha; d11 *= rootAlpha; d12 *= rootAlpha;
    d20 *= rootAlpha; d21 *= rootAlpha; d22 *= rootAlpha;
    
    // Bending-shear coupling terms
    ksi(0,1) += d01;
    ksi(0,2) += d02;
    ksi(5,1) -= y*d01;
    ksi(5,2) -= y*d02;
    ksi(4,1) += z*d01;
    ksi(4,2) += z*d02;
    ksi(1,0) += d10;
    ksi(2,0) += d20;
    ksi(1,5) -= y*d10;
    ksi(2,5) -= y*d20;
    ksi(1,4) += z*d10;
    ksi(2,4) += z*d20;
    
    // Torsion-shear coupling terms
    y2 =  y*d22;
    z2 = -z*d11;
    ksi(3,1) +=  z2 + y*d21;
    ksi(3,2) += -z*d12 + y2;
    ksi(1,3) +=  z2 + y*d12;
    ksi(2,3) += -z*d21 + y2;

    double sig0 = stress(0)*A;
    double sig1 = stress(1)*A;
    double sig2 = stress(2)*A;

    si(0) +=    sig0;
    si(1) += -y*sig0;
    si(2) +=  z*sig0;
    si(3) += rootAlpha*sig1;
    si(4) += rootAlpha*sig2;
    si(5) += -z*sig1 + y*sig2;
  }

  return res;
}

const Vector&
FrameSolidSection3d::getSectionDeformation()
{
  return e;
}

const Matrix&
FrameSolidSection3d::getInitialTangent()
{
  static double kInitial[36];
  static Matrix ki(kInitial, 6, 6);
  ki.Zero();

  double rootAlpha = 1.0;
  if (alpha != 1.0)
    rootAlpha = sqrt(alpha);

  for (int i = 0; i < fibers.size(); i++) {
    NDMaterial *theMat = fibers[i].material;
    double y = fibers[i].y - yBar;
    double z = fibers[i].z - zBar;
    double A = fibers[i].wgt;

    double y2 = y*y;
    double z2 = z*z;
    double yz = y*z;
    double tmp;

    const Matrix &tangent = theMat->getInitialTangent();

    double d00 = tangent(0,0)*A;
    double d01 = tangent(0,1)*A;
    double d02 = tangent(0,2)*A;
    double d10 = tangent(1,0)*A;
    double d11 = tangent(1,1)*A;
    double d12 = tangent(1,2)*A;
    double d20 = tangent(2,0)*A;
    double d21 = tangent(2,1)*A;
    double d22 = tangent(2,2)*A;

    // Bending terms
    ki(0,0) += d00;
    ki(1,1) += y2*d00;
    ki(2,2) += z2*d00;
    tmp = -y*d00;
    ki(0,1) += tmp;
    ki(1,0) += tmp;
    tmp = z*d00;
    ki(0,2) += tmp;
    ki(2,0) += tmp;
    tmp = -yz*d00;
    ki(1,2) += tmp;
    ki(2,1) += tmp;
    
    // Shear terms
    ki(3,3) += alpha*d11;
    ki(3,4) += alpha*d12;
    ki(4,3) += alpha*d21;
    ki(4,4) += alpha*d22;
    
    // Torsion term
    ki(5,5) += z2*d11 - yz*(d12+d21) + y2*d22;
    
    // Bending-torsion coupling terms
    tmp = -z*d01 + y*d02;
    ki(0,5) += tmp;
    ki(1,5) -= y*tmp;
    ki(2,5) += z*tmp;
    tmp = -z*d10 + y*d20;
    ki(5,0) += tmp;
    ki(5,1) -= y*tmp;
    ki(5,2) += z*tmp;
    
    // Hit tangent terms with rootAlpha
    d01 *= rootAlpha; d02 *= rootAlpha;
    d10 *= rootAlpha; d11 *= rootAlpha; d12 *= rootAlpha;
    d20 *= rootAlpha; d21 *= rootAlpha; d22 *= rootAlpha;
    
    // Bending-shear coupling terms
    ki(0,3) += d01;
    ki(0,4) += d02;
    ki(1,3) -= y*d01;
    ki(1,4) -= y*d02;
    ki(2,3) += z*d01;
    ki(2,4) += z*d02;
    ki(3,0) += d10;
    ki(4,0) += d20;
    ki(3,1) -= y*d10;
    ki(4,1) -= y*d20;
    ki(3,2) += z*d10;
    ki(4,2) += z*d20;
    
    // Torsion-shear coupling terms
    y2 =  y*d22;
    z2 = -z*d11;
    ki(5,3) +=  z2 + y*d21;
    ki(5,4) += -z*d12 + y2;
    ki(3,5) +=  z2 + y*d12;
    ki(4,5) += -z*d21 + y2;
  }

  return ki;
}

const Matrix&
FrameSolidSection3d::getSectionTangent()
{
  return *ks;
}

const Vector&
FrameSolidSection3d::getStressResultant()
{
  return *s;
}

FrameSection*
FrameSolidSection3d::getFrameCopy()
{
  FrameSolidSection3d *theCopy = new FrameSolidSection3d ();
  theCopy->setTag(this->getTag());


  if (fibers.size() != 0) {
    // TODO
  }

  theCopy->e = e;
  theCopy->QzBar = QzBar;
  theCopy->QyBar = QyBar;
  theCopy->Abar = Abar;
  theCopy->yBar = yBar;
  theCopy->zBar = zBar;
  theCopy->computeCentroid = computeCentroid;
  theCopy->alpha = alpha;
  theCopy->parameterID = parameterID;

  return theCopy;
}

const ID&
FrameSolidSection3d::getType()
{
  return code;
}

int
FrameSolidSection3d::getOrder () const
{
  return nsr;
}

int
FrameSolidSection3d::commitState()
{
  int err = 0;

  for (auto& fiber: fibers)
    err += fiber.material->commitState();

  return err;
}

int
FrameSolidSection3d::revertToLastCommit()
{
  int err = 0;

  ks->Zero();
  s->Zero();

  double rootAlpha = 1.0;
  if (alpha != 1.0)
    rootAlpha = sqrt(alpha);

  for (int i = 0; i < fibers.size(); i++) {
    auto& fiber = fibers[i];
    NDMaterial *theMat = fibers[i].material;
    const double y  = fiber.y - yBar;
    const double z  = fiber.z - zBar;
    const double A  = fiber.wgt;

    double y2 = y*y;
    double z2 = z*z;
    double yz = y*z;
    double tmp;

    // invoke revertToLast on the material
    err += theMat->revertToLastCommit();

    // get material stress & tangent for this strain and determine ks and fs
    const Matrix &tangent = theMat->getTangent();
    const Vector &stress = theMat->getStress();

    double d00 = tangent(0,0)*A;
    double d01 = tangent(0,1)*A;
    double d02 = tangent(0,2)*A;
    double d10 = tangent(1,0)*A;
    double d11 = tangent(1,1)*A;
    double d12 = tangent(1,2)*A;
    double d20 = tangent(2,0)*A;
    double d21 = tangent(2,1)*A;
    double d22 = tangent(2,2)*A;

    Matrix &ksi = *ks;
    Vector &si = *s;

    // Bending terms
    ksi(0,0) += d00;
    ksi(5,5) += y2*d00;
    ksi(4,4) += z2*d00;
    tmp = -y*d00;
    ksi(0,5) += tmp;
    ksi(5,0) += tmp;
    tmp = z*d00;
    ksi(0,4) += tmp;
    ksi(4,0) += tmp;
    tmp = -yz*d00;
    ksi(5,4) += tmp;
    ksi(4,5) += tmp;
    
    // Shear terms
    ksi(1,1) += alpha*d11;
    ksi(1,2) += alpha*d12;
    ksi(2,1) += alpha*d21;
    ksi(2,2) += alpha*d22;
    
    // Torsion term
    ksi(3,3) += z2*d11 - yz*(d12+d21) + y2*d22;
    
    // Bending-torsion coupling terms
    tmp       = -z*d01 + y*d02;
    ksi(0,3) += tmp;
    ksi(5,3) -= y*tmp;
    ksi(4,3) += z*tmp;
    tmp       = -z*d10 + y*d20;
    ksi(3,0) += tmp;
    ksi(3,5) -= y*tmp;
    ksi(3,4) += z*tmp;
    
    // Hit tangent terms with rootAlpha
    d01 *= rootAlpha; d02 *= rootAlpha;
    d10 *= rootAlpha; d11 *= rootAlpha; d12 *= rootAlpha;
    d20 *= rootAlpha; d21 *= rootAlpha; d22 *= rootAlpha;
    
    // Bending-shear coupling terms
    ksi(0,1) += d01;
    ksi(0,2) += d02;
    ksi(5,1) -= y*d01;
    ksi(5,2) -= y*d02;
    ksi(4,1) += z*d01;
    ksi(4,2) += z*d02;
    ksi(1,0) += d10;
    ksi(2,0) += d20;
    ksi(1,5) -= y*d10;
    ksi(2,5) -= y*d20;
    ksi(1,4) += z*d10;
    ksi(2,4) += z*d20;
    
    // Torsion-shear coupling terms
    y2 =  y*d22;
    z2 = -z*d11;
    ksi(3,1) +=  z2 + y*d21;
    ksi(3,2) += -z*d12 + y2;
    ksi(1,3) +=  z2 + y*d12;
    ksi(2,3) += -z*d21 + y2;

    double sig0 = stress(0)*A;
    double sig1 = stress(1)*A;
    double sig2 = stress(2)*A;

    si(0) +=    sig0;
    si(1) += -y*sig0;
    si(2) +=  z*sig0;
    si(3) += rootAlpha*sig1;
    si(4) += rootAlpha*sig2;
    si(5) += -z*sig1 + y*sig2;
  }

  if (alpha != 1.0) {

  }

  return err;
}

int
FrameSolidSection3d::revertToStart()
{
  // revert the fibers to start    
  int err = 0;

  ks->Zero();
  s->Zero();
  

  double rootAlpha = 1.0;
  if (alpha != 1.0)
    rootAlpha = sqrt(alpha);

  for (int i = 0; i < fibers.size(); i++) {
    auto & fiber = fibers[i];
    NDMaterial *theMat = fibers[i].material;

    double y  = fiber.y   - yBar;
    double z  = fiber.z - zBar;
    double A  = fiber.wgt;

    double y2 = y*y;
    double z2 = z*z;
    double yz = y*z;
    double tmp;

    // invoke revertToLast on the material
    err += theMat->revertToStart();

    // get material stress & tangent for this strain and determine ks and fs
    const Matrix &tangent = theMat->getTangent();
    const Vector &stress = theMat->getStress();

    double d00 = tangent(0,0)*A;
    double d01 = tangent(0,1)*A;
    double d02 = tangent(0,2)*A;
    double d10 = tangent(1,0)*A;
    double d11 = tangent(1,1)*A;
    double d12 = tangent(1,2)*A;
    double d20 = tangent(2,0)*A;
    double d21 = tangent(2,1)*A;
    double d22 = tangent(2,2)*A;

    Matrix &ksi = *ks;
    Vector &si = *s;

    // Bending terms
    ksi(0,0) += d00;
    ksi(5,5) += y2*d00;
    ksi(4,4) += z2*d00;
    tmp = -y*d00;
    ksi(0,5) += tmp;
    ksi(5,0) += tmp;
    tmp = z*d00;
    ksi(0,4) += tmp;
    ksi(4,0) += tmp;
    tmp = -yz*d00;
    ksi(5,4) += tmp;
    ksi(4,5) += tmp;
    
    // Shear terms
    ksi(1,1) += alpha*d11;
    ksi(1,2) += alpha*d12;
    ksi(2,1) += alpha*d21;
    ksi(2,2) += alpha*d22;
    
    // Torsion term
    ksi(3,3) += z2*d11 - yz*(d12+d21) + y2*d22;
    
    // Bending-torsion coupling terms
    tmp = -z*d01 + y*d02;
    ksi(0,3) += tmp;
    ksi(5,3) -= y*tmp;
    ksi(4,3) += z*tmp;
    tmp = -z*d10 + y*d20;
    ksi(3,0) += tmp;
    ksi(3,5) -= y*tmp;
    ksi(3,4) += z*tmp;
    
    // Hit tangent terms with rootAlpha
    d01 *= rootAlpha; d02 *= rootAlpha;
    d10 *= rootAlpha; d11 *= rootAlpha; d12 *= rootAlpha;
    d20 *= rootAlpha; d21 *= rootAlpha; d22 *= rootAlpha;
    
    // Bending-shear coupling terms
    ksi(0,1) += d01;
    ksi(0,2) += d02;
    ksi(5,1) -= y*d01;
    ksi(5,2) -= y*d02;
    ksi(4,1) += z*d01;
    ksi(4,2) += z*d02;
    ksi(1,0) += d10;
    ksi(2,0) += d20;
    ksi(1,5) -= y*d10;
    ksi(2,5) -= y*d20;
    ksi(1,4) += z*d10;
    ksi(2,4) += z*d20;
    
    // Torsion-shear coupling terms
    y2 =  y*d22;
    z2 = -z*d11;
    ksi(3,1) +=  z2 + y*d21;
    ksi(3,2) += -z*d12 + y2;
    ksi(1,3) +=  z2 + y*d12;
    ksi(2,3) += -z*d21 + y2;

    double sig0 = stress(0)*A;
    double sig1 = stress(1)*A;
    double sig2 = stress(2)*A;

    si[0] += sig0;
    si[1] += -y*sig0;
    si[2] += z*sig0;
    si[3] += rootAlpha*sig1;
    si[4] += rootAlpha*sig2;
    si[5] += -z*sig1 + y*sig2;
  }

  if (alpha != 1.0) {
  }

  return err;
}

int
FrameSolidSection3d::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;
  return res;
}

int
FrameSolidSection3d::recvSelf(int commitTag, Channel &theChannel,
                         FEM_ObjectBroker &theBroker)
{
  int res = 0;
  return res;
}

void
FrameSolidSection3d::Print(OPS_Stream &s, int flag)
{
  s << "\nFrameSolidSection3d, tag: " << this->getTag() << endln;
  s << "\tSection code: " << code;
  s << "\tNumber of Fibers: " << fibers.size() << endln;
  s << "\tCentroid (y,z): " << yBar << ' ' << zBar << endln;
  s << "\tShape factor, alpha = " << alpha << endln;

  if (flag == 1) {
    for (int i = 0; i < fibers.size(); i++) {
      auto & fiber = fibers[i];
      s << "\nLocation (y,z) = " << fiber.y << ' ' << fiber.z;
      s << "\nArea = " << fiber.wgt << endln;
      fibers[i].material->Print(s, flag);
    }
  }
}

Response*
FrameSolidSection3d::setResponse(const char **argv, int argc,
                              OPS_Stream &output)
{
  Response *theResponse = nullptr;

  if (argc > 2 && strcmp(argv[0],"fiber") == 0) {

    
    int key = fibers.size();
    int passarg = 2;
    
    if (argc <= 3) {                  // fiber number was input directly
      
      key = atoi(argv[1]);
      
    } else if (argc > 4) {  // find fiber closest to coord. with mat tag
      
      int matTag = atoi(argv[3]);
      double yCoord = atof(argv[1]);
      double zCoord = atof(argv[2]);
      double closestDist = 0;
      double ySearch, zSearch, dy, dz;
      double distance;
      int j;
      // Find first fiber with specified material tag
      for (j = 0; j < fibers.size(); j++) {
        auto& fiber = fibers[j];
        if (matTag == fiber.material->getTag()) {
          ySearch = fiber.y;
          zSearch = fiber.z;

          dy = ySearch-yCoord;
          dz = zSearch-zCoord;
          closestDist = dy*dy + dz*dz;
          key = j;
          break;
        }
      }
      // Search the remaining fibers
      for ( ; j < fibers.size(); j++) {
        auto& fiber = fibers[j];
        if (matTag == fiber.material->getTag()) {
          ySearch = fiber.y;
          zSearch = fiber.z;

          dy = ySearch - yCoord;
          dz = zSearch - zCoord;
          distance = dy*dy + dz*dz;
          if (distance < closestDist) {
            closestDist = distance;
            key = j;
          }
        }
      }
      passarg = 4;
    }
    
    else { // fiber near-to coordinate specified
      
      double yCoord = atof(argv[1]);
      double zCoord = atof(argv[2]);
      double closestDist;
      double ySearch, zSearch, dy, dz;
      double distance;
      
      ySearch = fibers[0].y;
      zSearch = fibers[1].z;

      dy = ySearch-yCoord;
      dz = zSearch-zCoord;
      closestDist = dy*dy + dz*dz;
      key = 0;
      for (int j = 1; j < fibers.size(); j++) {
        auto& fiber = fibers[j];
        ySearch = fiber.y;
        zSearch = fiber.z;
        dy = ySearch-yCoord;
        dz = zSearch-zCoord;
        distance = dy*dy + dz*dz;
        if (distance < closestDist) {
          closestDist = distance;
          key = j;
        }
      }
      passarg = 3;
    }
    
    if (key < fibers.size() && key >= 0) {
      output.tag("FiberOutput");
      output.attr("yLoc", fibers[key].y);
      output.attr("zLoc", fibers[key].z);
      output.attr("area", fibers[key].wgt);
      
      theResponse = fibers[key].material->setResponse(&argv[passarg], argc-passarg, output);
      
      output.endTag();
    }

  }

  if (theResponse == nullptr)
    return FrameSection::setResponse(argv, argc, output);

  return theResponse;
}


int 
FrameSolidSection3d::getResponse(int responseID, Information &sectInfo)
{
  // Just call the base class method ... don't need to define
  // this function, but keeping it here just for clarity
  return FrameSection::getResponse(responseID, sectInfo);
}



// AddingSensitivity:BEGIN ////////////////////////////////////
int
FrameSolidSection3d::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;

  int result = -1;

  if (strstr(argv[0],"alpha") != 0)
    return param.addObject(1, this);

  // Check if the parameter belongs to the material (only option for now)
  if (strstr(argv[0],"material") != 0) {
    
    if (argc < 3)
      return 0;

    // Get the tag of the material
    int materialTag = atoi(argv[1]);
    
    // Loop over fibers to find the right material
    for (int i = 0; i < fibers.size(); i++)
      if (materialTag == fibers[i].material->getTag()) {
        int ok = fibers[i].material->setParameter(&argv[2], argc-2, param);
        if (ok != -1)
          result = ok;
      }
    return result;
  }
  int ok = 0; 
  for (int i = 0; i < fibers.size(); i++) {
    ok = fibers[i].material->setParameter(argv, argc, param);
    if (ok != -1)
      result = ok;
  }
  return result;
}

int
FrameSolidSection3d::updateParameter(int paramID, Information &info)
{
  switch(paramID) {
  case 1:
    alpha = info.theDouble;
    return 0;
  default:
    return -1;
  }
}

int
FrameSolidSection3d::activateParameter(int paramID)
{
  parameterID = paramID;

  return 0;
}

const Vector &
FrameSolidSection3d::getSectionDeformationSensitivity(int gradIndex)
{
  return dedh;
}

const Vector &
FrameSolidSection3d::getStressResultantSensitivity(int gradIndex, bool conditional)
{
  static Vector ds(6);
  
  ds.Zero();
  
  static Vector stress(3);
  static Vector dsigdh(3);
  static Vector sig_dAdh(3);
  static Matrix tangent(3,3);

  static double dydh[10000];
  static double dzdh[10000];
  static double areaDeriv[10000];

  {
    for (int i = 0; i < fibers.size(); i++) {
      dydh[i] = 0.0;
      dzdh[i] = 0.0;
      areaDeriv[i] = 0.0;
    }
  }
  
  double rootAlpha = 1.0;
  if (alpha != 1.0)
    rootAlpha = sqrt(alpha);

  double drootAlphadh = 0.0;
  if (parameterID == 1)
    drootAlphadh = 0.5/rootAlpha;

  for (int i = 0; i < fibers.size(); i++) {
    double y = fibers[i].y - yBar;
    double z = fibers[i].z - zBar;
    double A = fibers[i].wgt;
    
    dsigdh = fibers[i].material->getStressSensitivity(gradIndex,true);

    ds[0] += dsigdh(0)*A;
    ds[1] += -y*dsigdh(0)*A;
    ds[2] +=  z*dsigdh(0)*A;
    ds[3] += rootAlpha*dsigdh(1)*A;
    ds[4] += rootAlpha*dsigdh(2)*A;
    ds[5] += (-z*dsigdh(1)+y*dsigdh(2))*A;

    if (areaDeriv[i] != 0.0 || dydh[i] != 0.0 ||  dzdh[i] != 0.0 || parameterID == 1)
      stress = fibers[i].material->getStress();

    if (dydh[i] != 0.0 || dzdh[i] != 0.0 || parameterID == 1)
      tangent = fibers[i].material->getTangent();

    if (areaDeriv[i] != 0.0) {
      sig_dAdh(0) = stress(0)*areaDeriv[i];
      sig_dAdh(1) = stress(1)*areaDeriv[i];
      sig_dAdh(2) = stress(2)*areaDeriv[i];
      
      ds[0] += sig_dAdh(0);
      ds[1] += -y*sig_dAdh(0);
      ds[2] +=  z*sig_dAdh(0);
      ds[3] += rootAlpha*sig_dAdh(1);
      ds[4] += rootAlpha*sig_dAdh(2);
      ds[5] += -z*sig_dAdh(1)+y*sig_dAdh(2);
    }

    if (dydh[i] != 0.0) {
      ds(1) += -dydh[i] * (stress(0)*A);
      ds(5) +=  dydh[i] * (stress(2)*A);
    }

    if (dzdh[i] != 0.0) {
      ds[2] +=  dzdh[i] * (stress(0)*A);
      ds[5] += -dzdh[i] * (stress(1)*A);
    }

    if (parameterID == 1) {
      ds[3] += drootAlphadh * (stress(1)*A);
      ds[4] += drootAlphadh * (stress(2)*A);
    }

    MatrixND<3,6> as;
    as(0,0) =  1;
    as(0,1) = -y;
    as(0,2) =  z;
    as(1,3) = rootAlpha;
    as(2,4) = rootAlpha;
    as(1,5) = -z;
    as(2,5) =  y;
    
    static Matrix dasdh(3,6);
    dasdh(0,1) = -dydh[i];
    dasdh(0,2) = dzdh[i];
    dasdh(1,3) = drootAlphadh;
    dasdh(2,4) = drootAlphadh;
    dasdh(1,5) = -dzdh[i];
    dasdh(2,5) = dydh[i];
    
    static Matrix tmpMatrix(6,6);
    tmpMatrix.addMatrixTripleProduct(0.0, as, tangent, dasdh, 1.0);
    
    ds.addMatrixVector(1.0, tmpMatrix, e, A);
  }

  return ds;
}

const Matrix &
FrameSolidSection3d::getInitialTangentSensitivity(int gradIndex)
{
  static Matrix dksdh(6,6);
  
  dksdh.Zero();
  /*
  double y, A, dydh, dAdh, tangent, dtangentdh;

  static double fiberLocs[10000];
  static double fiberArea[10000];

  if (sectionIntegr != 0) {
    sectionIntegr->getFiberLocations(fibers.size(), fiberLocs);
    sectionIntegr->getFiberWeights(fibers.size(), fiberArea);
  }  
  else {
    for (int i = 0; i < fibers.size(); i++) {
      fiberLocs[i] = matData[2*i];
      fibers[i].wgt = matData[2*i+1];
    }
  }

  static double locsDeriv[10000];
  static double areaDeriv[10000];

  if (sectionIntegr != 0) {
    sectionIntegr->getLocationsDeriv(fibers.size(), locsDeriv);  
    sectionIntegr->getWeightsDeriv(fibers.size(), areaDeriv);
  }
  else {
    for (int i = 0; i < fibers.size(); i++) {
      locsDeriv[i] = 0.0;
      areaDeriv[i] = 0.0;
    }
  }
  
  for (int i = 0; i < fibers.size(); i++) {
    y = fiberLocs[i] - yBar;
    A = fibers[i].wgt;
    dydh = locsDeriv[i];
    dAdh = areaDeriv[i];
    
    tangent = fibers[i].material->getInitialTangent();
    dtangentdh = fibers[i].material->getInitialTangentSensitivity(gradIndex);

    dksdh(0,0) += dtangentdh*A + tangent*dAdh;

    dksdh(0,1) += -y*(dtangentdh*A+tangent*dAdh) - dydh*(tangent*A);

    dksdh(1,1) += 2*(y*dydh*tangent*A) + y*y*(dtangentdh*A+tangent*dAdh);
  }

  dksdh(1,0) = dksdh(0,1);
  */
  return dksdh;
}

int
FrameSolidSection3d::commitSensitivity(const Vector& defSens,
                                       int gradIndex, int numGrads)
{
  double d0 = defSens(0);
  double d1 = defSens(1);
  double d2 = defSens(2);
  double d3 = defSens(3);
  double d4 = defSens(4);
  double d5 = defSens(5);

  dedh = defSens;

  static double dydh[10000];
  static double dzdh[10000];

  { // TODO
    for (int i = 0; i < fibers.size(); i++) {
      dydh[i] = 0.0;
      dzdh[i] = 0.0;
    }
  }

  static Vector depsdh(3);

  double rootAlpha = 1.0;
  if (alpha != 1.0)
    rootAlpha = sqrt(alpha);

  double drootAlphadh = 0.0;
  if (parameterID == 1)
    drootAlphadh = 0.5/rootAlpha;

  for (int i = 0; i < fibers.size(); i++) {
    auto& fiber = fibers[i];
    NDMaterial *theMat = fiber.material;
    const double y  = fiber.y - yBar;
    const double z  = fiber.z - zBar;

    // determine material strain and set it
    depsdh[0] = d0 - y*d1 + z*d2 - dydh[i]*e(1) + dzdh[i]*e(2);
    depsdh[1] = rootAlpha*d3 - z*d5 + drootAlphadh*e(3) - dzdh[i]*e(5);
    depsdh[2] = rootAlpha*d4 + y*d5 + drootAlphadh*e(4) + dydh[i]*e(5);

    theMat->commitSensitivity(depsdh,gradIndex,numGrads);
  }

  return 0;
}
// AddingSensitivity:END ///////////////////////////////////
