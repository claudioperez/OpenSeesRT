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
// Written: cmp
// Created: Spring 2025
//
#include <stdlib.h>
#include <string.h>
#include <cmath>
#include <array>
#include <Channel.h>
#include <Vector.h>
#include <VectorND.h>
#include <Matrix.h>
#include <MatrixND.h>
#include <classTags.h>
#include "FrameSolidSection3d.h"
#include <ID.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <SensitiveResponse.h>
typedef SensitiveResponse<FrameSection> SectionResponse;
#include <NDMaterial.h>
#include <Parameter.h>

#define SEC_TAG_FrameSolidSection3d 0

using OpenSees::VectorND;
using OpenSees::MatrixND;

ID FrameSolidSection3d::code(nsr);

FrameSolidSection3d::FrameSolidSection3d(int tag, int num, double a, bool compCentroid): 
    FrameSection(tag, SEC_TAG_FrameSolidSection3d),
    Abar(0.0), 
    QyBar(0.0), QzBar(0.0), 
    yBar(0.0), zBar(0.0), 
    computeCentroid(compCentroid),
    alpha(a),
    e(nsr), s(nullptr), ks(nullptr), 
    parameterID(0), dedh(nsr),
    fibers(new std::vector<FiberData>)
{
    s = new Vector(sData, nsr);
    ks = new Matrix(kData, nsr, nsr);

    code(inx) = SECTION_RESPONSE_P;
    code(imz) = SECTION_RESPONSE_MZ;
    code(imy) = SECTION_RESPONSE_MY;
    code(iny) = SECTION_RESPONSE_VY;
    code(inz) = SECTION_RESPONSE_VZ;
    code(imx) = SECTION_RESPONSE_T;
}

// constructor for blank object that recvSelf needs to be invoked upon
FrameSolidSection3d::FrameSolidSection3d():
  FrameSection(0, SEC_TAG_FrameSolidSection3d),
  Abar(0.0), QyBar(0.0), QzBar(0.0), 
  yBar(0.0), zBar(0.0), 
  computeCentroid(true),
  alpha(1.0), 
  e(nsr), s(nullptr), ks(nullptr),
  parameterID(0), dedh(nsr),
  fibers(new std::vector<FiberData>)
{
  s = new Vector(sData, nsr);
  ks = new Matrix(kData, nsr, nsr);

  code(inx) = SECTION_RESPONSE_P;
  code(imz) = SECTION_RESPONSE_MZ;
  code(imy) = SECTION_RESPONSE_MY;
  code(iny) = SECTION_RESPONSE_VY;
  code(inz) = SECTION_RESPONSE_VZ;
  code(imx) = SECTION_RESPONSE_T;
}

int
FrameSolidSection3d::addFiber(NDMaterial& theMat, 
                              double Area, 
                              double yLoc, 
                              double zLoc
                              )
{
  std::array<std::array<double,3>,3> warp{0};
  FiberData fiber {yLoc, zLoc, Area, warp, {0.0, yLoc, zLoc}};
  fibers->emplace_back(fiber);

  materials.emplace_back(theMat.getCopy("BeamFiber"));

  if (materials[materials.size()-1] == nullptr) {
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
  for (auto material : materials) {
    if (material != nullptr)
      delete material;
  }

  if (s != nullptr)
    delete s;

  if (ks != nullptr)
    delete ks;
}

int
FrameSolidSection3d::setTrialSectionDeformation(const Vector &e_trial)
{
  e = e_trial;
  s->Zero();
  ks->Zero();
  return stateDetermination(*ks, s, &e, CurrentTangent);
}

int
FrameSolidSection3d::stateDetermination(Matrix& ksi, Vector* s_trial, const Vector *e_trial, int tangentFlag)
{
  const double 
         // Position
         e0 = e_trial? (*e_trial)(inx) : 0.0, // N
         e3 = e_trial? (*e_trial)(iny) : 0.0, // Vy
         e4 = e_trial? (*e_trial)(inz) : 0.0, // Vz
         // Curvature
         kx = e_trial? (*e_trial)(imx) : 0.0, // T
         kz = e_trial? (*e_trial)(imz) : 0.0, // Mz
         ky = e_trial? (*e_trial)(imy) : 0.0; // My


  const double rootAlpha = std::sqrt(alpha);

  int res = 0;
  const int nf = fibers->size();
  for (int i = 0; i < nf; i++) {

    NDMaterial &theMat = *materials[i];
    auto & fiber = (*fibers)[i];
    const double y  = fiber.y - yBar;
    const double z  = fiber.z - zBar;
    const double A  = fiber.area;

    if (e_trial != nullptr) {
      // Form material strain
      //
      //      0  5 4       1       2  3
      // a = [1 -y z       0       0  0
      //      0  0 0 sqrt(a)       0 -z
      //      0  0 0       0 sqrt(a)  y]
      VectorND<3> eps;
      eps[0] = e0 - y*kz + z*ky;
      eps[1] = rootAlpha*e3 - z*kx;
      eps[2] = rootAlpha*e4 + y*kx;
      res += theMat.setTrialStrain(eps);
    }

    const Matrix &tangent = tangentFlag==CurrentTangent? 
                            theMat.getTangent()
                            : theMat.getInitialTangent();

    double d00 = tangent(0,0)*A;
    double d01 = tangent(0,1)*A;
    double d02 = tangent(0,2)*A;
    double d10 = tangent(1,0)*A;
    double d11 = tangent(1,1)*A;
    double d12 = tangent(1,2)*A;
    double d20 = tangent(2,0)*A;
    double d21 = tangent(2,1)*A;
    double d22 = tangent(2,2)*A;

    double y2 = y*y;
    double z2 = z*z;
    double yz = y*z;
    double tmp;
    // Bending terms
    ksi(inx,inx) +=    d00;
    ksi(imz,imz) += y2*d00;
    ksi(imy,imy) += z2*d00;
    tmp = -y*d00;
    ksi(inx,imz) += tmp;
    ksi(imz,inx) += tmp;
    tmp = z*d00;
    ksi(inx,imy) += tmp;
    ksi(imy,inx) += tmp;
    tmp = -yz*d00;
    ksi(imz,imy) += tmp;
    ksi(imy,imz) += tmp;
    
    // Shear terms
    ksi(iny,iny) += alpha*d11;
    ksi(iny,inz) += alpha*d12;
    ksi(inz,iny) += alpha*d21;
    ksi(inz,inz) += alpha*d22;
    
    // Torsion term
    ksi(imx,imx) += z2*d11 - yz*(d12+d21) + y2*d22;
    
    // Bending-torsion coupling terms
    tmp = -z*d01 + y*d02;
    ksi(inx,imx) += tmp;
    ksi(imz,imx) -= y*tmp;
    ksi(imy,imx) += z*tmp;
    tmp = -z*d10 + y*d20;
    ksi(imx,inx) += tmp;
    ksi(imx,imz) -= y*tmp;
    ksi(imx,imy) += z*tmp;
    
    // Scale tangent terms with rootAlpha
    d01 *= rootAlpha; d02 *= rootAlpha;
    d10 *= rootAlpha; d11 *= rootAlpha; d12 *= rootAlpha;
    d20 *= rootAlpha; d21 *= rootAlpha; d22 *= rootAlpha;
    
    // Bending-shear coupling terms
    ksi(inx,iny) += d01;
    ksi(inx,inz) += d02;
    ksi(imz,iny) -= y*d01;
    ksi(imz,inz) -= y*d02;
    ksi(imy,iny) += z*d01;
    ksi(imy,inz) += z*d02;
    ksi(iny,inx) += d10;
    ksi(inz,inx) += d20;
    ksi(iny,imz) -= y*d10;
    ksi(inz,imz) -= y*d20;
    ksi(iny,imy) += z*d10;
    ksi(inz,imy) += z*d20;
    
    // Torsion-shear coupling terms
    y2 =  y*d22;
    z2 = -z*d11;
    ksi(imx,iny) +=  z2 + y*d21;
    ksi(imx,inz) += -z*d12 + y2;
    ksi(iny,imx) +=  z2 + y*d12;
    ksi(inz,imx) += -z*d21 + y2;

    if (s_trial != nullptr) {
      const Vector &stress  = theMat.getStress();
      double sig0 = stress(0)*A;
      double sig1 = stress(1)*A;
      double sig2 = stress(2)*A;
      (*s_trial)(inx) +=    sig0;
      (*s_trial)(imz) += -y*sig0;
      (*s_trial)(imy) +=  z*sig0;
      (*s_trial)(iny) += rootAlpha*sig1;
      (*s_trial)(inz) += rootAlpha*sig2;
      (*s_trial)(imx) += -z*sig1 + y*sig2;
    }
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
  static double kInitial[nsr*nsr];
  static Matrix ksi(kInitial, nsr, nsr);
  ksi.Zero();

  this->stateDetermination(ksi, nullptr, nullptr, InitialTangent);

  return ksi;
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


  for (auto& material: materials)
    theCopy->materials.push_back(material->getCopy("BeamFiber"));

  theCopy->fibers = fibers;
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

  for (auto& material: materials)
    err += material->commitState();

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

  // invoke revertToLast on the material
  for (auto& material : materials)
    err += material->revertToLastCommit();

  err += this->stateDetermination(*ks, s, nullptr, CurrentTangent);

  return err;
}

int
FrameSolidSection3d::revertToStart()
{
  // revert the fibers to start    
  int err = 0;


  // invoke revertToLast on the material
  for (auto& material: materials)
    err += material->revertToStart();

  ks->Zero();
  s->Zero();
  err += stateDetermination(*ks, s, nullptr, CurrentTangent);

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
  const int nf = fibers->size();
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
      s << OPS_PRINT_JSON_MATE_INDENT << "{";
      s << "\"name\": \"" << this->getTag() << "\", ";
      s << "\"type\": \"" << this->getClassType() << "\", ";

      double mass;
      if (this->FrameSection::getIntegral(Field::Density, State::Init, mass) == 0)
        s << "\"mass\": " << mass;

      s << "\"fibers\": [\n";

      for (int i = 0; i < nf; i++) {
            s << OPS_PRINT_JSON_MATE_INDENT << "\t{\"coord\": [" << (*fibers)[i].y << ", " << (*fibers)[i].z << "], ";
            s << "\"area\": " << (*fibers)[i].area << ", ";
            s << "\"material\": " << materials[i]->getTag();
            if (i < nf - 1)
                s << "},\n";
            else
                s << "}\n";
      }
      s << OPS_PRINT_JSON_MATE_INDENT << "]}";
      return;
  } 
  else if (flag == 1) {
    for (int i = 0; i < nf; i++) {
      auto & fiber = (*fibers)[i];
      s << "\nLocation (y,z) = " << fiber.y << ' ' << fiber.z;
      s << "\nArea = " << fiber.area << endln;
      materials[i]->Print(s, flag);
    }
  } else {
    s << "\nFrameSolidSection3d, tag: " << this->getTag() << endln;
    s << "\tSection code: " << code;
    s << "\tNumber of Fibers: " << nf << endln;
    s << "\tCentroid (y,z): " << yBar << ' ' << zBar << endln;
    s << "\tShape factor, alpha = " << alpha << endln;
  }
}

Response*
FrameSolidSection3d::setResponse(const char **argv, int argc,
                              OPS_Stream &output)
{
  Response *theResponse = nullptr;

  if (argc > 2 && strcmp(argv[0],"fiber") == 0) {

    
    int key = fibers->size();
    int passarg = 2;
    
    if (argc <= 3) {                  // fiber number was input directly
      
      key = atoi(argv[1]);
      
    } else if (argc > 4) {  // find fiber closest to coord. with mat tag
      
      int matTag = atoi(argv[3]);
      double yCoord = atof(argv[1]);
      double zCoord = atof(argv[2]);
      double closestDist = 0;
      double y_search, z_search, dy, dz;
      double distance;
      int j;
      // Find first fiber with specified material tag

      const int nf = fibers->size();
      for (j = 0; j < nf; j++) {
        auto& fiber = (*fibers)[j];
        if (matTag == materials[j]->getTag()) {
          y_search = fiber.y;
          z_search = fiber.z;

          dy = y_search-yCoord;
          dz = z_search-zCoord;
          closestDist = dy*dy + dz*dz;
          key = j;
          break;
        }
      }

      // Search the remaining fibers
      for ( ; j < nf; j++) {
        auto& fiber = (*fibers)[j];
        if (matTag == materials[j]->getTag()) {
          y_search = fiber.y;
          z_search = fiber.z;

          double dy = y_search - yCoord;
          double dz = z_search - zCoord;
          distance = dy*dy + dz*dz;
          if (distance < closestDist) {
            closestDist = distance;
            key = j;
          }
        }
      }
      passarg = 4;
    }
    
    else {
      // fiber near-to coordinate specified
      
      double yCoord = atof(argv[1]);
      double zCoord = atof(argv[2]);
      double closestDist;
      double distance;
      
      double y_search = (*fibers)[0].y;
      double z_search = (*fibers)[0].z;

      double dy = y_search-yCoord;
      double dz = z_search-zCoord;
      closestDist = dy*dy + dz*dz;
      key = 0;

      const int nf = fibers->size();
      for (int j = 1; j < nf; j++) {
        auto& fiber = (*fibers)[j];
        double dy = fiber.y-yCoord;
        double dz = fiber.z-zCoord;
        distance = dy*dy + dz*dz;
        if (distance < closestDist) {
          closestDist = distance;
          key = j;
        }
      }
      passarg = 3;
    }

    if (key < fibers->size() && key >= 0) {
      output.tag("FiberOutput");
      output.attr("y",    (*fibers)[key].y);
      output.attr("z",    (*fibers)[key].z);
      output.attr("area", (*fibers)[key].area);
      
      theResponse = materials[key]->setResponse(&argv[passarg], argc-passarg, output);
      
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
  return FrameSection::getResponse(responseID, sectInfo);
}

int
FrameSolidSection3d::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;

  int result = -1;

  if (strcmp(argv[0],"alpha") == 0)
    return param.addObject(Param::alpha, this);


  if (strcmp(argv[0], "warp") == 0) {
    if (argc < 2) {
      opserr << "FrameSolidSection3d::setParameter - fiberID required\n";
      return -1;
    }
    int fiberID = atoi(argv[1]);
    if (fiberID < 0 || fiberID >= fibers->size()) {
      opserr << "FrameSolidSection3d::setParameter - fiberID out of range\n";
      return -1;
    }

    int field = 0;
    if (argc > 4) {
      field = atoi(argv[2]);
    }

    return param.addObject(Param::FiberFieldBase+fiberID*100+field, this);
  }

  // Check if the parameter belongs to the material (only option for now)
  if (strstr(argv[0],"material") != 0) {
    
    if (argc < 3)
      return 0;

    // Get the tag of the material
    int materialTag = atoi(argv[1]);
    
    // Loop over fibers to find the right material
    for (auto& material: materials)
      if (materialTag == material->getTag()) {
        int ok = material->setParameter(&argv[2], argc-2, param);
        if (ok != -1)
          result = ok;
      }
    return result;
  }

  int ok = 0; 
  for (auto& material: materials) {
    ok = material->setParameter(argv, argc, param);
    if (ok != -1)
      result = ok;
  }
  return result;
}

int
FrameSolidSection3d::updateParameter(int paramID, Information &info)
{

  if (paramID == Param::alpha){
    alpha = info.theDouble;
    return 0;
  }

  if (paramID >= Param::FiberFieldBase) {
    int fiberID = (paramID - Param::FiberFieldBase) / 100;
    int field   = (paramID - Param::FiberFieldBase) % 100;
  
    if (fiberID >= fibers->size()) 
      return -1;

    switch (field) {
      case Param::FiberArea:
        (*fibers)[fiberID].area = info.theDouble;
        break;
      case Param::FiberY:
        (*fibers)[fiberID].y = info.theDouble;
        break;
      case Param::FiberZ:
        (*fibers)[fiberID].z = info.theDouble;
        break;
      case Param::FiberWarpX:
        (*fibers)[fiberID].warp[0][0] = info.theDouble;
        break;
      case Param::FiberWarpXY:
        (*fibers)[fiberID].warp[0][1] = info.theDouble;
        break;
      case Param::FiberWarpXZ:
        (*fibers)[fiberID].warp[0][2] = info.theDouble;
        break;
      //
      case Param::FiberWarpY:
        (*fibers)[fiberID].warp[1][0] = info.theDouble;
        break;
      case Param::FiberWarpYY:
        (*fibers)[fiberID].warp[1][1] = info.theDouble;
        break;
      case Param::FiberWarpYZ:
        (*fibers)[fiberID].warp[1][2] = info.theDouble;
        break;

      case Param::FiberWarpZ:
        (*fibers)[fiberID].warp[2][0] = info.theDouble;
        break;
      case Param::FiberWarpZY:
        (*fibers)[fiberID].warp[2][1] = info.theDouble;
        break;
      case Param::FiberWarpZZ:
        (*fibers)[fiberID].warp[2][2] = info.theDouble;
        break;
      default:
        return -1;
    }
    return 0;
  }


  return -1;
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
  static Vector ds(nsr);
  
  ds.Zero();
  
  static Vector stress(3);
  static Vector dsigdh(3);
  static Vector sig_dAdh(3);
  static Matrix tangent(3,3);

  static double dydh[10000];
  static double dzdh[10000];
  static double areaDeriv[10000];
  const int nf = fibers->size();
  {
    for (int i = 0; i < nf; i++) {
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

  for (int i = 0; i < nf; i++) {
    double y = (*fibers)[i].y - yBar;
    double z = (*fibers)[i].z - zBar;
    double A = (*fibers)[i].area;
    
    dsigdh = materials[i]->getStressSensitivity(gradIndex,true);

    ds[0] += dsigdh(0)*A;
    ds[1] += -y*dsigdh(0)*A;
    ds[2] +=  z*dsigdh(0)*A;
    ds[3] += rootAlpha*dsigdh(1)*A;
    ds[4] += rootAlpha*dsigdh(2)*A;
    ds[5] += (-z*dsigdh(1)+y*dsigdh(2))*A;

    if (areaDeriv[i] != 0.0 || dydh[i] != 0.0 ||  dzdh[i] != 0.0 || parameterID == 1)
      stress = materials[i]->getStress();

    if (dydh[i] != 0.0 || dzdh[i] != 0.0 || parameterID == 1)
      tangent = materials[i]->getTangent();

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
      (*fibers)[i].area = matData[2*i+1];
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
    A = (*fibers)[i].area;
    dydh = locsDeriv[i];
    dAdh = areaDeriv[i];
    
    tangent = materials[i]->getInitialTangent();
    dtangentdh = materials[i]->getInitialTangentSensitivity(gradIndex);

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
  const int nf = fibers->size();
  
  { // TODO
    for (int i = 0; i < nf; i++) {
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

  for (int i = 0; i < nf; i++) {
    auto& fiber = (*fibers)[i];
    const double y  = fiber.y - yBar;
    const double z  = fiber.z - zBar;

    // determine material strain and set it
    depsdh[0] = d0 - y*d1 + z*d2 - dydh[i]*e(1) + dzdh[i]*e(2);
    depsdh[1] = rootAlpha*d3 - z*d5 + drootAlphadh*e(3) - dzdh[i]*e(5);
    depsdh[2] = rootAlpha*d4 + y*d5 + drootAlphadh*e(4) + dydh[i]*e(5);

    materials[i]->commitSensitivity(depsdh,gradIndex,numGrads);
  }

  return 0;
}
