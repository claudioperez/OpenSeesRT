//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
// Description: This file contains the class implementation of FrameSolidSection3d.
// FrameSolidSection3d provides the abstraction of a 3D beam section discretized by fibers.
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
#include <Matrix3D.h>
#include <Rotations.hpp>
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
using OpenSees::Matrix3D;

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
  code(iny) = SECTION_RESPONSE_VY;
  code(inz) = SECTION_RESPONSE_VZ;
  code(imx) = SECTION_RESPONSE_T;
  code(imy) = SECTION_RESPONSE_MY;
  code(imz) = SECTION_RESPONSE_MZ;
  code(iwx) = FrameStress::Bimoment;
  code(iwy) = FrameStress::By;
  code(iwz) = FrameStress::Bz;
  code(ivx) = FrameStress::Bishear;
  code(ivy) = FrameStress::Qy;
  code(ivz) = FrameStress::Qz;
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
  code(iny) = SECTION_RESPONSE_VY;
  code(inz) = SECTION_RESPONSE_VZ;
  code(imx) = SECTION_RESPONSE_T;
  code(imy) = SECTION_RESPONSE_MY;
  code(imz) = SECTION_RESPONSE_MZ;
  code(iwx) = FrameStress::Bimoment;
  code(iwy) = FrameStress::By;
  code(iwz) = FrameStress::Bz;
  code(ivx) = FrameStress::Bishear;
  code(ivy) = FrameStress::Qy;
  code(ivz) = FrameStress::Qz;
}

int
FrameSolidSection3d::getIntegral(Field field, State state, double& value) const
{
  value = 0.0;

  const int nf = fibers->size();
  switch (field) {
    case Field::Unit:
      for (int i=0; i<nf; i++) {
        const double A  = (*fibers)[i].area;
        value += A;
      }
      return 0;

    case Field::Density:
      // First check if density has been specified for the section
      if (this->FrameSection::getIntegral(field, state, value) == 0) 
        return 0;

      for (int i=0; i<nf; i++) {
        double density;
        const double A  = (*fibers)[i].area;
        if (materials[i]->getRho() != 0)
          value += A*density;
        else
          return -1;
      }
      return 0;

    case Field::UnitY: // TODO: Centroid
      for (int i=0; i<nf; i++) {
        const double A  = (*fibers)[i].area;
        const double y  = (*fibers)[i].y - yBar;
        value += A*y;
      }
      return 0;


    case Field::UnitYY:
    case Field::UnitCentroidYY:
      for (int i=0; i<nf; i++) {
        const double A  = (*fibers)[i].area;
        const double y  = (*fibers)[i].y
                        - yBar*(Field::UnitCentroidYY == field);
        value += A*y*y;
      }
      return 0;

    case Field::UnitZZ:
    case Field::UnitCentroidZZ:
      for (int i=0; i<nf; i++) {
        const double A  = (*fibers)[i].area;
        const double z  = (*fibers)[i].z
                        - zBar*(Field::UnitCentroidZZ == field);
        value += A*z*z;
      }

    default:
      return -1;
  }
  return -1;
}


int
FrameSolidSection3d::addFiber(NDMaterial& theMat, 
                              double Area, 
                              double yLoc, 
                              double zLoc)
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
  
  return materials.size()-1;
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
FrameSolidSection3d::stateDetermination(Matrix& ksi_, Vector* s_trial, const Vector * const e_trial, int tangentFlag)
{

  const Vector3D kappa {
    e_trial? (*e_trial)(imx) : 0.0,
    e_trial? (*e_trial)(imy) : 0.0,
    e_trial? (*e_trial)(imz) : 0.0 
  };
  const Vector3D gamma {
    e_trial? (*e_trial)(inx) : 0.0,
    e_trial? (*e_trial)(iny) : 0.0,
    e_trial? (*e_trial)(inz) : 0.0 
  };

  Kmm.zero();
  Kmn.zero();
  Kmw.zero();
  Knn.zero();
  Knw.zero();
  Kvv.zero();
  Knv.zero();
  Kww.zero();
  Kmv.zero();


  int res = 0;
  const int nf = fibers->size();
  for (int i = 0; i < nf; i++) {

    NDMaterial &theMat = *materials[i];
    auto & fiber = (*fibers)[i];
    auto & w = fiber.warp;

    if (e_trial != nullptr) {
      // Form material strain
      Vector3D eps = gamma + kappa.cross(fiber.r);
        for (int j=1; j<3; j++)
          for (int k=0; k<nwm; k++)
            eps[j] += w[k][j];
      res += theMat.setTrialStrain(eps);
    }

    const Matrix &tangent = tangentFlag==CurrentTangent
                            ? theMat.getTangent()
                            : theMat.getInitialTangent();


    Matrix3D C{};
    C.addMatrix(tangent, fiber.area);

    // NOTE: Matrix 3D is column major so these are transposed.
    const Matrix3D iow{{
      {w[0][0],     0.0,     0.0},
      {w[1][0],     0.0,     0.0},
      {w[2][0],     0.0,     0.0}
    }};

    const Matrix3D iodw{{
      {    0.0, w[0][1], w[0][2]},
      {    0.0, w[1][1], w[1][2]},
      {    0.0, w[2][1], w[2][2]}
    }};

    Knn.addMatrix(C, 1.0);
    {
      Matrix3D rxC{};
      rxC.addSpinMatrixProduct(fiber.r, C, 1.0);
      Kmm.addMatrixSpinProduct(rxC, fiber.r, -1.0);
      Kmn.addMatrix(rxC, 1.0);
      Kmw.addMatrixProduct(rxC, iow,  1.0);
      Kmv.addMatrixProduct(rxC, iodw, 1.0);
    }
    {
      Matrix3D Ciow{};
      Ciow.addMatrixProduct(C, iow, 1.0);
      Knw.addMatrix(Ciow,  1.0);
      Kww.addMatrixTransposeProduct(1.0, iow,  Ciow, 1.0);
    }
    {
      Matrix3D Ciodw{};
      Ciodw.addMatrixProduct(C, iodw, 1.0);
      Knv.addMatrix(Ciodw, 1.0);
      Kvv.addMatrixTransposeProduct(1.0, iodw,  Ciodw, 1.0);
    }

    if (s_trial != nullptr) {
      const double y = fiber.r[1];
      const double z = fiber.r[2];
      const Vector &stress  = theMat.getStress();
      double sig0 = stress(0)*fiber.area;
      double sig1 = stress(1)*fiber.area;
      double sig2 = stress(2)*fiber.area;
      // n += s da
      (*s_trial)(inx) +=    sig0;
      (*s_trial)(iny) +=    sig1;
      (*s_trial)(inz) +=    sig2;
      // m += r.cross(s) da
      (*s_trial)(imx) += -z*sig1 + y*sig2;
      (*s_trial)(imy) +=  z*sig0;
      (*s_trial)(imz) += -y*sig0;
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
  ks->Zero();
  ks->Assemble(Knn, 0, 0, 1.0);
  ks->Assemble(Knw, 0, 6, 1.0);
  ks->Assemble(Knv, 0, 9, 1.0);
  ks->Assemble(Kmn, 3, 0, 1.0);
  ks->Assemble(Kmm, 3, 3, 1.0);
  ks->Assemble(Kmw, 3, 6, 1.0);
  ks->Assemble(Kmv, 3, 9, 1.0);
  ks->Assemble(Kww, 6, 6, 1.0);
  ks->Assemble(Kvv, 9, 9, 1.0);

  ks->AssembleTranspose(Knw, 6, 0, 1.0);
  ks->AssembleTranspose(Knv, 9, 0, 1.0);
  ks->AssembleTranspose(Kmn, 0, 3, 1.0);
  ks->AssembleTranspose(Kmw, 6, 3, 1.0);
  ks->AssembleTranspose(Kmv, 9, 3, 1.0);
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
    theCopy->materials.push_back(material->getCopy()); // "BeamFiber"

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
FrameSolidSection3d::getOrder() const
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
  for (auto& material : materials)
    err += material->revertToLastCommit();


  ks->Zero();
  s->Zero();
  err += this->stateDetermination(*ks, s, nullptr, CurrentTangent);

  return err;
}

int
FrameSolidSection3d::revertToStart()
{
  int err = 0;
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
  return 0;
}

int
FrameSolidSection3d::recvSelf(int , Channel &,
                              FEM_ObjectBroker &)
{
  return 0;
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
      s << OPS_PRINT_JSON_MATE_INDENT << "\t{\"location\": [" 
        << (*fibers)[i].r[1] << ", " 
        << (*fibers)[i].r[2] << "], ";
      s << "\"area\": " << (*fibers)[i].area << ", ";
      s << "\"warp\": [";
      for (int j = 0; j < nwm; j++) {
        s << "[";
        for (int k = 0; k < 3; k++) {
          s << (*fibers)[i].warp[j][k];
          if (k < 2)
            s << ", ";
        }
        s << "]";
        if (j < nwm-1)
          s << ", ";
      }
      s << "], ";

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
      s << "\nLocation (y,z) = " << fiber.r[1] << ' ' << fiber.r[2];
      s << "\nArea = " << fiber.area << endln;
      materials[i]->Print(s, flag);
    }
  }
}

Response*
FrameSolidSection3d::setResponse(const char **argv, int argc,
                              OPS_Stream &output)
{
  Response *theResponse = nullptr;

  if (argc > 2 && strcmp(argv[0], "fiber") == 0) {

    
    int key = fibers->size();
    int passarg = 2;
    
    if (argc <= 3) {      // fiber number was input directly
      key = atoi(argv[1]);
    }

    else if (argc > 4) {  // find fiber closest to coord. with mat tag
      
      int matTag    = atoi(argv[3]);
      double yCoord = atof(argv[1]);
      double zCoord = atof(argv[2]);
      Vector3D r_given{{0, yCoord, zCoord}};
      double closestDist = 0;

      // Find first fiber with specified material tag
      const int nf = fibers->size();
      int j;
      for (j = 0; j < nf; j++) {
        auto& fiber = (*fibers)[j];
        if (matTag == materials[j]->getTag()) {
          Vector3D dr = fiber.r - r_given;
          closestDist = dr.dot(dr);
          key = j;
          break;
        }
      }

      // Search the remaining fibers
      double distance;
      for ( ; j < nf; j++) {
        auto& fiber = (*fibers)[j];
        if (matTag == materials[j]->getTag()) {
          Vector3D dr = fiber.r - r_given;
          distance = dr.dot(dr);
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
      Vector3D r_given{
        0.0, atof(argv[1]), atof(argv[2])
      };
      Vector3D dr = (*fibers)[0].r - r_given;
      double closestDist = dr.dot(dr);
      key = 0;
      double distance;

      const int nf = fibers->size();
      for (int j = 1; j < nf; j++) {
        auto& fiber = (*fibers)[j];
        Vector3D dr = fiber.r - r_given;
        distance = dr.dot(dr);
        if (distance < closestDist) {
          closestDist = distance;
          key = j;
        }
      }
      passarg = 3;
    }

    if (key < fibers->size() && key >= 0) {
      output.tag("FiberOutput");
      output.attr("y",    (*fibers)[key].r[1]);
      output.attr("z",    (*fibers)[key].r[2]);
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
    // ... warp $fiberID $warpField
    if (argc < 3) {
      opserr << "FrameSolidSection3d::setParameter - fiberID is required\n";
      return -1;
    }
    int fiberID = atoi(argv[1]);
    if (fiberID < 0 || fiberID >= fibers->size()) {
      opserr << "FrameSolidSection3d::setParameter - fiberID " << fiberID << " out of range\n";
      return -1;
    }

    int field = atoi(argv[2]);

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
        (*fibers)[fiberID].r[1] = info.theDouble;
        break;
      case Param::FiberZ:
        (*fibers)[fiberID].r[2] = info.theDouble;
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
    const double y = (*fibers)[i].r[1] - yBar;
    const double z = (*fibers)[i].r[2] - zBar;
    const double A = (*fibers)[i].area;
    
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
    const double y  = fiber.r[1] - yBar;
    const double z  = fiber.r[2] - zBar;

    // determine material strain sensitivity
    depsdh[0] = d0 - y*d1 + z*d2 - dydh[i]*e(1) + dzdh[i]*e(2);
    depsdh[1] = rootAlpha*d3 - z*d5 + drootAlphadh*e(3) - dzdh[i]*e(5);
    depsdh[2] = rootAlpha*d4 + y*d5 + drootAlphadh*e(4) + dydh[i]*e(5);

    materials[i]->commitSensitivity(depsdh,gradIndex,numGrads);
  }

  return 0;
}
