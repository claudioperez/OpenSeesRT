//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
// Claudio M. Perez
//
#pragma once
#include <vector>
#include <array>
#include <Flag.h>
#include <Element.h>
#include <VectorND.h>
#include <MatrixND.h>

#include <Frame/FrameMass.hpp>
#include <FrameSection.h>
#include <FrameTransform.h>
#include <Frame/FiniteElement.h>
#include <Rotations.hpp>
//#include <BeamIntegration.h>


namespace OpenSees {

template<int nen, int nip>
class ExactFrame3d: 
  public FiniteElement<nen, 3, 6>
{
public:
  enum Logarithm {
    None,
    LogLeft,
    LogRight
  };

  ExactFrame3d(int tag, std::array<int,nen>& nodes,
               FrameSection *section[nip], 
               FrameTransform3d& transf
//             int mass_type
  );

  ~ExactFrame3d();

  const char*
  getClassType() const final
  {
    return "ExactFrame3d";
  }

  virtual int setNodes();
  virtual const Vector &getResistingForce();
  virtual const Matrix &getTangentStiff();
  virtual const Matrix &getMass();
  virtual const Matrix &getInitialStiff();

  virtual int update() final;

  virtual int revertToStart();
  virtual int revertToLastCommit();
  virtual int commitState() final {

    if (this->Element::commitState() != 0)
      opserr << "ExactFrame3d::commitState () - failed in base class";

    past = points;

    for (GaussPoint& point : points) {
      if (point.material->commitState() != 0)
        return -1;
    }
    return OpenSees::Flag::Success;
  }

  // MovableObject
  int sendSelf(int cTag, Channel&);
  int recvSelf(int cTag, Channel&, FEM_ObjectBroker&);

  // TaggedObject
  void Print(OPS_Stream& s, int flag = 0);

  private:
    //
    // Constexpr
    //
    constexpr static int 
          nsr = 6,              // Number of section resultants
          ndm = 3,              // Dimension of the problem (3D)
          ndf = 6;              // Degrees of freedom per node


    // Layout of stress resultants
    static constexpr FrameStressLayout scheme = {
      FrameStress::N,
      FrameStress::Vy,
      FrameStress::Vz,
      FrameStress::T,
      FrameStress::My,
      FrameStress::Mz,
    };

    //
    //
    //
    struct GaussPoint {
      double point,
             weight;
      FrameSection* material;

      double   shape[2][nen];
      Matrix3D rotation;
      Vector3D curvature;
    };

    std::array<GaussPoint,nip> points;
    std::array<GaussPoint,nip> past;
    BeamIntegration*        stencil;
    FrameTransform3d*       transform;
    Logarithm               logarithm;


    //
    //
    VectorND<nen*ndf> p;
    MatrixND<nen*ndf,nen*ndf> K;
};

} // namespace OpenSees
#include "ExactFrame3d.tpp"
