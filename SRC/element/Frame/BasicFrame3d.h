#pragma once

#include <Frame/FiniteElement.h>
#include <CrdTransf.h>

#include <Matrix.h>
#include <MatrixND.h>
#include <Vector.h>
#include <VectorND.h>
#include <ElementalLoad.h>

using namespace OpenSees;


class BasicFrame3d : public FiniteElement<2, 3, 6> {
  constexpr static int ndm = 3;

  public:
    virtual ~BasicFrame3d();
    BasicFrame3d(int tag, int clstag, std::array<int, 2> &nodes, CrdTransf& tran,
                 int cMass_, int rz, int ry)
      : FiniteElement<2, 3, 6> (tag, clstag, nodes), theCoordTransf(tran.getCopy3d()),
        p_iner(12),
        releasez(rz), releasey(ry),
        cMass(cMass_), mass(0.0),
        wx(0.0), wy(0.0), wz(0.0),
        numEleLoads(0), sizeEleLoads(0), eleLoads(0), eleLoadFactors(0),
        parameterID(0)
    {
        zeroLoad();

        // Make no release if input not 0, 1, 2, or 3
        if (releasez < 0 || releasez > 3)
          releasez = 0;
        if (releasey < 0 || releasey > 3)
          releasey = 0;  

        if (theCoordTransf == nullptr) {
          opserr << "PrismFrame3d::PrismFrame3d -- failed to get copy of coordinate transformation\n";
        }
    }

    BasicFrame3d(int tag, int classtag)
      : FiniteElement<2, 3, 6> (tag, classtag), theCoordTransf(nullptr), p_iner(12), mass(0.0),
        releasez(0), releasey(0),
        wx(0.0), wy(0.0), wz(0.0),
        numEleLoads(0), sizeEleLoads(0), eleLoads(0), eleLoadFactors(0),
        parameterID(0)
    {
        zeroLoad();
    }

    //
    // For FiniteElement
    //
    virtual double          getTotalMass() final;
    virtual int             setNodes();
    virtual VectorND<12>    getForce(State state, int rate) override final;
    virtual MatrixND<12,12> getTangent(State state, int rate) override final {
      MatrixND<12,12> K;
      return K;
    }

    //
    // For Element
    //
    virtual void zeroLoad() final;
    virtual int addLoad(ElementalLoad *theLoad, double loadFactor) final;
    virtual int         addInertiaLoadToUnbalance(const Vector &accel) final;
    virtual const Vector &getResistingForceIncInertia() final;
    virtual const Matrix &getInitialStiff() final;
    virtual const Matrix &getMass() final;
    // Sensitivity
    const Matrix & getMassSensitivity(int gradNumber);
    virtual int setParameter(const char **argv, int argc, Parameter &param);
    virtual int            updateParameter(int parameterID, Information &info);
    virtual int            activateParameter(int parameterID);


    const Matrix &
    getTangentStiff() final
    {
      VectorND<6>   q  = this->getBasicForce();
      MatrixND<6,6> km = this->getBasicTangent(State::Pres, 0);
      
//    q += q0; // TODO!!! move this into PrismFrame and maybe DisplFrame

      // TODO
      const Matrix ktemp(km);
      const Vector qtemp(q);
      return theCoordTransf->getGlobalStiffMatrix(ktemp,qtemp);
    }


    const Vector &
    getResistingForce()
    {
      VectorND<6> q  = this->getBasicForce();

//    q += q0; // TODO!!! move this into PrismFrame and maybe DisplFrame

      const Vector p0Vec(p0);
      P = theCoordTransf->getGlobalResistingForce(q, p0Vec);

      // Subtract other external nodal loads ... P_res = P_int - P_ext
      if (mass != 0.0)
        P.addVector(1.0, p_iner, -1.0);

      return P;
    }


protected:

  // Implemented by children
  virtual VectorND<6>&   getBasicForce() = 0;
  virtual MatrixND<6,6>& getBasicTangent(State state, int rate) = 0;

  // Supplied to children
          double         getLength(State flag);
  const   VectorND<6>&   getBasicState(State flag);
  // Reactions of basic system due to element loads
  void computeReactionSensitivity(double *dp0dh, int gradNumber);
  void computeReactions(double *p0);

// to be made private
   int cMass;
   CrdTransf* theCoordTransf;
   OpenSees::VectorND<6>   q0;  // Fixed end forces in basic system (with torsion)
   OpenSees::VectorND<6>   p0;  // Reactions in basic system (with torsion)
                                // TODO(cmp): change to size 12
   static Matrix K;
   static Vector P;

   double rho;

   int parameterID;

   int releasez; // moment release for bending about z-axis 0=none, 1=I, 2=J, 3=I,J
   int releasey; // same for y-axis

   double wx;
   double wy;
   double wz;
   int numEleLoads;               // Number of element load objects
   int sizeEleLoads;
   ElementalLoad **eleLoads;
   double *eleLoadFactors;

  private:
   Vector p_iner;
   double mass;

};

