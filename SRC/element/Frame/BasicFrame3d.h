//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
#pragma  once
#include <Field.h>
#include <State.h>
#include <Frame/FiniteElement.h>
#include <FrameTransform.h>

#include <Matrix.h>
#include <MatrixND.h>
#include <Vector.h>
#include <VectorND.h>
#include <ElementalLoad.h>
#include <vector>
#include <utility>

using namespace OpenSees;

class BasicFrame3d : public FiniteElement<2, 3, 6> {
  constexpr static int ndm = 3;

  public:
    virtual ~BasicFrame3d();
    BasicFrame3d(int tag, int clstag, std::array<int, 2> &nodes, 
                 FrameTransform3d& tran)
      : FiniteElement<2, 3, 6> (tag, clstag, nodes), theCoordTransf(tran.getCopy()),
        p_iner(12),
        wx(0.0), wy(0.0), wz(0.0),
        // numEleLoads(0), // sizeEleLoads(0), eleLoads(0), eleLoadFactors(0),
        parameterID(0)
    {
        zeroLoad();

        if (theCoordTransf == nullptr) {
          opserr << "PrismFrame3d::PrismFrame3d -- failed to get copy of coordinate transformation\n";
        }
    }

    BasicFrame3d(int tag, int classtag)
      : FiniteElement<2, 3, 6> (tag, classtag), theCoordTransf(nullptr), 
        p_iner(12), 
        wx(0.0), wy(0.0), wz(0.0),
        // numEleLoads(0), // sizeEleLoads(0), eleLoads(0), eleLoadFactors(0),
        parameterID(0)
    {
        zeroLoad();
    }

    // FrameElement
    virtual int getIntegral(Field field, State state, double& total) {
      return -1;
    }

    //
    // For FiniteElement
    //
    virtual int             setNodes();
#ifdef FEFT
    virtual VectorND<12>    getForce(State state, int rate) override final {
      // TODO: Implement getForce?
      VectorND<12> p;
      return p;
    }
    virtual MatrixND<12,12> getTangent(State state, int rate) override final {
      MatrixND<12,12> K;
      return K;
    }
#endif

    //
    // For Element
    //
    virtual int   update();
    virtual const Matrix &getTangentStiff();
    virtual const Matrix &getMass();
#if 0
    virtual const Vector &getResistingForce();
#endif

    virtual void  zeroLoad() final;
    virtual int   addLoad(ElementalLoad *theLoad, double loadFactor) final;

    virtual int   addInertiaLoadToUnbalance(const Vector &accel) final;
    virtual const Vector &getResistingForceIncInertia() final;
    virtual const Matrix &getInitialStiff();

    // Sensitivity
    const Matrix & getMassSensitivity(int gradNumber);
    virtual int setParameter(const char **argv, int argc, Parameter &param);
    virtual int            updateParameter(int parameterID, Information &info);
    virtual int            activateParameter(int parameterID);


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
   FrameTransform3d* theCoordTransf;
   OpenSees::VectorND<6>   q0;  // Fixed end forces in basic system
   OpenSees::VectorND<6>   p0;  // Reactions in basic system
                                // { 
                                // TODO(cmp): change to size 12
   static Matrix K;
   static Vector P;


   int parameterID;

   double wx;
   double wy;
   double wz;

   std::vector<std::pair<ElementalLoad*,double>> eleLoads;


   Vector p_iner;
  private:
   int cMass;
   double rho;

   VectorND<12> pg;
   double total_mass,
          twist_mass;

   // TODO: Remove
    int releasez; // moment release for bending about z-axis 0=none, 1=I, 2=J, 3=I,J
    int releasey; // same for y-axis

};

