//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
// Description: This file contains the class definition for DisplEulerFrame3d.
// The element displacement field gives rise to constant axial strain,
// linear curvature, and constant twist angle.
//
#ifndef DisplEulerFrame3d_h
#define DisplEulerFrame3d_h

#include <array>
#include <vector>
#include <Frame/FiniteElement.h>
#include <Matrix.h>
#include <Vector.h>
#include <VectorND.h>
#include <ID.h>
#include <FrameSection.h>

class Node;
class FrameTransform3d;
class BeamIntegration;
class Response;

class DisplEulerFrame3d : public FiniteElement<2, 3, 6>
{
  public:
    DisplEulerFrame3d(int tag, std::array<int,2>& nodes,
             std::vector<FrameSection*> &secs,
             BeamIntegration &bi, FrameTransform3d &coordTransf,
             double rho, int mass_flag, bool use_mass);

    DisplEulerFrame3d();
    ~DisplEulerFrame3d();

    const char *getClassType() const {return "DisplEulerFrame3d";};
    static constexpr const char* class_name = "DisplEulerFrame3d";


    // public methods to set the state of the element    
    int commitState();
    int revertToLastCommit();
    int revertToStart();
    int update();

    // public methods to obtain stiffness, mass, damping and residual information    
    const Vector &getResistingForce();
    const Matrix &getTangentStiff();
    const Matrix &getInitialStiff();

    void zeroLoad();
    int addLoad(ElementalLoad *theLoad, double loadFactor);


    // public methods for element output
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker 
          &theBroker);

    void Print(OPS_Stream &s, int flag =0);

    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int responseID, Information &eleInfo);

    // Parameters
    int setParameter(const char **argv, int argc, Parameter &param);
    int            updateParameter(int parameterID, Information &info);
    int            activateParameter(int parameterID);

    virtual int getIntegral(Field field, State state, double& total);

  protected:
    int setNodes();
    
  private:
    constexpr static int
      nen = 2,
      ndf = 6,
      nsr = 4,
      maxNumSections = 20; // TODO: remove

    int mass_flag;
    bool use_density;

    struct ShapeBasis {
      double eval,
             grad,
             hess;
    };
    ShapeBasis getBasis(int, double);

    struct GaussPoint {
      double point,
             weight;
      FrameSection* material;

      MatrixND<nsr,nsr> fs;         // flexibility matrix
      VectorND<nsr>     es;         // deformations
      VectorND<nsr>     Ssr;        // stress resultants
      VectorND<nsr> es_save;        // committed section deformations
    };

    std::vector<GaussPoint> points;

    int numSections;
    FrameSection **sections;               // vector of Sections
    FrameTransform3d *theCoordTransf;      // pointer to coordinate transformation

    BeamIntegration *beamInt;


    OpenSees::VectorND<6*nen> p;   //
    OpenSees::VectorND<6*nen> q0;  // Fixed end forces in basic system
    OpenSees::VectorND<5> p0;      // Reactions in basic system

    double density;    // Mass density per unit length

    int parameterID;

    static constexpr FrameStressLayout scheme = {
      FrameStress::N,
      FrameStress::Mz,
      FrameStress::My,
      FrameStress::T,
    };

    constexpr static int max_nip = 20;
    double xi[max_nip];
    double wt[max_nip];
};

#endif

