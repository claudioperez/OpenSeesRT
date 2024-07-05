/* ------------------------------------------------------------------ *
      OpenSees - Open System for Earthquake Engineering Simulation      
            Pacific Earthquake Engineering Research Center              
 * ------------------------------------------------------------------ */
//
// Description: This file contains the class definition for CubicFrame3d.
// The element displacement field gives rise to constant axial strain,
// linear curvature, and constant twist angle.
//
#ifndef CubicFrame3d_h
#define CubicFrame3d_h

#include <array>
#include <Frame/BasicFrame3d.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <FrameSection.h>

class Node;
// class FrameSection;
class FrameTransform3d;
class BeamIntegration;
class Response;

class CubicFrame3d : public BasicFrame3d
{
  public:
    CubicFrame3d(int tag, std::array<int,2>& nodes,
                     int numSections, FrameSection **s,
                     BeamIntegration &bi, FrameTransform3d &coordTransf,
                     double rho = 0.0, int cMass=0, int rz=0, int ry=0);
    CubicFrame3d();
    ~CubicFrame3d();

    const char *getClassType() const {return "CubicFrame3d";};
    static constexpr const char* class_name = "CubicFrame3d";


    // public methods to set the state of the element    
    int commitState();
    int revertToLastCommit();
    int revertToStart();

    // public methods to obtain stiffness, mass, damping and residual information    
    int update();
//  const Matrix &getTangentStiff();
//  const Matrix &getInitialStiff();
//  const Matrix &getMass();

//  void zeroLoad();
//  int addLoad(ElementalLoad *theLoad, double loadFactor);
//  int addInertiaLoadToUnbalance(const Vector &accel);
//  const Vector &getResistingForceIncInertia();

    // public methods for element output
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker 
                 &theBroker);

    void Print(OPS_Stream &s, int flag =0);

    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int responseID, Information &eleInfo);

    // AddingSensitivity:BEGIN //////////////////////////////////////////
    virtual int setParameter(const char **argv, int argc, Parameter &param) final;
//  virtual int            updateParameter(int parameterID, Information &info);
//  virtual int            activateParameter(int parameterID) final;
    virtual const Vector & getResistingForceSensitivity(int gradNumber);
    virtual const Matrix & getInitialStiffSensitivity(int gradNumber);
//  virtual const Matrix & getMassSensitivity(int gradNumber);
    virtual int            commitSensitivity(int gradNumber, int numGrads);
    // AddingSensitivity:END ///////////////////////////////////////////

protected:
    // For BasicFrame3d
    virtual  OpenSees::VectorND<6>&   getBasicForce() final;
    virtual  OpenSees::MatrixND<6,6>& getBasicTangent(State state, int rate) final;
    // For FiniteElement
    virtual int setNodes() final;

private:
//  const Matrix &getInitialBasicStiff();

    int numSections;
    FrameSection **theSections; // pointer to the ND material objects
    BeamIntegration *beamInt;

//  Vector q;           // Basic force
    OpenSees::MatrixND<6,6> kb;
    OpenSees::VectorND<6>   q;

    double rho;         // Mass density per unit length

    constexpr static int max_nip = 20;
    double xi[max_nip];
    double wt[max_nip];

    static double workArea[];

//  static constexpr int scheme[] = {
    static constexpr FrameStressLayout scheme = {
      FrameStress::N,
      FrameStress::Mz,
      FrameStress::My,
      FrameStress::T,
    };
};

#endif

