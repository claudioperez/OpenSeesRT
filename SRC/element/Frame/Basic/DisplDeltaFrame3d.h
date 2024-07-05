/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */
//
// Description: This file contains the class definition for DisplDeltaFrame3d.
// The element displacement field gives rise to constant axial strain,
// linear curvature, and constant twist angle.

// Modified by: Xinlong Du and Jerome F. Hajjar, Northeastern University, USA; Year 2019
// Description: Adapted for analysis of asymmetric sections with introducing
// high-order axial terms for the basic element formulation
// References:
// Du, X., & Hajjar, J. (2021). Three-dimensional nonlinear displacement-based beam element
// for members with angle and tee sections. Engineering Structures, 239, 112239.
//
// Written: MHS
// Created: Feb 2001
//
#ifndef DisplDeltaFrame3d_h
#define DisplDeltaFrame3d_h

#include <array>
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

class DisplDeltaFrame3d : public FiniteElement<2, 3, 6>
{
  public:
    DisplDeltaFrame3d(int tag, std::array<int,2>& nodes,
             int numSections, FrameSection **s,
             BeamIntegration &bi, FrameTransform3d &coordTransf,
             double rho = 0.0, int cMass = 0);
    DisplDeltaFrame3d();
    ~DisplDeltaFrame3d();

    const char *getClassType() const {return "DisplDeltaFrame3d";};
    static constexpr const char* class_name = "DisplDeltaFrame3d";


    // public methods to set the state of the element    
    int commitState();
    int revertToLastCommit();
    int revertToStart();
    int update();

    // public methods to obtain stiffness, mass, damping and residual information    
    const Matrix &getTangentStiff();
    const Matrix &getInitialStiff();

    void zeroLoad();
    int addLoad(ElementalLoad *theLoad, double loadFactor);

    const Vector &getResistingForce();

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

  protected:
    int setNodes();
    
  private:
    constexpr static int
      nen = 2,
      ndf = 6,
      nsr = 4;

    int cMass;

    enum {maxNumSections = 20};


    int numSections;
    FrameSection **theSections; // pointer to the Sections
    FrameTransform3d *theCoordTransf;        // pointer to coordinate transformation

    BeamIntegration *beamInt;


    OpenSees::VectorND<6*nen> p;   //
    OpenSees::VectorND<6*nen> q0;  // Fixed end forces in basic system
    OpenSees::VectorND<5> p0;      // Reactions in basic system

    double rho;    // Mass density per unit length

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

