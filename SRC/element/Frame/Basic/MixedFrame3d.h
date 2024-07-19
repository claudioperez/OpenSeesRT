//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
#ifndef MixedFrame3d_h
#define MixedFrame3d_h

// Written: Mark D. Denavit, University of Illinois at Urbana-Champaign
//
// Description: This file contains the interface for the MixedFrame3d class.
// It defines the class interface and the class attributes.
//
#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <FrameTransform.h>

class Node;
class Channel;
class Response;
class BeamIntegration;
class FrameSection;

class MixedFrame3d : public Element {
public:
  // constructors
  MixedFrame3d(int tag, int nodeI, int nodeJ, int numSections,
               FrameSection** sectionPtrs, BeamIntegration& bi, FrameTransform3d& coordTransf,
               double massDensPerUnitLength, int doRayleigh, bool geomLinear = true);
  MixedFrame3d();

  // destructor
  ~MixedFrame3d();

  // public methods to obtain information about dof & connectivity
  int getNumExternalNodes() const;
  const ID& getExternalNodes();
  Node** getNodePtrs();
  int getNumDOF();
  void setDomain(Domain* theDomain);

  // public methods to set the state of the element
  int update();
  int commitState();
  int revertToLastCommit();
  int revertToStart();

  // public methods to obtain stiffness, mass, damping and residual information
  const Matrix& getTangentStiff();
  const Matrix& getInitialStiff();
  const Matrix& getDamp();

  void zeroLoad();
  int addLoad(ElementalLoad* theLoad, double loadFactor);

  const Vector& getResistingForce();
  const Vector& getResistingForceIncInertia();

  // public methods for output
  int sendSelf(int cTag, Channel& theChannel);
  int recvSelf(int cTag, Channel& theChannel, FEM_ObjectBroker& theBroker);
  void Print(OPS_Stream& s, int flag = 0);
  friend OPS_Stream& operator<<(OPS_Stream& s, MixedFrame3d& E);

  Response* setResponse(const char** argv, int argc, OPS_Stream& output);
  int getResponse(int responseID, Information& eleInfo);

  const char*
  getClassType() const
  {
    return "MixedFrame3d";
  };

protected:
private:
  // Private Functions - Shape Functions
  Matrix getNld_hat(int sec, const Vector& v, double L, bool geomLinear);
  Vector getd_hat(int sec, const Vector& v, double L, bool geomLinear);
  Matrix getNd1(int sec, const Vector& v, double L, bool geomLinear);
  Matrix getNd2(int sec, double P, double L);
  Matrix getKg(int sec, double P, double L);
  Matrix getMd(int sec, Vector dShapeFcn, Vector dFibers, double L);

  // Private Functions - Interaction With The Sections
  void getSectionTangent(int sec, int type, Matrix& kSection, double& GJ);
  void getSectionStress(int sec, Vector& fSection, double& torsion);
  void setSectionDeformation(int sec, Vector& defSection, double& twist);

  // Private Attributes - a copy for each object of the class
  BeamIntegration* beamIntegr;        //
  int numSections;                    //
  FrameSection** sections;            // array of pointers to sections
  FrameTransform3d* crdTransf;        // pointer to coordinate transformation object

  int doRayleigh;  // flag for whether or not rayleigh damping is active for this element
  bool geomLinear; // flag for whether or not the internal geometric nonlinearity is active
  double rho;      // mass density per unit length

  int itr; // Counts the number of iterations
  int initialFlag;

  // Attributes that do NOT change during the analysis
  double L0;
  Matrix* Ki;

  // Element Load Variables
  Matrix* sp;
  double p0[5]; // Reactions in the basic system due to element loads

  // Attributes that change during the analysis
  Vector V;
  Vector internalForceOpenSees;
  Vector naturalForce;
  Vector lastNaturalDisp;
  Matrix Hinv;
  Matrix GMH;
  Matrix kv; // stiffness matrix in the basic system
  Vector* sr_trial;
  Vector* es_trial;
  Matrix* fs_trial;

  // Committed versions
  Vector  committedV;
  Vector  qe_past;
  Vector  commitedNaturalForce;
  Vector  commitedLastNaturalDisp;
  Matrix  commitedHinv;
  Matrix  commitedGMH;
  Matrix  ke_past;
  Vector* sr_past;
  Vector* es_past;
  Matrix* commitedSectionFlexibility;

  // static data - single copy for all objects of the class
  static Matrix theMatrix;
  static Vector theVector;
  static Matrix transformNaturalCoords;
  static Matrix transformNaturalCoordsT;
  // matrix to transform the natural coordinates from what the coordinate transformation uses and what the element uses

  // These variable are always recomputed, so there is no need to store them for each instance of the element
  static Vector* ei_trial;
  static Vector* si_trial;
  static Matrix* nldhat;
  static Matrix* nd1;
  static Matrix* nd2;
  static Matrix* nd1T;
  static Matrix* nd2T;

  static constexpr FrameStressLayout scheme = {
    FrameStress::N,
    FrameStress::Mz,
    FrameStress::My,
    FrameStress::T,
  };
};

#endif
