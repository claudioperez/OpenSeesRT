//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
/*
 * References
 *

    State Determination Algorithm
    ---
    Neuenhofer, A. and F. C. Filippou (1997). "Evaluation of Nonlinear Frame Finite
    Element Models." Journal of Structural Engineering, 123(7):958-966.
    
    Spacone, E., V. Ciampi, and F. C. Filippou (1996). "Mixed Formulation of
    Nonlinear Beam Finite Element." Computers and Structures, 58(1):71-83.
    
    
    Plastic Hinge Integration
    ---
    Scott, M. H. and G. L. Fenves (2006). "Plastic Hinge Integration Methods for
    Force-Based Beam-Column Elements." Journal of Structural Engineering,
    132(2):244-252.
    
    
    Analytical Response Sensitivity (DDM)
    ---
    Scott, M. H., P. Franchin, G. L. Fenves, and F. C. Filippou (2004).
    "Response Sensitivity for Nonlinear Beam-Column Elements."
    Journal of Structural Engineering, 130(9):1281-1288.

 *
 */

#ifndef ForceDeltaFrame3d_h
#define ForceDeltaFrame3d_h

#include <element/Frame/BasicFrame3d.h>
#include <Node.h>
#include <Matrix.h>
#include <Vector.h>
#include <Channel.h>
#include <BeamIntegration.h>
#include <FrameSection.h>
#include <FrameTransform.h>

class Response;
class ElementalLoad;

class ForceDeltaFrame3d : public BasicFrame3d {
public:
  ForceDeltaFrame3d();
  ForceDeltaFrame3d(int tag, int nodeI, int nodeJ, int numSections, FrameSection** sec,
               BeamIntegration& beamIntegr, FrameTransform3d& coordTransf, double rho = 0.0,
               bool includeShear = false, int maxNumIters = 10, double tolerance = 1.0e-12);

  ~ForceDeltaFrame3d();

  const char*
  getClassType() const
  {
    return "ForceDeltaFrame3d";
  };

  int setNodes();

  int commitState();
  int revertToLastCommit();
  int revertToStart();
  int update();

  //const Matrix &getMass();

  //void zeroLoad();
  //int addLoad(ElementalLoad *theLoad, double loadFactor);
  //int addInertiaLoadToUnbalance(const Vector &accel);

  //const Vector &getResistingForce();
  //const Vector &getResistingForceIncInertia();

  int sendSelf(int cTag, Channel& theChannel);
  int recvSelf(int cTag, Channel& theChannel, FEM_ObjectBroker& theBroker);

  friend OPS_Stream& operator<<(OPS_Stream& s, ForceDeltaFrame3d& E);
  void Print(OPS_Stream& s, int flag = 0);

  Response* setResponse(const char** argv, int argc, OPS_Stream& s);
  int getResponse(int responseID, Information& eleInformation);
  int getResponseSensitivity(int responseID, int gradNumber, Information& eleInformation);

  // AddingSensitivity:BEGIN //////////////////////////////////////////
  int setParameter(const char** argv, int argc, Parameter& param);
  int updateParameter(int parameterID, Information& info);
  int activateParameter(int parameterID);
  const Vector& getResistingForceSensitivity(int gradNumber);
  const Matrix& getKiSensitivity(int gradNumber);
  const Matrix& getMassSensitivity(int gradNumber);
  int commitSensitivity(int gradNumber, int numGrads);
  // AddingSensitivity:END ///////////////////////////////////////////

protected:
  void setSectionPointers(int numSections, FrameSection** secPtrs);
  int getInitialFlexibility(Matrix& fe);
  int getInitialDeformations(Vector& v0);

private:
  constexpr static int
         nsr = 6,
         maxNumEleLoads = 100 ,
         NDM = 3 ,               // dimension of the problem (2d)
         NND = 6 ,               // number of nodal dof's
         NEGD = 12 ,             // number of element global dof's
         NEBD = 6 ;              // number of element dof's in the basic system

  enum { maxNumSections = 20 };
  enum { maxSectionOrder = 10 };

  double wt[maxNumSections];
  double xi[maxNumSections];


  static VectorND<nsr>      es_trial[maxNumSections]; //  strain
  static VectorND<nsr>      sr_trial[maxNumSections]; //  stress resultant
  static MatrixND<nsr,nsr>  Fs_trial[maxNumSections]; //  flexibility
  static constexpr FrameStressLayout scheme = {
    FrameStress::N,
    FrameStress::Vy,
    FrameStress::Vz,
    FrameStress::T,
    FrameStress::My,
    FrameStress::Mz,
  };

  void getForceInterpolatMatrix(double xi, Matrix& b, const ID& code);
  void getDistrLoadInterpolatMatrix(double xi, Matrix& bp, const ID& code);
  void computew(Vector& w, Vector& wp, double xi[], const Vector& kappa, const Vector& gamma);
  void computedwdq(Matrix& dwidq, const Vector& q, const Vector& w, const Vector& wp,
                   const Matrix& lsk, const Matrix& lsg, const Matrix& lskp, const Matrix& lsgp);
  void computedwzdq(Matrix& dwidzq, const Vector& q, const Vector& wz, const Vector& wpz,
                    const Matrix& lsk, const Matrix& lsg, const Matrix& lskp, const Matrix& lsgp);
  void compSectionDisplacements(Vector sectionCoords[], Vector sectionDispls[]) const;
  void initializeSectionHistoryVariables();

//void getG(int numSections, double xi[], Matrix& G);
//void getGinv(int numSections, double xi[], Matrix& Ginv);
  void getHk(int numSections, double xi[], Matrix& H);
  void getHg(int numSections, double xi[], Matrix& H);
  void getHkp(int numSections, double xi[], Matrix& H);
  void getHgp(int numSections, double xi[], Matrix& H);

  // Reactions of basic system due to element loads
  void computeReactions(double* p0);

  // Section forces due to element loads
  void computeSectionForces(VectorND<nsr>& sp, int isec);

  BeamIntegration* beamIntegr;
  int numSections;
  FrameSection** sections; // array of pointers to sections
  bool CSBDI;
  double rho;              // mass density per unit length
  int maxIters;            // maximum number of local iterations
  double tol;              // tolerance for relative energy norm for local iterations

  int initialFlag;         // indicates if the element has been initialized

  Matrix kv;               // stiffness matrix in the basic system
  Vector Se;               // element resisting forces in the basic system

  Matrix K_past; // committed stiffness matrix in the basic system
  Vector q_past; // committed element end forces in the basic system

  Matrix* fs;  // array of section flexibility matrices
  Vector* vs;  // array of section deformation vectors
  Vector* Ssr; // array of section resisting force vectors

  Vector* e_past; // array of committed section deformation vectors

  int numEleLoads; // Number of element load objects
  int sizeEleLoads;
  ElementalLoad** eleLoads;
  double* eleLoadFactors;

  Matrix* Ki;

  static Matrix theMatrix;
  static Vector theVector;
  static double workArea[];

  // following are added for subdivision of displacement increment
  int maxSubdivisions; // maximum number of subdivisons of dv for local iterations

  // AddingSensitivity:BEGIN //////////////////////////////////////////
  int parameterID;
  const Vector& computedqdh(int gradNumber);
  const Matrix& computedfedh(int gradNumber);
  void computedwdh(double dwidh[], int gradNumber, const Vector& q);
  void computeReactionSensitivity(double* dp0dh, int gradNumber);
  void computeSectionForceSensitivity(Vector& dspdh, int isec, int gradNumber);
  // AddingSensitivity:END ///////////////////////////////////////////
};

#endif
