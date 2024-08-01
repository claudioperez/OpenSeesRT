//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
//
// Adapted by CMP
//
/*
 * References
 *

    State Determination Algorithm
    ---

    de Souza, R.M. (2000) 
      "Force-based finite element for large displacement inelastic analysis of frames". University of California, Berkeley. Available at: https://www.proquest.com/docview/304624959/D8D738C3AC49427EPQ/1?accountid=14496.

    Neuenhofer, A. and Filippou, F.C. (1998) 
      "Geometrically Nonlinear Flexibility-Based Frame Finite Element", 
      Journal of Structural Engineering, 124(6), pp. 704–711. Available at: https://doi.org/10/d8jvb5.
 
    Spacone, E., V. Ciampi, and F. C. Filippou (1996). 
       "Mixed Formulation of Nonlinear Beam Finite Element."
       Computers and Structures, 58(1):71-83.
    

    Response Sensitivity
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
  ForceDeltaFrame3d(int tag, 
               std::array<int,2>& nodes,
               std::vector<FrameSection*>& sec,
               BeamIntegration& stencil, FrameTransform3d& coordTransf, 
               double rho, int mass_type, bool use_mass,
               int maxNumIters, double tolerance,
               bool includeShear
               );

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
  virtual const Matrix &getTangentStiff() final;

  //void zeroLoad();
  //int addLoad(ElementalLoad *theLoad, double loadFactor);
  //int addInertiaLoadToUnbalance(const Vector &accel);

  const Vector &getResistingForce();
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
  virtual VectorND<6>&   getBasicForce();
  virtual MatrixND<6,6>& getBasicTangent(State state, int rate);

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
         nq = 6 ;              // number of element dof's in the basic system

  enum { maxNumSections = 20 };

  static constexpr FrameStressLayout scheme = {
    FrameStress::N,
    FrameStress::Vy,
    FrameStress::Vz,
    FrameStress::T,
    FrameStress::My,
    FrameStress::Mz,
  };

  // TODO
  double wt[maxNumSections];
  double xi[maxNumSections];


  void getForceInterpolatMatrix(double xi, Matrix& b, const ID& code);
  void getDistrLoadInterpolatMatrix(double xi, Matrix& bp, const ID& code);
  void computew(Vector& w, Vector& wp, double xi[], const Vector& kappa, const Vector& gamma);
  void computedwdq(Matrix& dwidq, const Vector& q, const Vector& w, const Vector& wp,
                   const Matrix& lsk, const Matrix& lsg, const Matrix& lskp, const Matrix& lsgp);
  void computedwzdq(Matrix& dwidzq, const Vector& q, const Vector& wz, const Vector& wpz,
                    const Matrix& lsk, const Matrix& lsg, const Matrix& lskp, const Matrix& lsgp);
  void compSectionDisplacements(Vector sectionCoords[], Vector sectionDispls[]) const;
  void initializeSectionHistoryVariables();


  void getHk(int numSections, double xi[], Matrix& H);
  void getHg(int numSections, double xi[], Matrix& H);
  void getHkp(int numSections, double xi[], Matrix& H);
  void getHgp(int numSections, double xi[], Matrix& H);

  // Reactions of basic system due to element loads
//void computeReactions(double* p0);

  // Section forces due to element loads
  void computeSectionForces(VectorND<nsr>& sp, int isec);

  // Parameters
  int    shear_flag;
  double density;                // mass density per unit length
  double twist_mass;
  double total_mass;
  int  mass_flag;
  bool mass_initialized;
  bool use_density;
  int max_iter;            // maximum number of local iterations
  double tol;              // tolerance for relative energy norm for local iterations

  // Element State
  MatrixND<6,6> K_pres,          // stiffness matrix in the basic system 
                K_past;          // committed stiffness matrix in the basic system
  VectorND<6>   q_pres,          // element resisting forces in the basic system
                q_past;          // committed element end forces in the basic system
  
  int    state_flag;             // indicate if the element has been initialized

  //
  // Section State
  //
  struct GaussPoint {
    double point,
           weight;
    FrameSection* material;

    MatrixND<nsr,nsr> Fs;         // section flexibility
    VectorND<nsr>     es;         // section deformations
    VectorND<nsr>     sr;         // section stress resultants
    VectorND<nsr> es_save;        // committed section deformations
  };

  std::vector<GaussPoint> points;
  BeamIntegration* stencil;

  Matrix* Ki;

  static Matrix theMatrix;
  static Vector theVector;

  // following are added for subdivision of displacement increment
  int maxSubdivisions; // maximum number of subdivisons of dv for local iterations


  //
  // Sensitivity
  //
  int parameterID;
  const Vector& computedqdh(int gradNumber);
  const Matrix& computedfedh(int gradNumber);
  void computedwdh(double dwidh[], int gradNumber, const Vector& q);
//void computeReactionSensitivity(double* dp0dh, int gradNumber);
  void computeSectionForceSensitivity(Vector& dspdh, int isec, int gradNumber);

};

#endif
