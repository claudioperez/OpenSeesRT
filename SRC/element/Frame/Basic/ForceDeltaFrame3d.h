//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
//
#ifndef ForceDeltaFrame3d_h
#define ForceDeltaFrame3d_h

#include <array>
#include <vector>
#include <element/Frame/BasicFrame3d.h>
#include <Vector.h>
#include <Channel.h>
#include <FrameSection.h>
#include <FrameTransform.h>

class Response;
class ElementalLoad;
class BeamIntegration;

class ForceDeltaFrame3d : public BasicFrame3d {
public:
  ForceDeltaFrame3d(int tag, 
               std::array<int,2>& nodes,
               std::vector<FrameSection*>& sec,
               BeamIntegration&  stencil,
               FrameTransform3d& coordTransf, 
               double rho, int mass_type, bool use_mass,
               int maxNumIters, double tolerance,
               bool includeShear
  );

  ForceDeltaFrame3d();

  ~ForceDeltaFrame3d();

  const char*
  getClassType() const final {
    return "ForceDeltaFrame3d";
  }

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

  int sendSelf(int cTag, Channel&);
  int recvSelf(int cTag, Channel&, FEM_ObjectBroker&);

  void Print(OPS_Stream& s, int flag = 0);

  Response* setResponse(const char** argv, int argc, OPS_Stream& s);
  int getResponse(int responseID, Information& information);
  int getResponseSensitivity(int responseID, int gradNumber, Information& information);

  // Parameters
  int setParameter(const char** argv, int argc, Parameter&);
  int updateParameter(int parameterID, Information&);
  int activateParameter(int parameterID);

  // Sensitivity
  const Vector& getResistingForceSensitivity(int gradNumber);
  const Matrix& getKiSensitivity(int gradNumber);
  const Matrix& getMassSensitivity(int gradNumber);
  int commitSensitivity(int gradNumber, int numGrads);


protected:
  virtual VectorND<6>&   getBasicForce();
  virtual MatrixND<6,6>& getBasicTangent(State state, int rate);


private:
  constexpr static int
         nsr = 6,
         nen = 2,               // number of element nodes
         ndm = 3 ,              // dimension of the problem (2d)
         ndf = 6 ,              // dofs per node
         nq = 6 ,               // number of element dof's in the basic system
         maxNumSections = 20,
         maxNumEleLoads = 100;
//  enum { maxNumSections = 20 };

  static constexpr FrameStressLayout scheme = {
    FrameStress::N,
    FrameStress::Vy,
    FrameStress::Vz,
    FrameStress::T,
    FrameStress::My,
    FrameStress::Mz,
  };

  void setSectionPointers(int numSections, FrameSection** secPtrs);
  int getInitialFlexibility(Matrix& fe);
  int getInitialDeformations(Vector& v0);

//void getForceInterpolatMatrix(double xi, Matrix& b, const ID& code);
//void getDistrLoadInterpolatMatrix(double xi, Matrix& bp, const ID& code);
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
  // TODO
  double wt[maxNumSections];
  double xi[maxNumSections];



  // Parameters
  double density;                // mass density per unit length
  double twist_mass;
  double total_mass;
  int    mass_flag;
  bool   mass_initialized;
  bool   use_density;

  int max_iter;            // maximum number of local iterations
  double tol;              // tolerance for relative energy norm for local iterations
  int    shear_flag;

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
  BeamIntegration*        stencil;

  Matrix* Ki;

  static Vector theVector;



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
