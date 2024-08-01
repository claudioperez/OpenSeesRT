//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
/*
 * References
 *
 *
 *  State Determination Algorithm
 *  ---
 *  Neuenhofer, A. and F. C. Filippou (1997). "Evaluation of Nonlinear Frame Finite
 *  Element Models." Journal of Structural Engineering, 123(7):958-966.
 *
 *  Spacone, E., V. Ciampi, and F. C. Filippou (1996). "Mixed Formulation of
 *  Nonlinear Beam Finite Element." Computers and Structures, 58(1):71-83.
 *
 *  Analytical Response Sensitivity (DDM)
 *  ---
 *  Scott, M. H., P. Franchin, G. L. Fenves, and F. C. Filippou (2004).
 *  "Response Sensitivity for Nonlinear Beam-Column Elements."
 *  Journal of Structural Engineering, 130(9):1281-1288.
 *
 */
//
#ifndef ForceFrame3d_h
#define ForceFrame3d_h
//
#include <array>
#include <vector>
#include <element/Frame/BasicFrame3d.h>
#include <Vector.h>
#include <VectorND.h>
#include <MatrixND.h>
#include <FrameSection.h>

class Matrix;
class Channel;
class Response;
class ElementalLoad;
class BeamIntegration;

using namespace OpenSees;

class ForceFrame3d: public BasicFrame3d
{
 public:
  ForceFrame3d();
  ForceFrame3d(int tag, std::array<int,2>& nodes,
        std::vector<FrameSection*>& sections,
		    BeamIntegration &beamIntegr,
		    FrameTransform3d &coordTransf, 
        double density, int cMass, bool use_density,
		    int max_iter, double tolerance
  );
  
  ~ForceFrame3d();

  const char *getClassType() const {return "ForceFrame3d";};

  int setNodes();
  int commitState();
  int revertToLastCommit();        
  int revertToStart();
  int update();    

  virtual const Matrix &getMass() final;
  virtual const Matrix &getTangentStiff() final;
  const Matrix &getInitialStiff() final;
  const Vector &getResistingForce();

  /*
  void zeroLoad();	
  int addLoad(ElementalLoad *theLoad, double loadFactor);
  const Vector &getResistingForceIncInertia();
  int addInertiaLoadToUnbalance(const Vector &accel); 
  */
  
  Response *setResponse(const char **argv, int argc, OPS_Stream &s);
  int getResponse(int responseID, Information &eleInformation);
  

  int setParameter(const char **argv, int argc, Parameter &param);
  int updateParameter(int parameterID, Information &info);
  int activateParameter(int parameterID);
  const Vector &getResistingForceSensitivity(int gradNumber);
  const Matrix &getKiSensitivity(int gradNumber);
  const Matrix &getMassSensitivity(int gradNumber);
  int commitSensitivity(int gradNumber, int numGrads);
  int getResponseSensitivity(int responseID, int gradNumber, Information &info);


  virtual int getIntegral(Field field, State state, double& total);


  // MovableObject
  int sendSelf(int cTag, Channel &theChannel);
  int recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
  
  // TaggedObject
  void Print(OPS_Stream &s, int flag =0);    
  
 protected:

  // For BasicFrame3d
  virtual VectorND<6>&   getBasicForce();
  virtual MatrixND<6,6>& getBasicTangent(State state, int rate);
  
 private:
  //
  // Constexpr
  //
  constexpr static int 
        nsr = 6,              // number of section resultants
        ndm = 3,              // dimension of the problem (3D)
        nq  = 6,              // number of element dof's in the basic system
        maxNumSections = 20,
        maxSubdivisions= 10;

  static constexpr FrameStressLayout scheme = {
    FrameStress::N,
    FrameStress::Vy,
    FrameStress::Vz,
    FrameStress::T,
    FrameStress::My,
    FrameStress::Mz,
  };

  enum Respond: int {
    GlobalForce = 1,
    BasicPlasticDeformation = 4,
    LocalForce  = 2,
    BasicForce  = 7,
    BasicStiff  =19,
  };
  
  int setSectionPointers(std::vector<FrameSection*>&);
  int getInitialFlexibility(MatrixND<nq,nq> &fe);
  int getInitialDeformations(Vector &v0);
  
  // Add section forces due to element loads
  void addLoadAtSection(VectorND<nsr> &sp, double x);

  void compSectionDisplacements(Vector sectionCoords[], Vector sectionDispls[]) const;
  void initializeSectionHistoryVariables();

  // Sensitivity
  int parameterID;
  const Vector &computedqdh(int gradNumber);
  const Matrix &computedfedh(int gradNumber);
  void computeSectionForceSensitivity(Vector &dspdh, int isec, int gradNumber);

 
  //
  // Element State
  //
  // Parameters
  double density;                // mass density per unit length
  double twist_mass;
  double total_mass;
  int  mass_flag;
  bool mass_initialized;
  bool use_density;

  Matrix *Ki;

  int    max_iter;               // maximum number of local iterations
  double tol;	                   // tolerance for relative energy norm for local iterations

  // State
  MatrixND<12,12> tangent;
  VectorND<12>    residual,
                  inertia;

  MatrixND<6,6> K_pres,          // stiffness matrix in the basic system 
                K_save;          // committed stiffness matrix in the basic system
  VectorND<6>   q_pres,          // element resisting forces in the basic system
                q_save;          // committed element end forces in the basic system
  
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
  


//void getForceInterpolatMatrix(double xi, Matrix &b, const ID &code);
//void getDistrLoadInterpolatMatrix(double xi, Matrix &bp, const ID &code);
};

#endif
