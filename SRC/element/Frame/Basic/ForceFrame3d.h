/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */
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
#include <element/Frame/BasicFrame3d.h>
#include <Vector.h>
#include <VectorND.h>
#include <MatrixND.h>

class Matrix;
class Channel;
class Response;
class ElementalLoad;
class BeamIntegration;
#include <FrameSection.h>

using namespace OpenSees;

class ForceFrame3d: public BasicFrame3d
{
 public:
  ForceFrame3d();
  ForceFrame3d(int tag, std::array<int,2>& nodes,
		    int numSections, FrameSection **sec,
		    BeamIntegration &beamIntegr,
		    FrameTransform3d &coordTransf, double rho = 0.0, 
		    int maxNumIters = 10, double tolerance = 1.0e-12);
  
  ~ForceFrame3d();

  const char *getClassType() const {return "ForceFrame3d";};

  int setNodes();
  int commitState();
  int revertToLastCommit();        
  int revertToStart();
  int update();    
  
  /*
  const Matrix &getTangentStiff();
  const Matrix &getInitialStiff();
  const Matrix &getMass();    

  void zeroLoad();	
  int addLoad(ElementalLoad *theLoad, double loadFactor);
  int addInertiaLoadToUnbalance(const Vector &accel); 
  const Vector &getResistingForce();
  const Vector &getResistingForceIncInertia();
  */
  
  int sendSelf(int cTag, Channel &theChannel);
  int recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
  
  friend OPS_Stream &operator<<(OPS_Stream &s, ForceFrame3d &E);        
  void Print(OPS_Stream &s, int flag =0);    
  
  Response *setResponse(const char **argv, int argc, OPS_Stream &s);
  int getResponse(int responseID, Information &eleInformation);
  
 // AddingSensitivity:BEGIN //////////////////////////////////////////
  int setParameter(const char **argv, int argc, Parameter &param);
  int updateParameter(int parameterID, Information &info);
  int activateParameter(int parameterID);
  const Vector &getResistingForceSensitivity(int gradNumber);
  const Matrix &getKiSensitivity(int gradNumber);
  const Matrix &getMassSensitivity(int gradNumber);
  int commitSensitivity(int gradNumber, int numGrads);
  int getResponseSensitivity(int responseID, int gradNumber,
			     Information &eleInformation);
  // AddingSensitivity:END ///////////////////////////////////////////

 protected:

  // For BasicFrame3d
  virtual VectorND<6>&   getBasicForce();
  virtual MatrixND<6,6>& getBasicTangent(State state, int rate);
  
 private:
  constexpr static int 
//      ndf  = 6,             // number of nodal dofs
        nsr = 6,              // number of section resultants
        NDM  = 3,             // dimension of the problem (3D)
//      NEGD = 12,            // number of element global dof's
        NEBD = 6,             // number of element dof's in the basic system
        maxNumEleLoads = 100,
        maxNumSections = 20,
        maxSubdivisions= 10;
  

  void setSectionPointers(int numSections, FrameSection **secPtrs);
  int getInitialFlexibility(MatrixND<NEBD,NEBD> &fe);
  int getInitialDeformations(Vector &v0);

  void compSectionDisplacements(Vector sectionCoords[], Vector sectionDispls[]) const;
  void initializeSectionHistoryVariables ();
  

  // Add section forces due to element loads
  void addLoadAtSection(VectorND<nsr> &sp, int isec);

 
  //
  // Element State
  //
  // Parameters
  double rho;                    // mass density per unit length
  int    maxIters;               // maximum number of local iterations
  double tol;	                   // tolerance for relative energy norm for local iterations
  // State
  MatrixND<6,6> kv,              // stiffness matrix in the basic system 
                K_save;          // committed stiffness matrix in the basic system
  VectorND<6>   q_pres,          // element resisting forces in the basic system
                q_save;          // committed element end forces in the basic system
  
  int    initialFlag;            // indicates if the element has been initialized

  Matrix *Ki;

  bool isTorsion;
  static constexpr BasicForceLayout force_layout = {
  };


  //
  // Section State
  //
  double wt[maxNumSections];
  double xi[maxNumSections];
  BeamIntegration* beamIntegr;
  int numSections;
  FrameSection** sections;          // array of pointers to sections
  
  MatrixND<nsr,nsr> *fs;         // array of section flexibility matrices
  VectorND<nsr>     *es;         // array of section deformations
  VectorND<nsr>     *Ssr;        // array of section stress resultants
  VectorND<nsr> *es_save;        // array of committed section deformation vectors

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

  //
  // Other
  //
  static Matrix theMatrix;
  static Vector theVector;
  static double workArea[];

  // AddingSensitivity:BEGIN //////////////////////////////////////////
  int parameterID;
  const Vector &computedqdh(int gradNumber);
  const Matrix &computedfedh(int gradNumber);
  void computeSectionForceSensitivity(Vector &dspdh, int isec, int gradNumber);
  // AddingSensitivity:END ///////////////////////////////////////////

//void getForceInterpolatMatrix(double xi, Matrix &b, const ID &code);
//void getDistrLoadInterpolatMatrix(double xi, Matrix &bp, const ID &code);
};

#endif
