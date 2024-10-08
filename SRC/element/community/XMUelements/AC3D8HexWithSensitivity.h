
// Written: Quan Gu, Yichao Gao and Zhijian Qiu  
// Created: 2015/01/25 
// Acoustic fulid Material
//------------------------------------------

#ifndef AC3D8HexWithSensitivity_H
#define AC3D8HexWithSensitivity_H

#ifndef _bool_h
#include <stdbool.h>
#endif

#include <Element.h>
#include <Node.h>

#include <ID.h>
#include <Renderer.h>
#include <Domain.h>
#include <string.h>

#include <OPS_Globals.h>

#include <Matrix.h>
#include <Vector.h>

#include <NDMaterial.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <Renderer.h>
#include <Domain.h>
#include <string.h>
#include <Information.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ElementResponse.h>

#include <ElementalLoad.h>


#define ELE_TAG_AC3D8HexWithSensitivity 100001

class Node;

class AC3D8HexWithSensitivity: public Element
{

  public:
    AC3D8HexWithSensitivity(int element_number,
			    int node_numb_1,  int node_numb_2,  int node_numb_3,  int node_numb_4,
			    int node_numb_5,  int node_numb_6,  int node_numb_7,  int node_numb_8,
			    NDMaterial * Globalmmodel);

    AC3D8HexWithSensitivity(int element_number,
			    int node_numb_1,  int node_numb_2,  int node_numb_3,  int node_numb_4,
			    int node_numb_5,  int node_numb_6,  int node_numb_7,  int node_numb_8);

    AC3D8HexWithSensitivity();
    ~AC3D8HexWithSensitivity();

    const char* getClassType() const {
      return "AC3D8HexWithSensitivity";
    }

    int getNumExternalNodes() const;
    const ID &getExternalNodes();
    Node **getNodePtrs();
    int getNumDOF();
    
    void setDomain(Domain *theDomain);

    // public methods to set the state of the element
    int commitState();
    int revertToLastCommit();
    int revertToStart();
    
    int update();

    // public methods to obtain stiffness, mass, damping and residual information
    const Matrix &getTangentStiff();
    const Matrix &getInitialStiff();
    const Matrix &getMass();
    const Matrix &getConsMass();
    
    const Matrix &getDamp();

    void zeroLoad();
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);

    const Vector &getResistingForce();
    const Vector &getResistingForceIncInertia();

    // public methods for element output
    int sendSelf (int commitTag, Channel &theChannel);
    int recvSelf (int commitTag, Channel &theChannel, FEM_ObjectBroker
		  &theBroker);

    Response *setResponse (const char **argv, int argc, OPS_Stream &output);
    int getResponse (int responseID, Information &matInformation);
      
    int displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode);
    void Print(OPS_Stream &s, int flag =0);
    
    // const Vector &getExternalLoadIncInertia();
    
    double get_Gauss_p_c(short order, short point_numb);
    double get_Gauss_p_w(short order, short point_numb);
    	// AddingSensitivity:BEGIN //////////////////////////////////////////
    int            setParameter                (const char **argv, int argc, Parameter &param);
    int            updateParameter             (int parameterID, Information &info);
    int            activateParameter           (int parameterID);
    const Vector & getResistingForceSensitivity(int gradNumber);
    const Matrix & getKiSensitivity            (int gradNumber);
    const Matrix & getMassSensitivity          (int gradNumber);
    int            commitSensitivity           (int gradNumber, int numGrads);
    // AddingSensitivity:END ///////////////////////////////////////////

    
  protected:

  private:
    constexpr static int NEN = 8,          // number of nodes in element
                         NIP = 8;
    static const int numDOF;               // DOF number of element
    static const int r_integration_order;  // Gauss-Legendre integration order in r direction
    static const int s_integration_order;  // Gauss-Legendre integration order in s direction
    static const int t_integration_order;  // Gauss-Legendre integration order in t direction
    static const int dim;                  // spatial dimension

    // private attributes - a copy for each object of the class
    ID  connectedExternalNodes; // Tags of quad nodes

    Matrix *Ki;
    Node *theNodes[8];
    
    Matrix **L;       // global differential operator
    double *detJ;     // determinant of Jacobian matrix
    
    NDMaterial **theMaterial; // pointer to the ND material objects
    
    Vector Q;		           // Applied nodal loads
    Node *nodePointers[8] ;      //pointers to four nodes   
    double *impVals;
    
    unsigned char hasConstrained; // 0:not set, 1: constrained, 2: not constrained

  private:
    static Matrix K;		   // Element stiffness Matrix
    static Matrix C;		   // Element damping matrix
    static Matrix M;		   // Element mass matrix
    static Vector P;		   // Element resisting force vector
    // static Matrix B;		   // [B] matrix
    static Matrix mass ; 
    static ID actDOFs;     // activated element dofs, add Yichao Gao
    
    
    static Matrix **H;   // Matrix array holds h
    static Matrix **DH;  // Matrix array holds h
    static Matrix **HH;  // Matrix array holds h*h
    
    static Vector VecA;  // Vector(numDOF);
    static Vector VecV;  // Vector(numDOF);
    static Vector VecD;  // Vector(numDOF);

    // AddingSensitivity:BEGIN ///////////////////////////
    int parameterID;
    // AddingSensitivity:END /////////////////////////////
    //material information
    NDMaterial *materialPointers[8] ; //pointers to eight materials
	
  public:
    // shape function
    Matrix interp_fun(double r, double s, double t);
    Matrix diff_interp_fun(double r, double s, double t);
    
    // compute Jacobian matrix and its determinant
    // Matrix Jacobian(Matrix dh);
    // Matrix Jacobian(Matrix dh, Matrix h);
    double Jacobian_det(Matrix Jac);
    
    int computeH();
    int computeHH();
    int computeDiff();
    
    // get nodal coordinates
    Matrix getNodalCoords();
    
    // get nodal forces from stress
    Matrix getNodalForces();
    
    // get total displacement
    Matrix getTotalDisp();
    
    // get [B] Matrix
    // Matrix &getBMatrix(Matrix dhGlobal);
    
    // handling active dofs
    ID *getActiveDofs();
    
    // for acoustic elements
    int getIntegrateFlag();
    
    int setNDMaterial(NDMaterial *Globalmmodel);
    
    AC3D8HexWithSensitivity & operator[](int subscript);
    
    // for external forces
    Vector nodal_forces_from_displacement(const Vector &u);
    
    // face shape function
    Matrix interp_fun_face(double r, double s);
    Matrix diff_interp_fun_face(double r, double s);
    
    Matrix getFaceNodalCoords(int face_num);
    void localFaceMapping(int face_num, ID &face_nodes);
    
    int setImpedance(int face, double val);
    Matrix get_face_impedance(int face_num);

};


#endif





