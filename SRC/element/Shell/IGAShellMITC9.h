/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */
//
// https://github.com/li-ming-jiang/MyOpenSees/blob/master/SRC/element/shell/IGAShellMITC9.cpp
// Written: Leopoldo Tesser, Diego Talledo
// 9-node lagrandian shell element with membrane and drill
//
#include <ID.h> 
#include <Vector.h>
#include <Matrix.h>
#include <Element.h>
#include <Node.h>
#include <SectionForceDeformation.h>


class IGAShellMITC9 : public Element
{
  public:

    //null constructor
    IGAShellMITC9();
    //full constructor
    IGAShellMITC9(int tag, 
	     int node1,
	     int node2,
	     int node3,
	     int node4,
	     int node5,
	     int node6,
	     int node7,
	     int node8,
	     int node9,
	     SectionForceDeformation &theMaterial ) ;
    //destructor 
    virtual ~IGAShellMITC9( ) ;

    //get the number of external nodes
    int getNumExternalNodes( ) const ;
    //return connected external nodes
    const ID &getExternalNodes( ) ;
    Node **getNodePtrs();
    //return number of dofs
    int getNumDOF( ) ;

    //set domain 
    void setDomain( Domain *theDomain ) ;
    //commit state
    int commitState( ) ;
    //revert to last commit 
    int revertToLastCommit( ) ;
    //revert to start 
    int revertToStart( ) ;

    //print out element data
    void Print( OPS_Stream &s, int flag ) ;
	
    //return stiffness matrix 
    const Matrix &getTangentStiff( ) ;
    const Matrix &getInitialStiff();
    const Matrix &getMass();

    // methods for applying loads
    void zeroLoad(void);	
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);

    //get residual
    const Vector &getResistingForce( ) ;
    
    //get residual with inertia terms
    const Vector &getResistingForceIncInertia( ) ;

    // public methods for element output
    int sendSelf (int commitTag, Channel &theChannel);
    int recvSelf (int commitTag, Channel &theChannel, FEM_ObjectBroker 
		  &theBroker);


    Response* setResponse(const char **argv, int argc, OPS_Stream &output);
    int getResponse(int responseID, Information &eleInfo);

  private : 

    //static data
    Matrix stiff ;
    Vector resid ;
    Matrix mass ;
    Matrix damping ;

    //quadrature data
    static const double root3 ;
    static const double root3_over_root5 ;
    double* sg;
    double* tg;
    double* wg;

    int numCPs;
    int ndof;
    int obfs[2];
    Vector KnotVect_x;
    Vector KnotVect_y;
    ID eleIdInfo;// stores ex, nex, ey, ney
    Matrix* N0nx; // stores B-Spline Basis Function and 1st order derivatives
    Matrix* N0ny;

    //node information
    ID connectedExternalNodes ;  //nine node numbers
    //pointers to nine nodes
    Node **nodePointers;

    //drilling stiffness
    double Ktt ;

    //material information: pointers to four materials
    SectionForceDeformation **materialPointers;
					  
    //local nodal coordinates, two coordinates for each of nine nodes
    //static double xl[][9] ; 
    double** xl;

    //shell basis vectors
    double g1[3] ;
    double g2[3] ;
    double g3[3] ;

    //compute local coordinates and basis
    void computeBasis( ) ;
        
    //inertia terms
    void formInertiaTerms( int tangFlag ) ;

    //form residual and tangent					  
    void formResidAndTangent( int tang_flag ) ;

    //compute Jacobian matrix and inverse at point {L1,L2}
    //void  computeJacobian( double L1, double L2,const double x[2][9], 
    //                       Matrix &JJ,Matrix &JJinv ) ;

    //compute Bdrill matrix
    double* computeBdrill( int node, double** shp ) ;

    //assemble a B matrix 
    const Matrix& assembleB( const Matrix &Bmembrane,
			     const Matrix &Bbend, 
			     const Matrix &Bshear ) ;
    
    //compute Bmembrane matrix
    const Matrix& computeBmembrane( int node, double** shp ) ;
    
    //compute Bbend matrix
    const Matrix& computeBbend( int node, double** shp ) ;
    
    //compute standard Bshear matrix
    const Matrix&  computeBshear( int node, double** shp ) ;
    
    //Matrix transpose
    Matrix transpose( int dim1, int dim2, const Matrix &M ) ;

    //shape function routine for four node quads
    double shape2d(int qx, int qy, ID eleIdInfo, Vector KnotVect_x, Vector KnotVect_y, Matrix x);
    int findSpan(int order, int ncp, int eleId, int nele, Vector KnotVect);
    void DerBasisFuns(double Idx, Vector Pts, int obf, int n, Vector KnotVect, Matrix** N0n);
    void calcDersBasisFunsAtGPs(int obf, int ncp, Vector KnotVect, int d, int NGPs, int Idx, double* J2, Vector* Xg, Matrix** N0n);
    void Rationalize(Vector WeightsCP, Vector N0, Matrix N1, Vector* R0, Matrix* R1);
    void GaussRule(int NGPs1, Vector* Xg, Vector* WXg);
    Vector dotProduct(Vector VecA, Vector VecB);
    Vector dotDivide(Vector VecA, Vector VecB);
    // vector for applying loads
    Vector *load;
    Matrix *Ki;
    double **shp;
} ; 

