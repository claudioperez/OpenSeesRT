/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */
//
// Ed "C++" Love
//
// Mixed Presssure/Volume Nine Node Quadrilateral
// Plane Strain (NOT PLANE STRESS)
//
#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

#include <ID.h> 
#include <Vector.h>
#include <Matrix.h>
#include <Element.h>
#include <Node.h>
#include <NDMaterial.h>

class NineNodeMixedQuad : public Element {

  public :

    NineNodeMixedQuad();

    NineNodeMixedQuad( int tag, 
		       int node1,
		       int node2,
		       int node3,
		       int node4,
		       int node5,
		       int node6,
		       int node7,
		       int node8,
		       int node9,
		       NDMaterial &theMaterial ) ;

    ~NineNodeMixedQuad( ) ;

    const char *getClassType(void) const {return "NineNodeMixedQuad";}
    static constexpr const char* class_name = "NineNodeMixedQuad";

    //set domain 
    void setDomain( Domain *theDomain ) ;

    //get the number of external nodes
    int getNumExternalNodes() const;
    int getNumDOF();
 
    //return connected external nodes
    const ID &getExternalNodes( ) ;
    Node **getNodePtrs();


    //
    int commitState( ) ;
    int revertToLastCommit( ) ;
    int revertToStart( ) ;

    void Print( OPS_Stream &s, int flag ) ;
	
    // 
    const Matrix &getTangentStiff();
    const Matrix &getInitialStiff();     
    const Matrix &getMass();

    void zeroLoad( ) ;
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);

    //get residual
    const Vector &getResistingForce( ) ;
    
    //get residual with inertia terms
    const Vector &getResistingForceIncInertia( ) ;

    // public methods for element output
    Response *setResponse(const char **argv, int argc, OPS_Stream &s);

    int getResponse(int responseID, Information &eleInformation);

    int sendSelf (int commitTag, Channel &);
    int recvSelf (int commitTag, Channel &, FEM_ObjectBroker &);

  private : 

    //static data
    static Matrix stiff ;
    static Vector resid ;
    static Matrix mass ;
    static Matrix damping ;
    
    
    // quadrature data
    static double root06;
    static double sg[3] ;
    static double wg[3] ;

    //node information
    ID connectedExternalNodes ;  //nine node numbers
    Node *nodePointers[9] ;      //pointers to nine nodes
					
    //material information
    NDMaterial *materialPointers[9] ; //pointers to nine materials
					  
    //nodal coordinates, two coordinates for each of nine nodes
    static double xl[][9] ; 
    
    //form residual and tangent					  
    void formResidAndTangent( int tang_flag ) ;

    //inertia terms
    void formInertiaTerms( int tangFlag ) ;

    const Matrix& computeBbar( int node, 
			       const double natCoor[2], 
			       const double shp[3][9], 
			       double shpBar[3][9][3] ) ;

    //shape function routine for four node quads
    void shape2dNine( double coor[2], 
		  const double x[2][9], 
		  double shp[3][9], 
		  double &xsj ) ;

    //nodal coordinates
    void computeBasis( ) ;

    //1d quadratic shape functions
    double shape1d( int code, int node, double xi ) ;

    Vector *load;
    Matrix *Ki;
} ; 
