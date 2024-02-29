/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
//
#include <stdio.h> 
#include <string.h>
#include <stdlib.h> 
#include <math.h> 

#include <ID.h> 
#include <Vector.h>
#include <Matrix.h>
#include <MatrixND.h>
#include <Element.h>
#include <Node.h>
#include <NDMaterial.h>
#include <Domain.h>
#include <ErrorHandler.h>
#include <EnhancedQuad.h>
#include <Renderer.h>
#include <ElementResponse.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
using namespace OpenSees;
#if 0
void * OPS_ADD_RUNTIME_VPV(OPS_EnhancedQuad)
{
    if (OPS_GetNDM() != 2 || OPS_GetNDF() != 2) {
        opserr << "WARNING -- model dimensions and/or nodal DOF not compatible with quad element\n";
        return 0;
    }
    
    if (OPS_GetNumRemainingInputArgs() < 8) {
        opserr << "WARNING insufficient arguments\n";
        opserr << "Want: element ConstantPressureVolumeQuad eleTag? iNode? jNode? kNode? lNode? thk? type? matTag?\n";
        return 0;
    }

    // EnhancedQuadId, iNode, jNode, kNode, lNode
    int data[5];
    int num = 5;
    if (OPS_GetIntInput(&num,data) < 0) {
        opserr<<"WARNING: invalid integer input\n";
        return 0;
    }

    double thk = 1.0;
    num = 1;
    if (OPS_GetDoubleInput(&num, &thk) < 0) {
        opserr << "WARNING: invalid double inputs\n";
        return 0;
    }

    const char* type = OPS_GetString();

    int matTag;
    num = 1;
    if (OPS_GetIntInput(&num,&matTag) < 0) {
        opserr<<"WARNING: invalid matTag\n";
        return 0;
    }
    NDMaterial* mat = OPS_getNDMaterial(matTag);
    if (mat == 0) {
        opserr << "WARNING material not found\n";
        opserr << "Material: " << matTag;
        opserr << "\nConstantPressureVolumeQuad element: " << data[0] << endln;
        return 0;
    }

    return new EnhancedQuad(data[0],data[1],data[2],data[3],data[4],
                            *mat,type,thk);
}
#endif


// static data
double  EnhancedQuad::xl[2][4] ;
Matrix  EnhancedQuad::stiff(8,8) ;
Vector  EnhancedQuad::resid(8) ;
Matrix  EnhancedQuad::mass(8,8) ;

double  EnhancedQuad::stressData[3][4] ;
double  EnhancedQuad::tangentData[3][3][4] ;
  

// null constructor
EnhancedQuad::EnhancedQuad( ) :
Element( 0, ELE_TAG_EnhancedQuad ),
connectedExternalNodes(4),
alpha(4), thickness(0.0), load(0), Ki(0)
{ 
  for ( int i = 0 ;  i < 4; i++ ) {
    materialPointers[i] = nullptr;
  }

  // zero enhanced parameters
  alpha.Zero( ) ;
}

// full constructor
EnhancedQuad::EnhancedQuad(int tag, 
                           int node1,
                           int node2,
                           int node3,
                           int node4,
                           NDMaterial &theMaterial,
                           const char *type,
                           double t) 
 :
  Element( tag, ELE_TAG_EnhancedQuad ),
  connectedExternalNodes(4),
  alpha(4), thickness(t), load(0), Ki(0)
{

  connectedExternalNodes(0) = node1 ;
  connectedExternalNodes(1) = node2 ;
  connectedExternalNodes(2) = node3 ;
  connectedExternalNodes(3) = node4 ;

  if (strcmp(type,"PlaneStrain") != 0 && strcmp(type,"PlaneStress") != 0
      && strcmp(type,"PlaneStrain2D") != 0 && strcmp(type,"PlaneStress2D") != 0) {
    opserr << "EnhancedQuad::EnhancedQuad -- improper material type " << type << " for EnhancedQuad\n";
    exit(-1);
  }

  for (int i = 0 ;  i < 4; i++ ) {

      materialPointers[i] = theMaterial.getCopy(type) ;

      if (materialPointers[i] == 0) {
        opserr << "EnhancedQuad::EnhancedQuad -- failed to get a material of type " << type << endln;
        exit(-1);
      }
      
  } // end for i 

  // zero enhanced parameters
  alpha.Zero( ) ;

}

// destructor 
EnhancedQuad::~EnhancedQuad( )
{
  for (int i = 0 ;  i < 4; i++ ) 
    if (materialPointers[i] != nullptr)
      delete materialPointers[i] ;

  if (load != nullptr)
    delete load;

  if (Ki != nullptr)
    delete Ki;
}

// set domain
void  EnhancedQuad::setDomain( Domain *theDomain ) 
{  
  // node pointers
  for (int i = 0; i < 4; i++ ) 
    nodePointers[i] = theDomain->getNode( connectedExternalNodes(i) ) ;

  this->DomainComponent::setDomain(theDomain) ;
}

// get the number of external nodes
int  EnhancedQuad::getNumExternalNodes( ) const
{
  return 4 ;
} 

// return connected external nodes
const ID&  EnhancedQuad::getExternalNodes( ) 
{
  return connectedExternalNodes ;
} 

Node **
EnhancedQuad::getNodePtrs(void) 
{
  return nodePointers;
} 

// return number of dofs
int  EnhancedQuad::getNumDOF( ) 
{
  return 8 ;
}

// commit state
int  EnhancedQuad::commitState( )
{
  int success = 0 ;

  // call element commitState to do any base class stuff
  if ((success = this->Element::commitState()) != 0) {
    opserr << "EnhancedQuad::commitState () - failed in base class";
  }    

  for (int i = 0; i < 4; i++ ) 
    success += materialPointers[i]->commitState( ) ;
  
  return success ;
}

// revert to last commit 
int  EnhancedQuad::revertToLastCommit( ) 
{
  int success = 0 ;

  for (int i = 0; i < nip; i++ ) 
    success += materialPointers[i]->revertToLastCommit( ) ;
  
  return success ;
}

// revert to start 
int  EnhancedQuad::revertToStart( ) 
{
  int success = 0 ;

  // zero enhanced parameters
  this->alpha.Zero( ) ;

  for (int i = 0; i < 4; i++ ) 
    success += materialPointers[i]->revertToStart( ) ;
  
  return success ;
}

// print out element data
void  EnhancedQuad::Print( OPS_Stream &s, int flag )
{
    if (flag == OPS_PRINT_CURRENTSTATE) {
        s << endln;
        s << "Enhanced Strain Four Node Quad \n";
        s << "Element Number: " << this->getTag() << endln;
        s << "Node 1 : " << connectedExternalNodes(0) << endln;
        s << "Node 2 : " << connectedExternalNodes(1) << endln;
        s << "Node 3 : " << connectedExternalNodes(2) << endln;
        s << "Node 4 : " << connectedExternalNodes(3) << endln;
        s << "thickness : " << thickness << endln;
        s << "Material Information : \n ";
        materialPointers[0]->Print(s, flag);
        s << endln;
    }
    
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"EnhancedQuad\", ";
        s << "\"nodes\": [" << connectedExternalNodes(0) << ", ";
        s << connectedExternalNodes(1) << ", ";
        s << connectedExternalNodes(2) << ", ";
        s << connectedExternalNodes(3) << "], ";
        s << "\"thickness\": " << thickness << ", ";
        s << "\"material\": \"" << materialPointers[0]->getTag() << "\"}";
    }
}

// return stiffness matrix 
const Matrix&  EnhancedQuad::getTangentStiff( ) 
{
  int tang_flag = 1 ; // get the tangent 

  // do tangent and residual here
  formResidAndTangent( tang_flag ) ;  

  return stiff ;
}    

// return secant matrix 
const Matrix&  EnhancedQuad::getInitialStiff( ) 
{

  if (Ki != 0)
    return *Ki;

  static const int ndm = 2 ;
  static const int ndf = 2 ; 
  static const int nstress = 3 ; 
  static const int numberNodes = 4 ;
  static const int nip = 4 ;
  static const int nShape = 3 ;
  static const int nEnhanced = 4 ; 
  static const int nModes = 2 ;
  static const int numberDOF = 8 ;


  int i, j, k, p, q ;
  int jj, kk ;

  static double xsj[nip] ;  // determinant jacaobian matrix 
  static double dvol[nip] ; // volume element
  static Vector strain(nstress) ;  // strain
  static double shp[nShape][numberNodes] ;  // shape functions at a gauss point

  static double Shape[nShape][numberNodes][nip] ; // all the shape functions

  static Vector residJ(ndf) ; // nodeJ residual 
  static Matrix stiffJK(ndf,ndf) ; // nodeJK stiffness 
  static Matrix stiffKJ(ndf,ndf) ; // nodeKJ stiffness
  static Vector stress(nstress) ;  // stress
  static Matrix dd(nstress,nstress) ;  // material tangent
  static Matrix J0(ndm,ndm) ; //Jacobian matrix at center of element
  static Matrix J0inv(ndm,ndm) ; // inverse of above

  static Matrix Kee(nEnhanced,nEnhanced) ;

  static Vector residE(nEnhanced) ;

  static Vector Umode(ndf) ;

  static Vector dalpha(nEnhanced) ;

  static Matrix Kue(numberDOF,nEnhanced) ;

  static Matrix Keu(nEnhanced,numberDOF) ;

  static Matrix KeeInvKeu(nEnhanced,numberDOF) ;


  // ---------B-matrices------------------------------------
    MatrixND<nstress, ndf> B[numberNodes];
    static Matrix BJ(nstress,ndf) ;      // B matrix node J
//  static Matrix BJtran(ndf,nstress) ;
    static Matrix BK(nstress,ndf) ;      // B matrix node k
    static Matrix BKtran(ndf,nstress) ;
    static Matrix BJtranD(ndf,nstress) ;
    static Matrix BKtranD(ndf,nstress) ;
  // -------------------------------------------------------

  
  // zero stiffness and residual 
  stiff.Zero( ) ;

  Kee.Zero( ) ;
  residE.Zero( ) ;

  Kue.Zero( ) ;
  Keu.Zero( ) ;


  // compute Jacobian and inverse at center
  double L1 = 0.0 ;
  double L2 = 0.0 ;
  computeJacobian( L1, L2, xl, J0, J0inv ) ; 

  // gauss loop to compute and save shape functions 
  double det ;
  for (int i = 0; i < nip; i++ ) {

    // get shape functions    
    shape2d( pts[i][0], pts[i][1], xl, shp, det ) ;

    // save shape functions
    for (int p = 0; p < nShape; p++ ) {
       for ( q = 0; q < numberNodes; q++ )
          Shape[p][q][i] = shp[p][q] ;
    }

    // save jacobian determinant
    xsj[i] = det ;

    // volume element to also be saved
    dvol[i] = wts[i] * det * thickness;  

  } // end for i gauss loop

  Kee.Zero( ) ;
  // gauss loop
  for (int i = 0; i < nip; i++ ) {
    // tangent 
    dd = materialPointers[i]->getInitialTangent( ) ;
    
    // multiply by volume element
    dd *= dvol[i] ;
    saveData( i, stress, dd ) ;   
  }
  Kee.Zero( ) ;
  for (int i = 0; i < nip; i++ ) {
    // extract shape functions from saved array
    for (int p = 0; p < nShape; p++ ) {
      for (int q = 0; q < numberNodes; q++ )
        shp[p][q]  = Shape[p][q][i] ;
    }
    
    // enhanced residual and tangent calculations loops
    jj = 0 ;
    for (int j = 0; j < nModes; j++ ) {
      
      // compute B matrix 
      BJ = computeBenhanced( j, pts[i][0], pts[i][1], xsj[i], J0inv ) ; 
      
      // transpose 
//    BJtran = transpose( BJ ) ;

      //BJtranD = BJtran * dd ;
      BJtranD.addMatrixTransposeProduct(0.0, BJ, dd, 1.0) ;
      
      kk = 0 ;
      for (int k = 0; k < nModes; k++ ) {

        BK = computeBenhanced( k, pts[i][0], pts[i][1], xsj[i], J0inv ) ;

        // stiffJK =  BJtranD * BK  ;
        stiffJK.addMatrixProduct(0.0, BJtranD,BK,1.0) ;

        for (int p = 0; p < ndf; p++ )  {
          for (int q = 0; q < ndf; q++ )
            Kee( jj+p, kk+q ) += stiffJK( p, q ) ;
        }
        
        kk += ndf ;
      }
      
      jj += ndf ;
    } // end for j loop
  } // end for i gauss loop 
  
  /// HERE
  for (int i = 0; i < nip; i++ ) {
    // extract shape functions from saved array
    for ( p = 0; p < nShape; p++ ) {
       for ( q = 0; q < numberNodes; q++ )
          shp[p][q]  = Shape[p][q][i] ;
    }


    // recover stress and tangent from saved data
    getData( i, stress, dd ) ;
    jj = 0 ;
    for (int j = 0; j < numberNodes; j++ ) {

      BJ = computeB( j, shp, B[j]) ;
   
      // transpose 
//    BJtran = B[j].transpose(); // transpose( BJ ) ;

      //BJtranD = BJtran * dd ;
      BJtranD.addMatrixTransposeProduct(0.0, BJ, dd, 1.0) ;
      
      // node-node stiffness
      kk = 0 ;
      for (int k = 0; k < numberNodes; k++ ) {
        
        BK = computeB( k, shp, B[k]) ;
        
        // stiffJK =  BJtranD * BK  ;
        stiffJK.addMatrixProduct(0.0, BJtranD, B[k],1.0) ;
        
        for ( p = 0; p < ndf; p++ )  {
          for ( q = 0; q < ndf; q++ )
            stiff( jj+p, kk+q ) += stiffJK( p, q ) ;
        } // end for p
        
        kk += ndf ;
      } // end for k loop
      

      // node-enhanced stiffness Kue 
      kk = 0 ;
      for ( k = 0; k < nModes; k++ ) {
        
        BK = computeBenhanced( k, pts[i][0], pts[i][1], xsj[i], J0inv ) ;
        
        // stiffJK =  BJtranD * BK  ;
        stiffJK.addMatrixProduct(0.0, BJtranD,BK,1.0) ;
        
        for ( p = 0; p < ndf; p++ )  {
          for ( q = 0; q < ndf; q++ )
            Kue( jj+p, kk+q ) += stiffJK( p, q ) ;
        } // end for p 
        
        kk += ndf ;
      } // end for k loop

      // enhanced-node stiffness Keu 
      kk = 0 ;
      for ( k = 0; k < nModes; k++ ) {
        
        BK = computeBenhanced( k, pts[i][0], pts[i][1], xsj[i], J0inv ) ;
        
        // transpose 
        BKtran = transpose( BK ) ;
        
        //BKtranD = BKtran * dd ;
        BKtranD.addMatrixProduct(0.0, BKtran,dd,1.0 ) ;
        
        // stiffKJ = (BKtran*dd)*BJ ;
        stiffKJ.addMatrixProduct(0.0, BKtranD,BJ,1.0) ;
        
        for (int p = 0; p < ndf; p++ )  {
          for (int q = 0; q < ndf; q++ )
            Keu( kk+p, jj+q ) += stiffKJ( p, q ) ;
        } // end for p  
        
        kk += ndf ;
      } // end for k loop
    
      jj += ndf ;
    } // end for j loop
  
  
  } // end for i gauss loop 


  Kee.Solve( Keu, KeeInvKeu ) ;

  // stiff -= ( Kue * KeeInvKeu ) ;
  stiff.addMatrixProduct(1.0,  Kue, KeeInvKeu, -1.0);

  Ki = new Matrix(stiff);

  return stiff;
}

// return mass matrix
const Matrix&  EnhancedQuad::getMass( ) 
{

  int tangFlag = 1 ;

  formInertiaTerms( tangFlag ) ;

  return mass ;

} 

void  EnhancedQuad::zeroLoad( )
{
  if (load != 0)
    load->Zero();

  return ;
}

int 
EnhancedQuad::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  opserr << "EnhancedQuad::addLoad - load type unknown for ele with tag: " << this->getTag() << endln;
  return -1;
}

int
EnhancedQuad::addInertiaLoadToUnbalance(const Vector &accel)
{
  static const int nip = 4 ;
  static const int numberNodes = 4 ;
  static const int ndf = 2 ; 

  int i;

  // check to see if have mass
  int haveRho = 0;
  for (i = 0; i < nip; i++) {
    if (materialPointers[i]->getRho() != 0.0)
      haveRho = 1;
  }

  if (haveRho == 0)
    return 0;

  // Compute mass matrix
  int tangFlag = 1 ;
  formInertiaTerms( tangFlag ) ;

  // store computed RV for nodes in resid vector
  int count = 0;
  for (i=0; i<numberNodes; i++) {
    const Vector &Raccel = nodePointers[i]->getRV(accel);
    for (int j=0; j<ndf; j++)
      resid(count++) = Raccel(i);
  }

  // create the load vector if one does not exist
  if (load == 0) 
    load = new Vector(numberNodes*ndf);

  // add -M * RV(accel) to the load vector
  load->addMatrixVector(1.0, mass, resid, -1.0);
  
  return 0;
}

// get residual
const Vector&  EnhancedQuad::getResistingForce( ) 
{
  int tang_flag = 0 ; // don't get the tangent

  formResidAndTangent( tang_flag ) ;

  // subtract external loads 
  if (load != 0)
    resid -= *load;

  return resid ;   
}

// get residual with inertia terms
const Vector&  EnhancedQuad::getResistingForceIncInertia( )
{
  int tang_flag = 0 ; // don't get the tangent

  static Vector res(8);

  // do tangent and residual here 
  formResidAndTangent( tang_flag ) ;

  // inertia terms
  formInertiaTerms( tang_flag ) ;

  res = resid;

  // add the damping forces if rayleigh damping
  if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
    res += this->getRayleighDampingForces();

  // subtract external loads 
  if (load != 0)
    res -= *load;

  return res;
}

//*********************************************************************
// form inertia terms

void   EnhancedQuad::formInertiaTerms( int tangFlag ) 
{

  static const int ndf = 2 ; 
  static const int numberNodes = 4 ;
  static const int nip = 4 ;
  static const int nShape = 3 ;
  static const int massIndex = nShape - 1 ;

  double xsj ;  // determinant jacaobian matrix 

  double dvol ; // volume element

  static double shp[nShape][numberNodes] ;  // shape functions at a gauss point

  static Vector momentum(ndf) ;

  int i, j, k, p ;
  int jj, kk ;

  double temp, rho, massJK ;

  // zero mass 
  mass.Zero( ) ;

  // gauss loop 
  for (int i = 0; i < nip; i++ ) {

    // get shape functions    
    shape2d( pts[i][0], pts[i][1], xl, shp, xsj ) ;
    
    // volume element
    dvol = wts[i] * xsj * thickness;

    // node loop to compute acceleration
    momentum.Zero( ) ;
    for (int j = 0; j < numberNodes; j++ ) 
      // momentum += shp[massIndex][j] * ( nodePointers[j]->getTrialAccel()  ) ; 
      momentum.addVector( 1.0,
                          nodePointers[j]->getTrialAccel(),
                          shp[massIndex][j] ) ;

    // density
    rho = materialPointers[i]->getRho() ;

    // multiply acceleration by density to form momentum
    momentum *= rho ;


    // residual and tangent calculations node loops
    jj = 0 ;
    for (int j = 0; j < numberNodes; j++ ) {

      temp = shp[massIndex][j] * dvol ;

      for (int p = 0; p < ndf; p++ )
        resid( jj+p ) += ( temp * momentum(p) )  ;

      
      if ( tangFlag == 1 ) {

         // multiply by density
         temp *= rho ;

         // node-node mass
         kk = 0 ;
         for (int k = 0; k < numberNodes; k++ ) {

            massJK = temp * shp[massIndex][k] ;

            for (int p = 0; p < ndf; p++ )  
              mass( jj+p, kk+p ) += massJK ;
            
            kk += ndf ;
          }

      } // end if tang_flag 

      jj += ndf ;
    } // end for j loop
  } // end for i gauss loop 
}

//*********************************************************************
// form residual and tangent
void  EnhancedQuad::formResidAndTangent( int tang_flag ) 
{

  static constexpr double tolerance = 1.0e-08 ;
  static constexpr int nIterations = 10 ;

  static constexpr int ndm = 2 ;
  static constexpr int ndf     = 2 ; 
  static constexpr int nstress = 3 ; 
  static constexpr int numberNodes = 4 ;
  static constexpr int nip    = 4 ;
  static constexpr int nShape = 3 ;
  static constexpr int nEnhanced = 4 ; 
  static constexpr int nModes = 2 ;
  static constexpr int numberDOF = 8 ;

  int i, j, k, p, q ;
  int jj, kk ;
  int success ;

  OPS_STATIC double xsj[nip] ;  // determinant jacaobian matrix 
  OPS_STATIC double dvol[nip] ; // volume element
  OPS_STATIC double Shape[nip][nShape][numberNodes]; // [nip] ; // all the shape functions

  static VectorND<nstress> strain ;  // strain
  static Vector residJ(ndf) ; // nodeJ residual 
  static Matrix stiffJK(ndf,ndf) ; // nodeJK stiffness 
  static Matrix stiffKJ(ndf,ndf) ; // nodeKJ stiffness
  static Vector stress(nstress) ;  // stress
  static Matrix dd(nstress,nstress) ;  // material tangent

  static Matrix J0(ndm,ndm) ; //Jacobian matrix at center of element
  static Matrix J0inv(ndm,ndm) ; // inverse of above

  static Matrix Kee(nEnhanced,nEnhanced) ;

  static Vector residE(nEnhanced);
  static Vector Umode(ndf);
  static Vector dalpha(nEnhanced);
  static Matrix Kue(numberDOF,nEnhanced) ;
  static Matrix Keu(nEnhanced,numberDOF) ;
  static Matrix KeeInvKeu(nEnhanced,numberDOF) ;

  // ---------B-matrices------------------------------------
  static MatrixND<nstress, ndf> B[numberNodes];
  static Matrix BJ(nstress,ndf) ;      // B matrix node J
  static Matrix BK(nstress,ndf) ;      // B matrix node k
//static Matrix BJtran(ndf,nstress) ;
//static Matrix BKtran(ndf,nstress) ;
  static Matrix BJtranD(ndf,nstress) ;
  static Matrix BKtranD(ndf,nstress) ;
  // -------------------------------------------------------
  
  // zero stiffness and residual 
  stiff.Zero( );
  resid.Zero( );

  Kee.Zero( );
  residE.Zero( );

  Kue.Zero( );
  Keu.Zero( );

  // compute Jacobian and inverse at center
  double L1 = 0.0 ;
  double L2 = 0.0 ;
  computeJacobian( L1, L2, xl, J0, J0inv ) ; 

  // gauss loop to compute and save shape functions 
  double det ;
  for (int i = 0; i < nip; i++ ) {

    // get shape functions    
    shape2d( pts[i][0], pts[i][1], xl, Shape[i], det ) ;

    // save jacobian determinant
    xsj[i] = det ;

    // volume element to also be saved
    dvol[i] = wts[i] * det * thickness;  

  }

  // -------------------------------------------------------------------
  // Newton loop to solve for enhanced strain parameters

  int count = 0 ;
  do {

    residE.Zero( ) ;
    Kee.Zero( ) ;

    // gauss loop
    for (int i = 0; i < nip; i++ ) {

      // zero the strain
      strain.zero( ) ;

      // j-node loop to compute nodal strain contributions
      for (int j = 0; j < numberNodes; j++ )  {

        // compute B matrix 
        computeB( j, Shape[i], B[j]) ;
      
        // nodal displacements 
        const Vector &ul = nodePointers[j]->getTrialDisp( ) ;

        // compute the strain
        // strain += (BJ*ul) ; 
        strain.addMatrixVector(1.0, B[j], ul, 1.0) ;

      }

      // j-node loop to compute enhanced strain contributions
      for (int j = 0; j < nModes; j++ )  {

        // compute B matrix 
        BJ = computeBenhanced( j, pts[i][0], pts[i][1], xsj[i], J0inv ) ; 
      
        // enhanced "displacements" 
        Umode(0) = this->alpha( 2*j     ) ;
        Umode(1) = this->alpha( 2*j + 1 ) ;

        // compute the strain
        // strain += (BJ*Umode) ; 
        strain.addMatrixVector(1.0, BJ, Umode, 1.0) ;

      }

      success = materialPointers[i]->setTrialStrain( strain ) ;

      // compute the stress
      stress = materialPointers[i]->getStress( ) ;

      // multiply by volume element
      stress  *= dvol[i] ;

      // tangent 
      dd = materialPointers[i]->getTangent( ) ;

      // multiply by volume element
      dd *= dvol[i] ;

      // save stress and tangent (already multiplied by volume element)
      saveData( i, stress, dd ) ; 

      // enhanced residual and tangent calculations loops

      jj = 0 ;
      for (int j = 0; j < nModes; j++ ) {

        // compute B matrix 
        BJ = computeBenhanced( j, pts[i][0], pts[i][1], xsj[i], J0inv ) ; 

        // residual
        // residJ = -BJ^T * stress) ;
        residJ.addMatrixTransposeVector(0.0, BJ, stress,-1.0) ;

        // residual 
        for (int p = 0; p < ndf; p++ )
          residE( jj+p ) += residJ(p)  ;


        //BJtranD = BJtran * dd ;
        BJtranD.addMatrixTransposeProduct(0.0, BJ, dd, 1.0) ;

        kk = 0 ;
        for (int k = 0; k < nModes; k++ ) {

          BK = computeBenhanced( k, pts[i][0], pts[i][1], xsj[i], J0inv ) ;
  
          // stiffJK =  BJtranD * BK  ;
          stiffJK.addMatrixProduct(0.0, BJtranD, BK, 1.0) ;

          for (int p = 0; p < ndf; p++ )  {
             for (int q = 0; q < ndf; q++ )
                Kee( jj+p, kk+q ) += stiffJK( p, q ) ;
          }

          kk += ndf ;
        }

        jj += ndf ;
      } // end for j loop

    } // end for i gauss loop 


    // solve for enhanced strain parameters
    Kee.Solve( residE, dalpha ) ;

    if (dalpha(0) > 1.0e10)  opserr << "dalpha: " << residE << dalpha;

    this->alpha += dalpha ;

    count++ ;
    if ( count > nIterations ) {
      opserr << "Exceeded " << nIterations
             << " iterations solving for enhanced strain parameters "
             << endln ;
      break ;
    }

    // do at least 2 iterations so saved data is good
  } while ( residE.Norm() > tolerance  ||  count < 2 ) ;


  // end enhanced strain parameters newton loop
  // -------------------------------------------------------------------

  // Gauss loop 
  for (int i = 0; i < nip; i++ ) {

    // recover stress and tangent from saved data
    getData( i, stress, dd ) ;

    // residual and tangent calculations node loops

    jj = 0 ;
    for (int j = 0; j < numberNodes; j++ ) {

      // TODO(cmp): check this
      BJ = computeB( j, Shape[i], B[j]) ;

      // residual
      // residJ = BJtran * stress ;
      residJ.addMatrixTransposeVector(0.0, B[j], stress, 1.0) ;

      for (int p = 0; p < ndf; p++ )
        resid( jj+p ) += residJ(p)  ;


      if ( tang_flag == 1 ) {

        //BJtranD = BJtran * dd ;
        BJtranD.addMatrixTransposeProduct(0.0, B[j], dd, 1.0) ;

         // node-node stiffness
         kk = 0 ;
         for (int k = 0; k < numberNodes; k++ ) {

            computeB( k, Shape[i], B[k]) ;
  
            // stiffJK =  BJtranD * BK  ;
            stiffJK.addMatrixProduct(0.0, BJtranD, B[k], 1.0) ;

            for (int p = 0; p < ndf; p++ )  {
               for (int q = 0; q < ndf; q++ )
                  stiff( jj+p, kk+q ) += stiffJK( p, q ) ;
            }

            kk += ndf ;
          }

         // node-enhanced stiffness Kue 
         kk = 0 ;
         for (int k = 0; k < nModes; k++ ) {

            BK = computeBenhanced( k, pts[i][0], pts[i][1], xsj[i], J0inv ) ;
  
            // stiffJK =  BJtranD * BK  ;
            stiffJK.addMatrixProduct(0.0, BJtranD,BK,1.0) ;
           
            for (int p = 0; p < ndf; p++ ) {
               for (int q = 0; q < ndf; q++ )
                  Kue( jj+p, kk+q ) += stiffJK( p, q ) ;
            }

            kk += ndf ;
          }

         // enhanced-node stiffness Keu 
         kk = 0 ;
         for (int k = 0; k < nModes; k++ ) {

            BK = computeBenhanced( k, pts[i][0], pts[i][1], xsj[i], J0inv ) ;
            
            //BKtranD = BKtran * dd ;
            BKtranD.addMatrixTransposeProduct(0.0, BK, dd, 1.0 ) ;

            // stiffKJ = (BKtran*dd)*BJ ;
            stiffKJ.addMatrixProduct(0.0, BKtranD,BJ,1.0) ;

            for (int p = 0; p < ndf; p++ )  {
               for (int q = 0; q < ndf; q++ )
                  Keu( kk+p, jj+q ) += stiffKJ( p, q ) ;
            }  

            kk += ndf ;
          }

      } // end if tang_flag 

      jj += ndf ;
    } // end for j loop


  } // end for i gauss loop 


  // static condensation of enhanced parameters

  if ( tang_flag == 1 ) {  
     Kee.Solve( Keu, KeeInvKeu ) ;

     // stiff -= ( Kue * KeeInvKeu ) ;
     stiff.addMatrixProduct(1.0,  Kue,KeeInvKeu, -1.0 ) ;
  }

  return ;
}

int  
EnhancedQuad::update(void) 
{
  // compute basis vectors and local nodal coordinates
  computeBasis( ) ;

  return 0;
}
  
//************************************************************************
void   
EnhancedQuad::saveData(int gp, 
                       const Vector &stress,
                       const Matrix &tangent ) 
{
  // save stress
  for (int i=0; i<3; i++ )
    stressData[i][gp] = stress(i) ;

  // save tangent
  for (int i=0; i<3; i++ ) {
    for (int j=0; j<3; j++ ) 
      tangentData[i][j][gp] = tangent(i,j) ;
  }

  return ;
}

//************************************************************************
// recover stress and tangent data

void  EnhancedQuad::getData( int gp,
                             Vector &stress,
                             Matrix &tangent ) 
{

  // get stress
  for (int i=0; i<3; i++ )
    stress(i) = stressData[i][gp] ;

  // get tangent
  for (int i=0; i<3; i++ ) {
    for (int j=0; j<3; j++ ) 
      tangent(i,j) = tangentData[i][j][gp] ;
  }

  return ;
}

//************************************************************************
// compute local coordinates and basis

void EnhancedQuad::computeBasis( ) 
{
  // nodal coordinates 

  for (int i = 0; i < 4; i++ ) {
     const Vector &coorI = nodePointers[i]->getCrds( ) ;

     xl[0][i] = coorI(0);
     xl[1][i] = coorI(1);
  }
}

//*************************************************************************
// compute B matrix

const Matrix
EnhancedQuad::computeB( int node, const double shp[3][4] , MatrixND<3,2> & B)
{
// ---B Matrix in standard {1,2,3} mechanics notation---------------
//
//                -             -
//               | +N,1      0   | 
// B    =        |   0     +N,2  |    (3x2)
//               | +N,2    +N,1  |
//                -             -  
//
// -------------------------------------------------------------------

// static Matrix B(3,2) ;
  B.zero( ) ;

  B(0,0) = shp[0][node];
  B(1,1) = shp[1][node];
  B(2,0) = shp[1][node];
  B(2,1) = shp[0][node];

  return B ;

}

//***********************************************************************
// compute enhanced strain B-matrices

const Matrix&   
EnhancedQuad::computeBenhanced(int node, 
                               double L1,
                               double L2,
                               double j, 
                               const Matrix &Jinv )
{
  static Matrix B(3,2) ;
  static double JinvTran[2][2] ;
  static double shape[2] ;
  static double parameter ;


  // compute JinvTran
  JinvTran[0][0] = Jinv(0,0) ;
  JinvTran[1][1] = Jinv(1,1) ;
  JinvTran[0][1] = Jinv(1,0) ;
  JinvTran[1][0] = Jinv(0,1) ;      // residual 

  if ( node == 0 ) {

    // first column of JinvTran 
    shape[0] = JinvTran[0][0] ;
    shape[1] = JinvTran[1][0] ;

    parameter = L1 / j ;

  }
  else if ( node == 1 ) {
    // second column of JinvTran 
    shape[0] = JinvTran[0][1] ;
    shape[1] = JinvTran[1][1] ;

    parameter = L2 / j ;
  } 

  shape[0] *= parameter ;
  shape[1] *= parameter ; 


  B.Zero( ) ;

  B(0,0) = shape[0] ;
  B(1,1) = shape[1] ;
  B(2,0) = shape[1] ;
  B(2,1) = shape[0] ;
  
  return B ;

}

//***********************************************************************
// compute Jacobian matrix and inverse at point {L1,L2} 

void EnhancedQuad::computeJacobian(double L1, double L2,
                                   const double x[2][4], 
                                   Matrix &JJ, 
                                   Matrix &JJinv )
{
  int i, j, k ;
     
  static constexpr double s[] = { -0.5,  0.5, 0.5, -0.5 } ;
  static constexpr double t[] = { -0.5, -0.5, 0.5,  0.5 } ;

  static double shp[2][4] ;

  double ss = L1;
  double tt = L2;

  for (int i = 0; i < 4; i++ ) {
      shp[0][i] = s[i] * ( 0.5 + t[i]*tt ) ;
      shp[1][i] = t[i] * ( 0.5 + s[i]*ss ) ;
  }
  
  // Construct jacobian and its inverse
  
  JJ.Zero( ) ;
  for (int i = 0; i < 2; i++ ) {
    for (int j = 0; j < 2; j++ ) {

      for (int k = 0; k < 4; k++ )
          JJ(i,j) +=  x[i][k] * shp[j][k] ;

    }
  }

  double xsj = JJ(0,0)*JJ(1,1) - JJ(0,1)*JJ(1,0) ;

  // inverse jacobian
  double jinv = 1.0 / xsj ;
  JJinv(0,0) =  JJ(1,1) * jinv ;
  JJinv(1,1) =  JJ(0,0) * jinv ;
  JJinv(0,1) = -JJ(0,1) * jinv ;
  JJinv(1,0) = -JJ(1,0) * jinv ;

  return ;

}

//************************************************************************
// shape function routine for four node quads

void  EnhancedQuad::shape2d(double ss, double tt, 
                            const double x[2][4], 
                            double shp[3][4], 
                            double &xsj            )
{ 

  int i, j, k ;

  double temp ;
     
  static constexpr double s[] = { -0.5,  0.5, 0.5, -0.5 } ;
  static constexpr double t[] = { -0.5, -0.5, 0.5,  0.5 } ;

  static Matrix xs(2,2) ;
  static Matrix sx(2,2) ;

  for (int i = 0; i < 4; i++ ) {
      shp[2][i] = ( 0.5 + s[i]*ss )*( 0.5 + t[i]*tt ) ;
      shp[0][i] = s[i] * ( 0.5 + t[i]*tt ) ;
      shp[1][i] = t[i] * ( 0.5 + s[i]*ss ) ;
  }

  
  // Construct jacobian and its inverse
  
  for (int i = 0; i < 2; i++ ) {
    for (int j = 0; j < 2; j++ ) {
      xs(i,j) = 0.0;
      for ( k = 0; k < 4; k++ )
          xs(i,j) +=  x[i][k] * shp[j][k] ;

    }
  }

  xsj = xs(0,0)*xs(1,1) - xs(0,1)*xs(1,0) ;

  // inverse jacobian
  // inverse jacobian
  double jinv = 1.0 / xsj ;
  sx(0,0) =  xs(1,1) * jinv ;
  sx(1,1) =  xs(0,0) * jinv ;
  sx(0,1) = -xs(0,1) * jinv ;
  sx(1,0) = -xs(1,0) * jinv ;


  // form global derivatives 
  
  for ( i = 0; i < 4; i++ ) {
    temp      = shp[0][i]*sx(0,0) + shp[1][i]*sx(1,0) ;
    shp[1][i] = shp[0][i]*sx(0,1) + shp[1][i]*sx(1,1) ;
    shp[0][i] = temp ;
  }

  return ;
}

const Matrix&
EnhancedQuad::transpose( const Matrix &M ) 
{
  // we're always transposing 3x2 matrices for this element,
  // so always return a 2x3 .

  static constexpr int dim1 = 2 ;
  static constexpr int dim2 = 3 ;
  static Matrix Mtran(dim1,dim2) ;

  for (int i = 0; i < dim1; i++ ) {
     for (int j = 0; j < dim2; j++ ) 
         Mtran(i,j) = M(j,i) ;
  }

  return Mtran ;
}

Response*
EnhancedQuad::setResponse(const char **argv, int argc, 
                          OPS_Stream &output)
{
  Response *theResponse = nullptr;

  output.tag("ElementOutput");
  output.attr("eleType","EnhancedQuad");
  output.attr("eleTag",this->getTag());
  output.attr("node1",connectedExternalNodes[0]);
  output.attr("node2",connectedExternalNodes[1]);
  output.attr("node3",connectedExternalNodes[2]);
  output.attr("node4",connectedExternalNodes[3]);

  char dataOut[10];
  if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0) {
    
    for (int i=1; i<=4; i++) {
      sprintf(dataOut,"P1_%d",i);
      output.tag("ResponseType",dataOut);
      sprintf(dataOut,"P2_%d",i);
      output.tag("ResponseType",dataOut);
    }
    
    theResponse =  new ElementResponse(this, 1, resid);
  }  else if (strcmp(argv[0],"material") == 0 || strcmp(argv[0],"integrPoint") == 0) {
    int pointNum = atoi(argv[1]);
    if (pointNum > 0 && pointNum <= 4) {

      output.tag("GaussPoint");
      output.attr("number",pointNum);
      output.attr("eta",   pts[pointNum-1][0]);
      output.attr("neta",  pts[pointNum-1][1]);

      theResponse =  materialPointers[pointNum-1]->setResponse(&argv[2], argc-2, output);
      
      output.endTag();

    } 
  }
  else if ((strcmp(argv[0],"stress") == 0) || (strcmp(argv[0],"stresses") == 0)) {

      for (int i=0; i<4; i++) {
        output.tag("GaussPoint");
        output.attr("number",i+1);
        output.attr("eta",pts[i][0]);
        output.attr("neta",pts[i][1]);

        output.tag("NdMaterialOutput");
        output.attr("classType", materialPointers[i]->getClassTag());
        output.attr("tag", materialPointers[i]->getTag());

        output.tag("ResponseType","sigma11");
        output.tag("ResponseType","sigma22");
        output.tag("ResponseType","sigma12");

        output.endTag(); // GaussPoint
        output.endTag(); // NdMaterialOutput
      }

      theResponse =  new ElementResponse(this, 3, Vector(12));
  }
  
  else if ((strcmp(argv[0],"strain") == 0) ||(strcmp(argv[0],"strains") == 0)) {

      for (int i=0; i<4; i++) {
        output.tag("GaussPoint");
        output.attr("number",i+1);
        output.attr("eta",pts[i][0]);
        output.attr("neta",pts[i][1]);

        output.tag("NdMaterialOutput");
        output.attr("classType", materialPointers[i]->getClassTag());
        output.attr("tag", materialPointers[i]->getTag());

        output.tag("ResponseType","eta11");
        output.tag("ResponseType","eta22");
        output.tag("ResponseType","eta12");

        output.endTag(); // GaussPoint
        output.endTag(); // NdMaterialOutput
      }

      theResponse =  new ElementResponse(this, 4, Vector(12));
  }
        
  output.endTag(); // ElementOutput

  return theResponse;
}

int 
EnhancedQuad::getResponse(int responseID, Information &eleInfo)
{
  if (responseID == 1) {

    return eleInfo.setVector(this->getResistingForce());

  } else if (responseID == 3) {

    // Loop over the integration points
    static Vector stresses(12);
    int cnt = 0;
    for (int i = 0; i < 4; i++) {

      // Get material stress response
      const Vector &sigma = materialPointers[i]->getStress();
      stresses(cnt) = sigma(0);
      stresses(cnt+1) = sigma(1);
      stresses(cnt+2) = sigma(2);
      cnt += 3;
    }
    return eleInfo.setVector(stresses);

  } else if (responseID == 4) {

    // Loop over the integration points
    static Vector stresses(12);
    int cnt = 0;
    for (int i = 0; i < 4; i++) {

      // Get material stress response
      const Vector &sigma = materialPointers[i]->getStrain();
      stresses(cnt) = sigma(0);
      stresses(cnt+1) = sigma(1);
      stresses(cnt+2) = sigma(2);
      cnt += 3;
    }
    return eleInfo.setVector(stresses);
        
  } else

    return -1;
}

int
EnhancedQuad::sendSelf (int commitTag, Channel &theChannel)
{
  int res = 0;
  
  // note: we don't check for dataTag == 0 for Element
  // objects as that is taken care of in a commit by the Domain
  // object - don't want to have to do the check if sending data
  int dataTag = this->getDbTag();
  
  // Quad packs its data into a Vector and sends this to theChannel
  // along with its dbTag and the commitTag passed in the arguments
  static Vector data(6);
  data(0) = this->getTag();
  data(1) = thickness;

  data(2) = alphaM;
  data(3) = betaK;
  data(4) = betaK0;
  data(5) = betaKc;
  
  res += theChannel.sendVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING EnhancedQuad::sendSelf() - " << this->getTag() << " failed to send Vector\n";
    return res;
  }              
  

  // Now quad sends the ids of its materials
  int matDbTag;
  
  static ID idData(12);
  
  int i;
  for (i = 0; i < 4; i++) {
    idData(i) = materialPointers[i]->getClassTag();
    matDbTag = materialPointers[i]->getDbTag();
    // NOTE: we do have to ensure that the material has a database
    // tag if we are sending to a database channel.
    if (matDbTag == 0) {
      matDbTag = theChannel.getDbTag();
                        if (matDbTag != 0)
                materialPointers[i]->setDbTag(matDbTag);
    }
    idData(i+4) = matDbTag;
  }
  
  idData(8) = connectedExternalNodes(0);
  idData(9) = connectedExternalNodes(1);
  idData(10) = connectedExternalNodes(2);
  idData(11) = connectedExternalNodes(3);

  res += theChannel.sendID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING EnhancedQuad::sendSelf() - " << this->getTag() << " failed to send ID\n";
    return res;
  }

  // Finally, quad asks its material objects to send themselves
  for (i = 0; i < 4; i++) {
    res += materialPointers[i]->sendSelf(commitTag, theChannel);
    if (res < 0) {
      opserr << "WARNING EnhancedQuad::sendSelf() - " << this->getTag() << " failed to send its Material\n";
      return res;
    }
  }
  
  return res;
}

int
EnhancedQuad::recvSelf(int commitTag, Channel &theChannel, 
                       FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  int dataTag = this->getDbTag();

  // Quad creates a Vector, receives the Vector and then sets the 
  // internal data with the data in the Vector
  static Vector data(6);
  res += theChannel.recvVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING EnhancedQuad::recvSelf() - failed to receive Vector\n";
    return res;
  }
  
  this->setTag((int)data(0));
  thickness = data(1);

  alphaM = data(2);
  betaK = data(3);
  betaK0 = data(4);
  betaKc = data(5);

  static ID idData(12);
  // Quad now receives the tags of its four external nodes
  res += theChannel.recvID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING EnhancedQuad::recvSelf() - " << this->getTag() << " failed to receive ID\n";
    return res;
  }

  connectedExternalNodes(0) = idData(8);
  connectedExternalNodes(1) = idData(9);
  connectedExternalNodes(2) = idData(10);
  connectedExternalNodes(3) = idData(11);
  
  if (materialPointers[0] == 0) {
    for (int i = 0; i < 4; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i+4);
      // Allocate new material with the sent class tag
      materialPointers[i] = theBroker.getNewNDMaterial(matClassTag);
      if (materialPointers[i] == 0) {
            opserr << "EnhancedQuad::recvSelf() - Broker could not create NDMaterial of class type " << matClassTag << endln;
            return -1;
      }
      // Now receive materials into the newly allocated space
      materialPointers[i]->setDbTag(matDbTag);
      res += materialPointers[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
        opserr << "EnhancedQuad::recvSelf() - material " << i << "failed to recv itself\n";
            return res;
      }
    }
  }

  // materials exist , ensure materials of correct type and recvSelf on them
  else {
    for (int i = 0; i < 4; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i+4);
      // Check that material is of the right type; if not,
      // delete it and create a new one of the right type
      if (materialPointers[i]->getClassTag() != matClassTag) {
        delete materialPointers[i];
        materialPointers[i] = theBroker.getNewNDMaterial(matClassTag);
        if (materialPointers[i] == nullptr) {
          opserr << "EnhancedQuad::recvSelf() - material " << i << "failed to create\n";
          return -1;
        }
      }
      // Receive the material
      materialPointers[i]->setDbTag(matDbTag);
      res += materialPointers[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
        opserr << "EnhancedQuad::recvSelf() - material " << i << "failed to recv itself\n";
        return res;
      }
    }
  }
  
  return res;
}

int
EnhancedQuad::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
  // get the end point display coords
  static Vector v1(3);
  static Vector v2(3);
  static Vector v3(3);
  static Vector v4(3);
  nodePointers[0]->getDisplayCrds(v1, fact, displayMode);
  nodePointers[1]->getDisplayCrds(v2, fact, displayMode);
  nodePointers[2]->getDisplayCrds(v3, fact, displayMode);
  nodePointers[3]->getDisplayCrds(v4, fact, displayMode);

  // place values in coords matrix
  static Matrix coords(4, 3);
  for (int i = 0; i < 3; i++) {
      coords(0, i) = v1(i);
      coords(1, i) = v2(i);
      coords(2, i) = v3(i);
      coords(3, i) = v4(i);
  }

  // set the quantity to be displayed at the nodes;
  // if displayMode is 1 through 3 we will plot material stresses otherwise 0.0
  // until someone projects the stress to the nodes will display the stress 
  // at the guass points at the nodes .. could also just display the average!
  static Vector values(4);
  if (displayMode < 4 && displayMode > 0) {
      for (int i = 0; i < 4; i++) {
          const Vector& stress = materialPointers[i]->getStress();
          values(i) = stress(displayMode - 1);
      }
  }
  else {
      for (int i = 0; i < 4; i++)
          values(i) = 0.0;
  }

  // draw the polygon
  return theViewer.drawPolygon(coords, values, this->getTag());
}
