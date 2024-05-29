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
// Written: cmp 2024
//
// Purpose: This file contains the class definition for PrismFrame3d.
// PrismFrame3d is a plane frame member.
//
#ifndef PrismFrame3d_h
#define PrismFrame3d_h

#include <Element.h>
#include <Node.h>
#include <Matrix.h>
#include <Vector.h>
#include <MatrixND.h>
#include <VectorND.h>

class Channel;
class Information;
class CrdTransf;
class Response;
class SectionForceDeformation;

class PrismFrame3d : public Element
{
  public:
    PrismFrame3d();        
    PrismFrame3d(int tag, double A, double E, double G, 
		  double Jx, double Iy, double Iz,
          int Nd1, int Nd2, CrdTransf &theTransf,
          double rho = 0.0, int cMass = 0,
		  int releasez = 0, int releasey = 0);

    PrismFrame3d(int tag, int Nd1, int Nd2, SectionForceDeformation *section, 
		  CrdTransf &theTransf, double rho = 0.0, int cMass = 0,
		  int releasez = 0, int releasey = 0);

    ~PrismFrame3d();

    const char *getClassType() const {return "PrismFrame3d";};

    int getNumExternalNodes() const;
    const ID &getExternalNodes();
    Node **getNodePtrs();

    int getNumDOF();
    void setDomain(Domain *theDomain);
    
    int commitState();
    int revertToLastCommit();        
    int revertToStart();
    
    int update();
    const Matrix &getTangentStiff();
    const Matrix &getInitialStiff();
    const Matrix &getMass();    

    void zeroLoad();	
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);

    const Vector &getResistingForce();
    const Vector &getResistingForceIncInertia();            
    
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    
    void Print(OPS_Stream &s, int flag =0);    

    Response *setResponse (const char **argv, int argc, OPS_Stream &s);
    int getResponse (int responseID, Information &info);
 
    int setParameter (const char **argv, int argc, Parameter &param);
    int updateParameter (int parameterID, Information &info);

  private:
    double A,E,G,Jx,Iy,Iz;

    double rho;
    int cMass;

    int releasez; // moment release for bending about z-axis 0=none, 1=I, 2=J, 3=I,J
    int releasey; // same for y-axis
    
    Vector Q;
    
    OpenSees::MatrixND<6,6> kb;
    OpenSees::VectorND<6>   q;
    OpenSees::VectorND<6>   q0;  // Fixed end forces in basic system (no torsion)
    OpenSees::VectorND<6>   p0;  // Reactions in basic system (no torsion)

    double wx;
    double wy;
    double wz;
    
    Node *theNodes[2];

    ID  connectedExternalNodes;    

    CrdTransf *theCoordTransf;

    static Matrix K;
    static Vector P;
};

#endif
