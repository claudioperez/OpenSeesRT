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
                                                                        
// $Revision: 1.13 $
// $Date: 2008-08-19 22:50:49 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/node/Node.h,v $
// Written: fmk 
// Created: 11/96
// Revision: A
//
// Purpose: This file contains the class interface for Node.
// A Node provides the abstraction for a node in the FEM.
// Nodes have original position, trial displacement, velocity and 
// acceleration and committed displacement, velocity and acceleration 
// (the last committed() trial quantities).
//
// What: "@(#) Node.h, revA"
//
#ifndef Node_h
#define Node_h

#define VIRTUAL

//#include <DomainComponent.h>
#include <TaggedObject.h>
#include <MovableObject.h>

// TODO: Remove include of NodeData
#include "NodeData.h"

class Element;
class Vector;
class Matrix;
class Channel;

class DOF_Group;
class NodalThermalAction; //L.Jiang [ SIF ]

class Node :
#if 0
  public DomainComponent
#else
  public TaggedObject,
  public MovableObject
#endif
{
  public:
    // constructors
    Node(int classTag);
    Node(int tag, int classTag);
    Node(int tag, int ndof, double Crd1, Vector *displayLoc = 0);
    Node(int tag, int ndof, double Crd1, double Crd2, Vector *displayLoc = 0);
    Node(int tag, int ndof, double Crd1, double Crd2, double Crd3, Vector *displayLoc = 0);
    Node(const Node &theCopy, bool copyMass = true);
    
    // destructor
    VIRTUAL ~Node();

    // public methods dealing with the DOF at the node
    VIRTUAL int  getNumberDOF(void) const;
    VIRTUAL void setDOF_GroupPtr(DOF_Group *theDOF_Grp);
    VIRTUAL DOF_Group *getDOF_GroupPtr(void);

    // public methods for obtaining the nodal coordinates
    VIRTUAL const Vector &getCrds(void) const;
    VIRTUAL void  setCrds(double Crd1);
    VIRTUAL void  setCrds(double Crd1, double Crd2);
    VIRTUAL void  setCrds(double Crd1, double Crd2, double Crd3);
    VIRTUAL void  setCrds(const Vector &);

    //
    // State
    //
    // public methods dealing with the committed state of the node
    VIRTUAL int commitState();
    VIRTUAL int revertToLastCommit();
    VIRTUAL int revertToStart();

    //
    // Response
    //
    // public methods for obtaining committed and trial 
    // response quantities of the node
    VIRTUAL const Vector &getDisp(void);
    VIRTUAL const Vector &getVel(void);
    VIRTUAL const Vector &getAccel(void);
    VIRTUAL const Vector &getIncrDisp(void);
    VIRTUAL const Vector &getIncrDeltaDisp(void);
    VIRTUAL const Vector &getTrialDisp(void);
    VIRTUAL const Vector &getTrialVel(void);
    VIRTUAL const Vector &getTrialAccel(void);

    // public methods for updating the trial response quantities
    VIRTUAL int setTrialDisp(double value, int dof);
    VIRTUAL int setTrialDisp(const Vector &);
    VIRTUAL int setTrialVel(const Vector &);
    VIRTUAL int setTrialAccel(const Vector &);
    VIRTUAL int incrTrialDisp(const Vector &);
    VIRTUAL int incrTrialVel(const Vector &);
    VIRTUAL int incrTrialAccel(const Vector &);

    // Thermodynamics
    virtual NodalThermalAction* getNodalThermalActionPtr(void);
    virtual void setNodalThermalActionPtr(NodalThermalAction* theAction);

    // Dynamics
    VIRTUAL const Matrix &getMass(void);
    VIRTUAL int setMass(const Matrix &theMass);
    VIRTUAL int setNumColR(int numCol);
    VIRTUAL int setR(int row, int col, double Value);
    VIRTUAL const Vector &getRV(const Vector &V);
    VIRTUAL int setRayleighDampingFactor(double alphaM);
    VIRTUAL const Matrix &getDamp(void);

    // Eigen vectors
    VIRTUAL int setNumEigenvectors(int numVectorsToStore);
    VIRTUAL int setEigenvector(int mode, const Vector &eigenVector);
    VIRTUAL const Matrix &getEigenvectors(void);

    // Misc. Responses
    VIRTUAL const Vector &getReaction();
    VIRTUAL int   addReactionForce(const Vector &, double factor);
    VIRTUAL int   resetReactionForce(int flag);

    VIRTUAL const Vector *getResponse(NodeData);
    int fillResponse(NodeData responseType, Vector& result, int offset=0);
    
    //
    // Load information
    //
    VIRTUAL void zeroUnbalancedLoad(void);
    VIRTUAL int addUnbalancedLoad(const Vector &load, double fact = 1.0);
    VIRTUAL int addInertiaLoadToUnbalance(const Vector &accel, double fact = 1.0);
    VIRTUAL const Vector &getUnbalancedLoad(void);
    VIRTUAL const Vector &getUnbalancedLoadIncInertia(void);

    
    //
    // Parallel
    //
    VIRTUAL int sendSelf(int commitTag, Channel &theChannel);
    VIRTUAL int recvSelf(int commitTag, Channel &theChannel, 
                         FEM_ObjectBroker &theBroker);

    //
    // Misc
    //
    VIRTUAL void Print(OPS_Stream &s, int flag = 0);

    // AddingSensitivity:BEGIN /////////////////////////////////////////
    int addInertiaLoadSensitivityToUnbalance(const Vector &accel, 
                                             double fact = 1.0, 
                                             bool tag=false);
    Matrix getMassSensitivity(void);
    VIRTUAL const Matrix &getDampSensitivity(void);
    int    getCrdsSensitivity(void);
    int    saveDispSensitivity(const Vector &v, int gradNum, int numGrads);
    int    saveVelSensitivity(const Vector &vdot, int gradNum, int numGrads);
    int    saveAccelSensitivity(const Vector &vdot, int gradNum, int numGrads);
    double getDispSensitivity(int dof, int gradNum);
    double getVelSensitivity(int dof, int gradNum);
    double getAccSensitivity(int dof, int gradNum);
    int    setParameter(const char **argv, int argc, Parameter &param);
    int    updateParameter(int parameterID, Information &info);
    int    activateParameter(int parameterID);
    // AddingSensitivity:END ///////////////////////////////////////////


    //
    // Display
    //
    VIRTUAL int getDisplayCrds(Vector &results, double fact, int displayMode=0);
    VIRTUAL int getDisplayRots(Vector& results, double fact, int displayMode=0);
    VIRTUAL int setDisplayCrds(const Vector &theCrds);


    Domain *getDomain(void) {return theDomain;};
    void setDomain(Domain *model) {theDomain = model;};

  protected:

  private:
#if 1
    Domain* theDomain;
#endif

    // Private global state
    int setGlobalMatrices();
    static Matrix **theMatrices;
    static int numMatrices;
    static Matrix **theVectors;
    static int numVectors;
    int index;


    // priavte methods used to create the Vector objects 
    // for the committed and trial response quantaties.
    int createDisp(void);
    int createVel(void);
    int createAccel(void);

    // private data associated with each node object
    int numberDOF;                    // number of dof at Node
    DOF_Group *theDOF_GroupPtr;       // pointer to associated DOF_Group
    Vector *Crd;                      // original nodal coords
    Vector *commitDisp, *commitVel, *commitAccel;  // committed quantities
    Vector *trialDisp, *trialVel, *trialAccel;     // trial quantities
    Vector *unbalLoad;                             // unbalanced load
    Vector *incrDisp;
    Vector *incrDeltaDisp;
    NodalThermalAction *theNodalThermalActionPtr; //Added by Liming Jiang for pointer to nodalThermalAction, [SIF]
    
    double *disp, *vel, *accel; // double arrays holding the displ, 
                                // vel and accel values

    Matrix *R;                          // nodal participation matrix
    Matrix *mass;                       // pointer to mass matrix
    Vector *unbalLoadWithInertia;
    double alphaM;                    // rayleigh damping factor 
    Matrix *theEigenvectors;


    int dbTag1, dbTag2, dbTag3, dbTag4; // needed for database

    // AddingSensitivity:BEGIN /////////////////////////////////////////
    Matrix *dispSensitivity;
    Matrix *velSensitivity;
    Matrix *accSensitivity;
    int parameterID;
    // AddingSensitivity:END ///////////////////////////////////////////


    Vector *reaction;

    Vector *displayLocation;

};

#endif

