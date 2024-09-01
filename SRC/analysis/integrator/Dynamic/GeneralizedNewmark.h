//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
// Description: This file contains the class definition for GeneralizedNewmark.
// GeneralizedNewmark is an algorithmic class for performing a transient analysis
// using the GeneralizedNewmark integration scheme.
//
// Written : cmp
// Created : 06/2024
// Adapted from Newmark.cpp
//
#ifndef GeneralizedNewmark_h
#define GeneralizedNewmark_h

#include <TransientIntegrator.h>
#include <Vector.h>

class DOF_Group;
class FE_Element;

class GeneralizedNewmark : public TransientIntegrator
{
public:
    // constructors
    GeneralizedNewmark(double gamma, double beta, 
            double alphaF, double alpahM,
            int uflag=1,                   // choose which "u"nknown is solved for: d, v or a
            int iflag=3,                   // choose how to "i"nitialize the unknown: Dd=0, Dv=0 or Da=0
            bool aflag=false);

    // destructor
    ~GeneralizedNewmark();
    
    // methods which define what the FE_Element and DOF_Groups add
    // to the system of equation object.
    virtual int formEleTangent(FE_Element *theEle)  final;
    virtual int formNodTangent(DOF_Group *theDof)   final;
    virtual int formEleResidual(FE_Element* theEle) final;
    virtual int formNodUnbalance(DOF_Group* theDof) final;
    
    int domainChanged();    
    int newStep(double deltaT);
    int revertToLastStep();
    virtual int update(const Vector &deltaU);

    double getCFactor();

    const Vector &getVel();
    
    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    
    void Print(OPS_Stream &s, int flag = 0);        
    
    // AddingSensitivity:BEGIN //////////////////////////////////
    int revertToStart();
    int formSensitivityRHS(int gradNum);
    int formIndependentSensitivityRHS();
    int saveSensitivity   (const Vector &v, int gradNum, int numGrads);
    int commitSensitivity (int gradNum, int numGrads);  
    int computeSensitivities( );
    // AddingSensitivity:END ////////////////////////////////////

protected:

private:
    enum Unknown {
      Displacement=1,
      Velocity=2,
      Acceleration=3
    };
    int unknown;                    // flag indicating whether displ(1), vel(2) or accel(3) increments
    int unknown_initialize = 1;     //

    double gamma;
    double beta;
    double alphaF;
    double alphaM;

    int step;                       // track step number to initialize accelerations
    double dt;                      // store time step to determine step number
    double c1, c2, c3;              // some constants we need to keep
    Vector *Ut, *Utdot, *Utdotdot;  // response quantities at time t
    Vector *Ua, *Uadot, *Uadotdot;  // response quantities at time t+alpha
    Vector *U, *Udot, *Udotdot;     // response quantities at time t+deltaT
    bool determiningMass;           // flag to check if just want the mass contribution

    // Adding Sensitivity
    int sensitivityFlag;
    int gradNumber;
    Vector *massMatrixMultiplicator;
    Vector *dampingMatrixMultiplicator;
    int assemblyFlag;
    Vector independentRHS;
    Vector dUn, dVn, dAn;
    //////////////////////
};

#endif
