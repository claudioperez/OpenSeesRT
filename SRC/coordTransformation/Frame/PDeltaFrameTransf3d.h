//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
// Description: This file contains the class definition for
// PDeltaFrameTransf3d.h. PDeltaFrameTransf3d provides the
// abstraction of a linear transformation for a spatial frame
// between the global and basic coordinate systems
//
// Written: Remo Magalhaes de Souza (rmsouza@ce.berkeley.edu)
// Created: 04/2000
//
#ifndef PDeltaFrameTransf3d_h
#define PDeltaFrameTransf3d_h

#include <FrameTransform.h>
#include <Vector.h>
#include <Matrix.h>

class PDeltaFrameTransf3d: public FrameTransform3d
{
public:
    PDeltaFrameTransf3d(int tag, const Vector &vecInLocXZPlane);
    PDeltaFrameTransf3d(int tag, const Vector &vecInLocXZPlane,
                        const Vector &rigJntOffsetI,
                        const Vector &rigJntOffsetJ);
    
    PDeltaFrameTransf3d();
    ~PDeltaFrameTransf3d();
    
    const char *getClassType() const {return "PDeltaFrameTransf3d";};
    
    double getInitialLength();
    double getDeformedLength();
    
    virtual int initialize(Node *node1Pointer, Node *node2Pointer);
    virtual int update();
    virtual int commitState();
    virtual int revertToLastCommit();        
    virtual int revertToStart();
    
    const Vector &getBasicTrialDisp();
    const Vector &getBasicIncrDisp();
    const Vector &getBasicIncrDeltaDisp();
    const Vector &getBasicTrialVel();
    const Vector &getBasicTrialAccel();
    
    
    virtual VectorND<12>    pushResponse(VectorND<12>&pl) override final;
    virtual VectorND<12>    pushConstant(const VectorND<12>&pl) const override final;
    virtual MatrixND<12,12> pushResponse(MatrixND<12,12>& kl, const VectorND<12>& pl) override final;
    virtual MatrixND<12,12> pushConstant(const MatrixND<12,12>& kl) override final;

    const Vector &getGlobalResistingForce(const Vector &basicForce, const Vector &p0);
    const Matrix &getGlobalStiffMatrix(const Matrix &basicStiff, const Vector &basicForce);
    const Matrix &getInitialGlobalStiffMatrix(const Matrix &basicStiff);
    
    FrameTransform3d *getCopy();
    
    int sendSelf(int cTag, Channel &theChannel);
    int recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    
    void Print(OPS_Stream &s, int flag = 0);
    
    // method used to rotate consistent mass matrix
    const Matrix &getGlobalMatrixFromLocal(const Matrix &local);

    // methods used in post-processing only
    const Vector &getPointGlobalCoordFromLocal(const Vector &localCoords);
    const Vector &getPointGlobalDisplFromBasic(double xi, const Vector &basicDisps);
    const Vector &getPointLocalDisplFromBasic(double xi, const Vector &basicDisps);    
    
    int getLocalAxes(Vector &xAxis, Vector &yAxis, Vector &zAxis);

private:
    int computeElemtLengthAndOrient();
    void compTransfMatrixLocalGlobal(Matrix &Tlg);
    
    // internal data
    Node *nodeIPtr, *nodeJPtr;  // pointers to the element two endnodes
    
    double *nodeIOffset, *nodeJOffset;    // rigid joint offsets
    
    double R[3][3];      // rotation matrix
    double L;            // undeformed element length
    double ul17;         // Transverse local displacement offsets of P-Delta
    double ul28;


    double *nodeIInitialDisp,
           *nodeJInitialDisp;
    bool initialDispChecked;
};

#endif
