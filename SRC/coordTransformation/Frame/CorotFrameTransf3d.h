//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
// Description: This file contains the class definition for
// CorotFrameTransf3d.h. CorotFrameTransf3d provides the
// abstraction of a corotation transformation for a spatial frame element
//
// Written by : cmp, March 2024
//
// Adapted from work by: Remo Magalhaes de Souza
//
#ifndef CorotFrameTransf3d_h
#define CorotFrameTransf3d_h

#include <array>
#include "FrameTransform.h"
#include <Vector.h>
#include <Matrix.h>
#include <Matrix3D.h>

class Triad;

class CorotFrameTransf3d: public FrameTransform3d
{
public:
    CorotFrameTransf3d(int tag, const Vector &vecInLocXZPlane,
                       const Vector &rigJntOffsetI, 
                       const Vector &rigJntOffsetJ);

    CorotFrameTransf3d();
    ~CorotFrameTransf3d();

    const char *getClassType() const {return "CorotFrameTransf3d";};

    FrameTransform3d *getCopy();

    int initialize(Node *nodeIPointer, Node *nodeJPointer);
    int update();       // Set RI,RJ,Rbar, Ln, e and ul
    int commitState();
    int revertToLastCommit();        
    int revertToStart();

    double getInitialLength();
//  double getPresentLength();
//  double getCurrentLength();
    double getDeformedLength();
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

    // method used to rotate consistent mass matrix
    const Matrix &getGlobalMatrixFromLocal(const Matrix &local);

    // methods used in post-processing only
    const Vector &getPointGlobalCoordFromLocal(const Vector &localCoords);
    const Vector &getPointGlobalDisplFromBasic(double xi, const Vector &basicDisps);
    const Vector &getPointLocalDisplFromBasic(double xi, const Vector &basicDisps);    
    
    int  getLocalAxes(Vector &xAxis, Vector &yAxis, Vector &zAxis);

    // Movable Object
    int sendSelf(int cTag, Channel &theChannel);
    int recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

    // Tagged Object
    void Print(OPS_Stream &s, int flag = 0);

protected:
    int addTangent(MatrixND<12,12>& M, const VectorND<12>& pl);

    virtual const Layout& getNodeLayout() {
      static std::vector<int> l {
         1,
         1,
         0,
         0,
         1,
         1,
      };
      return l;
    }
    virtual const Layout& getForceLayout() {
      static std::vector<int> l {
        FrameTransform3d::N,  // 1
        FrameTransform3d::T,  // 1
        FrameTransform3d::My, // 0
        FrameTransform3d::Mz, // 0
        FrameTransform3d::My, // 1
        FrameTransform3d::Mz, // 1
      };
      return l;
    }

protected:

private:
    // compute the transformation matrix
    void compTransfMatrixBasicGlobal(const Triad&  r, const Triad&  E, const Triad&  rI, const Triad&  rJ);
    void compTransfMatrixLocalGlobal(Matrix &Tlg);


    //
    // Internal data
    //
    std::array<Node*, 2> nodes;                // pointers to the element two endnodes

    Vector xAxis;                              // local x axis
    Vector vAxis;                              // Vector that lies in local plane xz

    // Rigid joint offsets
    enum {
      end_i = 0<<1,
      end_j = 0<<2,
    };
    int joint_offsets;

    Vector nodeIOffset, 
           nodeJOffset;                        

    double L;                                  // initial element length
    double Ln;                                 // current element length (at trial state)
    
                                               // (the columns of which are the element local axes)
    OpenSees::VectorND<4> alphaIq;             // quaternion for node I
    OpenSees::VectorND<4> alphaJq;             // quaternion for node I 
    OpenSees::VectorND<4> alphaIqcommit;       // commited quaternion for node I
    OpenSees::VectorND<4> alphaJqcommit;       // commited quaternion for node J
    Vector alphaI;                             // last trial rotations end i
    Vector alphaJ;                             // last trial rotatations end j

    Vector ul;                                 // local displacements
    Vector ulcommit;                           // commited local displacements
    Vector ulpr;                               // previous local displacements

    OpenSees::MatrixND<12,12> ag;

    OpenSees::MatrixND<7,12> T;    // transformation matrix from local to global system
    OpenSees::Matrix3D A;
    OpenSees::Matrix3D R0;         // rotation matrix from local to global coordinates
    OpenSees::Matrix3D e, RI, RJ, Rbar;
    
    double *nodeIInitialDisp,
           *nodeJInitialDisp;
    bool initialDispChecked;

    // Static workspace variables
    static Matrix Tp;                 // transformation matrix to renumber dofs
    static MatrixND<12,3> Lr2, Lr3;   // auxiliary matrices
};

#endif
