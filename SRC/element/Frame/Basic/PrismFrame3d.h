//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
// Purpose: This file contains the class definition for PrismFrame3d.
// PrismFrame3d is a plane frame member.
//
// RELEASES
// [ ] 
//
// Written: cmp 2024
//
#ifndef PrismFrame3d_h
#define PrismFrame3d_h

#include <array>
#include <element/Frame/BasicFrame3d.h>
#include <Node.h>
#include <Matrix.h>
#include <Vector.h>
#include <MatrixND.h>
#include <VectorND.h>

class Channel;
class Information;
class FrameTransform3d;
class Response;
class FrameSection;

class PrismFrame3d : public BasicFrame3d
{
  public:
    PrismFrame3d();
    PrismFrame3d(int tag, std::array<int, 2>& nodes,
                 double A, double E, double G, 
		             double Jx, double Iy, double Iz,
                 FrameTransform3d &theTransf,
                 double density, int mass_flag,
		             int releasez, int releasey);

    PrismFrame3d(int tag, std::array<int,2>& nodes,
      FrameSection &section, 
		  FrameTransform3d &theTransf,
      double density, int mass_flag, bool use_mass,
		  int releasez, int releasey);

//  ~PrismFrame3d();

    const char *getClassType() const {return "PrismFrame3d";};
/*
//  void zeroLoad();	
//  int addLoad(ElementalLoad *theLoad, double loadFactor);
//  int addInertiaLoadToUnbalance(const Vector &accel);
*/
    
    int update();
    virtual const Matrix &getMass() final;
    int commitState();
    int revertToLastCommit();        
    int revertToStart();

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    
    void Print(OPS_Stream &s, int flag =0);

    Response *setResponse(const char **argv, int argc, OPS_Stream &s) final;
    int getResponse(int responseID, Information &info) final;
 
    virtual int setParameter(const char **argv, int argc, Parameter &param) final;
    virtual int updateParameter(int parameterID, Information &info) final;

  protected:
    // For BasicFrame3d
    virtual  OpenSees::VectorND<6>&   getBasicForce() final;
    virtual  OpenSees::MatrixND<6,6>& getBasicTangent(State state, int rate) final;
    // For FiniteElement
    virtual int setNodes() final;

  private:
    void formBasicStiffness(OpenSees::MatrixND<6,6>& kb) const;

    double A,E,G,Jx,Iy,Iz;
    double L;

    int mass_flag;
    double total_mass,
           twist_mass,
           density;

    int section_tag;
    int geom_flag = 0; 
    int releasez; // moment release for bending about z-axis 0=none, 1=I, 2=J, 3=I,J
    int releasey; // same for y-axis

    OpenSees::MatrixND<6,6> ke;
    OpenSees::MatrixND<6,6> km;
    OpenSees::MatrixND<6,6> kg;
    OpenSees::VectorND<6>   q;


};

#endif
