/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */
//
// Purpose: This file contains the class definition for PrismFrame3d.
// PrismFrame3d is a plane frame member.
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
                 double rho = 0.0, int cMass = 0,
		             int releasez = 0, int releasey = 0);

    PrismFrame3d(int tag, std::array<int,2>&, FrameSection *section, 
		  FrameTransform3d &theTransf, double rho = 0.0, int cMass = 0,
		  int releasez = 0, int releasey = 0);

//  ~PrismFrame3d();

    const char *getClassType() const {return "PrismFrame3d";};
/*
//  void zeroLoad();	
//  int addLoad(ElementalLoad *theLoad, double loadFactor);
//  int addInertiaLoadToUnbalance(const Vector &accel);
*/
    
    int update();
    int commitState();
    int revertToLastCommit();        
    int revertToStart();

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    
    void Print(OPS_Stream &s, int flag =0);    

    Response *setResponse (const char **argv, int argc, OPS_Stream &s);
    int getResponse (int responseID, Information &info);
 
    virtual int setParameter (const char **argv, int argc, Parameter &param) final;
    virtual int updateParameter (int parameterID, Information &info) final;

  protected:
    // For BasicFrame3d
    virtual  OpenSees::VectorND<6>&   getBasicForce() final;
    virtual  OpenSees::MatrixND<6,6>& getBasicTangent(State state, int rate) final;
    // For FiniteElement
    virtual int setNodes() final;

  private:
    void formBasicStiffness(OpenSees::MatrixND<6,6>& kb) const;

    double A,E,G,Jx,Iy,Iz,rho;
    double L;

    int geom_flag = 0; 

    OpenSees::MatrixND<6,6> ke;
    OpenSees::MatrixND<6,6> km;
    OpenSees::MatrixND<6,6> kg;
    OpenSees::VectorND<6>   q;

};

#endif
