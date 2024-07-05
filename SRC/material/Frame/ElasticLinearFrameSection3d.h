/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */
//
//
#ifndef ElasticLinearFrameSection3d_h
#define ElasticLinearFrameSection3d_h

#include <memory>
#include <FrameSection.h>
#include <MatrixND.h>

class Matrix;
class Vector;
class Channel;
class FEM_ObjectBroker;
class Information;
class Parameter;

class ElasticLinearFrameSection3d : public FrameSection
{

 public:
  ElasticLinearFrameSection3d(int tag, 
      double E, 
      double A, 
      double Az,
      double Ay,
      double Qz,
      double Qy,
      double Iz,
      double Iy,
      double Iyz,
      double G, 
      double J,
      double yc, // Centroid location
      double zc, 
      double ys, // Shear center location
      double zs
  );

  ElasticLinearFrameSection3d();
  ~ElasticLinearFrameSection3d();

  const char *getClassType() const {return "ElasticLinearFrameSection3d";};
  
  int commitState();
  int revertToLastCommit();
  int revertToStart();
  
  int setTrialSectionDeformation(const Vector&);
  const Vector &getSectionDeformation();
  
  const Vector &getStressResultant();
  const Matrix &getSectionTangent();
  const Matrix &getInitialTangent();
  const Matrix &getSectionFlexibility();
  const Matrix &getInitialFlexibility();
  
  FrameSection *getFrameCopy();
  const ID &getType();
  int getOrder() const;
  
  int sendSelf(int commitTag, Channel &theChannel);
  int recvSelf(int commitTag, Channel &theChannel,
               FEM_ObjectBroker &theBroker);
  
  void Print(OPS_Stream &s, int flag = 0);

  int setParameter(const char **argv, int argc, Parameter &param);
  int updateParameter(int parameterID, Information &info);
  int activateParameter(int parameterID);
  const Vector& getStressResultantSensitivity(int gradIndex,
                                              bool conditional);
  const Matrix& getInitialTangentSensitivity(int gradIndex);

 protected:
  
 private:

  double E, A, Iz, Iy, I0, G, J;
  
  constexpr static int nr = 6;
  Vector  v;
  Matrix  M,         // Generic matrix for returning Ks or Fs (nr x nr)
         *Ksen;      // Tangent sensitivity (nr x nr)

  OpenSees::VectorND<nr> s;
  OpenSees::VectorND<nr> e;                      // section trial deformations
  std::shared_ptr<OpenSees::MatrixND<nr,nr>> Ks;
  std::shared_ptr<OpenSees::MatrixND<nr,nr>> Fs;

  int parameterID;

//static int laydat[nr];
  static ID layout;

};

#endif
