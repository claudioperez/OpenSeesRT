/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */
//
#include <ElasticLinearFrameSection3d.h>
#include <Matrix.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Parameter.h>
#include <classTags.h>

#include <stdlib.h>
#include <string.h>

using OpenSees::MatrixND;

constexpr int SEC_TAG_ElasticLinearFrame3d = 0;

static int layout_array[] = {
    FrameStress::N,
    FrameStress::Vy,
    FrameStress::Vz,
    FrameStress::T,
    FrameStress::My,
    FrameStress::Mz
};

ID ElasticLinearFrameSection3d::layout(layout_array, 6);


ElasticLinearFrameSection3d::ElasticLinearFrameSection3d()
: FrameSection(0, SEC_TAG_ElasticLinearFrame3d),
  E(0.0), A(0.0), Iz(0.0), Iy(0.0), G(0.0), J(0.0),// e(4)
  Ksen(nullptr)
{
}

ElasticLinearFrameSection3d::ElasticLinearFrameSection3d(int tag, 
    double E_in, 
    double A_in, 
    double Qz,
    double Qy,
    double Az,
    double Ay,
    double Iz_in,
    double Iy_in,
    double Iyz,
    double G_in, 
    double J_in,
    double yc,     // Centroid location
    double zc,     
    double ys,     // Shear center location
    double zs
    )
: FrameSection(tag, SEC_TAG_ElasticLinearFrame3d),
  E(E_in), A(A_in), 
  Iz(Iz_in), 
  Iy(Iy_in), 
  G(G_in),
  J(J_in), //e(4),
  I0(Iy+Iz),
  Ksen(nullptr),
  Fs(new MatrixND<6,6> {}),
  Ks(new MatrixND<6,6> {
         {{  E*A,     0.,    0.,    0.,   E*Qy,  -E*Qz}, //   0.,},
          {    0.,  G*Ay,    0., -G*Qy,     0.,     0.}, // G*Qy, },
          {    0.,    0.,  G*Az,  G*Qz,     0.,     0.},
          {    0., -G*Qy,  G*Qz,  G*I0,     0.,     0.},
          {  E*Qy,    0.,    0.,    0.,  E*Iy ,  E*Iyz},
          { -E*Qz,    0.,    0.,    0., -E*Iyz,  E*Iz }}}
  )
{
  //
  Ks->invert(*Fs);

  if (E <= 0.0)  {
    //opserr << "ElasticLinearFrameSection3d::ElasticLinearFrameSection3d -- Input E <= 0.0\n";
  }
  
  if (A <= 0.0)  {
    //opserr << "ElasticLinearFrameSection3d::ElasticLinearFrameSection3d -- Input A <= 0.0\n";
  }
  
  if (Iz <= 0.0)  {
    //opserr << "ElasticLinearFrameSection3d::ElasticLinearFrameSection3d -- Input Iz <= 0.0\n";
  }
  
  if (Iy <= 0.0)  {
    //opserr << "ElasticLinearFrameSection3d::ElasticLinearFrameSection3d -- Input Iy <= 0.0\n";
  }
  
  if (G <= 0.0)  {
    //opserr << "ElasticLinearFrameSection3d::ElasticLinearFrameSection3d -- Input G <= 0.0\n";
  }
  
  if (J <= 0.0)  {
    //opserr << "ElasticLinearFrameSection3d::ElasticLinearFrameSection3d -- Input J <= 0.0\n";
  }
  
}


ElasticLinearFrameSection3d::~ElasticLinearFrameSection3d()
{
  if (Ksen != nullptr)
    delete Ksen;

  return;
}


FrameSection*
ElasticLinearFrameSection3d::getFrameCopy()
{
  ElasticLinearFrameSection3d *theCopy = new ElasticLinearFrameSection3d();
  *theCopy = *this;

  return theCopy;
}


int 
ElasticLinearFrameSection3d::commitState()
{
  return 0;
}


int 
ElasticLinearFrameSection3d::revertToLastCommit()
{
  return 0;
}


int 
ElasticLinearFrameSection3d::revertToStart()
{
  return 0;
}


int
ElasticLinearFrameSection3d::setTrialSectionDeformation(const Vector &def)
{
  e = def;
    return 0;
}


const Vector &
ElasticLinearFrameSection3d::getSectionDeformation() // needed ?
{
  v.setData(e);
  return v;
}


const Vector &
ElasticLinearFrameSection3d::getStressResultant()
{

  s = (*Ks)*e;

  v.setData(s);
  
  return v;
}


const Matrix &
ElasticLinearFrameSection3d::getSectionTangent()
{
  M.setData(*Ks);
  return M;
}


const Matrix &
ElasticLinearFrameSection3d::getInitialTangent()
{
  M.setData(*Ks);
  return M;
}


const Matrix &
ElasticLinearFrameSection3d::getSectionFlexibility ()
{
  M.setData(*Fs);
  return M;
}


const Matrix &
ElasticLinearFrameSection3d::getInitialFlexibility ()
{
  M.setData(*Fs);
  return M;
}


const ID&
ElasticLinearFrameSection3d::getType()
{
  return layout;
}


int
ElasticLinearFrameSection3d::getOrder() const
{
  return 6;
}

int
ElasticLinearFrameSection3d::sendSelf(int commitTag, Channel &theChannel)
{
    int res = 0;

    static Vector data(7);

    int dataTag = this->getDbTag();
    
      data(0) = this->getTag();
    data(1) = E;
    data(2) = A;    
    data(3) = Iz;
    data(4) = Iy;
    data(5) = G;
    data(6) = J;
    
    res += theChannel.sendVector(dataTag, commitTag, data);
    if(res < 0) {
      opserr << "ElasticLinearFrameSection3d::sendSelf -- failed to send data\n";
      return res;
    }
    
    return res;
}


int
ElasticLinearFrameSection3d::recvSelf(int commitTag, Channel &theChannel,
                                          FEM_ObjectBroker &theBroker)
{
    int res = 0;
    
      static Vector data(7);

    int dataTag = this->getDbTag();

    res += theChannel.recvVector(dataTag, commitTag, data);
    if(res < 0) {
      opserr << "ElasticLinearFrameSection3d::recvSelf -- failed to receive data\n";
      return res;
    }

    this->setTag((int)data(0));
    E  = data(1);
    A  = data(2);    
    Iz = data(3);
    Iy = data(4);
    G  = data(5);
    J  = data(6);    

    return res;
}
 
void
ElasticLinearFrameSection3d::Print(OPS_Stream &s, int flag)
{
    if (flag == OPS_PRINT_PRINTMODEL_SECTION) {
        s << "ElasticLinearFrameSection3d, tag: " << this->getTag() << endln;
        s << "\t E: " << E << endln;
        s << "\t A: " << A << endln;
        s << "\tIz: " << Iz << endln;
        s << "\tIy: " << Iy << endln;
        s << "\t G: " << G << endln;
        s << "\t J: " << J << endln;
    }

    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": \"" << this->getTag() << "\", ";
        s << "\"type\": \"ElasticLinearFrameSection3d\", ";
        s << "\"E\": " << E << ", ";
        s << "\"G\": " << G << ", ";
        s << "\"A\": " << A << ", ";
        s << "\"Jx\": " << J << ", ";
        s << "\"Iy\": " << Iy << ", ";
        s << "\"Iz\": " << Iz << "}";
    }
}

int
ElasticLinearFrameSection3d::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;

  if (strcmp(argv[0],"E") == 0) {
    param.setValue(E);
    return param.addObject(1, this);
  }
  if (strcmp(argv[0],"A") == 0) {
    param.setValue(A);
    return param.addObject(2, this);
  }
  if (strcmp(argv[0],"Iz") == 0) {
    param.setValue(Iz);
    return param.addObject(3, this);
  }
  if (strcmp(argv[0],"Iy") == 0) {
    param.setValue(Iy);
    return param.addObject(4, this);
  }
  if (strcmp(argv[0],"G") == 0) {
    param.setValue(G);
    return param.addObject(5, this);
  }
  if (strcmp(argv[0],"J") == 0) {
    param.setValue(J);
    return param.addObject(6, this);
  }
  return -1;
}

int
ElasticLinearFrameSection3d::updateParameter(int paramID, Information &info)
{
  if (paramID == 1)
    E = info.theDouble;
  if (paramID == 2)
    A = info.theDouble;
  if (paramID == 3)
    Iz = info.theDouble;
  if (paramID == 4)
    Iy = info.theDouble;
  if (paramID == 5)
    G = info.theDouble;
  if (paramID == 6)
    J = info.theDouble;

  return 0;
}

int
ElasticLinearFrameSection3d::activateParameter(int paramID)
{
  parameterID = paramID;

  return 0;
}

const Vector&
ElasticLinearFrameSection3d::getStressResultantSensitivity(int gradIndex,
                                                           bool conditional)
{
  s.zero();

  if (parameterID == 1) { // E
    s[0] = A*e[0];
    s[1] = Iz*e[1];
    s[2] = Iy*e[2];
  }
  if (parameterID == 2) // A
    s[0] = E*e[0];
  if (parameterID == 3) // Iz
    s[1] = E*e[1];
  if (parameterID == 4) // Iy
    s[2] = E*e[2];
  if (parameterID == 5) // G
    s[3] = J*e[3];
  if (parameterID == 6) // J
    s[3] = G*e[3];

  v.setData(s);
  return v;
}

const Matrix&
ElasticLinearFrameSection3d::getInitialTangentSensitivity(int gradIndex)
{
  if (Ksen == nullptr)
    Ksen = new Matrix(nr,nr);

  return *Ksen;
}

