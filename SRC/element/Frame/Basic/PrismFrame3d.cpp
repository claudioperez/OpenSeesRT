//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
// 3D prismatic frame element. 
//
// Written: cmp 2024
//
#include <Frame/BasicFrame3d.h>
#include <PrismFrame3d.h>
#include <Domain.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <Information.h>
#include <Parameter.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>
#include <FrameSection.h>
#include <FrameTransform.h>
#include <ID.h>
#include <math.h>
#include <stdlib.h>
#include <Exponential.h>

PrismFrame3d::PrismFrame3d()
  :BasicFrame3d(0,ELE_TAG_ElasticBeam3d),
   A(0.0), E(0.0), G(0.0), Jx(0.0), Iy(0.0), Iz(0.0),
   total_mass(0.0), twist_mass(0.0)
{
  q.zero();
}

PrismFrame3d::PrismFrame3d(int tag, std::array<int, 2>& nodes,
                           double  a, double  e, double  g, 
                           double jx, double iy, double iz,
                           FrameTransform3d &coordTransf, 
                           double r, int cm, int rz, int ry)
  :BasicFrame3d(tag,ELE_TAG_ElasticBeam3d, nodes, coordTransf),
   A(a), E(e), G(g), Jx(jx), Iy(iy), Iz(iz), 
   mass_flag(cm), density(r),
   releasez(rz), releasey(ry),
   section_tag(-1)
{
  q.zero();
}

PrismFrame3d::PrismFrame3d(int tag, 
                           std::array<int,2>& nodes,
                           FrameSection &section,  
                           FrameTransform3d &coordTransf, 
                           double rho_, int cMass, bool use_mass, 
                           int rz, int ry
                           )
  :BasicFrame3d(tag,ELE_TAG_ElasticBeam3d, nodes, coordTransf),
   mass_flag(cMass), density(rho_),
   releasez(rz), releasey(ry)
{
  q.zero();

  Jx = 1.0;

  section.getIntegral(Field::Unit,   State::Init, A);
  section.getIntegral(Field::UnitZZ, State::Init, Iy);
  section.getIntegral(Field::UnitYY, State::Init, Iz);

  section_tag = section.getTag();

  const Matrix &sectTangent = section.getInitialTangent();
  const ID &sectCode = section.getType();
  for (int i=0; i<sectCode.Size(); i++) {
    int code = sectCode(i);
    switch(code) {
    case SECTION_RESPONSE_P:
      E = sectTangent(i,i)/A;
      break;
    case SECTION_RESPONSE_T:
      G  = sectTangent(i,i)/Jx;
      break;
    default:
      break;
    }
  }

  // TODO
  if (!use_mass) {
    if (section.getIntegral(Field::Density, State::Init, density) == 0) {
      ;
    }
  }

}

int
PrismFrame3d::setNodes()
{

  int status = this->BasicFrame3d::setNodes();

  if (status != 0)
    return status;

  L = this->getLength(State::Init);

  if (L == 0.0) {
    opserr << "PrismFrame3d::setDomain  tag: " << this->getTag() << " -- Element has zero length\n";
    return -1;
  }

  //
  formBasicStiffness(km);

  total_mass = density*L; 
  twist_mass = (density/A)*Jx*L;

  return 0;
}

inline void
PrismFrame3d::formBasicStiffness(OpenSees::MatrixND<6,6>& kb) const
{
//  const double L = this->getLength(State::Init);
    const double oneOverL = 1.0/L;
    const double EoverL   = E*oneOverL;
    const double EAoverL  = A*EoverL;              // EA/L
    const double GJoverL  = G*Jx*oneOverL;         // GJ/L
    
    kb.zero();
    kb(0,0) = EAoverL;
    kb(5,5) = GJoverL;
    if (releasez == 0) {
      double EIzoverL2 = 2.0*Iz*EoverL;            // 2EIz/L
      double EIzoverL4 = 2.0*EIzoverL2;            // 4EIz/L
      kb(1,1) = kb(2,2) = EIzoverL4;
      kb(2,1) = kb(1,2) = EIzoverL2;
    }
    if (releasez == 1) { // release I
      kb(2,2) = 3.0*Iz*EoverL;
    }
    if (releasez == 2) { // release J
      kb(1,1) = 3.0*Iz*EoverL;
    }

    if (releasey == 0) {
      double EIyoverL2 = 2.0*Iy*EoverL;            // 2EIy/L
      double EIyoverL4 = 2.0*EIyoverL2;            // 4EIy/L
      kb(3,3) = kb(4,4) = EIyoverL4;
      kb(4,3) = kb(3,4) = EIyoverL2;
    }
    if (releasey == 1) { // release I
      double EIoverL3 = 3.0*Iy*EoverL;
      kb(4,4) = EIoverL3;
    }
    if (releasey == 2) { // release J
      double EIoverL3 = 3.0*Iy*EoverL;
      kb(3,3) = EIoverL3;
    }
}

int
PrismFrame3d::commitState()
{
  int retVal = 0;
  if ((retVal = this->Element::commitState()) != 0) {
    opserr << "PrismFrame3d::commitState () - failed in base class";
  }    
  retVal += theCoordTransf->commitState();
  return retVal;
}

int
PrismFrame3d::revertToLastCommit()
{
  return theCoordTransf->revertToLastCommit();
}

int
PrismFrame3d::revertToStart()
{
  return theCoordTransf->revertToStart();
}

int
PrismFrame3d::update()
{
  int ok = theCoordTransf->update();

  const Vector &v = theCoordTransf->getBasicTrialDisp();

  double N = E*A/L*v(0);

  switch (geom_flag) {
    case 0:
      ke = km ;

      break;

    case 1:
      kg.zero();
      kg(1,1) = kg(2,2) =  4.0*N/L;
      kg(1,2) = kg(2,1) = -1.0*N/L;
      ke = km + kg;
      break;

    case 2:
      kg.zero();
      { // zz
        double psi = L*std::sqrt(std::fabs(N)/(E*Iz));

        if (psi > 1e-12) {
          double cs = std::cos(psi);
          double sn = std::sin(psi);
          double C  = E*Iz/(L*(2.0 - 2.0*cs - psi*sn));
          ke(1,1) = ke(2,2) =  C*psi*(sn - psi*cs);
          ke(1,2) = ke(2,1) =  C*psi*(psi - sn);
        } else {
          ke = km;
        }
      }
      { // yy
        double psi = L*std::sqrt(std::fabs(N)/(E*Iy));

        if (psi > 1e-12) {
          double cs = std::cos(psi);
          double sn = std::sin(psi);
          double C  = E*Iy/(L*(2.0 - 2.0*cs - psi*sn));
          ke(3,3) = ke(4,4) =  C*psi*(sn - psi*cs);
          ke(4,3) = ke(3,4) =  C*psi*(psi - sn);
        } else {
          ke = km;
        }
      }
      break;
  }

  q = ke*v;
  q += q0;

  return ok;
}


OpenSees::VectorND<6>&
PrismFrame3d::getBasicForce()
{
  return q;
}

const Vector &
PrismFrame3d::getResistingForce()
{
  double q0 = q[0];
  double q1 = q[1];
  double q2 = q[2];
  double q3 = q[3];
  double q4 = q[4];
  double q5 = q[5];

  double oneOverL = 1.0 / theCoordTransf->getInitialLength();

  thread_local VectorND<12> pl;
  pl[0]  = -q0;                    // Ni
  pl[1]  =  oneOverL * (q1 + q2);  // Viy
  pl[2]  = -oneOverL * (q3 + q4);  // Viz
  pl[3]  = -q5;                    // Ti
  pl[4]  =  q3;
  pl[5]  =  q1;
  pl[6]  =  q0;                    // Nj
  pl[7]  = -pl[1];                 // Vjy
  pl[8]  = -pl[2];                 // Vjz
  pl[9]  = q5;                     // Tj
  pl[10] = q4;
  pl[11] = q2;

  thread_local VectorND<12> pf{0.0};
  pf[0] = p0[0];
  pf[1] = p0[1];
  pf[7] = p0[2];
  pf[2] = p0[3];
  pf[8] = p0[4];


  thread_local VectorND<12> pg;
  thread_local Vector wrapper(pg);
//    const Vector p0Vec(p0);
//    P = theCoordTransf->getGlobalResistingForce(q, p0Vec);
  pg  = theCoordTransf->pushResponse(pl);
  pg += theCoordTransf->pushConstant(pf);

  // Subtract other external nodal loads ... P_res = P_int - P_ext
  if (total_mass != 0.0)
    wrapper.addVector(1.0, p_iner, -1.0);

  return wrapper;
}

OpenSees::MatrixND<6,6>&
PrismFrame3d::getBasicTangent(State flag, int rate)
{
  return ke;
}

const Matrix &
PrismFrame3d::getMass()
{
    if (total_mass == 0.0) {

        thread_local MatrixND<12,12> M{0.0};
        thread_local Matrix Wrapper{M};
        return Wrapper;

    } else if (mass_flag == 0)  {

        thread_local MatrixND<12,12> M{0.0};
        thread_local Matrix Wrapper{M};
        // lumped mass matrix
        double m = 0.5*total_mass;
        M(0,0) = m;
        M(1,1) = m;
        M(2,2) = m;
        M(6,6) = m;
        M(7,7) = m;
        M(8,8) = m;
        return Wrapper;

    } else {
        // consistent (cubic, prismatic) mass matrix

        double L  = this->getLength(State::Init);
        double m  = total_mass/420.0;
        double mx = twist_mass;
        thread_local MatrixND<12,12> ml{0};

        ml(0,0) = ml(6,6) = m*140.0;
        ml(0,6) = ml(6,0) = m*70.0;

        ml(3,3) = ml(9,9) = mx/3.0; // Twisting
        ml(3,9) = ml(9,3) = mx/6.0;

        ml( 2, 2) = ml( 8, 8) =  m*156.0;
        ml( 2, 8) = ml( 8, 2) =  m*54.0;
        ml( 4, 4) = ml(10,10) =  m*4.0*L*L;
        ml( 4,10) = ml(10, 4) = -m*3.0*L*L;
        ml( 2, 4) = ml( 4, 2) = -m*22.0*L;
        ml( 8,10) = ml(10, 8) = -ml(2,4);
        ml( 2,10) = ml(10, 2) =  m*13.0*L;
        ml( 4, 8) = ml( 8, 4) = -ml(2,10);

        ml( 1, 1) = ml( 7, 7) =  m*156.0;
        ml( 1, 7) = ml( 7, 1) =  m*54.0;
        ml( 5, 5) = ml(11,11) =  m*4.0*L*L;
        ml( 5,11) = ml(11, 5) = -m*3.0*L*L;
        ml( 1, 5) = ml( 5, 1) =  m*22.0*L;
        ml( 7,11) = ml(11, 7) = -ml(1,5);
        ml( 1,11) = ml(11, 1) = -m*13.0*L;
        ml( 5, 7) = ml( 7, 5) = -ml(1,11);

        // Transform local mass matrix to global system
        return theCoordTransf->getGlobalMatrixFromLocal(ml);
    }
}


int
PrismFrame3d::sendSelf(int cTag, Channel &theChannel)
{
    int res = 0;

    static Vector data(19);
    
    data(0) = A;
    data(1) = E;
    data(2) = G;
    data(3) = Jx;
    data(4) = Iy;
    data(5) = Iz;

// 
    data( 6) = total_mass; // TODO
    data( 7) = mass_flag;
    data( 8) = this->getTag();
    data( 9) = connectedExternalNodes(0);
    data(10) = connectedExternalNodes(1);
    data(11) = theCoordTransf->getClassTag();            

    int dbTag = theCoordTransf->getDbTag();
    
    if (dbTag == 0) {
      dbTag = theChannel.getDbTag();
      if (dbTag != 0)
        theCoordTransf->setDbTag(dbTag);
    }

    data(12) = dbTag;
    
    data(13) = alphaM;
    data(14) = betaK;
    data(15) = betaK0;
    data(16) = betaKc;
    data(17) = releasez;
    data(18) = releasey;    
    
    // Send the data vector
    res += theChannel.sendVector(this->getDbTag(), cTag, data);
    if (res < 0) {
      opserr << "PrismFrame3d::sendSelf -- could not send data Vector\n";
      return res;
    }

    // Ask the CoordTransf to send itself
    res += theCoordTransf->sendSelf(cTag, theChannel);
    if (res < 0) {
      opserr << "PrismFrame3d::sendSelf -- could not send CoordTransf\n";
      return res;
    }

    return res;
}

int
PrismFrame3d::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = 0;
  static Vector data(19);

  res += theChannel.recvVector(this->getDbTag(), cTag, data);
  if (res < 0) {
    opserr << "PrismFrame3d::recvSelf -- could not receive data Vector\n";
    return res;
  }
  
  A  = data(0);
  E  = data(1); 
  G  = data(2); 
  Jx = data(3); 
  Iy = data(4); 
  Iz = data(5);     

  total_mass = data(6);
  mass_flag  = (int)data(7);
  this->setTag((int)data(8));
  connectedExternalNodes(0) = (int)data(9);
  connectedExternalNodes(1) = (int)data(10);
  
  alphaM   = data(13);
  betaK    = data(14);
  betaK0   = data(15);
  betaKc   = data(16);
  releasez = (int)data(17);
  releasey = (int)data(18);
  
  // Check if the CoordTransf is null; if so, get a new one
  int crdTag = (int)data(11);
  if (theCoordTransf == nullptr) {
    // TODO(cmp)
    theCoordTransf = nullptr; // theBroker.getNewFrameTransform3d(crdTag);
    if (theCoordTransf == 0) {
      opserr << "PrismFrame3d::recvSelf -- could not get a FrameTransform3d\n";
      return -1;
    }
  }

  // Check that the CoordTransf is of the right type; if not, delete
  // the current one and get a new one of the right type
  if (theCoordTransf->getClassTag() != crdTag) {
    delete theCoordTransf;
    // TODO(cmp)
    theCoordTransf = nullptr; // theBroker.getNewFrameTransform3d(crdTag);
    if (theCoordTransf == 0) {
      opserr << "PrismFrame3d::recvSelf -- could not get a FrameTransform3d\n";
      return -1;
    }
  }

  // Now, receive the CoordTransf
  theCoordTransf->setDbTag((int)data(12));
  res += theCoordTransf->recvSelf(cTag, theChannel, theBroker);
  if (res < 0) {
    opserr << "PrismFrame3d::recvSelf -- could not receive CoordTransf\n";
    return res;
  }

  return res;
}

void
PrismFrame3d::Print(OPS_Stream &s, int flag)
{

    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << OPS_PRINT_JSON_ELEM_INDENT << "{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"PrismFrame3d\", ";
        s << "\"nodes\": [" << connectedExternalNodes(0) << ", " 
                            << connectedExternalNodes(1) << "], ";
        s << "\"massperlength\": " << total_mass/L << ", ";
        s << "\"releasez\": "<< releasez << ", ";
        s << "\"releasey\": "<< releasey << ", ";                
        s << "\"crdTransformation\": " << theCoordTransf->getTag()  << ", ";
        // 
        if (section_tag > 0) {
          s << "\"section\": " << section_tag ; // << ", ";
        } else {
          s << "\"E\": "  << E  << ", ";
          s << "\"G\": "  << G  << ", ";
          s << "\"A\": "  << A  << ", ";
          s << "\"Jx\": " << Jx << ", ";
          s << "\"Iy\": " << Iy << ", ";
          s << "\"Iz\": " << Iz;
        }
        // 
        s << "}";
    }
    
    this->getResistingForce(); 

    if (flag == -1) {
        int eleTag = this->getTag();
        s << "EL_BEAM\t" << eleTag << "\t";
        s << "\t" << connectedExternalNodes(0) << "\t" << connectedExternalNodes(1);
        s << "\t0\t0.0000000\n";
    }

    else if (flag < -1) {
      int counter = (flag + 1) * -1;
      int eleTag = this->getTag();
      const Vector &force = this->getResistingForce();

      double P, MZ1, MZ2, VY, MY1, MY2, VZ, T;
      double L = theCoordTransf->getInitialLength();
      double oneOverL = 1.0 / L;

      P   = q[0];
      MZ1 = q[1];
      MZ2 = q[2];
      VY  = (MZ1 + MZ2)*oneOverL;
      T   = q[5];
      MY1 = q[3];
      MY2 = q[4];
      VZ  = (MY1 + MY2)*oneOverL;

      s << "FORCE\t" << eleTag << "\t" << counter << "\t0";
      s << "\t" << -P + p0[0] << "\t" << VY + p0[1] << "\t" << -VZ + p0[3] << "\n";
      s << "FORCE\t" << eleTag << "\t" << counter << "\t1";
      s << "\t" << P << ' ' << -VY + p0[2] << ' ' << VZ + p0[4] << "\n";
      s << "MOMENT\t" << eleTag << "\t" << counter << "\t0";
      s << "\t" << -T << "\t" << MY1 << "\t" << MZ1 << "\n";
      s << "MOMENT\t" << eleTag << "\t" << counter << "\t1";
      s << "\t" << T << ' ' << MY2 << ' ' << MZ2 << "\n";
    }

    else if (flag == 2) {
        this->getResistingForce(); // in case linear algo

        static Vector xAxis(3);
        static Vector yAxis(3);
        static Vector zAxis(3);

        theCoordTransf->getLocalAxes(xAxis, yAxis, zAxis);

        s << "#ElasticBeamColumn3D\n";
        s << "#LocalAxis " << xAxis(0) << " " << xAxis(1) << " " << xAxis(2);
        s << " " << yAxis(0) << " " << yAxis(1) << " " << yAxis(2) << " ";
        s << zAxis(0) << " " << zAxis(1) << " " << zAxis(2) << "\n";

        const Vector &xi = theNodes[0]->getCrds();
        const Vector &node2Crd = theNodes[1]->getCrds();
        const Vector &node1Disp = theNodes[0]->getDisp();
        const Vector &node2Disp = theNodes[1]->getDisp();

        s << "#NODE " << xi(0) << " " << xi(1) << " " << xi(2)
                << " " << node1Disp(0) << " " << node1Disp(1) << " " << node1Disp(2)
                << " " << node1Disp(3) << " " << node1Disp(4) << " " << node1Disp(5) << "\n";

        s << "#NODE " << node2Crd(0) << " " << node2Crd(1) << " " << node2Crd(2)
                << " " << node2Disp(0) << " " << node2Disp(1) << " " << node2Disp(2)
                << " " << node2Disp(3) << " " << node2Disp(4) << " " << node2Disp(5) << "\n";

        double N, Mz1, Mz2, Vy, My1, My2, Vz, T;
        double L = theCoordTransf->getInitialLength();
        double oneOverL = 1.0 / L;

        N = q[0];
        Mz1 = q[1];
        Mz2 = q[2];
        Vy = (Mz1 + Mz2)*oneOverL;
        My1 = q[3];
        My2 = q[4];
        Vz = -(My1 + My2)*oneOverL;
        T = q[5];

        s << "#END_FORCES " << -N + p0[0] << ' ' << Vy + p0[1] << ' ' << Vz + p0[3] << ' '
                            << -T << ' ' 
                            << My1 << ' ' 
                            << Mz1 << "\n";
        s << "#END_FORCES " << N << ' ' 
                            << -Vy + p0[2] << ' ' 
                            << -Vz + p0[4] << ' '
                            << T << ' ' 
                            << My2 << ' ' 
                            << Mz2 << "\n";
    }
    
    if (flag == OPS_PRINT_CURRENTSTATE) {
        this->getResistingForce(); // in case linear algo

        s << "\n  PrismFrame3d: " << this->getTag() << "\n";
        s << "\tConnected Nodes: " << connectedExternalNodes;
        s << "\tCoordTransf: " << theCoordTransf->getTag() << "\n";
        s << "\tmass density:  " << total_mass/L << ", mass_type: " << mass_flag << "\n";
        s << "\trelease about z:  " << releasez << "\n";
        s << "\trelease about y:  " << releasey << "\n";                
        double N, Mz1, Mz2, Vy, My1, My2, Vz, T;
        double L = theCoordTransf->getInitialLength();
        double oneOverL = 1.0 / L;

        N = q[0];
        Mz1 = q[1];
        Mz2 = q[2];
        Vy = (Mz1 + Mz2)*oneOverL;
        My1 = q[3];
        My2 = q[4];
        Vz = -(My1 + My2)*oneOverL;
        T = q[5];

        s << "\tEnd 1 Forces (P Mz Vy My Vz T): "
                << -N + p0[0] << ' ' << Mz1 << ' ' << Vy + p0[1] << ' ' << My1 << ' ' << Vz + p0[3] << ' ' << -T << "\n";
        s << "\tEnd 2 Forces (P Mz Vy My Vz T): "
                << N << ' ' << Mz2 << ' ' << -Vy + p0[2] << ' ' << My2 << ' ' << -Vz + p0[4] << ' ' << T << "\n";
    }

}


Response*
PrismFrame3d::setResponse(const char **argv, int argc, OPS_Stream &output)
{

  Response *theResponse = nullptr;

  output.tag("ElementOutput");
  output.attr("eleType","PrismFrame3d");
  output.attr("eleTag",this->getTag());
  output.attr("node1",connectedExternalNodes[0]);
  output.attr("node2",connectedExternalNodes[1]);
  
  // global forces
  if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0 ||
      strcmp(argv[0],"globalForce") == 0 || strcmp(argv[0],"globalForces") == 0) {


    output.tag("ResponseType","Px_1");
    output.tag("ResponseType","Py_1");
    output.tag("ResponseType","Pz_1");
    output.tag("ResponseType","Mx_1");
    output.tag("ResponseType","My_1");
    output.tag("ResponseType","Mz_1");
    output.tag("ResponseType","Px_2");
    output.tag("ResponseType","Py_2");
    output.tag("ResponseType","Pz_2");
    output.tag("ResponseType","Mx_2");
    output.tag("ResponseType","My_2");
    output.tag("ResponseType","Mz_2");

    theResponse =  new ElementResponse(this, 2, Vector(12));

  // local forces
  } else if (strcmp(argv[0],"localForce") == 0 || 
             strcmp(argv[0],"localForces") == 0) {

    output.tag("ResponseType","N_1");
    output.tag("ResponseType","Vy_1");
    output.tag("ResponseType","Vz_1");
    output.tag("ResponseType","T_1");
    output.tag("ResponseType","My_1");
    output.tag("ResponseType","Mz_1");
    output.tag("ResponseType","N_2");
    output.tag("ResponseType","Vy_2");
    output.tag("ResponseType","Vz_2");
    output.tag("ResponseType","T_2");
    output.tag("ResponseType","My_2");
    output.tag("ResponseType","Mz_2");

    theResponse =  new ElementResponse(this, 3, Vector(12));

  // basic forces
  } else if (strcmp(argv[0],"basicForce") == 0 || 
             strcmp(argv[0],"basicForces") == 0) {

    output.tag("ResponseType","N");
    output.tag("ResponseType","Mz_1");
    output.tag("ResponseType","Mz_2");
    output.tag("ResponseType","My_1");
    output.tag("ResponseType","My_2");
    output.tag("ResponseType","T");
    
    theResponse = new ElementResponse(this, 4, Vector(6));

  }  else if (strcmp(argv[0],"deformations") == 0 || 
              strcmp(argv[0],"basicDeformations") == 0) {
    
    output.tag("ResponseType","eps");
    output.tag("ResponseType","theta11");
    output.tag("ResponseType","theta12");
    output.tag("ResponseType","theta21");
    output.tag("ResponseType","theta22");
    output.tag("ResponseType","phi");
    theResponse = new ElementResponse(this, 5, Vector(6));
  }

  else if (strcmp(argv[0],"sectionX") == 0) {
    if (argc > 2) {
      float xL = atof(argv[1]);
      if (xL < 0.0)
        xL = 0.0;
      if (xL > 1.0)
        xL = 1.0;
      if (strcmp(argv[2],"forces") == 0) {
        theResponse = new ElementResponse(this,6,Vector(6));
        Information &info = theResponse->getInformation();
        info.theDouble = xL;
      }
    }   
  }

  if (theResponse == nullptr)
    theResponse = theCoordTransf->setResponse(argv, argc, output);

  output.endTag(); // ElementOutput

  return theResponse;
}

int
PrismFrame3d::getResponse(int responseID, Information &info)
{
  double L = theCoordTransf->getInitialLength();
  double oneOverL = 1.0/L;
  static Vector Res(12);
  Res = this->getResistingForce();
  opserr << Res << "\n";
  static Vector s(6);
  
  switch (responseID) {
  case 1: // stiffness
    return info.setMatrix(this->getTangentStiff());
    
  case 2: // global forces
    return info.setVector(Res);
    
  case 3: { // local forces
    double V, Mi, Mj, T;
    // Axial
    double N = q[0];
    P(6) =  N;
    P(0) = -N + p0[0];
    
    // Torsion
    T = q[5];
    P(9) =  T;
    P(3) = -T;
    
    // Moments about z and shear along y
    Mi    = q[1];
    Mj    = q[2];
    P(5)  = Mi;
    P(11) = Mj;
    V     = (Mi + Mj)*oneOverL;
    P(1)  =  V+p0[1];
    P(7)  = -V+p0[2];
    
    // Moments about y and shear along z
    Mi    = q[3];
    Mj    = q[4];
    P(4)  = Mi;
    P(10) = Mj;
    V     = (Mi + Mj)*oneOverL;
    P(2)  = -V + p0[3];
    P(8)  =  V + p0[4];

    return info.setVector(P);
  } 
  case 4: // basic forces
    return info.setVector(q);

  case 5:
    return info.setVector(theCoordTransf->getBasicTrialDisp());

  case 6: {
    double xL = info.theDouble;
    double x = xL*L;
    
    s(0) = q[0] + wx*(L-x);
    s(1) = q[1]*(xL-1.0) + q[2]*xL + 0.5*wy*x*(x-L);
    s(2) = (q[1] + q[2])/L + wy*(x-0.5*L);
    s(3) = q[3]*(xL-1.0) + q[4]*xL - 0.5*wz*x*(x-L);
    s(4) = (q[3] + q[4])/L - wz*(x-0.5*L);
    s(5) = q[5];

    return info.setVector(s);
  }
  default:
    break;
  }
  return -1;
}


int
PrismFrame3d::setParameter(const char **argv, int argc, Parameter &param)
{

  int status = this->BasicFrame3d::setParameter(argv, argc, param);
  if (status != -1)
    return status;

  if (argc < 1)
    return -1;

  if (strcmp(argv[0],"E") == 0) {
    param.setValue(E);
    return param.addObject(11, this);
  }

  if (strcmp(argv[0],"A") == 0) {
    param.setValue(A);
    return param.addObject(12, this);
  }

  if (strcmp(argv[0],"Iz") == 0) {
    param.setValue(Iz);
    return param.addObject(13, this);
  }

  if (strcmp(argv[0],"Iy") == 0) {
    param.setValue(Iy);
    return param.addObject(14, this);
  }

  if (strcmp(argv[0],"G") == 0) {
    param.setValue(G);
    return param.addObject(15, this);
  }

  if (strcmp(argv[0],"J") == 0) {
    param.setValue(Jx);
    return param.addObject(16, this);
  }

  return -1;
}

int
PrismFrame3d::updateParameter(int parameterID, Information &info)
{
    int status = this->BasicFrame3d::updateParameter(parameterID, info);
    if (status != -1)
      return status;

    switch (parameterID) {
      case -1:
        return -1;
      case 11:
        E = info.theDouble;
        return 0;
      case 12:
        A = info.theDouble;
        return 0;
      case 13:
        Iz = info.theDouble;
        return 0;
      case 14:
        Iy = info.theDouble;
        return 0;
      case 15:
        G = info.theDouble;
        return 0;
      case 16:
        Jx = info.theDouble;
        return 0;

      default:
        return -1;

    }

    // Update the element state

    formBasicStiffness(km);
}

