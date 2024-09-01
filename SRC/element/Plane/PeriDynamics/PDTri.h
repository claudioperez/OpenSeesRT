/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **

** ****************************************************************** */
//
// Description: This file contains the class definition for PDTri.
//
// Written: Xuan Hu
//
// Adapted from: Roozbeh Geraili Mikola (roozbehg@berkeley.edu)
//
//
#ifndef PDTri_h
#define PDTri_h

#ifndef _bool_h
#  include <stdbool.h>
#endif

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

class Node;
class NDMaterial;
class Response;

class PDTri : public Element {
public:
  PDTri(int tag, int nd1, int nd2, int nd3, NDMaterial& m, const char* type, double t,
        double pressure = 0.0, double rho = 0.0, double b1 = 0.0, double b2 = 0.0);
  PDTri();
  ~PDTri();

  const char*
  getClassType() const
  {
    return "PDTri";
  };
  static constexpr const char* class_name = "PDTri";

  int getNumExternalNodes() const;
  const ID& getExternalNodes();
  Node** getNodePtrs();

  int getNumDOF();
  void setDomain(Domain* theDomain);

  // public methods to set the state of the element
  int commitState();
  int revertToLastCommit();
  int revertToStart();
  int update();

  // public methods to obtain stiffness, mass, damping and residual information
  const Matrix& getTangentStiff();
  const Matrix& getInitialStiff();
  const Matrix& getMass();

  void zeroLoad();
  int addLoad(ElementalLoad* theLoad, double loadFactor);
  int addInertiaLoadToUnbalance(const Vector& accel);

  const Vector& getResistingForce();
  const Vector& getResistingForceIncInertia();

  // public methods for element output
  int sendSelf(int commitTag, Channel& theChannel);
  int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);

  void Print(OPS_Stream& s, int flag = 0);

  Response* setResponse(const char** argv, int argc, OPS_Stream& s);

  int getResponse(int responseID, Information& eleInformation);

  int setParameter(const char** argv, int argc, Parameter& param);
  int updateParameter(int parameterID, Information& info);

  // RWB; PyLiq1 & TzLiq1 need to see the excess pore pressure and initial stresses.
  friend class PyLiq1;
  friend class TzLiq1;
  friend class QzLiq1; // Sumeet

protected:
private:
  static constexpr int numgp    = 1; // number of gauss points
  static constexpr int numnodes = 3; // number of nodes


  NDMaterial** theMaterial; // pointer to the ND material objects

  ID connectedExternalNodes; // Tags of PDTri nodes

  Node* theNodes[3];

  static double matrixData[36]; // array data for matrix
  static Matrix K;              // Element stiffness, damping, and mass Matrix
  static Vector P;              // Element resisting force vector
  Vector Q;                     // Applied nodal loads
  double b[2];                  // Body forces

  double appliedB[2]; // Body forces applied with load pattern
  int applyLoad;      // flag for body force in load

  Vector pressureLoad; // Pressure load at nodes

  double thickness; // Element thickness
  double pressure;  // Normal surface traction (pressure) over entire element
                    // Note: positive for outward normal
  double rho;
  static double shp[3][3]; // Stores shape functions and derivatives (overwritten)
  static double pts[1][2]; // Stores quadrature points
  static double wts[1];    // Stores quadrature weights

  // private member functions - only objects of this class can call these
  double shapeFunction(double xi, double eta);
  void setPressureLoadAtNodes();

  Matrix* Ki;
};

#endif