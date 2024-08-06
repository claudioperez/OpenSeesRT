#include <Flag.h>
#include <Element.h>
#include <Domain.h>
#include <Node.h>
#include <Lagrange1D.cpp>
#include <SectionForceDeformation.h>
#include <quadrature/GaussLegendre1D.hpp>
#include <VectorND.h>

#include <BeamIntegration.h>
#include <Rotations.hpp>
#if 0

class ID;
class Domain;
using OpenSees::VectorND;

template<int nen, int nip>
class SimoFrame: 
  public Element 
{
  SimoFrame(int tag, int nodes[nen], 
            SectionForceDeformation *section[nip])
    : Element(tag, 0) {

    for (int i=0; i<nip; i++)
      sections[i] = section[i]->getCopy();

    for (int i=0; i<nen; i++)
      conn[i] = nodes[i];

//  beamIntegr->getSectionLocations(numSections, L, xi);
//  beamIntegr->getSectionWeights(numSections, L, wt);

  }

  virtual int update() override final {
    constexpr int ndf = 6;

    // Gauss loop
    for (int i=0; i<nip; i++) {

      Vector3D dx     {{0.0}};
      Vector3D theta  {{0.0}};
      Vector3D dtheta {{0.0}};

#if 0
      for (int j=0; j < nen; j++) {
        dx     += shape[i][1][j]*nodes[i]->getIncrDeltaDispl();
        theta  += shape[i][0][j]*nodes[i]->getIncrDeltaDispl();
        dtheta += shape[i][1][j]*nodes[i]->getIncrDeltaDispl();
      }

      auto gib = GibSO3(theta.norm());

      MatrixND<3,3> dR = ExpSO3(theta, gib);
      rotation[j] = dR*rotation[j];

   // curvature[j]->incr(rotations[j]->
#endif
    }
    return OpenSees::Flag::Success;
  }

  virtual int commitState() override final {
    return OpenSees::Flag::Success;
  }

  virtual int revertToLastCommit() override final {
    return OpenSees::Flag::Success;
  }

  virtual void setDomain(Domain* theDomain) override final {
    double x[nen];
    if (theDomain == nullptr) {
      for (int i=0; i < nen; i++)
        nodes[i] = nullptr;
      return;
    }

    for (int i=0; i < nen; i++) {
      nodes[i] = theDomain->getNode(conn[i]);
      if (!nodes[i] || nodes[i]->getNumberDOF() != 6) {
        opserr << "Bad Node\n";
        return;
      }
      lagrange<nen>(gauss.pts[i], x, shape[i]);
    }

    return this->Element::setDomain(theDomain);
  }


  private:
    int                      conn[nen];
    SectionForceDeformation *sections[nip];
    MatrixND<3,3>            rotation[nip];

    Node                    *nodes[nen];
    GaussLegendre<1, nip>    gauss;
    double                   shape[nip][2][nen];
    double wt[nip];
    double xi[nip];
};


typedef SimoFrame<2, 1> A;
#endif
