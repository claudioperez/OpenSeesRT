#pragma once
#include <material/section/SectionForceDeformation.h>

enum FrameStress : int {
  End   =     0,
  N     =     2, //
  Vy    =     3, // 0b00000010
  Vz    =     5, // 0b00000100
  T     =     6, // 0b00000000
  My    =     4, // 0b00000000
  Mz    =     1, // 0b00000000
  R     =     7, // 0b00000000
  Q     =     8, // 0b00000000
  B     =     9, // 0b00000000
  W     =    10, // 0b00000000
  Max   =    11,
};

typedef       int             FrameStressLayout[10];
static constexpr FrameStressLayout UnknownScheme {0};



class FrameSection : public SectionForceDeformation {

public:
  FrameSection(int tag, int clstag): SectionForceDeformation(tag, clstag) {};

  virtual FrameSection* getFrameCopy() =0;

  virtual SectionForceDeformation* getCopy() {return this->getFrameCopy();};

  template <int n, const FrameStressLayout& scheme>
  int setTrialState(OpenSees::VectorND<n, double> e) {
    double data[MaxResultants]{};

    int m = this->getOrder();
    Vector trial(data, m);

    const ID& layout = this->getType();

    for (int i=0; i<n; i++) {
      for (int j=0; j<m; j++)
        if (layout(j) == scheme[i])
          trial[j] = e[i];
    }

    return this->setTrialSectionDeformation(trial);
  }

  template <int n, const FrameStressLayout& scheme>
  OpenSees::VectorND<n, double> getDeformation() {

    OpenSees::VectorND<n,double> sout;

    const ID& layout = this->getType();

    int m = this->getOrder();

    const Vector& es = this->getSectionDeformation();
    for (int i=0; i<n; i++) {
      sout[i] = 0.0;
      for (int j=0; j<m; j++)
        if (layout(j) == scheme[i])
          sout[i] = es(j);
    }

    return sout;
  }

  template <int n, const FrameStressLayout& scheme>
  OpenSees::VectorND<n, double> getResultant() {

    OpenSees::VectorND<n,double> sout;

    const ID& layout = this->getType();

    int m = this->getOrder();

    const Vector& s = this->getStressResultant();
    for (int i=0; i<n; i++) {
      sout[i] = 0.0;
      for (int j=0; j<m; j++)
        if (layout(j) == scheme[i])
          sout[i] = s(j);
    }

    return sout;
  }

  template <int n, const FrameStressLayout& scheme>
  OpenSees::MatrixND<n,n, double> getTangent(State state) {

    OpenSees::MatrixND<n,n,double> kout;

    const ID& layout = this->getType();

    int m = this->getOrder();

//  const Matrix& ks = this->getSectionTangent();
    const Matrix& ks = (state == State::Init)
                      ? this->getInitialTangent()
                      : this->getSectionTangent();

    for (int i=0; i<n; i++) {
      for (int j=0; j<n; j++) {  
        kout(i,j) = 0.0;
        for (int k=0; k<m; k++) {
          if (layout(k) == scheme[i]) {
            for (int l=0; l<m; l++)
              if (layout(l) == scheme[j]) {
                kout(i,j) = ks(k,l);
              }
          }
        }
      }
    }

    return kout;
  }

  template <int n, const FrameStressLayout& scheme>
  OpenSees::MatrixND<n,n, double> getFlexibility(State state=State::Pres) {

    OpenSees::MatrixND<n,n,double> fout;

    const ID& layout = this->getType();

    int m = this->getOrder();

//  const Matrix& ks = this->getSectionFlexibility();
    const Matrix& ks = (state == State::Init)
                      ? this->getInitialFlexibility()
                      : this->getSectionFlexibility();

    for (int i=0; i<n; i++) {
      for (int j=0; j<n; j++) {  
        fout(i,j) = 0.0;
        for (int k=0; k<m; k++) {
          if (layout(k) == scheme[i]) {
            for (int l=0; l<m; l++)
              if (layout(l) == scheme[j]) {
                fout(i,j) = ks(k,l);
              }
          }
        }
      }
    }

    return fout;
  }

private:
  static constexpr int MaxResultants = 20;
};
