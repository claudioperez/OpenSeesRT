#include <vector>
#include <Element.h>

class Domain;

template<int ndim>
class PeriElement : public Element {
  PeriElement(int tag, PeriDomain<ndim>& domain);


  virtual void setDomain(Domain*);

  virtual int getNumExternalNodes() const;
  virtual const ID &getExternalNodes()  =0;    
  virtual Node **getNodePtrs()  =0;    
  virtual int    getNumDOF()    =0;

  virtual int  commitState() final;
  virtual int  revertToLastCommit() final;
  virtual int  revertToStart() final;
  virtual int  update() final;


  virtual const Vector &getResistingForce() =0;
  virtual const Matrix &getTangentStiff();
  virtual const Matrix &getInitialStiff();
  virtual const Matrix &getDamp();
  virtual const Matrix &getMass();

private:
  Domain*           edomain;
  PeriDomain<ndim>& pdomain;

  Vector* pforce;

};

#include "PeriElement.tpp"

