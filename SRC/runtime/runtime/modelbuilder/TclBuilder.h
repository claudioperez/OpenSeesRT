
// Description: This file contains the class definition for TclBuilder.
// A TclBuilder adds the commands to create the model for the standard
// models that can be generated using the elements released with the g3
// framework. currently these elements include:
//
// What: "@(#) TclBuilder.h, revA"

#ifndef TclBuilder_h
#define TclBuilder_h

#include <ModelBuilder.h>
#include <string>
#include <tcl.h>

class LoadPattern;
class SectionForceDeformation;
class SectionRepres;
class NDMaterial;
class CrdTrasnf;
class TaggedObjectStorage;

class TclBuilder : public ModelBuilder {
public:
  TclBuilder(Domain &theDomain, int ndm, int ndf);
  virtual ~TclBuilder();

  // eventually make private
  int currentSectionTag = -1;

  int buildFE_Model(void);

  int getNDM(void) const;
  int getNDF(void) const;

  LoadPattern *getCurrentLoadPattern(void);

private:
  int ndm; // space dimension of the mesh
  int ndf; // number of degrees of freedom per node

  LoadPattern* m_current_load_pattern = nullptr;

protected:
  Tcl_Interp *theInterp;
};

#endif

