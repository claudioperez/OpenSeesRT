//
// Description: This file contains the class definition for TclBuilder.
//
// written: cmp
//
#include <stdlib.h>

#include <ID.h>

#include <Domain.h>

#include <TclBuilder.h>
#include <MultiSupportPattern.h>

extern MultiSupportPattern *theTclMultiSupportPattern;

TclBuilder::TclBuilder(Domain &theDomain, int NDM, int NDF)
    : ModelBuilder(theDomain), ndm(NDM), ndf(NDF)
{

}

TclBuilder::~TclBuilder()
{

}

int
TclBuilder::buildFE_Model(void) {return 0;}

int
TclBuilder::getNDM(void) const {return ndm;}

int
TclBuilder::getNDF(void) const {return ndf;}

LoadPattern*
TclBuilder::getCurrentLoadPattern(void) 
{
  return m_current_load_pattern;
}

