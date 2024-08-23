/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
//
// Purpose: This file contains the implementation of class MP_Constraint.
//
// Written: fmk 
// Created: 11/96
// Revision: A
//
#include <MP_Constraint.h>
#include <string.h>
#include <Matrix.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Domain.h>
#include <assert.h>

static int numMPs = 0;
static int nextTag = 0;

 
// constructor for FEM_ObjectBroker			// Arash
MP_Constraint::MP_Constraint(int clasTag )		
:DomainComponent(nextTag++, clasTag),
 nodeRetained(0),nodeConstrained(0),constraint(nullptr),constrDOF(0),retainDOF(0),
 dbTag1(0), dbTag2(0)
{
  numMPs++;
}

// constructor for Subclass
MP_Constraint::MP_Constraint(int nodeRetain, int nodeConstr, 
			     ID &constrainedDOF, 
			     ID &retainedDOF, int clasTag)
:DomainComponent(nextTag++, clasTag),
 nodeRetained(nodeRetain), nodeConstrained(nodeConstr), 
 constraint(nullptr), constrDOF(0), retainDOF(0),  dbTag1(0), dbTag2(0)
{
  numMPs++;
  
  constrDOF = new ID(constrainedDOF);
  retainDOF = new ID(retainedDOF);
}


// general constructor for ModelBuilder
MP_Constraint::MP_Constraint(int nodeRetain, int nodeConstr, Matrix &constr,
			     ID &constrainedDOF, ID &retainedDOF)
:DomainComponent(nextTag++, CNSTRNT_TAG_MP_Constraint), 
 nodeRetained(nodeRetain), nodeConstrained(nodeConstr), 
 constraint(0), constrDOF(0), retainDOF(0), dbTag1(0), dbTag2(0)
{
  numMPs++;    
  constrDOF = new ID(constrainedDOF);
  retainDOF = new ID(retainedDOF);    
  
  constraint = new Matrix(constr);
}



MP_Constraint::~MP_Constraint()
{
    // invoke the destructor on the matrix and the two ID objects
    if (constraint != 0)
	delete constraint;
    if (constrDOF != 0)
	delete constrDOF;
    if (retainDOF != 0)
	delete retainDOF;    
    
    numMPs--;
    if (numMPs == 0)
      nextTag = 0;
}


int
MP_Constraint::getNodeRetained(void) const
{
    // return id of retained node
    return nodeRetained;
}

int
MP_Constraint::getNodeConstrained(void) const
{
    // return id of constrained node    
    return nodeConstrained;
}


const ID &
MP_Constraint::getConstrainedDOFs(void) const
{
  assert(constrDOF != nullptr);
  // return the ID corresponding to constrained DOF of Ccr
  return *constrDOF;    
}


// return the ID corresponding to retained DOF of Ccr
const ID &
MP_Constraint::getRetainedDOFs(void) const
{
    assert(retainDOF != nullptr);
    return *retainDOF;    
}

// does nothing MP_Constraint objects are time invariant
int 
MP_Constraint::applyConstraint(double timeStamp)
{
    return 0;
}

bool
MP_Constraint::isTimeVarying(void) const
{
    return false;
}


// return the constraint matrix Ccr
const Matrix &
MP_Constraint::getConstraint(void)
{
    assert(constraint != nullptr);
    return *constraint;    
}

int 
MP_Constraint::sendSelf(int cTag, Channel &theChannel)
{
    static ID data(10);
    int dataTag = this->getDbTag();

    data(0) = this->getTag(); 
    data(1) = nodeRetained;
    data(2) = nodeConstrained;
    if (constraint == 0) data(3) = 0; else data(3) = constraint->noRows();
    if (constraint == 0) data(4) = 0; else data(4) = constraint->noCols();    
    if (constrDOF == 0) data(5) = 0; else data(5) = constrDOF->Size();    
    if (retainDOF == 0) data(6) = 0; else data(6) = retainDOF->Size();        
    
    // need two database tags for ID objects
    if (constrDOF != 0 && dbTag1 == 0) 
      dbTag1 = theChannel.getDbTag();
    if (retainDOF != 0 && dbTag2 == 0) 
      dbTag2 = theChannel.getDbTag();

    data(7) = dbTag1;
    data(8) = dbTag2;
    data(9) = nextTag;

    int result = theChannel.sendID(dataTag, cTag, data);
    if (result < 0) {
	// opserr << "WARNING MP_Constraint::sendSelf - error sending ID data\n";
	return result;  
    }    
    
    if (constraint != 0 && constraint->noRows() != 0) {
	int result = theChannel.sendMatrix(dataTag, cTag, *constraint);
	if (result < 0) {
	    // opserr << "WARNING MP_Constraint::sendSelf ";
	    // opserr << "- error sending Matrix data\n"; 
	    return result;  
	}
    }

    if (constrDOF != 0 && constrDOF->Size() != 0) {
	int result = theChannel.sendID(dbTag1, cTag, *constrDOF);
	if (result < 0) {
	    // opserr << "WARNING MP_Constraint::sendSelf ";
	    // opserr << "- error sending constrained data\n"; 
	    return result;  
	}
    }

    if (retainDOF != 0 && retainDOF->Size() != 0) {
	int result = theChannel.sendID(dbTag2, cTag, *retainDOF);
	if (result < 0) {
	    // opserr << "WARNING MP_Constraint::sendSelf ";
	    // opserr << "- error sending retained data\n"; 
	    return result;  
	}
    }
    
    return 0;
}


int 
MP_Constraint::recvSelf(int cTag, Channel &theChannel, 
			FEM_ObjectBroker &theBroker)
{
    int dataTag = this->getDbTag();
    static ID data(10);
    int result = theChannel.recvID(dataTag, cTag, data);
    if (result < 0) {
	// opserr << "WARNING MP_Constraint::recvSelf - error receiving ID data\n";
	return result;  
    }    

    this->setTag(data(0));
    nodeRetained = data(1);
    nodeConstrained = data(2);
    int numRows = data(3); 
    int numCols = data(4);
    dbTag1 = data(7);
    dbTag2 = data(8);
    nextTag = data(9);

    if (numRows != 0 && numCols != 0) {
	constraint = new Matrix(numRows,numCols);
	
	int result = theChannel.recvMatrix(dataTag, cTag, *constraint);
	if (result < 0) {
	    // opserr << "WARNING MP_Constraint::recvSelf ";
	    // opserr << "- error receiving Matrix data\n"; 
	    return result;  
	}
    }    
    int size = data(5);
    if (size != 0) {
	constrDOF = new ID(size);
	int result = theChannel.recvID(dbTag1, cTag, *constrDOF);
	if (result < 0) {
	    // opserr << "WARNING MP_Constraint::recvSelf ";
	    // opserr << "- error receiving constrained data\n"; 
	    return result;  
	}	
    }
    
    size = data(6);
    if (size != 0) {
	retainDOF = new ID(size);
	int result = theChannel.recvID(dbTag2, cTag, *retainDOF);
	if (result < 0) {
	    // opserr << "WARNING MP_Retainaint::recvSelf ";
	    // opserr << "- error receiving retained data\n"; 
	    return result;  
	}	
    }    
    
    return 0;
}



void
MP_Constraint::Print(OPS_Stream &s, int flag)
{     

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
      // const char *indent = "            ";
      const char *indent = " ";
      const char *newln = " ";
      s << "            {";
      s << indent << "\"name\": \"" << this->getTag() << "\"," << newln;
      s << indent << "\"node_constrained\": " << nodeConstrained << "," << newln;
      s << indent << "\"node_retained\": " << nodeRetained << "," << newln;
      if (constrDOF != 0 && retainDOF != 0) {
        s << indent << "\"constrained_dof\": [";
        const int nc = (*constrDOF).Size();
        for (int i=0; i < nc; i++)
          s << (*constrDOF)(i)+1 << (i < nc-1? ", " : ""); 
        s << "]," << newln;

        const int nr = (*retainDOF).Size();
        s << indent << "\"retained_dof\": [";
        for (int i=0; i < nr; i++)
          s << (*retainDOF)(i)+1 << (i < nr-1? ", " : "");
        s << "]," << newln;

        if (constraint != 0) {
          s << indent << "\"constraint_matrix\": [";
          // TODO print constraint matrix
          // : s << *constraint ;
          s << "]";
        }
      }
      s << "}";
      return;

    } else {
      s << "MP_Constraint: " << this->getTag() << "\n";
      s << "\tNode Constrained: " << nodeConstrained;
      s << " node Retained: " << nodeRetained << "\n";
      if (constrDOF != 0 && retainDOF != 0) {
        s << " constrained dof: ";
        for (int i=0; i<(*constrDOF).Size(); i++)
          s << (*constrDOF)(i)+1 << " ";
        s << "\n";
          s << " retained dof: ";        
        for (int i=0; i<(*retainDOF).Size(); i++)
          s << (*retainDOF)(i)+1 << " ";
        s << "\n";
        if (constraint != 0)
          s << " constraint matrix: " << *constraint << "\n";
      }
    }
}


