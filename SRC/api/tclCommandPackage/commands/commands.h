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

// $Revision: 1.30 $
// $Date: 2010-09-13 21:33:06 $
// $Source: /usr/local/cvs/OpenSees/SRC/tcl/commands.h,v $

// Written: fmk
// Created: 04/98
// Revision: A
//
// Description: This file contains the functions that will be called by
// the interpreter when the appropriate command name is specified,
// see tkAppInit.C for command names.
//
// What: "@(#) commands.C, revA"

#include <OPS_Globals.h>

// by SAJalali
int OPS_recorderValue(ClientData clientData, G3_Runtime *rt, int argc,
                      TCL_Char **argv);

int OpenSeesAppInit(G3_Runtime *rt);

int OPS_SetObjCmd(ClientData clientData, G3_Runtime *rt, int argc,
                  Tcl_Obj *const *argv);

int OPS_SourceCmd(ClientData clientData, G3_Runtime *rt, int argc,
                  Tcl_Obj *const *argv);

int getNDM(ClientData clientData, G3_Runtime *rt, int argc,
           TCL_Char **argv);

int getNDF(ClientData clientData, G3_Runtime *rt, int argc,
           TCL_Char **argv);

int wipeModel(ClientData clientData, G3_Runtime *rt, int argc,
              TCL_Char **argv);

int wipeAnalysis(ClientData clientData, G3_Runtime *rt, int argc,
                 TCL_Char **argv);

int resetModel(ClientData clientData, G3_Runtime *rt, int argc,
               TCL_Char **argv);

int initializeAnalysis(ClientData clientData, G3_Runtime *rt, int argc,
                       TCL_Char **argv);

int setLoadConst(ClientData clientData, G3_Runtime *rt, int argc,
                 TCL_Char **argv);

int setCreep(ClientData clientData, G3_Runtime *rt, int argc,
             TCL_Char **argv);

int setTime(ClientData clientData, G3_Runtime *rt, int argc,
            TCL_Char **argv);

int getTime(ClientData clientData, G3_Runtime *rt, int argc,
            TCL_Char **argv);

int getLoadFactor(ClientData clientData, G3_Runtime *rt, int argc,
                  TCL_Char **argv);

int buildModel(ClientData clientData, G3_Runtime *rt, int argc,
               TCL_Char **argv);

int analyzeModel(ClientData clientData, G3_Runtime *rt, int argc,
                 TCL_Char **argv);

int printModel(ClientData clientData, G3_Runtime *rt, int argc,
               TCL_Char **argv);
int specifyAnalysis(ClientData clientData, G3_Runtime *rt, int argc,
                    TCL_Char **argv);
int specifySOE(ClientData clientData, G3_Runtime *rt, int argc,
               TCL_Char **argv);

int specifyNumberer(ClientData clientData, G3_Runtime *rt, int argc,
                    TCL_Char **argv);
int specifyConstraintHandler(ClientData clientData, G3_Runtime *rt,
                             int argc, TCL_Char **argv);
int specifyAlgorithm(ClientData clientData, G3_Runtime *rt, int argc,
                     TCL_Char **argv);

int specifyCTest(ClientData clientData, G3_Runtime *rt, int argc,
                 TCL_Char **argv);
int getCTestNorms(ClientData clientData, G3_Runtime *rt, int argc,
                  TCL_Char **argv);
int getCTestIter(ClientData clientData, G3_Runtime *rt, int argc,
                 TCL_Char **argv);

int specifyIntegrator(ClientData clientData, G3_Runtime *rt, int argc,
                      TCL_Char **argv);
int addRecorder(ClientData clientData, G3_Runtime *rt, int argc,
                TCL_Char **argv);
int addAlgoRecorder(ClientData clientData, G3_Runtime *rt, int argc,
                    TCL_Char **argv);

int addDatabase(ClientData clientData, G3_Runtime *rt, int argc,
                TCL_Char **argv);

int playbackRecorders(ClientData clientData, G3_Runtime *rt, int argc,
                      TCL_Char **argv);

int playbackAlgorithmRecorders(ClientData clientData, G3_Runtime *rt,
                               int argc, TCL_Char **argv);

int groundExcitation(ClientData clientData, G3_Runtime *rt, int argc,
                     TCL_Char **argv);

int eigenAnalysis(ClientData clientData, G3_Runtime *rt, int argc,
                  TCL_Char **argv);

int modalProperties(ClientData clientData, G3_Runtime *rt, int argc,
                    TCL_Char **argv);

int responseSpectrum(ClientData clientData, G3_Runtime *rt, int argc,
                     TCL_Char **argv);

int videoPlayer(ClientData clientData, G3_Runtime *rt, int argc,
                TCL_Char **argv);

int removeObject(ClientData clientData, G3_Runtime *rt, int argc,
                 TCL_Char **argv);

int eleForce(ClientData clientData, G3_Runtime *rt, int argc,
             TCL_Char **argv);

int localForce(ClientData clientData, G3_Runtime *rt, int argc,
               TCL_Char **argv);

int eleDynamicalForce(ClientData clientData, G3_Runtime *rt, int argc,
                      TCL_Char **argv);

int eleResponse(ClientData clientData, G3_Runtime *rt, int argc,
                TCL_Char **argv);

int findID(ClientData clientData, G3_Runtime *rt, int argc,
           TCL_Char **argv);

int nodeDisp(ClientData clientData, G3_Runtime *rt, int argc,
             TCL_Char **argv);

int nodeReaction(ClientData clientData, G3_Runtime *rt, int argc,
                 TCL_Char **argv);

int nodeUnbalance(ClientData clientData, G3_Runtime *rt, int argc,
                  TCL_Char **argv);

int nodeEigenvector(ClientData clientData, G3_Runtime *rt, int argc,
                    TCL_Char **argv);

int nodeCoord(ClientData clientData, G3_Runtime *rt, int argc,
              TCL_Char **argv);

int setNodeCoord(ClientData clientData, G3_Runtime *rt, int argc,
                 TCL_Char **argv);

int updateElementDomain(ClientData clientData, G3_Runtime *rt, int argc,
                        TCL_Char **argv);

int eleType(ClientData clientData, G3_Runtime *rt, int argc,
            TCL_Char **argv);

int eleNodes(ClientData clientData, G3_Runtime *rt, int argc,
             TCL_Char **argv);

int nodeBounds(ClientData clientData, G3_Runtime *rt, int argc,
               TCL_Char **argv);

int nodeVel(ClientData clientData, G3_Runtime *rt, int argc,
            TCL_Char **argv);

int setNodeVel(ClientData clientData, G3_Runtime *rt, int argc,
               TCL_Char **argv);

int setNodeDisp(ClientData clientData, G3_Runtime *rt, int argc,
                TCL_Char **argv);

int setNodeAccel(ClientData clientData, G3_Runtime *rt, int argc,
                 TCL_Char **argv);

int nodeAccel(ClientData clientData, G3_Runtime *rt, int argc,
              TCL_Char **argv);

int nodeResponse(ClientData clientData, G3_Runtime *rt, int argc,
                 TCL_Char **argv);

int calculateNodalReactions(ClientData clientData, G3_Runtime *rt, int argc,
                            TCL_Char **argv);

int getNodeTags(ClientData clientData, G3_Runtime *rt, int argc,
                TCL_Char **argv);

int getEleTags(ClientData clientData, G3_Runtime *rt, int argc,
               TCL_Char **argv);

int fixedNodes(ClientData clientData, G3_Runtime *rt, int argc,
               TCL_Char **argv);

int fixedDOFs(ClientData clientData, G3_Runtime *rt, int argc,
              TCL_Char **argv);

int constrainedNodes(ClientData clientData, G3_Runtime *rt, int argc,
                     TCL_Char **argv);

int constrainedDOFs(ClientData clientData, G3_Runtime *rt, int argc,
                    TCL_Char **argv);

int retainedNodes(ClientData clientData, G3_Runtime *rt, int argc,
                  TCL_Char **argv);

int retainedDOFs(ClientData clientData, G3_Runtime *rt, int argc,
                 TCL_Char **argv);

int nodeDOFs(ClientData clientData, G3_Runtime *rt, int argc,
             TCL_Char **argv);

int nodeMass(ClientData clientData, G3_Runtime *rt, int argc,
             TCL_Char **argv);

int nodePressure(ClientData clientData, G3_Runtime *rt, int argc,
                 TCL_Char **argv);

int getParamTags(ClientData clientData, G3_Runtime *rt, int argc,
                 TCL_Char **argv);

int getParamValue(ClientData clientData, G3_Runtime *rt, int argc,
                  TCL_Char **argv);

int sdfResponse(ClientData clientData, G3_Runtime *rt, int argc,
                TCL_Char **argv);

// AddingSensitivity:BEGIN /////////////////////////////////////////////////

int computeGradients(ClientData clientData, G3_Runtime *rt, int argc,
                     TCL_Char **argv);

int sensNodeDisp(ClientData clientData, G3_Runtime *rt, int argc,
                 TCL_Char **argv);

int sensLambda(ClientData clientData, G3_Runtime *rt, int argc,
               TCL_Char **argv); // Abbas

int sensNodeVel(ClientData clientData, G3_Runtime *rt, int argc,
                TCL_Char **argv);

int sensNodeAccel(ClientData clientData, G3_Runtime *rt, int argc,
                  TCL_Char **argv);

int sensNodePressure(ClientData clientData, G3_Runtime *rt, int argc,
                     TCL_Char **argv);

int sensSectionForce(ClientData clientData, G3_Runtime *rt, int argc,
                     TCL_Char **argv);

int sensitivityAlgorithm(ClientData clientData, G3_Runtime *rt, int argc,
                         TCL_Char **argv);

int sensitivityIntegrator(ClientData clientData, G3_Runtime *rt, int argc,
                          TCL_Char **argv);
// AddingSensitivity:END ///////////////////////////////////////////////////

int getNumElements(ClientData clientData, G3_Runtime *rt, int argc,
                   TCL_Char **argv);

int getEleClassTags(ClientData clientData, G3_Runtime *rt, int argc,
                    TCL_Char **argv);

int getEleLoadClassTags(ClientData clientData, G3_Runtime *rt, int argc,
                        TCL_Char **argv);

int getEleLoadTags(ClientData clientData, G3_Runtime *rt, int argc,
                   TCL_Char **argv);

int getEleLoadData(ClientData clientData, G3_Runtime *rt, int argc,
                   TCL_Char **argv);

int startTimer(ClientData clientData, G3_Runtime *rt, int argc,
               TCL_Char **argv);

int stopTimer(ClientData clientData, G3_Runtime *rt, int argc,
              TCL_Char **argv);

int rayleighDamping(ClientData clientData, G3_Runtime *rt, int argc,
                    TCL_Char **argv);

int modalDamping(ClientData clientData, G3_Runtime *rt, int argc,
                 TCL_Char **argv);

int modalDampingQ(ClientData clientData, G3_Runtime *rt, int argc,
                  TCL_Char **argv);

int setElementRayleighDampingFactors(ClientData clientData, G3_Runtime *rt,
                                     int argc, TCL_Char **argv);

int addRegion(ClientData clientData, G3_Runtime *rt, int argc,
              TCL_Char **argv);

int sectionForce(ClientData clientData, G3_Runtime *rt, int argc,
                 TCL_Char **argv);

int sectionDeformation(ClientData clientData, G3_Runtime *rt, int argc,
                       TCL_Char **argv);

int sectionStiffness(ClientData clientData, G3_Runtime *rt, int argc,
                     TCL_Char **argv);

int sectionFlexibility(ClientData clientData, G3_Runtime *rt, int argc,
                       TCL_Char **argv);

int sectionLocation(ClientData clientData, G3_Runtime *rt, int argc,
                    TCL_Char **argv);

int sectionWeight(ClientData clientData, G3_Runtime *rt, int argc,
                  TCL_Char **argv);

int basicDeformation(ClientData clientData, G3_Runtime *rt, int argc,
                     TCL_Char **argv);

int basicForce(ClientData clientData, G3_Runtime *rt, int argc,
               TCL_Char **argv);

int basicStiffness(ClientData clientData, G3_Runtime *rt, int argc,
                   TCL_Char **argv);

// added: Chris McGann, U.Washington for initial state analysis of nDMaterials
int InitialStateAnalysis(ClientData clientData, G3_Runtime *rt, int argc,
                         TCL_Char **argv);

int totalCPU(ClientData clientData, G3_Runtime *rt, int argc,
             TCL_Char **argv);

int solveCPU(ClientData clientData, G3_Runtime *rt, int argc,
             TCL_Char **argv);

int accelCPU(ClientData clientData, G3_Runtime *rt, int argc,
             TCL_Char **argv);

int numFact(ClientData clientData, G3_Runtime *rt, int argc,
            TCL_Char **argv);

int numIter(ClientData clientData, G3_Runtime *rt, int argc,
            TCL_Char **argv);

int systemSize(ClientData clientData, G3_Runtime *rt, int argc,
               TCL_Char **argv);

int elementActivate(ClientData clientData, G3_Runtime *rt, int argc,
                    TCL_Char **argv);
int elementDeactivate(ClientData clientData, G3_Runtime *rt, int argc,
                      TCL_Char **argv);
