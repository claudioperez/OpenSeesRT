/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

// Description: This file contains the functions that will be called by
// the interpreter when the appropriate command name is specified.



Tcl_CmdProc setLoadConst;

Tcl_CmdProc setCreep;

Tcl_CmdProc getLoadFactor;

Tcl_CmdProc buildModel;


Tcl_CmdProc printModel;

// Tcl_CmdProc addRecorder;
Tcl_CmdProc TclAddRecorder;

Tcl_CmdProc addAlgoRecorder;

Tcl_CmdProc addDatabase;

Tcl_CmdProc playbackRecorders;

Tcl_CmdProc playbackAlgorithmRecorders;

Tcl_CmdProc groundExcitation;

Tcl_CmdProc removeObject;

Tcl_CmdProc eleForce;

Tcl_CmdProc localForce;

Tcl_CmdProc eleDynamicalForce;

Tcl_CmdProc eleResponse;

Tcl_CmdProc findID;

Tcl_CmdProc nodeDisp;

Tcl_CmdProc nodeReaction;

Tcl_CmdProc nodeUnbalance;

Tcl_CmdProc nodeEigenvector;

Tcl_CmdProc nodeCoord;

Tcl_CmdProc setNodeCoord;

Tcl_CmdProc updateElementDomain;

Tcl_CmdProc eleType;

Tcl_CmdProc eleNodes;

Tcl_CmdProc nodeBounds;

Tcl_CmdProc nodeVel;

Tcl_CmdProc setNodeVel;

Tcl_CmdProc setNodeDisp;

Tcl_CmdProc setNodeAccel;

Tcl_CmdProc nodeAccel;

Tcl_CmdProc nodeResponse;

Tcl_CmdProc calculateNodalReactions;

Tcl_CmdProc getNodeTags;

Tcl_CmdProc getEleTags;

Tcl_CmdProc fixedNodes;

Tcl_CmdProc fixedDOFs;

Tcl_CmdProc constrainedNodes;

Tcl_CmdProc constrainedDOFs;

Tcl_CmdProc retainedNodes;

Tcl_CmdProc retainedDOFs;

Tcl_CmdProc nodeDOFs;

Tcl_CmdProc nodeMass;

Tcl_CmdProc nodePressure;

Tcl_CmdProc getParamTags;

Tcl_CmdProc getParamValue;

// AddingSensitivity:BEGIN /////////////////////////////////////////////////
Tcl_CmdProc computeGradients;
Tcl_CmdProc sensNodeDisp;
Tcl_CmdProc sensLambda; // Abbas
Tcl_CmdProc sensNodeVel;
Tcl_CmdProc sensNodeAccel;
Tcl_CmdProc sensNodePressure;
Tcl_CmdProc sensSectionForce;
Tcl_CmdProc sensitivityAlgorithm;
// Tcl_CmdProc sensitivityIntegrator;
// AddingSensitivity:END ///////////////////////////////////////////////////

Tcl_CmdProc getNumElements;

Tcl_CmdProc getEleClassTags;
Tcl_CmdProc getEleLoadClassTags;
Tcl_CmdProc getEleLoadTags;
Tcl_CmdProc getEleLoadData;

// Tcl_CmdProc startTimer;
// Tcl_CmdProc stopTimer;
// Tcl_CmdProc setTime;
// Tcl_CmdProc getTime;

Tcl_CmdProc rayleighDamping;

Tcl_CmdProc modalDamping;

Tcl_CmdProc modalDampingQ;

Tcl_CmdProc setElementRayleighDampingFactors;

Tcl_CmdProc addRegion;

Tcl_CmdProc sectionForce;

Tcl_CmdProc sectionDeformation;

Tcl_CmdProc sectionStiffness;

Tcl_CmdProc sectionFlexibility;

Tcl_CmdProc sectionLocation;

Tcl_CmdProc sectionWeight;

Tcl_CmdProc basicDeformation;

Tcl_CmdProc basicForce;

Tcl_CmdProc basicStiffness;

// added: Chris McGann, U.Washington for initial state analysis of nDMaterials
Tcl_CmdProc InitialStateAnalysis;

Tcl_CmdProc elementActivate;
Tcl_CmdProc elementDeactivate;

// by SAJalali
Tcl_CmdProc OPS_recorderValue;


