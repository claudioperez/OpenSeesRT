#include <ID.h>
#include <tcl.h>
#include <assert.h>
#include <elementAPI.h>
#include <G3_Logging.h>
#include <BeamIntegration.h>
#include <BasicModelBuilder.h>

extern int OPS_ResetInput(ClientData clientData, Tcl_Interp *interp, int cArg,
                          int mArg, TCL_Char ** const argv, Domain *domain,
                          TclBuilder *builder);

extern void *OPS_LobattoBeamIntegration(int &integrationTag, ID &secTags);
extern void *OPS_LegendreBeamIntegration(int &integrationTag, ID &secTags);
extern void *OPS_NewtonCotesBeamIntegration(int &integrationTag, ID &secTags);
extern void *OPS_RadauBeamIntegration(int &integrationTag, ID &secTags);
extern void *OPS_TrapezoidalBeamIntegration(int &integrationTag, ID &secTags);
extern void *OPS_CompositeSimpsonBeamIntegration(int &integrationTag,
                                                 ID &secTags);
extern void *OPS_UserDefinedBeamIntegration(int &integrationTag, ID &secTags);
extern void *OPS_FixedLocationBeamIntegration(int &integrationTag, ID &secTags);
extern void *OPS_LowOrderBeamIntegration(int &integrationTag, ID &secTags);
extern void *OPS_MidDistanceBeamIntegration(int &integrationTag, ID &secTags);
extern void *OPS_UserHingeBeamIntegration(int &integrationTag, ID &secTags);
extern void *OPS_HingeMidpointBeamIntegration(int &integrationTag, ID &secTags);
extern void *OPS_HingeRadauBeamIntegration(int &integrationTag, ID &secTags);
extern void *OPS_HingeRadauTwoBeamIntegration(int &integrationTag, ID &secTags);
extern void *OPS_HingeEndpointBeamIntegration(int &integrationTag, ID &secTags);

extern int
TclCommand_addBeamIntegration(ClientData clientData, Tcl_Interp *interp,
                              int argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;

  if (argc < 2) {
    opserr << "WARNING: want beamIntegration type tag...\n";
    return TCL_ERROR;
  }

  OPS_ResetInput(clientData, interp, 2, argc, argv, nullptr, nullptr);

  int iTag;
  ID secTags;
  BeamIntegration *bi = 0;
  if (strcmp(argv[1], "Lobatto") == 0) {
    bi = (BeamIntegration *)OPS_LobattoBeamIntegration(iTag, secTags);
  } else if (strcmp(argv[1], "Legendre") == 0) {
    bi = (BeamIntegration *)OPS_LegendreBeamIntegration(iTag, secTags);
  } else if (strcmp(argv[1], "NewtonCotes") == 0) {
    bi = (BeamIntegration *)OPS_NewtonCotesBeamIntegration(iTag, secTags);
  } else if (strcmp(argv[1], "Radau") == 0) {
    bi = (BeamIntegration *)OPS_RadauBeamIntegration(iTag, secTags);
  } else if (strcmp(argv[1], "Trapezoidal") == 0) {
    bi = (BeamIntegration *)OPS_TrapezoidalBeamIntegration(iTag, secTags);
  } else if (strcmp(argv[1], "CompositeSimpson") == 0) {
    bi = (BeamIntegration *)OPS_CompositeSimpsonBeamIntegration(iTag, secTags);
  } else if (strcmp(argv[1], "UserDefined") == 0) {
    bi = (BeamIntegration *)OPS_UserDefinedBeamIntegration(iTag, secTags);
  } else if (strcmp(argv[1], "FixedLocation") == 0) {
    bi = (BeamIntegration *)OPS_FixedLocationBeamIntegration(iTag, secTags);
  } else if (strcmp(argv[1], "LowOrder") == 0) {
    bi = (BeamIntegration *)OPS_LowOrderBeamIntegration(iTag, secTags);
  } else if (strcmp(argv[1], "MidDistance") == 0) {
    bi = (BeamIntegration *)OPS_MidDistanceBeamIntegration(iTag, secTags);
  } else if (strcmp(argv[1], "UserHinge") == 0) {
    bi = (BeamIntegration *)OPS_UserHingeBeamIntegration(iTag, secTags);
  } else if (strcmp(argv[1], "HingeMidpoint") == 0) {
    bi = (BeamIntegration *)OPS_HingeMidpointBeamIntegration(iTag, secTags);
  } else if (strcmp(argv[1], "HingeRadau") == 0) {
    bi = (BeamIntegration *)OPS_HingeRadauBeamIntegration(iTag, secTags);
  } else if (strcmp(argv[1], "HingeRadauTwo") == 0) {
    bi = (BeamIntegration *)OPS_HingeRadauTwoBeamIntegration(iTag, secTags);
  } else if (strcmp(argv[1], "HingeEndpoint") == 0) {
    bi = (BeamIntegration *)OPS_HingeEndpointBeamIntegration(iTag, secTags);
  } else {
    opserr << "WARNING: integration type " << argv[1] << " is unknown\n";
    return TCL_ERROR;
  }

  if (bi == 0) {
    opserr << "WARNING: failed to create beam integration\n";
    return TCL_ERROR;
  }

  BeamIntegrationRule *rule = new BeamIntegrationRule(iTag, bi, secTags);
  if (rule == 0) {
    opserr << "WARNING: failed to create beam integration\n";
    delete bi;
    return TCL_ERROR;
  }

//if (OPS_addBeamIntegrationRule(rule) == false) {
  if (!builder->addRegistryObject("BeamIntegrationRule", iTag, rule)) {
    opserr << G3_ERROR_PROMPT << "could not add BeamIntegrationRule.";
    delete rule;
    return TCL_ERROR;
  }

  return TCL_OK;
}

