
int C_PathIndependent(G3_Runtime* rt, int argc, const char **argv)
{
   Tcl_Interp *interp = G3_getInterpreter(rt);

   if (strcmp(argv[1], "PathIndependent") == 0) {
      if (argc < 4) {
        opserr << "WARNING insufficient arguments\n";
        // printCommand(argc, argv);
        opserr << "Want: uniaxialMaterial PathIndependent tag? matTag?"
               << endln;
        return TCL_ERROR;
      }
      int tag, matTag;

      if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
        opserr << "WARNING invalid tag\n";
        opserr << "PathIndependent material: " << tag << endln;
        return TCL_ERROR;
      }

      if (Tcl_GetInt(interp, argv[3], &matTag) != TCL_OK) {
        opserr << "WARNING invalid matTag\n";
        opserr << "PathIndependent material: " << tag << endln;
        return TCL_ERROR;
      }

      UniaxialMaterial *material = G3_getUniaxialMaterialInstance(rt, matTag);

      if (material == 0) {
        opserr << "WARNING material does not exist\n";
        opserr << "material: " << matTag;
        opserr << "\nuniaxialMaterial PathIndependent: " << tag << endln;
        return TCL_ERROR;
      }

      theMaterial = new PathIndependentMaterial(tag, *material);
    }
}

