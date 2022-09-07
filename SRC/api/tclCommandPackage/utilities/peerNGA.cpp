
extern int peerSearchNGA(const char *eq, const char *soilType,
                         const char *fault, const char *magLo,
                         const char *magHi, const char *distLo,
                         const char *distHi, const char *vsLo, const char *vsHi,
                         const char *pgaLo, const char *pgaHi,
                         const char *latSW, const char *latNE,
                         const char *lngSW, const char *lngNW,
                         StringContainer &recordNames);

int
peerNGA(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  StringContainer ngaRecordNames;
  const char *eq = 0;
  const char *soilType = 0;
  const char *fault = 0;
  const char *magLo = 0;
  const char *magHi = 0;
  const char *distLo = 0;
  const char *distHi = 0;
  const char *vsLo = 0;
  const char *vsHi = 0;
  const char *pgaLo = 0;
  const char *pgaHi = 0;
  const char *latSW = 0;
  const char *latNE = 0;
  const char *lngSW = 0;
  const char *lngNW = 0;

  int currentArg = 1;
  while (currentArg + 1 < argc) {
    if (strcmp(argv[currentArg], "-eq") == 0) {
      eq = argv[currentArg + 1];
    } else if (strcmp(argv[currentArg], "-fault") == 0) {
      fault = argv[currentArg + 1];
    } else if (strcmp(argv[currentArg], "-soil") == 0) {
      soilType = argv[currentArg + 1];
    } else if (strcmp(argv[currentArg], "-magLo") == 0) {
      magLo = argv[currentArg + 1];
    } else if (strcmp(argv[currentArg], "-magHi") == 0) {
      magHi = argv[currentArg + 1];
    } else if (strcmp(argv[currentArg], "-distLo") == 0) {
      distLo = argv[currentArg + 1];
    } else if (strcmp(argv[currentArg], "-distHi") == 0) {
      distHi = argv[currentArg + 1];
    } else if (strcmp(argv[currentArg], "-vsLo") == 0) {
      vsLo = argv[currentArg + 1];
    } else if (strcmp(argv[currentArg], "-vsHi") == 0) {
      vsHi = argv[currentArg + 1];
    } else if (strcmp(argv[currentArg], "-pgaLo") == 0) {
      pgaLo = argv[currentArg + 1];
    } else if (strcmp(argv[currentArg], "-pgaHi") == 0) {
      pgaHi = argv[currentArg + 1];
    } else if (strcmp(argv[currentArg], "-latSW") == 0) {
      latSW = argv[currentArg + 1];
    } else if (strcmp(argv[currentArg], "-latNE") == 0) {
      latNE = argv[currentArg + 1];
    } else if (strcmp(argv[currentArg], "-lngSW") == 0) {
      lngSW = argv[currentArg + 1];
    } else if (strcmp(argv[currentArg], "-lngNW") == 0) {
      lngNW = argv[currentArg + 1];
    }
    // unrecognized
    currentArg += 2;
  }

  peerSearchNGA(eq, soilType, fault, magLo, magHi, distLo, distHi, vsLo, vsHi,
                pgaLo, pgaHi, latSW, latNE, lngSW, lngNW, ngaRecordNames);

  int numStrings = ngaRecordNames.getNumStrings();
  for (int i = 0; i < numStrings; i++) {
    Tcl_AppendResult(interp, ngaRecordNames.getString(i), NULL);
    Tcl_AppendResult(interp, " ", NULL);
  }

  return TCL_OK;
}
