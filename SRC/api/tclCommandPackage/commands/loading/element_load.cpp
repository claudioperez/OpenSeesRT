
int
TclCommand_addElementalLoad(ClientData clientData, Tcl_Interp *interp, int argc,
                         TCL_Char **argv)
{
  // ensure the destructor has not been called -
  if (theTclBuilder == 0) {
    opserr << "WARNING current builder has been destroyed - eleLoad\n";
    return TCL_ERROR;
  }

  if (theTclLoadPattern == 0) {
    opserr << "WARNING no active load pattern - eleLoad\n";
    return TCL_ERROR;
  }

  int ndm = theTclBuilder->getNDM();
  ElementalLoad *theLoad = 0;

  ID theEleTags(0,16);

  // we first create an ID containing the ele tags of all elements
  // for which the load applies.
  int count = 1;
  int doneEle = 0;
  int eleCount = 0;
  while (doneEle == 0 && count < argc) {
    if (strcmp(argv[count],"-ele") == 0) {
      count ++;
      int eleStart = count;
      int eleEnd = 0;
      int eleID;
      while (count < argc && eleEnd == 0) {
        if (Tcl_GetInt(interp, argv[count], &eleID) != TCL_OK)
          eleEnd = count;
        else
          count++;
      }
      if (eleStart != eleEnd) {
        for (int i=eleStart; i<eleEnd; i++) {
          Tcl_GetInt(interp, argv[i], &eleID);
          theEleTags[eleCount++] = eleID;
        }
      }
    }
    else if (strcmp(argv[count],"-range") == 0) {
      count ++;
      int eleStart, eleEnd;
      if (Tcl_GetInt(interp, argv[count], &eleStart) != TCL_OK) {
        opserr << "WARNING eleLoad -range invalid eleStart " << argv[count] << "\n"; return TCL_ERROR;
      }
      count++;
      if (Tcl_GetInt(interp, argv[count], &eleEnd) != TCL_OK) {
        opserr << "WARNING eleLoad -range invalid eleEnd " << argv[count] << "\n"; return TCL_ERROR;
      }
      count++;
      for (int i=eleStart; i<=eleEnd; i++)
        theEleTags[eleCount++] = i;
    } else
      doneEle = 1;
  }
  // we then create the load
  if (strcmp(argv[count],"-type") != 0) {
    opserr << "WARNING eleLoad - expecting -type option but got "
         << argv[count] << endln;
    return TCL_ERROR;
  }
  count++;
  if (strcmp(argv[count],"-beamUniform") == 0 ||
      strcmp(argv[count],"beamUniform") == 0){
    count++;
    if (ndm == 2) {
      double wta;
      double waa = 0.0;
      if (count >= argc || Tcl_GetDouble(interp, argv[count], &wta) != TCL_OK) {
        opserr << "WARNING eleLoad - invalid wt for beamUniform \n";
        return TCL_ERROR;
      }
      count++;
      if (count < argc && Tcl_GetDouble(interp, argv[count], &waa) != TCL_OK) {
        opserr << "WARNING eleLoad - invalid wa for beamUniform \n";
        return TCL_ERROR;
      }
      double aL = 0.0;
      double bL = 1.0;
      count++;
      if (count < argc && Tcl_GetDouble(interp, argv[count], &aL) != TCL_OK) {
        opserr << "WARNING eleLoad - invalid aOverL for beamUniform \n";
        return TCL_ERROR;
      }
      count++;
      if (count < argc && Tcl_GetDouble(interp, argv[count], &bL) != TCL_OK) {
        opserr << "WARNING eleLoad - invalid bOverL for beamUniform \n";
        return TCL_ERROR;
      }
      double wab = waa;
      double wtb = wta;
      count++;
      if (count < argc && Tcl_GetDouble(interp, argv[count], &wtb) != TCL_OK) {
        opserr << "WARNING eleLoad - invalid wt for beamUniform \n";
        return TCL_ERROR;
      }
      count++;
      if (count < argc && Tcl_GetDouble(interp, argv[count], &wab) != TCL_OK) {
        opserr << "WARNING eleLoad - invalid wa for beamUniform \n";
        return TCL_ERROR;
      }
      for (int i=0; i<theEleTags.Size(); i++) {
        if (aL > 0.0 || bL < 1.0 || wta != wtb || waa != wab)
          theLoad = new Beam2dPartialUniformLoad(eleLoadTag, wta, wtb, waa, wab, aL, bL, theEleTags(i)); else theLoad = new Beam2dUniformLoad(eleLoadTag, wta, waa, theEleTags(i));

        if (theLoad == 0) {
          opserr << "WARNING eleLoad - out of memory creating load of type " << argv[count] ;
          return TCL_ERROR;
        }

        // get the current pattern tag if no tag given in i/p
        int loadPatternTag = theTclLoadPattern->getTag();

        // add the load to the domain
        if (theTclDomain->addElementalLoad(theLoad, loadPatternTag) == false) {
          opserr << "WARNING eleLoad - could not add following load to domain:\n"; 
          opserr << theLoad; 
          delete theLoad; 
          return TCL_ERROR;
        }
        eleLoadTag++;
      }

      return 0;
    }
    else if (ndm == 3) {
      double wy, wz;
      double wx = 0.0;
      if (count >= argc || Tcl_GetDouble(interp, argv[count], &wy) != TCL_OK) {
        opserr << "WARNING eleLoad - invalid wy for beamUniform \n";
        return TCL_ERROR;
      }
      count++;
      if (count >= argc || Tcl_GetDouble(interp, argv[count], &wz) != TCL_OK) {
        opserr << "WARNING eleLoad - invalid wz for beamUniform \n";
        return TCL_ERROR;
      }
      count++;
      if (count < argc && Tcl_GetDouble(interp, argv[count], &wx) != TCL_OK) {
        opserr << "WARNING eleLoad - invalid wx for beamUniform \n";
        return TCL_ERROR;
      }
      double aL = 0.0;
      double bL = 1.0;
      count++;
      if (count < argc && Tcl_GetDouble(interp, argv[count], &aL) != TCL_OK) {
        opserr << "WARNING eleLoad - invalid aOverL for beamUniform \n";
        return TCL_ERROR;
      }
      count++;
      if (count < argc && Tcl_GetDouble(interp, argv[count], &bL) != TCL_OK) {
        opserr << "WARNING eleLoad - invalid bOverL for beamUniform \n";
        return TCL_ERROR;
      }

      for (int i=0; i<theEleTags.Size(); i++) {
        if (aL > 0.0 || bL < 1.0)
          theLoad = new Beam3dPartialUniformLoad(eleLoadTag, wy, wz, wx, aL, bL, theEleTags(i)); else theLoad = new Beam3dUniformLoad(eleLoadTag, wy, wz, wx, theEleTags(i));

        if (theLoad == 0) {
          opserr << "WARNING eleLoad - out of memory creating load of type " << argv[count] ; return TCL_ERROR;
        }

        // get the current pattern tag if no tag given in i/p
        int loadPatternTag = theTclLoadPattern->getTag();

        // add the load to the domain
        if (theTclDomain->addElementalLoad(theLoad, loadPatternTag) == false) {
          opserr << "WARNING eleLoad - could not add following load to domain:\n"; opserr << theLoad; delete theLoad; return TCL_ERROR;
        }
        eleLoadTag++;
      }

      return 0;

    }
    else {
      opserr << "WARNING eleLoad beamUniform currently only valid only for ndm=or 3\n"; return TCL_ERROR;
    }

  } else if (strcmp(argv[count],"-beamPoint") == 0 ||
             strcmp(argv[count],"beamPoint") == 0 ) {
    count++;
    if (ndm == 2) {
      double P, x;
      double N = 0.0;
      if (count >= argc || Tcl_GetDouble(interp, argv[count], &P) != TCL_OK) {
        opserr << "WARNING eleLoad - invalid P for beamPoint\n";
        return TCL_ERROR;
      }
      if (count+1 >= argc || Tcl_GetDouble(interp, argv[count+1], &x) != TCL_OK){ 
        opserr << "WARNING eleLoad - invalid xDivL for beamPoint\n"; 
        return TCL_ERROR;
      }
      if (count+2 < argc && Tcl_GetDouble(interp, argv[count+2], &N) != TCL_OK) { 
        opserr << "WARNING eleLoad - invalid N for beamPoint\n"; 
        return TCL_ERROR;
      }

      if (x < 0.0 || x > 1.0) {
        opserr << "WARNING eleLoad - invalid xDivL of " << x;
        opserr << " for beamPoint (valid range [0.0, 1.0]\n";
        return TCL_ERROR;
      }


      for (int i=0; i<theEleTags.Size(); i++) {
        theLoad = new Beam2dPointLoad(eleLoadTag, P, x, theEleTags(i), N);

        if (theLoad == 0) {
          opserr << "WARNING eleLoad - out of memory creating load of type " << argv[count] ; return TCL_ERROR;
        }

        // get the current pattern tag if no tag given in i/p
        int loadPatternTag = theTclLoadPattern->getTag();

        // add the load to the domain
        if (theTclDomain->addElementalLoad(theLoad, loadPatternTag) == false) {
          opserr << "WARNING eleLoad - could not add following load to domain:\n"; opserr << theLoad; delete theLoad; return TCL_ERROR;
        }
        eleLoadTag++;
      }

      return 0;

    }
    else if (ndm == 3) {
      double Py, Pz, x;
      double N = 0.0;
      if (count >= argc || Tcl_GetDouble(interp, argv[count], &Py) != TCL_OK) {
        opserr << "WARNING eleLoad - invalid Py for beamPoint\n";
        return TCL_ERROR;
      }
      if (count+1 >= argc || Tcl_GetDouble(interp, argv[count+1], &Pz) !=TCL_OK) { opserr << "WARNING eleLoad - invalid Pz  for beamPoint\n";
        return TCL_ERROR;
      }
      if (count+2 >= argc || Tcl_GetDouble(interp, argv[count+2], &x) != TCL_OK){ 
        opserr << "WARNING eleLoad - invalid xDivL for beamPoint\n";
        return TCL_ERROR;
      }
      if (count+3 < argc && Tcl_GetDouble(interp, argv[count+3], &N) != TCL_OK){ 
        opserr << "WARNING eleLoad - invalid N for beamPoint\n";
        return TCL_ERROR;
      }

      if (x < 0.0 || x > 1.0) {
        opserr << "WARNING eleLoad - invalid xDivL of " << x;
        opserr << " for beamPoint (valid range [0.0, 1.0]\n";
        return TCL_ERROR;
      }

      for (int i=0; i<theEleTags.Size(); i++) {
        theLoad = new Beam3dPointLoad(eleLoadTag, Py, Pz, x, theEleTags(i), N);

        if (theLoad == 0) {
          opserr << "WARNING eleLoad - out of memory creating load of type " << argv[count] ; return TCL_ERROR;
        }

        // get the current pattern tag if no tag given in i/p
        int loadPatternTag = theTclLoadPattern->getTag();

        // add the load to the domain
        if (theTclDomain->addElementalLoad(theLoad, loadPatternTag) == false) {
          opserr << "WARNING eleLoad - could not add following load to domain:\n"; 
          opserr << theLoad; 
          delete theLoad; 
          return TCL_ERROR;
        }
        eleLoadTag++;
      }
      return 0;
    }
    else {
      opserr << "WARNING eleLoad beamPoint type currently only valid only for ndm=2 or 3\n"; return TCL_ERROR;
    }
  }
  // Added Joey Yang UC Davis
  else if (strcmp(argv[count],"-BrickW") == 0) {

      for (int i=0; i<theEleTags.Size(); i++) {
        theLoad = new BrickSelfWeight(eleLoadTag, theEleTags(i));

        if (theLoad == 0) {
          opserr << "WARNING eleLoad - out of memory creating load of type " << argv[count] ; return TCL_ERROR;
        }

        // get the current pattern tag if no tag given in i/p
        int loadPatternTag = theTclLoadPattern->getTag();

        // add the load to the domain
        if (theTclDomain->addElementalLoad(theLoad, loadPatternTag) == false) {
          opserr << "WARNING eleLoad - could not add following load to domain:\n"; opserr << theLoad; delete theLoad; return TCL_ERROR;
        }
        eleLoadTag++;
      }
      return 0;
  }
  // Added: C.McGann, U.Washington
  else if ((strcmp(argv[count],"-surfaceLoad") == 0) ||
           (strcmp(argv[count],"-SurfaceLoad") == 0)) { 
    count++; 
    for (int i=0; i<theEleTags.Size(); i++) { 
      theLoad = new SurfaceLoader(eleLoadTag, theEleTags(i));

          if (theLoad == 0) {
              opserr << "WARNING eleLoad - out of memory creating load of type " << argv[count] ; return TCL_ERROR;
          }

              // get the current pattern tag if no tag given in i/p
              int loadPatternTag = theTclLoadPattern->getTag();

              // add the load to the domain
              if (theTclDomain->addElementalLoad(theLoad, loadPatternTag) == false) { opserr << "WARNING eleLoad - could not add following load to domain:\n "; 
                opserr << theLoad; 
                delete theLoad; 
                return TCL_ERROR;
              }
              eleLoadTag++;
      }
      return 0;
  }
  // Added: C.McGann, U.Washington
  else if ((strcmp(argv[count],"-selfWeight") == 0) || 
           (strcmp(argv[count],"-SelfWeight") == 0)) { 
    count++;
    double xf = 0.0, yf = 0.0, zf = 0.0;
    if (Tcl_GetDouble(interp, argv[count], &xf) != TCL_OK) {
          opserr << "WARNING eleLoad - invalid xFactor " << argv[count] << " for -selfWeight\n"; return TCL_ERROR;
        }
    if (Tcl_GetDouble(interp, argv[count+1], &yf) != TCL_OK) {
          opserr << "WARNING eleLoad - invalid yFactor " << argv[count+1] << " for -selfWeight\n"; return TCL_ERROR;
        }
    if (count+2 < argc) { // adding to stop seg faults
      if (Tcl_GetDouble(interp, argv[count+2], &zf) != TCL_OK) {
            opserr << "WARNING eleLoad - invalid zFactor " << argv[count+2] << " for -selfWeight\n"; return TCL_ERROR;
      }
    }

    for (int i=0; i<theEleTags.Size(); i++) {
      theLoad = new SelfWeight(eleLoadTag, xf, yf, zf, theEleTags(i));

      if (theLoad == 0) {
        opserr << "WARNING eleLoad - out of memory creating load of type " << argv[count] ; return TCL_ERROR;
      }

      // get the current pattern tag if no tag given in i/p
      int loadPatternTag = theTclLoadPattern->getTag();

      // add the load to the domain
      if (theTclDomain->addElementalLoad(theLoad, loadPatternTag) == false) {
        opserr << "WARNING eleLoad - could not add following load to domain:\n"; opserr << theLoad; delete theLoad; return TCL_ERROR;
      }
      eleLoadTag++;
    }
    return 0;
  }

 ///---------------------- Adding identifier for ThermalAction : [END] [SIF]----------------------------------------------//

  //-----------------Adding tcl command for shell thermal action, 2013..[Begin]--------------------- 
  else if (strcmp(argv[count], "-shellThermal") == 0) { 
    count++;
          //so far three kinds of temperature distribution
          //(1) 9 temperature points, i.e. 8 layers
          //(2) 5 temperature points, i.e. 4 layers
          //(3) 2 temperature points, i.e. 1 layers: linear or uniform

          double t1, locY1, t2, locY2; //t3, locY3, t4, locY4, t5, locY5, t6,
locY6, t7, locY7, t8, locY8, t9, locY9;
          // 9 temperature points are given,i.e. 8 layers are defined; Also the
9 corresponding vertical coordinate is given.
          // the temperature at each fiber is obtained by interpolating of
temperatures at the nearby temperature points.
          //Start to add source file
          if (strcmp(argv[count], "-source") == 0) {
                  if (strcmp(argv[count + 1], "-node") != 0) {
                          count++;

                          const char *pwd = getInterpPWD(interp);
                          // simulationInfo.addInputFile(argv[count], pwd);
                          TimeSeries* theSeries = new
PathTimeSeriesThermal(eleLoadTag, argv[count]);

                          count++;

                          double RcvLoc1, RcvLoc2;
                          if (argc - count == 2) {

                                  if (Tcl_GetDouble(interp, argv[count],
&RcvLoc1) != TCL_OK) { opserr << "WARNING eleLoad - invalid single loc  " <<
argv[count] << " for -beamThermal\n"; return TCL_ERROR;
                                  }
                                  if (Tcl_GetDouble(interp, argv[count + 1],
&RcvLoc2) != TCL_OK) { opserr << "WARNING eleLoad - invalid single loc  " <<
argv[count + 1] << " for -beamThermal\n"; return TCL_ERROR;
                                  }

                          }
                          else {
                                  opserr << "WARNING eleLoad - invalid input for
-shellThermal\n";
                          }

                          for (int i = 0; i<theEleTags.Size(); i++) {
                                  theLoad = new ShellThermalAction(eleLoadTag,
RcvLoc1, RcvLoc2, theSeries, theEleTags(i)); if (theLoad == 0) { opserr <<
"WARNING eleLoad - out of memory creating load of type " << argv[count]; return
TCL_ERROR;
                                  }

                                  // get the current pattern tag if no tag given
in i/p int loadPatternTag = theTclLoadPattern->getTag();

                                  // add the load to the domain
                                  if (theTclDomain->addElementalLoad(theLoad,
loadPatternTag) == false) { opserr << "WARNING eleLoad - could not add following
load to domain:\n "; opserr << theLoad; delete theLoad; return TCL_ERROR;
                                  }
                                  eleLoadTag++;
                          }
                  }
                  //if not using nodal thermal action input
                  else {
                          for (int i = 0; i<theEleTags.Size(); i++) {
                                  theLoad = new ShellThermalAction(eleLoadTag,
theEleTags(i)); if (theLoad == 0) { opserr << "WARNING eleLoad - out of memory
creating load of type " << argv[count]; return TCL_ERROR;
                                  }

                                  // get the current pattern tag if no tag given
in i/p int loadPatternTag = theTclLoadPattern->getTag();

                                  // add the load to the domain
                                  if (theTclDomain->addElementalLoad(theLoad,
loadPatternTag) == false) { opserr << "WARNING eleLoad - could not add following
load to domain:\n "; opserr << theLoad; delete theLoad; return TCL_ERROR;
                                  }
                                  eleLoadTag++;
                          }//end of for loop
                          return 0;
                  }//end of <if(strcmp(argv[count+1],"-node") = 0)>
          }
          //end of the interface for importing temperature data from external
file else
          {
                  if (argc - count == 18) {
                          double indata[18];
                          double BufferData;

                          for (int i = 0; i<18; i++) {
                                  if (Tcl_GetDouble(interp, argv[count],
&BufferData) != TCL_OK) { opserr << "WARNING eleLoad - invalid data " <<
argv[count] << " for -beamThermal 3D\n"; return TCL_ERROR;
                                  }
                                  indata[i] = BufferData;
                                  count++;
                          }

                          //temp1,loc1,temp2,loc2...temp9,loc9

                          for (int i = 0; i<theEleTags.Size(); i++) {
                                  theLoad = new ShellThermalAction(eleLoadTag,
                                          indata[0], indata[1], indata[2],
indata[3], indata[4], indata[5], indata[6], indata[7], indata[8], indata[9],
indata[10], indata[11], indata[12], indata[13], indata[14], indata[15],
                                          indata[16], indata[17],
theEleTags(i));


                                  if (theLoad == 0) {
                                          opserr << "WARNING eleLoad - out of
memory creating load of type " << argv[count]; return TCL_ERROR;
                                  }

                                  // get the current pattern tag if no tag given
in i/p int loadPatternTag = theTclLoadPattern->getTag();

                                  // add the load to the domain
                                  if (theTclDomain->addElementalLoad(theLoad,
loadPatternTag) == false) { opserr << "WARNING eleLoad - could not add following
load to domain:\n "; opserr << theLoad; delete theLoad; return TCL_ERROR;
                                  }
                                  eleLoadTag++;
                          }
                          return 0;
                  }
                  // 5 temperatures are given, i.e. 4 layers are defined.
                  else if (argc - count == 10) {
                          double indata[10];
                          double BufferData;

                          for (int i = 0; i<10; i++) {
                                  if (Tcl_GetDouble(interp, argv[count],
&BufferData) != TCL_OK) { opserr << "WARNING eleLoad - invalid data " <<
argv[count] << " for -beamThermal 3D\n"; return TCL_ERROR;
                                  }
                                  indata[i] = BufferData;
                                  count++;
                          }

                          //temp1,loc1,temp2,loc2...temp5,loc5

                          for (int i = 0; i<theEleTags.Size(); i++) {
                                  theLoad = new ShellThermalAction(eleLoadTag,
                                          indata[0], indata[1], indata[2],
indata[3], indata[4], indata[5], indata[6], indata[7], indata[8], indata[9],
theEleTags(i));


                                  if (theLoad == 0) {
                                          opserr << "WARNING eleLoad - out of
memory creating load of type " << argv[count]; return TCL_ERROR;
                                  }

                                  // get the current pattern tag if no tag given
in i/p int loadPatternTag = theTclLoadPattern->getTag();

                                  // add the load to the domain
                                  if (theTclDomain->addElementalLoad(theLoad,
loadPatternTag) == false) { opserr << "WARNING eleLoad - could not add following
load to domain:\n "; opserr << theLoad; delete theLoad; return TCL_ERROR;
                                  }
                                  eleLoadTag++;
                          }
                          return 0;
                  }
                  // two temperature is given,
                  //if the two temperatures are equal,i.e. uniform Temperature
change in element
                  //if the two temperatures are different,i.e. linear
Temperature change in element else if (argc - count == 4) { if
(Tcl_GetDouble(interp, argv[count], &t1) != TCL_OK) { opserr << "WARNING eleLoad
- invalid T1 " << argv[count] << " for -shellThermal\n"; return TCL_ERROR;
                          }

                          if (Tcl_GetDouble(interp, argv[count + 1], &locY1) !=
TCL_OK) { opserr << "WARNING eleLoad - invalid LocY1 " << argv[count + 1] << "
for -shellThermal\n"; return TCL_ERROR;
                          }
                          if (Tcl_GetDouble(interp, argv[count + 2], &t2) !=
TCL_OK) { opserr << "WARNING eleLoad - invalid T2 " << argv[count] << " for
-shellThermal\n"; return TCL_ERROR;
                          }

                          if (Tcl_GetDouble(interp, argv[count + 3], &locY2) !=
TCL_OK) { opserr << "WARNING eleLoad - invalid LocY2 " << argv[count + 1] << "
for -shellThermal\n"; return TCL_ERROR;
                          }

                          for (int i = 0; i<theEleTags.Size(); i++) {
                                  theLoad = new ShellThermalAction(eleLoadTag,
                                          t1, locY1, t2, locY2, theEleTags(i));


                                  if (theLoad == 0) {
                                          opserr << "WARNING eleLoad - out of
memory creating load of type " << argv[count]; return TCL_ERROR;
                                  }

                                  // get the current pattern tag if no tag given
in i/p int loadPatternTag = theTclLoadPattern->getTag();

                                  // add the load to the domain
                                  if (theTclDomain->addElementalLoad(theLoad,
loadPatternTag) == false) { opserr << "WARNING eleLoad - could not add following
load to domain:\n "; opserr << theLoad; delete theLoad; return TCL_ERROR;
                                  }
                                  eleLoadTag++;
                          }
                          return 0;
                  }
                  //finish the temperature arguments
                  else {
                          opserr << "WARNING eleLoad -shellThermalLoad invalid
number of temperature aguments,/n looking for 0, 1, 2 or 4 arguments.\n";
                  }
          }
          //end of if(recieved argument is not "source" or direct temperature
input)//Liming,2014
  }
  //-----------------Adding tcl command for shell thermal action,
2013..[End]-----------------------

  else if (strcmp(argv[count], "-ThermalWrapper") == 0 || strcmp(argv[count],
"-thermalWrapper") == 0) {

          count++;
          Vector loc = 0;;
          ID NodalThermal = 0;
          int numNodal = 0;

          if (strcmp(argv[count], "-nodeLoc") == 0 || strcmp(argv[count],
"nodeLoc") == 0)
          {
                  count++;
                  numNodal = (argc - count) / 2;
                  loc.resize(numNodal);
                  NodalThermal.resize(numNodal);

                  for (int i = 0; i<numNodal; i++) {

                          double Dblloc;
                          int NodalTtag;

                          if (Tcl_GetDouble(interp, argv[count + i * 2 + 1],
&Dblloc) != TCL_OK) { opserr << "WARNING NodalLoad - invalid loc  " <<
argv[count] << " for NodalThermalAction\n"; return TCL_ERROR;
                          }

                          if (Tcl_GetInt(interp, argv[count + 2 * i],
&NodalTtag) != TCL_OK) { opserr << "WARNING invalid nodeId: " << argv[1]; return
TCL_ERROR;
                          }

                          loc(i) = Dblloc; NodalThermal(i) = NodalTtag;
                  }
                  //end of for loop over numNodal
          }
          //end of nodelLoc
          else if (strcmp(argv[count], "-node") == 0 || strcmp(argv[count],
"node") == 0)
          {
                  count++;
                  numNodal = argc - count;
                  NodalThermal.resize(numNodal);

                  for (int i = 0; i<numNodal; i++) {

                          int NodalTtag;

                          if (Tcl_GetInt(interp, argv[count + i], &NodalTtag) !=
TCL_OK) { opserr << "WARNING invalid nodeId: " << argv[1]; return TCL_ERROR;
                          }

                          NodalThermal(i) = NodalTtag;
                  }
                  //end of for loop over numNodal
          }
          //end of node tag

          //Obtain Pointers to NodalThermalAction;
          Node* theNode = 0;
          NodalThermalAction** theNodalThermals = 0;
          theNodalThermals = new NodalThermalAction*[numNodal];

          for (int i = 0; i<numNodal; i++) {
            theNode = theTclDomain->getNode(NodalThermal(i));
            theNodalThermals[i] = theNode->getNodalThermalActionPtr();
            if (theNodalThermals[i] == 0) {
              opserr << "WARNING:: An empty nodalThermalAction detected for
ThermalActionWrapper" << endln; return TCL_ERROR;
            }
          }

          for (int i = 0; i<theEleTags.Size(); i++) {
            if (numNodal == 2) {
              theLoad = new ThermalActionWrapper(eleLoadTag, theEleTags(i),
theNodalThermals[0], theNodalThermals[1]);
            }
            else if (numNodal == 3) {
              theLoad = new ThermalActionWrapper(eleLoadTag, theEleTags(i),
theNodalThermals[0], theNodalThermals[1], theNodalThermals[2]);
            }
            else if (numNodal == 4) {
              theLoad = new ThermalActionWrapper(eleLoadTag, theEleTags(i),
theNodalThermals[0], theNodalThermals[1], theNodalThermals[2],
theNodalThermals[3]);
            }
            else if (numNodal == 5) {
              theLoad = new ThermalActionWrapper(eleLoadTag, theEleTags(i),
theNodalThermals[0], theNodalThermals[1], theNodalThermals[2],
theNodalThermals[3], theNodalThermals[4]);
            }
            else if (numNodal == 6) {
              theLoad = new ThermalActionWrapper(eleLoadTag, theEleTags(i),
theNodalThermals[0], theNodalThermals[1], theNodalThermals[2],
theNodalThermals[3], theNodalThermals[4], theNodalThermals[5]);
            }


            if (theLoad == 0) {
              opserr << "WARNING eleLoad - out of memory creating load of type "
<< argv[count]; return TCL_ERROR;
                  }

                  if (loc != 0)
                          ((ThermalActionWrapper*)theLoad)->setRatios(loc);

                  // get the current pattern tag if no tag given in i/p
                  int loadPatternTag = theTclLoadPattern->getTag();

                  // add the load to the domain
                  if (theTclDomain->addElementalLoad(theLoad, loadPatternTag) ==
false) { opserr << "WARNING eleLoad - could not add following load to domain:\n
"; opserr << theLoad; delete theLoad; return TCL_ERROR;
                  }
                  eleLoadTag++;
          }//end of for loop
          return 0;
  }
  //------------------------end  of using
ThermalActionWrapper--------------------------
  //-----------------Adding tcl command for beam thermal action(2D&3D),
2013..[Begin]---------------

  else if (strcmp(argv[count], "-beamThermal") == 0) {
          count++;
          //For two dimensional model
          if (ndm == 2) {
                  // The thermal action can be defined with external data file,
which is identified with "-source"
                  // The external file itself can either be elemental data or
nodal data, the latter is identified with "-node" if (strcmp(argv[count],
"-source") == 0) {

                          if (strcmp(argv[count + 1], "-node") == 0) {
                                  for (int i = 0; i<theEleTags.Size(); i++) {
                                          theLoad = new
Beam2dThermalAction(eleLoadTag, theEleTags(i)); if (theLoad == 0) { opserr <<
"WARNING eleLoad - out of memory creating load of type " << argv[count]; return
TCL_ERROR;
                                          }

                                          // get the current pattern tag if no
tag given in i/p int loadPatternTag = theTclLoadPattern->getTag();

                                          // add the load to the domain
                                          if
(theTclDomain->addElementalLoad(theLoad, loadPatternTag) == false) { opserr <<
"WARNING eleLoad - could not add following load to domain:\n "; opserr <<
theLoad; delete theLoad; return TCL_ERROR;
                                          }
                                          eleLoadTag++;
                                  }//end of for loop
                                  return 0;

                          }
                          //end of <if(strcmp(argv[count+1],"-node") != 0)>
                          else {
                                  count++;
                                  const char *pwd = getInterpPWD(interp);
                                  // simulationInfo.addInputFile(argv[count],
pwd); TimeSeries* theSeries = new PathTimeSeriesThermal(eleLoadTag,
argv[count]);

                                  count++;
                                  Vector locs(9);
                                  //---------------for recieving 2 arguments
                                  if (argc - count == 2) {
                                          double RcvLoc1, RcvLoc2;
                                          if (Tcl_GetDouble(interp, argv[count],
&RcvLoc1) != TCL_OK) { opserr << "WARNING eleLoad - invalid single loc  " <<
argv[count] << " for -beamThermal\n"; return TCL_ERROR;
                                          }
                                          if (Tcl_GetDouble(interp, argv[count +
1], &RcvLoc2) != TCL_OK) { opserr << "WARNING eleLoad - invalid single loc  " <<
argv[count + 1] << " for -beamThermal\n"; return TCL_ERROR;
                                          }

                                          locs(0) = RcvLoc1; locs(8) = RcvLoc2;
                                          for (int i = 1; i<8; i++) {
                                                  locs(i) = locs(0) - i*(locs(0)
- locs(8)) / 8;
                                          }
                                  }
                                  //----------------for recieving 9 arguments
                                  else if (argc - count == 9) {

                                          int ArgStart = count;
                                          int ArgEnd = argc;
                                          double data;

                                          if (ArgStart != ArgEnd) {
                                                  for (int i = ArgStart;
i<ArgEnd; i++) { Tcl_GetDouble(interp, argv[i], &data); locs(i - ArgStart) =
data;
                                                  }
                                          }

                                  }
                                  //end of recieving 9 arguments
                                  else {
                                          opserr << "WARNING eleLoad - invalid
input for -beamThermal\n";
                                  }
#ifdef _DEBUG
                                  opserr << "TclSafeBuilder:: locs" <<
locs << endln; #endif for (int i = 0; i<theEleTags.Size(); i++) { theLoad = new
Beam2dThermalAction(eleLoadTag, locs, theSeries, theEleTags(i)); if (theLoad ==
0) { opserr << "WARNING eleLoad - out of memory creating load of type " <<
argv[count]; return TCL_ERROR;
                                          }

                                          // get the current pattern tag if no
tag given in i/p int loadPatternTag = theTclLoadPattern->getTag();

                                          // add the load to the domain
                                          if
(theTclDomain->addElementalLoad(theLoad, loadPatternTag) == false) { opserr <<
"WARNING eleLoad - could not add following load to domain:\n "; opserr <<
theLoad; delete theLoad; return TCL_ERROR;
                                          }
                                          eleLoadTag++;
                                  }//end of for loop
                                  return 0;
                          }//end of <if(strcmp(argv[count+1],"-node") = 0)>
                           //--------------------------end for
beam2DThermalAction with time series ----------------------------------------
                  }
                  else
                  {
                          //(1) 9 temperature points, i.e. 8 layers
                          //(2) 5 temperature points, i.e. 4 layers
                          //(3) 2 temperature points, i.e. 1 layers: linear or
uniform
                          //double t1, locY1, t2, locY2, t3, locY3, t4, locY4,
t5, locY5, t6, locY6, t7, locY7, t8, locY8, t9, locY9; double Temp[9]; double
Loc[9];
                          // 9 temperature points are given,i.e. 8 layers are
defined; Also the 9 corresponding vertical coordinate is given.
                          // the temperature at each fiber is obtained by
interpolating of temperatures at the nearby temperature points. if (argc - count
== 18) { double indata[18]; double BufferData;

                                  for (int i = 0; i<18; i++) {
                                          if (Tcl_GetDouble(interp, argv[count],
&BufferData) != TCL_OK) { opserr << "WARNING eleLoad - invalid data " <<
argv[count] << " for -beamThermal 3D\n"; return TCL_ERROR;
                                          }
                                          indata[i] = BufferData;
                                          count++;
                                  }

                                  for (int i = 0; i<9; i++) {
                                          Temp[i] = indata[2 * i];
                                          Loc[i]  = indata[2 * i + 1];
                                  }

                          }

                          // 5 temperatures are given, i.e. 4 layers are
defined. else if (argc - count == 10) {

                                  double indata[10];
                                  double BufferData;

                                  for (int i = 0; i<10; i++) {
                                          if (Tcl_GetDouble(interp, argv[count],
&BufferData) != TCL_OK) { opserr << "WARNING eleLoad - invalid data " <<
argv[count] << " for -beamThermal 3D\n"; return TCL_ERROR;
                                          }
                                          indata[i] = BufferData;
                                          count++;
                                  }

                                  Temp[0] = indata[0]; Temp[2] = indata[2];
Temp[4] = indata[4]; Temp[6] = indata[6]; Temp[8] = indata[8]; Loc[0]  =
indata[1]; Loc[2]  = indata[3]; Loc[4]  = indata[5]; Loc[6]  = indata[7]; Loc[8]
= indata[9];

                                  for (int i = 1; i<5; i++) {
                                          Temp[2 * i - 1] = (Temp[2 * i - 2] +
Temp[2 * i]) / 2; Loc[2 * i - 1]  = (Loc[2 * i - 2] + Loc[2 * i]) / 2;
                                  }
                          }
                          //End for 5 inputs
                          // two temperature is given,
                          //if the two temperatures are equal,i.e. uniform
Temperature change in element
                          //if the two temperatures are different,i.e. linear
Temperature change in element else if (argc - count == 4) { double indata[4];
                                  double BufferData;

                                  for (int i = 0; i<4; i++) {
                                          if (Tcl_GetDouble(interp, argv[count],
&BufferData) != TCL_OK) { opserr << "WARNING eleLoad - invalid data " <<
argv[count] << " for -beamThermal 3D\n"; return TCL_ERROR;
                                          }
                                          indata[i] = BufferData;
                                          count++;
                                  }

                                  Temp[0] = indata[0]; Temp[8] = indata[2];
                                  Loc[0]  = indata[1]; Loc[8]  = indata[3];
                                  for (int i = 1; i<8; i++) {
                                          Temp[i] = Temp[0] - i*(Temp[0] -
Temp[8]) / 8; Loc[i] = Loc[0] - i*(Loc[0] - Loc[8]) / 8;
                                  }

                          }
                          //end for 2 inputs

                          for (int i = 0; i<theEleTags.Size(); i++) {
                                  theLoad = new Beam2dThermalAction(eleLoadTag,
                                          Temp[0], Loc[0], Temp[1], Loc[1],
                                          Temp[2], Loc[2], Temp[3], Loc[3],
                                          Temp[4], Loc[4], Temp[5], Loc[5],
                                          Temp[6], Loc[6], Temp[7], Loc[7],
                      Temp[8], Loc[8], theEleTags(i));

                                  if (theLoad == 0) {
                                          opserr << "WARNING eleLoad - out of
memory creating load of type beamThermal"; return TCL_ERROR;
                                  }

                                  // get the current pattern tag if no tag given
in i/p int loadPatternTag = theTclLoadPattern->getTag();

                                  // add the load to the domain
                                  if (theTclDomain->addElementalLoad(theLoad,
loadPatternTag) == false) { opserr << "WARNING eleLoad - could not add following
load to domain:\n "; opserr << theLoad; delete theLoad; return TCL_ERROR;
                                  }
                                  eleLoadTag++;
                          }

                          return 0;
                          //-End of Adding BeamThermalAction defined with direct
input.
                  }
                  //end for no sourced pattern
          }
          // End of for the if (ndm==2)
          else if (ndm == 3)
          {
                  //so far three kinds of temperature distribution
                  double t1, locY1, t2, locY2, t3, locY3, t4, locY4, t5, locY5,
                          t6, t7, locZ1, t8, t9, locZ2, t10, t11, locZ3, t12,
t13, locZ4, t14, t15, locZ5;

                  // the temperature at each fiber is obtained by interpolating
of temperatures at the nearby temperature points. if (strcmp(argv[count],
"-source") == 0) { count++; if (strcmp(argv[count], "-node") == 0) { for (int i
= 0; i<theEleTags.Size(); i++) { theLoad = new Beam3dThermalAction(eleLoadTag,
theEleTags(i)); if (theLoad == 0) { opserr << "WARNING eleLoad - out of memory
creating load of type " << argv[count]; return TCL_ERROR;
                                          }
                                          // get the current pattern tag if no
tag given in i/p int loadPatternTag = theTclLoadPattern->getTag();
                                          // add the load to the domain
                                          if
(theTclDomain->addElementalLoad(theLoad, loadPatternTag) == false) { opserr <<
"WARNING eleLoad - could not add following load to domain:\n "; opserr <<
theLoad; delete theLoad; return TCL_ERROR;
                                          }
                                          eleLoadTag++;
                                  }	//end of loop tf all elements defined
                                  return 0;
                          }//end for defing thermal action with nodal input
                          else {
                                  const char *pwd = getInterpPWD(interp);
                                  // simulationInfo.addInputFile(argv[count],
pwd); count++; bool using2Ddata = false;

                                  double RcvLoc1, RcvLoc2, RcvLoc3, RcvLoc4;
                                  TimeSeries* theSeries;

                                  if (argc - count == 4) {
                                          theSeries = new
PathTimeSeriesThermal(eleLoadTag, argv[count - 1], 15); using2Ddata = false;

                                          if (Tcl_GetDouble(interp, argv[count],
&RcvLoc1) != TCL_OK) { opserr << "WARNING eleLoad - invalid single loc  " <<
argv[count] << " for -beamThermal\n"; return TCL_ERROR;
                                          }
                                          if (Tcl_GetDouble(interp, argv[count +
1], &RcvLoc2) != TCL_OK) { opserr << "WARNING eleLoad - invalid single loc  " <<
argv[count + 1] << " for -beamThermal\n"; return TCL_ERROR;
                                          }
                                          if (Tcl_GetDouble(interp, argv[count +
2], &RcvLoc3) != TCL_OK) { opserr << "WARNING eleLoad - invalid single loc  " <<
argv[count + 2] << " for -beamThermal\n"; return TCL_ERROR;
                                          }
                                          if (Tcl_GetDouble(interp, argv[count +
3], &RcvLoc4) != TCL_OK) { opserr << "WARNING eleLoad - invalid single loc  " <<
argv[count + 3] << " for -beamThermal\n"; return TCL_ERROR;
                                          }

                                          //end for recieving input
                                          for (int i = 0; i<theEleTags.Size();
i++) {

                                                  theLoad = new
Beam3dThermalAction(eleLoadTag, RcvLoc1, RcvLoc2, RcvLoc3, RcvLoc4, theSeries,
theEleTags(i));

                                                  if (theLoad == 0) {
                                                          opserr << "WARNING
eleLoad - out of memory creating load of type " << argv[count]; return
TCL_ERROR;
                                                  }

                                                  // get the current pattern tag
if no tag given in i/p int loadPatternTag = theTclLoadPattern->getTag();

                                                  // add the load to the domain
                                                  if
(theTclDomain->addElementalLoad(theLoad, loadPatternTag) == false) { opserr <<
"WARNING eleLoad - could not add following load to domain:\n "; opserr <<
theLoad; delete theLoad; return TCL_ERROR;
                                                  }
                                                  eleLoadTag++;
                                          }
                                          return 0;

                                  }
                                  //end of defining 15 data points with external
file else if (argc - count == 2 || argc - count == 9) {
                                          // for receiving data which has the
similiar structure as 2D beam section Vector locs(9); using2Ddata = true;
                                          TimeSeries* theSeries = new
PathTimeSeriesThermal(eleLoadTag, argv[count - 1], 9); if (argc - count == 2) {

                                                  double RcvLoc1, RcvLoc2;
                                                  if (Tcl_GetDouble(interp,
argv[count], &RcvLoc1) != TCL_OK) { opserr << "WARNING eleLoad - invalid single
loc  " << argv[count] << " for -beamThermal\n"; return TCL_ERROR;
                                                  }
                                                  if (Tcl_GetDouble(interp,
argv[count + 1], &RcvLoc2) != TCL_OK) { opserr << "WARNING eleLoad - invalid
single loc  " << argv[count + 1] << " for -beamThermal\n"; return TCL_ERROR;
                                                  }

                                                  locs(0) = RcvLoc1; locs(8) =
RcvLoc2; for (int i = 1; i<8; i++) { locs(i) = locs(0) - i*(locs(0) - locs(8)) /
8;
                                                  }
                                          }
                                          //end of receiving 2 data points
                                          else {
                                                  double indata[9];
                                                  double BufferData;
                                                  for (int i = 0; i<9; i++) {
                                                          if
(Tcl_GetDouble(interp, argv[count], &BufferData) != TCL_OK) { opserr << "WARNING
eleLoad - invalid data " << argv[count] << " for -beamThermal 3D\n"; return
TCL_ERROR;
                                                          }
                                                          indata[i] =
BufferData; count++;
                                                  }

                                                  for (int i = 0; i<9; i++) {
                                                          locs(i) = indata[i];
                                                  }

                                          }
                                          //end of receiving 9data points

                                          for (int i = 0; i<theEleTags.Size();
i++) { theLoad = new Beam3dThermalAction(eleLoadTag, locs, theSeries,
theEleTags(i)); if (theLoad == 0) { opserr << "WARNING eleLoad - out of memory
creating load of type " << argv[count]; return TCL_ERROR;
                                                  }

                                                  // get the current pattern tag
if no tag given in i/p int loadPatternTag = theTclLoadPattern->getTag();

                                                  // add the load to the domain
                                                  if
(theTclDomain->addElementalLoad(theLoad, loadPatternTag) == false) { opserr <<
"WARNING eleLoad - could not add following load to domain:\n "; opserr <<
theLoad; delete theLoad; return TCL_ERROR;
                                                  }
                                                  eleLoadTag++;
                                          }//end of for loop
                                          return 0;
                                  }

                                  else {
                                          opserr << "WARNING eleLoad - invalid
input for -beamThermal\n";
                                  }

                          }//end for source beam element temperature data
-source -filePath

                  }
                  //-------------------------end for importing temp data from
external files --source ------------------------- else {

                          //double t1, locY1, t2, locY2, t3, locY3, t4, locY4,
t5, locY5,
                          //t6, t7, locZ1, t8, t9, locZ2, t10,t11, locZ3, t12,
t13, locZ4, t14,t15, locZ5;

                          if (argc - count == 25) {
                                  double indata[25];
                                  double BufferData;

                                  for (int i = 0; i<25; i++) {
                                          if (Tcl_GetDouble(interp, argv[count],
&BufferData) != TCL_OK) { opserr << "WARNING eleLoad - invalid data " <<
argv[count] << " for -beamThermal 3D\n"; return TCL_ERROR;
                                          }
                                          indata[i] = BufferData;
                                          count++;
                                  }

                                  for (int i = 0; i<theEleTags.Size(); i++) {
                                          theLoad = new
Beam3dThermalAction(eleLoadTag, indata[0], indata[1], indata[2], indata[3],
indata[4], indata[5], indata[6], indata[7], indata[8], indata[9], indata[10],
indata[11], indata[12], indata[13], indata[14], indata[15], indata[16],
indata[17], indata[18], indata[19], indata[20], indata[21], indata[22],
indata[23], indata[24], theEleTags(i)); if (theLoad == 0) { opserr << "WARNING
eleLoad - out of memory creating load of type " << argv[count]; return
TCL_ERROR;
                                          }
                                          // get the current pattern tag if no
tag given in i/p int loadPatternTag = theTclLoadPattern->getTag();

                                          // add the load to the domain
                                          if
(theTclDomain->addElementalLoad(theLoad, loadPatternTag) == false) { opserr <<
"WARNING eleLoad - could not add following load to domain:\n "; opserr <<
theLoad; delete theLoad; return TCL_ERROR;
                                          }
                                          eleLoadTag++;
                                  }
                                  return 0;
                          }
                          //end of  if (argc-count == 25){
                          else if (argc - count == 4) {

                                  if (Tcl_GetDouble(interp, argv[count], &t1) !=
TCL_OK) { opserr << "WARNING eleLoad - invalid T1 " << argv[count] << " for
-beamThermal\n"; return TCL_ERROR;
                                  }
                                  if (Tcl_GetDouble(interp, argv[count + 1],
&locY1) != TCL_OK) { opserr << "WARNING eleLoad - invalid LocY1 " << argv[count
+ 1] << " for -beamThermal\n"; return TCL_ERROR;
                                  }
                                  if (Tcl_GetDouble(interp, argv[count + 2],
&t5) != TCL_OK) { opserr << "WARNING eleLoad - invalid T1 " << argv[count] << "
for -beamThermal\n"; return TCL_ERROR;
                                  }
                                  if (Tcl_GetDouble(interp, argv[count + 3],
&locY5) != TCL_OK) { opserr << "WARNING eleLoad - invalid LocY1 " << argv[count
+ 1] << " for -beamThermal\n"; return TCL_ERROR;
                                  }

                                  locY2 = locY1 + (locY5 - locY1) / 4;
                                  locY3 = locY1 + (locY5 - locY1) / 2;
                                  locY4 = locY1 + 3 * (locY5 - locY1) / 4;
                                  t2 = t1 + (t5 - t1) / 4;
                                  t3 = t1 + (t5 - t1) / 2;
                                  t4 = t1 + 3 * (t5 - t1) / 4;
                                  locZ1 = locZ2 = locZ3 = locZ4 = locZ5 = 0;
                                  t6 = t7 = t8 = t9 = t10 = 0;
                                  t11 = t12 = t13 = t14 = t15 = 0;

                                  for (int i = 0; i<theEleTags.Size(); i++) {
                                          theLoad = new
Beam3dThermalAction(eleLoadTag, t1, locY1, t2, locY2, t3, locY3, t4, locY4, t5,
locY5, t6, t7, locZ1, t8, t9, locZ2, t10, t11, locZ3, t12, t13, locZ4, t14, t15,
locZ5, theEleTags(i)); if (theLoad == 0) { opserr << "WARNING eleLoad - out of
memory creating load of type " << argv[count]; return TCL_ERROR;
                                          }
                                          // get the current pattern tag if no
tag given in i/p int loadPatternTag = theTclLoadPattern->getTag();

                                          // add the load to the domain
                                          if
(theTclDomain->addElementalLoad(theLoad, loadPatternTag) == false) { opserr <<
"WARNING eleLoad - could not add following load to domain:\n "; opserr <<
theLoad; delete theLoad; return TCL_ERROR;
                                          }
                                          eleLoadTag++;
                                  }
                                  return 0;
                          }
                          //end of  if (argc-count == 4){
                          else {
                                  opserr << "WARNING eleLoad
Beam3dThermalAction: invalid number of temperature aguments,/n looking for
arguments for Temperatures and cordinates.\n";
                          }
                  } //end for not sourced pattern
          }//else if ndm=3
  }//else if '-beamThermal'

   //--Adding identifier for ThermalAction:[END] [SIF]--//



  // Added by Scott R. Hamilton   - Stanford
  else if (strcmp(argv[count],"-beamTemp") == 0) {
    count++;
    if (ndm == 2) {
      double temp1, temp2, temp3, temp4;

      // Four temps given, Temp change at top node 1, bottom node 1, top node 2,
bottom node 2. if (argc-count == 4){ if (Tcl_GetDouble(interp, argv[count],
&temp1) != TCL_OK) { opserr << "WARNING eleLoad - invalid Ttop1 " << argv[count]
<< " for -beamTemp\n"; return TCL_ERROR;
        }

        if (Tcl_GetDouble(interp, argv[count+1],&temp2 ) != TCL_OK) {
          opserr << "WARNING eleLoad - invalid Tbot1 " << argv[count+1] << " for
-beamTemp\n"; return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[count+2], &temp3) != TCL_OK) {
          opserr << "WARNING eleLoad - invalid Ttop2 " << argv[count+1] << " for
-beamTemp\n"; return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[count+3], &temp4) != TCL_OK) {
          opserr << "WARNING eleLoad - invalid Tbot2 " << argv[count+1] << " for
-beamTemp\n"; return TCL_ERROR;
        }

        for (int i=0; i<theEleTags.Size(); i++) {
          theLoad = new Beam2dTempLoad(eleLoadTag, temp1, temp2, temp3, temp4, theEleTags(i));

          if (theLoad == 0) {
            opserr << "WARNING eleLoad - out of memory creating load of type " << argv[count] ; return TCL_ERROR;
          }

          // get the current pattern tag if no tag given in i/p
          int loadPatternTag = theTclLoadPattern->getTag();

          // add the load to the domain
          if (theTclDomain->addElementalLoad(theLoad, loadPatternTag) == false){
            opserr << "WARNING eleLoad - could not add following load to domain:\n ";
            opserr << theLoad;
            delete theLoad;
            return TCL_ERROR;
          }
          eleLoadTag++;
        }

        return 0;

      }
      // Two temps given, temp change at top, temp at bottom of element
      else if (argc-count == 2) {
        if (Tcl_GetDouble(interp, argv[count], &temp1) != TCL_OK) {
          opserr << "WARNING eleLoad - invalid Ttop " << argv[count] << " for -beamTemp\n"; return TCL_ERROR;
        }

        if (Tcl_GetDouble(interp, argv[count+1],&temp2 ) != TCL_OK) {
          opserr << "WARNING eleLoad - invalid Tbot " << argv[count+1] << " for -beamTemp\n"; return TCL_ERROR;
        }

        for (int i=0; i<theEleTags.Size(); i++) {
          theLoad = new Beam2dTempLoad(eleLoadTag, temp1, temp2, theEleTags(i));

          if (theLoad == 0) {
            opserr << "WARNING eleLoad - out of memory creating load of type " << argv[count] ; return TCL_ERROR;
          }

          // get the current pattern tag if no tag given in i/p
          int loadPatternTag = theTclLoadPattern->getTag();

          // add the load to the domain
          if (theTclDomain->addElementalLoad(theLoad, loadPatternTag) == false){
            opserr << "WARNING eleLoad - could not add following load to domain:\n ";
            opserr << theLoad;
            delete theLoad;
            return TCL_ERROR;
          }
          eleLoadTag++;
        }
      }
      // One twmp change give, uniform temp change in element
      else if (argc-count == 1) {
        if (Tcl_GetDouble(interp, argv[count],&temp1 ) != TCL_OK) {
          opserr << "WARNING eleLoad - invalid Tbot " << argv[count+1] << " for -beamTemp\n"; return TCL_ERROR;
        }
        theLoad=0;

        for (int i=0; i<theEleTags.Size(); i++) {
          theLoad = new Beam2dTempLoad(eleLoadTag, temp1, theEleTags(i));

          if (theLoad == 0) {
            opserr << "WARNING eleLoad - out of memory creating load of type " << argv[count] ; return TCL_ERROR;
          }

          // get the current pattern tag if no tag given in i/p
          int loadPatternTag = theTclLoadPattern->getTag();

          // add the load to the domain
          if (theTclDomain->addElementalLoad(theLoad, loadPatternTag) == false) { opserr << "WARNING eleLoad - could not add following load to domain:\n ";
            opserr << theLoad;
            delete theLoad;
            return TCL_ERROR;
          }
          eleLoadTag++;
        }

        return 0;
      }

      else {
        opserr << "WARNING eleLoad -beamTempLoad invalid number of temperature aguments,/n looking for 0, 1, 2 or 4 arguments.\n";
      }
    } else {
      opserr << "WARNING eleLoad -beamTempLoad type currently only valid only for ndm=2\n"; return TCL_ERROR;
    }
  }
  // if get here we have sucessfully created the load and added it to the domain
  return TCL_OK;
}
