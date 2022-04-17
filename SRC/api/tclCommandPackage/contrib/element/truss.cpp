OPS_Export void *
OPS_TrussElement(G3_Runtime *rt=0)
{
  Element *theElement = 0;

  int numRemainingArgs = OPS_GetNumRemainingInputArgs(rt);

  if (numRemainingArgs < 4) {
    opserr << "Invalid Args want: element Truss $tag $iNode $jNode $sectTag <-rho $rho> <-cMass $flag> <-doRayleigh $flag>\n";
    opserr << " or: element Truss $tag $iNode $jNode $A $matTag <-rho $rho> <-cMass $flag> <-doRayleigh $flag>\n";
    return 0;	
  }

  if (numRemainingArgs == 4 || numRemainingArgs == 6 || numRemainingArgs == 8 || numRemainingArgs == 10)
    return 0; // it's a TrussSection

  int iData[3];
  double A = 0.0;
  double rho = 0.0;
  int matTag = 0;
  int doRayleigh = 0; // by default rayleigh not done
  int cMass = 0;      // by default use lumped mass matrix
  int ndm = G3_getNDM(rt);

  int numData = 3;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer (tag, iNode, jNode) in element Truss " << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetDouble(&numData, &A) != 0) {
    opserr << "WARNING: Invalid A: element Truss " << iData[0] << 
      " $iNode $jNode $A $matTag <-rho $rho> <-cMass $flag> <-doRayleigh $flag>\n";
    return 0;	
  }

  numData = 1;
  if (OPS_GetInt(&numData, &matTag) != 0) {
    opserr << "WARNING: Invalid matTag: element Truss " << iData[0] << 
      " $iNode $jNode $A $matTag <-rho $rho> <-cMass $flag> <-doRayleigh $flag>\n";
    return 0;
  }

  UniaxialMaterial *theUniaxialMaterial = G3_getUniaxialMaterialInstance(rt,matTag);
    
  if (theUniaxialMaterial == 0) {
    opserr << "WARNING: Invalid material not found element Truss " << iData[0] << " $iNode $jNode $A " << 
      matTag << " <-rho $rho> <-cMass $flag> <-doRayleigh $flag>\n";
    return 0;
  }
  
  numRemainingArgs -= 5;
  while (numRemainingArgs > 1) {
    const char *argvS = OPS_GetString();
  
    if (strcmp(argvS,"-rho") == 0) {
      numData = 1;
      if (OPS_GetDouble(&numData, &rho) != 0) {
	opserr << "WARNING Invalid rho in element Truss " << iData[0] << 
	  " $iNode $jNode $A $matTag <-rho $rho> <-cMass $flag> <-doRayleigh $flag>\n";
	return 0;
      }
    } else if (strcmp(argvS,"-cMass") == 0) {
      numData = 1;
      if (OPS_GetInt(&numData, &cMass) != 0) {
	opserr << "WARNING: Invalid cMass in element Truss " << iData[0] << 
	  " $iNode $jNode $A $matTag <-rho $rho> <-cMass $flag> <-doRayleigh $flag>\n";
	return 0;
      }
    } else if (strcmp(argvS,"-doRayleigh") == 0) {
      numData = 1;
      if (OPS_GetInt(&numData, &doRayleigh) != 0) {
	opserr << "WARNING: Invalid doRayleigh in element Truss " << iData[0] << 
	  " $iNode $jNode $A $matTag <-rho $rho> <-cMass $flag> <-doRayleigh $flag>\n";
	return 0;
      }
    } else {
      opserr << "WARNING: Invalid option " << argvS << "  in: element Truss " << iData[0] << 
	" $iNode $jNode $A $matTag <-rho $rho> <-cMass $flag> <-doRayleigh $flag>\n";
      return 0;
    }      
    numRemainingArgs -= 2;
  }

  // now create the Truss
  theElement = new Truss(iData[0], ndm, iData[1], iData[2], *theUniaxialMaterial, A, rho, doRayleigh, cMass);

  if (theElement == 0) {
    opserr << "WARNING: out of memory: element Truss " << iData[0] << 
      " $iNode $jNode $A $matTag <-rho $rho> <-cMass $flag> <-doRayleigh $flag>\n";
  }

  return theElement;
}

