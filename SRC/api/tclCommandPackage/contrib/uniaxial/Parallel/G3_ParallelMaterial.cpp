#include <ParallelMaterial.h>
// #include <TclSafeBuilder.h>
#include <G3Parse.h>

UniaxialMaterial*
G3Parse_newParallelMaterial(G3_Runtime* rt, int argc, G3_Char** argv)
{
// else if (strcmp(argv[1],"Parallel") == 0) {
    if (argc < 4) {
        opserr << "WARNING insufficient arguments\n";
        printCommand(argc,argv);
        opserr << "Want: uniaxialMaterial Parallel tag? tag1? tag2? ...";
        opserr << " <-min min?> <-max max?>" << endln;
        return nullptr;
    }

    int tag;
    UniaxialMaterial* theMaterial = nullptr;

    if (G3Parse_getInt(rt, argv[2], &tag) != TCL_OK) {
        opserr << "WARNING invalid uniaxialMaterial Parallel tag" << endln;
        return nullptr;
    }

    int numMaterials = argc-3;
    
    if (numMaterials == 0) {
        opserr << "WARNING no component material(s) provided\n";
        opserr << "uniaxialMaterial Parallel: " << tag << endln;
        return nullptr;
    }

    // Create an array to hold pointers to component materials
    UniaxialMaterial **theMats = new UniaxialMaterial *[numMaterials];
    
    // For each material get the tag and ensure it exists in model already
    for (int i=0; i<numMaterials; i++) {
      int tagI;
      if (G3Parse_getInt(rt, argv[i+3], &tagI) != TCL_OK) {
        opserr << "WARNING invalid component tag\n";
        opserr << "uniaxialMaterial Parallel: " << tag << endln;
        return nullptr;
      }
      
      UniaxialMaterial *theMat = G3_getUniaxialMaterialInstance(rt, tagI);
      
      if (theMat == 0) {
        opserr << "WARNING component material does not exist\n";
        opserr << "Component material: " << argv[i+3]; 
        opserr << "\nuniaxialMaterial Parallel: " << tag << endln;
        delete [] theMats;
        return nullptr;
      } else
        theMats[i] = theMat;
    }
    
    // Parsing was successful, allocate the material
    theMaterial = new ParallelMaterial(tag, numMaterials, theMats);
    
    // Deallocate the temporary pointers
    delete [] theMats;
    return theMaterial;
}
