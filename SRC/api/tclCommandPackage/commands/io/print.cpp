
void OPS_printNDMaterial(OPS_Stream &s, int flag) {
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\"ndMaterials\": [\n";
        MapOfTaggedObjectsIter theObjects = theNDMaterialObjects.getIter();
        theObjects.reset();
        TaggedObject *theObject;
        int count = 0;
        int numComponents = theNDMaterialObjects.getNumComponents();
        while ((theObject = theObjects()) != 0) {
            NDMaterial *theMaterial = (NDMaterial *)theObject;
            theMaterial->Print(s, flag);
            if (count < numComponents - 1)
                s << ",\n";
            count++;
        }
        s << "\n\t\t]";
    }
}


