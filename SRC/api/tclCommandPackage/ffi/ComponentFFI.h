
#ifdef __cplusplus

extern "C" int         OPS_GetNodeCrd(int* nodeTag, int* sizeData, double* data);
extern "C" int         OPS_GetNodeDisp(int* nodeTag, int* sizeData, double* data);
extern "C" int         OPS_GetNodeVel(int* nodeTag, int* sizeData, double* data);
extern "C" int         OPS_GetNodeAccel(int* nodeTag, int* sizeData, double* data);
extern "C" int         OPS_GetNodeIncrDisp(int* nodeTag, int* sizeData, double* data);
extern "C" int         OPS_GetNodeIncrDeltaDisp(int* nodeTag, int* sizeData, double* data);
extern "C" matObj*     OPS_GetMaterial(int* matTag, int* matType);
//extern "C" void      OPS_GetMaterialPtr(int *matTag, matObj *theRes);
extern "C" eleObj*     OPS_GetElement(int* eleTag);
extern "C" matObj*     OPS_GetMaterialType(char* type, int sizeType);
extern "C" eleObj*     OPS_GetElementType(char* type, int sizeType);
extern "C" int         OPS_AllocateElement(eleObject * theEle, int* matTags, int* matType);
extern "C" int         OPS_AllocateMaterial(matObject * theMat);

extern "C" int         OPS_InvokeMaterial(eleObject*, int*, modelState*, double*, double*, double*, int*);
extern "C" int         OPS_InvokeMaterialDirectly(matObject**, modelState*, double*, double*, double*, int*);
extern "C" int         OPS_InvokeMaterialDirectly2(matObject*, modelState*, double*, double*, double*, int*);

#else // __cplusplus

matObj* OPS_GetMaterial(int* matTag, int* matType);
void    OPS_GetMaterialPtr(int*, matObj*);
eleObj* OPS_GetElement(int*);
matObj* OPS_GetMaterialType(char* type, int sizeType);
eleObj* OPS_GetElementType(char*, int);
int     OPS_AllocateElement(eleObj*, int* matTags, int* maType);
int     OPS_AllocateMaterial(matObj*);

int    OPS_InvokeMaterial(struct eleObj*, int*, modelState*, double*, double*, double*, int*);
int    OPS_InvokeMaterialDirectly(matObj**, modelState*, double*, double*, double*, int*);
int    OPS_InvokeMaterialDirectly2(matObj*, modelState*, double*, double*, double*, int*);

int    OPS_GetNodeCrd(int* nodeTag, int* sizeData, double* data);
int    OPS_GetNodeDisp(int* nodeTag, int* sizeData, double* data);
int    OPS_GetNodeVel(int* nodeTag, int* sizeData, double* data);
int    OPS_GetNodeAcc(int* nodeTag, int* sizeData, double* data);
int    OPS_GetNodeIncrDisp(int* nodeTag, int* sizeData, double* data);
int    OPS_GetNodeIncrDeltaDisp(int* nodeTag, int* sizeData, double* data);
#endif
