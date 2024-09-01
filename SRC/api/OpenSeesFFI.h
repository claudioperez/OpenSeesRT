#define ISW_INIT 0
#define ISW_COMMIT 1
#define ISW_REVERT 2
#define ISW_FORM_TANG_AND_RESID 3
#define ISW_FORM_MASS 4
#define ISW_REVERT_TO_START 5
#define ISW_DELETE 6

#define ISW_SET_RESPONSE 7
#define ISW_GET_RESPONSE 8

#define OPS_UNIAXIAL_MATERIAL_TYPE 1
#define OPS_SECTION2D_TYPE 2
#define OPS_SECTION3D_TYPE 3
#define OPS_PLANESTRESS_TYPE 4
#define OPS_PLANESTRAIN_TYPE 5
#define OPS_THREEDIMENSIONAL_TYPE 6
#define OPS_SECTION_TYPE 7

struct modState {
    double time;
    double dt;
};

typedef struct modState modelState;

typedef void (*matFunct)(struct matObject*, modelState*, double* strain, double* tang, double* stress, int* isw, int* error);

struct matObject {
    int tag;
    int matType;
    int nParam;
    int nState;
    double* theParam;
    double* cState;
    double* tState;
    matFunct matFunctPtr;
    void* matObjectPtr;
};

typedef struct matObject matObj;

typedef void (*eleFunct)(struct eleObject*, modelState*, double* tang, double* resid, int* isw, int* error);

struct eleObject {
    int tag;
    int nNode;
    int nDOF;
    int nParam;
    int nState;
    int nMat;
    int* node;
    double* param;
    double* cState;
    double* tState;
    matObj** mats;
    eleFunct eleFunctPtr;
};

typedef struct eleObject eleObj;


#define OPS_InvokeMaterial ops_invokematerial_
#define OPS_InvokeMaterialDirectly ops_invokematerialdirectly_
