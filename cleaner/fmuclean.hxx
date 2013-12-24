/*

*/



#ifndef FMUCLEAN_HXX
#define FMUCLEAN_HXX

#include "cleaner/nodarray.hxx"
#include "cleaner/fmuvector.hxx"

#include "interface_sam/ISizeMap.h"
#include "samcommon/sam_ptrtofunc.h"

#include <cstdio>

#define fmu_kCLEAN_BAD_INPUT     -1
#define fmu_kCLEAN_NO_ACTION     0
#define fmu_kCLEAN_MESH_CHANGED  1
#define fmu_kCLEAN_ABORT         2
#define fmu_kCLEAN_BACKSTEP      3
#define fmu_kCLEAN_NO_MEMORY     6401

class fmuMaster;

class CleanTool
{
public:
    virtual ~CleanTool();
    CleanTool();

    void SetSizeMap(PtrToSMGradeFunc pSMGradeFunc);

protected:
    
    void smooth_mesh( nodClassDynArray* );
      
    nodClass* one_middle_node( nodClass*,nodClass* );
      
    void draw_mesh();
    void draw_open();
    void draw_close();

protected:

    fmuMaster       *pCMaster;
    fmuVector       CNormal;
    PtrToSMGradeFunc pSMGradeFunc;
    bool            lDrawCleaning;

private:
    CleanTool( const CleanTool& );
    CleanTool& operator=( const CleanTool& );
};

#endif
