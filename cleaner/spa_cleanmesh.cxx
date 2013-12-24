/*

*/




#define DEBUG_WriteMeshToUnv    0

#include "cleaner/spa_cleanmesh.x"
#include "errorhandler/ErrorReporter.hxx"

#if DEBUG_WriteMeshToUnv
#include "wksp/CFilePrint.h"
#include "FileImportExport/CspaUnvWrite.h"
#endif


// Argument                   I/O               Description
// iNumLoops                  I                 No of loops (outer +
//                                              all inner) the mesh
//                                              area has
// azContourLoops             I                 Array of loop data structs
// fAllowTrias                I                 Quad Dominant Flag
// poSizeMap                  I                 Ptr to SizeMap
// pzINodeList                I                 Interior node list stored
//                                              as a vector
// pzIElemList                I                 Interior element list
//                                              stored as a vector
// pzONodeList                O                 Ptr to the node list of
//                                              the cleaned mesh
// pzOElemList                O                 Ptr to the elem list of
//                                              the cleaned mesh
//
// RETURNS   Error
//=======================================================================
int spa_CleanMesh
(
    int                 idFace,
    int                 iNumLoops,
    fmu_TzContourLoop   *azContourLoops,
    int                 fAllowTrias,
    PtrToSMGradeFunc    pSMGradeFunc,
    spa_TzNodeList      *pzINodeList,
    spa_TzElemList      *pzIElemList,
    spa_TzNodeList      *pzONodeList,
    spa_TzElemList      *pzOElemList
)
{
    static const char scFuncName[] = "cleaner/spa_cleanmesh";
    int nElm = 0, iRet = 0;

    try
    {

#if DEBUG_WriteMeshToUnv
        static const size_t sNameLen = 1024;
        static char sFileName[sNameLen];
        static size_t countMesh = 0;
        ++countMesh;
        SNPRINTF( sFileName, sNameLen, "debug_CleanMesh_input_%u.unv", (unsigned int)idFace );
        CspaUnvWrite unvI( sFileName );
        unvI.Write( azContourLoops, *pzINodeList, *pzIElemList );
#endif

        // =====================================================
        // Call the Cleaner
        // =====================================================
        iRet = fmuCleanInterface(   idFace,
                                    iNumLoops,
                                    azContourLoops,
                                    (fAllowTrias==0?false:true),
                                    pSMGradeFunc,
                                    pzINodeList,
                                    pzIElemList,
                                    pzONodeList,
                                    pzOElemList,
                                    &nElm
                                );

#if DEBUG_WriteMeshToUnv
        SNPRINTF( sFileName, sNameLen, "DEBUG_CleanMesh_output_%u.unv", (unsigned int)idFace );
        CspaUnvWrite unvO( sFileName );
        unvO.Write( *pzONodeList, *pzOElemList );
#endif

    }
    catch( int iRc )
    {
        iRet = iRc;
    }
    catch( ... )
    {
        iRet = (int)spa_eFAILED_TO_CLEAN_MESH;
    }

    if (iRet == (int)spa_eNOT_ENOUGH_MEMORY)
    {
        ErrorReporter::report_warning(scFuncName, spa_eNOT_ENOUGH_MEMORY,0,0);
        throw (iRet);
    }
    else if (iRet)
        ErrorReporter::report_warning(scFuncName, spa_eFAILED_TO_CLEAN_MESH,0,0);

    return(iRet);
}
