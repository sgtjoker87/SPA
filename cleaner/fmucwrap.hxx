/*

*/

#ifndef FMUCWRAP_HXX
#define FMUCWRAP_HXX

#include "samcommon/sammesh.h"
#include "samcommon/sam_ptrtofunc.h"
#include "cleaner/fmuvector.hxx"
#include "cleaner/fmu.x"
#include "cleaner/fmufhs.hxx"
#include <cstddef>


class fmuMaster;

int fmuCleanInterface( int idFace,
                       int              iNumLoops,
                       fmu_TzContourLoop *azContourLoops,
                       bool             fAllowTrias,
                       PtrToSMGradeFunc pSMGradeFunc,
                       sam_TzNodeList   *pzINodeList,
                       sam_TzElemList   *pzIElemList,
                       sam_TzNodeList   *pzONodeList,
                       sam_TzElemList   *pzOElemList,
                       int              *nElem);

int fmuCleanUnpack( int                 iNumLoops,
                    fmu_TzContourLoop   *azContourLoops,
                    sam_TzNodeList      *pzINodeList,
                    sam_TzElemList      *pzElemList,
                    int                 *iFirstElem,
                    fmuMaster           *pCMaster );

int fmuProcessBoundary( int                 iNumLoops,
                        fmu_TzContourLoop   *azContourLoops,
                        NodeTable           *pzNodeLink,
                        fmuMaster           *pCMaster );

int fmuPlaneNormal( fmuMaster&, fmuVector& );

int fmuCleanPack( fmuMaster         &CMaster, 
                  int               iFirstElem,
                  sam_TzNodeList    *pzINodeList,
                  sam_TzNodeList    *pzONodeList,
                  sam_TzElemList    *pzOElemList,
                  int               *nElem );

void fmuCopyBoundaryNodes( int                  iNumLoops,
                           fmu_TzContourLoop    *azContourLoops,
                           sam_TzNodeList       *pzONodeList,
                           int                  *piRet);

#endif
