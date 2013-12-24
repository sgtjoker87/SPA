
/*
*/

#ifndef SPA_CLEANMESH_HXX
#define SPA_CLEANMESH_HXX

#include "spacommon/spamesh.h"

#include "cleaner/fmucwrap.hxx"
#include "cleaner/fmufhs.hxx"

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
);

#endif
