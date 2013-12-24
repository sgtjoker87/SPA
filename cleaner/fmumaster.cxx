/*

*/

#define DEBUG_Master_dump 0

#include "cleaner/fmumaster.hxx"
#include "cleaner/nodclass.hxx"
#include "cleaner/fmuedge.hxx"
#include "cleaner/elmqclass.hxx"
#include "cleaner/elmtclass.hxx"

#include "wksp/CSamEnvVar.hxx"
#include "FileImportExport/samFileNameUtility.h"
#include "FileImportExport/snprintf_macro.h"

fmuMaster::fmuMaster() :
    iHighestNodeId( 0 ),
    lAbortFlag( false ),
    fAllowTrias( 0 ),
    m_idFace( 0 )
{
}

void
fmuMaster::dump( int iteration ) const
{
#if DEBUG_Master_dump

    const int idFace = this->FaceId();

    const size_t len = 128;
    char fn[len];
    SNPRINTF( fn, len, "DEBUG_Cleaner_masterDump_%u_%d.txt", (unsigned int)idFace, iteration );
    samFileNameUtility sfn( fn, CSamEnvVar::getSamReportPath() );
    FILE* ps = fopen( sfn.PathName(), "wbc" );
    if( ps == NULL )
        return;

    const nodClassDynArray& dynNodeInt = global_interior();
    int nNI = dynNodeInt.length();
    fprintf( ps, "------------------\n" );
    fprintf( ps, "interior nodes: %5d nodes\n", nNI );
    for( int in=0; in<nNI; ++in )
    {
        nodClass* pN = dynNodeInt[in];
        if( pN!=NULL )
            fprintf( ps, "%8d", pN->id() );
        else
            fprintf( ps, "    NULL" );
        if( (in+1)%10 == 0 || in==(nNI-1) )
            fprintf( ps, "\n" );
    }
    fprintf( ps, "\n" );
    fprintf( ps, "node details...\n" );
    for( int in=0; in<nNI; ++in )
    {
        nodClass* pN = dynNodeInt[in];
        if( pN==NULL )
            continue;
        pN->dump( ps );
    }

    const nodClassDynArrayDynArray& dynNodeBL = global_boundary();
    int nBL = dynNodeBL.length();

    fprintf( ps, "------------------\n" );
    fprintf( ps, "CBoundaryList: %d loops\n", nBL );
    for( int inl=0; inl<nBL; ++inl )
    {
        nodClassDynArray& dynNodeLoop = ( *(dynNodeBL[inl]) );
        int nN = dynNodeLoop.length();
        fprintf( ps, "..------------------\n" );
        fprintf( ps, "..loop: %4d... %5d nodes\n", inl, nN );
        for( int in=0; in<nN; ++in )
        {
            nodClass* pN = dynNodeLoop[in];
            if( pN!=NULL )
                fprintf( ps, "%8d", pN->id() );
            else
                fprintf( ps, "    NULL" );
            fprintf( ps, "%8d", pN->id() );
            if( (in+1)%10 == 0 || in==(nN-1) )
                fprintf( ps, "\n" );
        }
        fprintf( ps, "\n" );
        fprintf( ps, "node details...\n" );
        for( int in=0; in<nN; ++in )
        {
            nodClass* pN = dynNodeLoop[in];
            if( pN==NULL )
                continue;
            pN->dump( ps );
        }
    }

    const fmuEdgeDynArray& dynEdge = global_edges();
    size_t nE = dynEdge.length();
    fprintf( ps, "------------------\n" );
    fprintf( ps, "edge list: %8d\n", nE );
    for( size_t ie=0; ie<nE; ++ie )
    {
        fmuEdge* pE = dynEdge[ie];
        if( pE==NULL )
            continue;
        pE->dump( ps );
    }

    const elmTClassDynArray& dynTria = global_trias();
    size_t nT = dynTria.length();
    fprintf( ps, "------------------\n" );
    fprintf( ps, "tria list: %8d\n", nT );
    for( size_t it=0; it<nT; ++it )
    {
        elmTClass* pT = dynTria[it];
        if( pT==NULL )
            continue;
        pT->dump( ps );
    }

    const elmQClassDynArray& dynQuad = global_quads();
    size_t nQ = dynQuad.length();
    fprintf( ps, "------------------\n" );
    fprintf( ps, "quad list: %8d\n", nQ );
    for( size_t iq=0; iq<nQ; ++iq )
    {
        elmQClass* pQ = dynQuad[iq];
        if( pQ==NULL )
            continue;
        pQ->dump( ps );
    }

    fflush( ps );
    fclose( ps );

#endif
}
