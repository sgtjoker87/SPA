/*

*/


#include "cleaner/elmbase.hxx"
#include "cleaner/nodclass.hxx"
#include "cleaner/fmuedge.hxx"
#include "cleaner/elmqarray.hxx"
#include "cleaner/elmtarray.hxx"

#include <SamTypes.h>

#include <cstddef>
#include <cstdio>

//   
//   HISTORY
//     
//     10/07/98 14:09 <               st qa> Rel:ms7      SCCA:39240  Delta:7.3
//   
//     09/12/97 14:11 <        Paul Kinney > Rel:ms6      SCCA:22677  Delta:6.2

// SDRC FE quad Member Function Descriptor Block -----
//
// elmBase::elmBase
// class - elmBase
// date 04-22-97  SDRC/P. Kinney
// ------------------------------------------

elmBase::~elmBase()
{
}

elmBase::elmBase() :
   mId              ( 0 ),
   mindexCurrentNode( 0 ),
   mIsFrozen        ( false )
{
}


void elmBase::dump( FILE* ps ) const
{
    static const bool printPointerValue = false;

    if( ps == NULL )
        return;

    // matches SAM_ElementShape
    static const char* saiShapeNames[9] =
    {
        "SAM_UNKNOWN",
        "SAM_POINT",
        "SAM_BEAM",
        "SAM_TRIA",
        "SAM_QUAD",
        "SAM_TETRA",
        "SAM_PYRAMID",
        "SAM_PENTA",
        "SAM_HEXA"
    };

    fprintf( ps, "element (%8d,%s)\n", id(), saiShapeNames[type()] );
    int nN = ( type() == SAM_QUAD ) ? 
        4:
        3;
    nodClass* aNode[4];
    nodes( aNode );
    fprintf( ps, "nodes...\n" );
    for( int in=0; in<nN; ++in )
    {
        if( printPointerValue )
            fprintf( ps, "%p, ", aNode[in] );
        fprintf( ps, "%8d\n", aNode[in]->id() );
    }

    int nE = ( type() == SAM_QUAD ) ? 
        4:
        3;
    fmuEdge* aEdge[4];
    edges( aEdge );
    fprintf( ps, "edges...\n" );
    for( int ie=0; ie<nE; ++ie )
    {
        if( printPointerValue )
            fprintf( ps, "%p, ", aEdge[ie] );
        fprintf( ps, "%8d\n", aEdge[ie]->id() );
    }
}
