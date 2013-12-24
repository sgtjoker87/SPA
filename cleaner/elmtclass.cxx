/*

*/

#include "cleaner/elmtclass.hxx"
#include "cleaner/elmtarray.hxx"
#include "cleaner/nodclass.hxx"
#include "cleaner/fmuedge.hxx"
#include "cleaner/fmumaster.hxx"
#include "cleaner/fmuvector.hxx"

#include "wksp/CSamEnvVar.hxx"

#include "FileImportExport/samFileNameUtility.h"
#include "FileImportExport/snprintf_macro.h"

#include <cstdio>
#include <cstddef>

static const double sam_CosineMinChopAngle = 0.707;

// SDRC FE tria Member Function Descriptor Block -----
//
// elmTClass::elmTClass
// class - elmTClass
// date 04-22-97  SDRC/P. Kinney
//                                            defaults to NULL
// Additional comments:
//
// If no edge currently exists between a pair of nodes, it is created.
// Even though the tria uses only 3 nodes & edges, the base element
// class declares 4, so initialize 4. The nodes input are assumed to be
// in order counterclockwise around the element normal.
//
// Trias are given labels that match their index position in the tria
// master global array. This is done for speed and is possible as trias
// are relabeled on output.
// ------------------------------------------

elmTClass::~elmTClass()
{
    apCEdges[0]->delete_elem( this );
    apCEdges[1]->delete_elem( this );
    apCEdges[2]->delete_elem( this );
}

// elmTClass( pCFirstNode, pCSecondNode, pCThirdNode, pCMaster )
//
// Input Parameters:
//
// pCFirstNode              nodClass*         the first node of the tria
// pCSecondNode             nodClass*         the second node of the tria
// pCThirdNode              nodClass*         the third node of the tria
// pCMaster                 fmuMaster*        global list of entities,

elmTClass::elmTClass
(
    nodClass* pCFirstNode,
    nodClass* pCSecondNode,
    nodClass* pCThirdNode,
    fmuMaster *pCMaster
) :
    elmBase()
{
    apCNodes[0] = pCFirstNode;
    apCNodes[1] = pCSecondNode;
    apCNodes[2] = pCThirdNode;

    for( int ii=0, jj=0; ii < 3; ++ii )
    {
        jj = ii + 1;
        if( jj == 3 )
            jj = 0;

        apCEdges[ii] = apCNodes[ii]->shared_edge( apCNodes[jj] );
        if( apCEdges[ii] == NULL )
            apCEdges[ii] = new fmuEdge( apCNodes[ii], apCNodes[jj], this, pCMaster );
        else
            apCEdges[ii]->add_elem( this );
    }

    if( pCMaster )
    {
        const int id = pCMaster->global_trias()->length();
        pCMaster->global_trias()->append( this );
        elmBase::id( id );
    }

    check_freeze();
}


// SDRC FE tria Member Function Descriptor Block -----
//
// elmTClass::check_freeze
// class - elmTClass
// date 04-22-97  SDRC/P. Kinney
//
// Description:
//
// Check whether this tria should be frozen (unchanged) due to its position
// in the geometry (two sides to boundaries) and the angle tolerance that
// was specified by the user
//
// Access:
//
// check_freeze()
//
// ------------------------------------------

void elmTClass::check_freeze()
{
    fmuVector CVec01( apCNodes[0], apCNodes[1] );
    CVec01.unitize();
    fmuVector CVec12( apCNodes[1], apCNodes[2] );
    CVec12.unitize();
    fmuVector CVec20( apCNodes[2], apCNodes[0] );
    CVec20.unitize();

    double dCos;
    if( apCEdges[0]->frozen() && apCEdges[1]->frozen() )
    {
        dCos = CVec12 % -CVec01;
        if( dCos > sam_CosineMinChopAngle )
        {
            mIsFrozen = true;
            return;
        }
    }
    if( apCEdges[1]->frozen() && apCEdges[2]->frozen() )
    {
        dCos = CVec20 % -CVec12;
        if( dCos > sam_CosineMinChopAngle )
        {
            mIsFrozen = true;
            return;
        }
    }
    if( apCEdges[2]->frozen() && apCEdges[0]->frozen() )
    {
        dCos = CVec01 % -CVec20;
        if( dCos > sam_CosineMinChopAngle )
        {
            mIsFrozen = true;
            return;
        }
    }
}


// SDRC FE tria Member Function Descriptor Block -----
//
// elmTClass::nodes
// class - elmTClass
// date 04-22-97  SDRC/P. Kinney
//
// Description:
//
// Retrieve an array of nodes for the tria. The nodes are in order
// counterclockwise around the element normal.
//
// Access:
//
// nodes( apCNodeList )
//
// Input Parameters:
//
// apCNodeList          nodclass*[3]          Array of pointers to nodes
//
// Output Parameters:
//
// none
//
// Return Code:
//
// none
// ------------------------------------------

void  elmTClass::nodes( nodClass *apCNodeList[3] ) const
{
    apCNodeList[0]  = apCNodes[0];
    apCNodeList[1]  = apCNodes[1];
    apCNodeList[2]  = apCNodes[2];
}


// SDRC FE tria Member Function Descriptor Block -----
//
// elmTClass::prev_node
// class - elmTClass
// date 04-22-97  SDRC/P. Kinney
//
// Description:
//
// get the node previous in order counterclockwise around the normal to the
// node passed in.
//
// Access:
//
// prev_node( pCThisNode )
//
// Input Parameters:
//
// pCThisNode                    nodClass*         the reference node
//
// Output Parameters:
//
// none
//
// Return Code:
//
// output_node                  nodClass*         the previous node
//
// ------------------------------------------

nodClass* elmTClass::prev_node( nodClass* pCThisNode ) const
{
    if( apCNodes[ mindexCurrentNode ] != pCThisNode )
        if( !find_index( pCThisNode ) )
            return NULL;

    if( mindexCurrentNode == 0 )
        return apCNodes[2];

    return apCNodes[ mindexCurrentNode-1 ];
}


// SDRC FE tria Member Function Descriptor Block -----
//
// elmTClass::next_node
// class - elmTClass
// date 04-22-97  SDRC/P. Kinney
//
// Description:
//
// get the node next in order counterclockwise around the normal to the
// node passed in.
//
// Access:
//
// next_node( pCThisNode )
//
// Input Parameters:
//
// pCThisNode                    nodClass*         the reference node
//
// Output Parameters:
//
// none
//
// Return Code:
//
// output_node                  nodClass*         the next node
// ------------------------------------------

nodClass* elmTClass::next_node( nodClass* pCThisNode ) const
{
    if( apCNodes[ mindexCurrentNode ] != pCThisNode )
        if( !find_index( pCThisNode ) )
            return NULL;

    if( mindexCurrentNode == 2 )
        return apCNodes[0];

    return apCNodes[ mindexCurrentNode+1 ];
}


// SDRC FE tria Member Function Descriptor Block -----
//
// elmTClass::draw_me
// class - elmTClass
// date 04-22-97  SDRC/P. Kinney
//
// Description:
//
// display this tria on the screen
//
// Access:
//
// draw_me( iColor )
//
// Input Parameters:
//
// iColor          int         The color to be used for the tria. These
//                                 colors match dp_kCLR_xxx as defined in
//                                 dpiinc/dpityp.h
//
// ------------------------------------------

void elmTClass::draw_me( int iColor ) const
{
    apCEdges[0]->draw_me( iColor );
    apCEdges[1]->draw_me( iColor );
    apCEdges[2]->draw_me( iColor );
}


// SDRC FE tria Member Function Descriptor Block -----
//
// elmTClass::find_index
// class - elmTClass
// date 04-22-97  SDRC/P. Kinney
//
// Description:
//
// set the internal node index for use in the prev, next, and opposite node
// retrieval routines
//
// Access:
//
// find_index( pCTestNode )
//
// Input Parameters:
//
// pCTestNode                    nodClass*         the new reference node
//
// Output Parameters:
//
// none
//
// Return Code:
//
// true if the pCTestNode is in the tria
// false if it is not
// ------------------------------------------

bool elmTClass::find_index( nodClass* pCTestNode ) const
{
    mindexCurrentNode = -1;
    for( int ii=0; ii < 3 && mindexCurrentNode == -1; ++ii )
        if( pCTestNode == apCNodes[ii] )
            mindexCurrentNode = ii;

    if( mindexCurrentNode > -1 )
        return true;

    mindexCurrentNode = 0;
    return false;
}

void elmTClass::edges( fmuEdge *apCEdgeList[4] ) const
{
    apCEdgeList[0] = apCEdges[0];
    apCEdgeList[1] = apCEdges[1];
    apCEdgeList[2] = apCEdges[2];
}

SAM_ElementShape elmTClass::type() const
{
    return SAM_TRIA;
}

// ----------------------------------------------------------------------------
namespace elemDynArray
{
void dump( FILE* ps, const elmTClassDynArray& dynTria )
{
    size_t nT = dynTria.length();
    for( size_t ie=0; ie<nT; ++ie )
    {
        elmBase* pT = dynTria[ie];
        if( pT != NULL )
            pT->dump( ps );
    }
}

void dump( int context, const elmTClassDynArray& dynTria )
{
    static const size_t len = 128;
    static char buf[len];
    SNPRINTF( buf, len, "DEBUG_triaDynArray_%u.txt", context );
    samFileNameUtility sfn( buf, CSamEnvVar::getSamReportPath() );
    FILE* ps = fopen( sfn.PathName(), "wbc" );
    dump( ps, dynTria );
}
} // end namespace elemDynArray
