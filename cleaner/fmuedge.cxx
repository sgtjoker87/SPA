/*

*/

#include "cleaner/fmuedge.hxx"
#include "cleaner/nodclass.hxx"
#include "cleaner/elmqclass.hxx"
#include "cleaner/fmumaster.hxx"
#include "wksp/CSamEnvVar.hxx"

#include "FileImportExport/samFileNameUtility.h"
#include "FileImportExport/snprintf_macro.h"

#include <cstdio>


fmuEdge::~fmuEdge()
{
    if( iNumberElems != 0 )
    {
        printf("Destructor called for referenced edge\n");
    }
    pCStartNode->delete_edge( this );
    pCEndNode->delete_edge( this );
}

// SDRC FE edge Member Function Descriptor Block -----
//
// fmuEdge::fmuEdge
// class - fmuEdge
// date 04-22-97  SDRC/P. Kinney
//
// Description:
//
// constructor with nodes supplied
//
// Access:
//
// fmuEdge( first_node, second_node ))
//
// Input Parameters:
//
// pCFirstNode                    nodClass*         first node of the edge
// pCSecondNode                   nodClass*         first node of the edge
//
// Output Parameters:
//
// none
//
// Return Code:
//
// address of new entity
//
// Error conditions:
//
// insufficient memory
//
// Additional comments:
//
// Edges are given labels that match their index position in the edge
// master global array. This is done for speed and is possible as edges
// are not output.
//
// Global Functions Referenced:
//
// new
//
// Global Data Referenced:
//
// none
//
// History:
//
// Created 4-22-97 as part of mesh valence analyzer and mesh cleaning
//
// ------------------------------------------
///
//   
//   HISTORY
//     
//     09/08/00 10:25 <          Mingwu Lai> Rel:ms9      SCCA:78384  Delta:9.2
//     
//     08/24/98 10:46 <               st qa> Rel:ms7      SCCA:36991  Delta:7.3
//     
//     03/26/98 16:40 <         Paul Kinney> Rel:ms6a     SCCA:32091  Delta:61.2
//     
//     12/02/97 16:32 <        Paul Kinney > Rel:ms6      SCCA:28618  Delta:6.4
//   
//     11/06/97 11:55 <        Paul Kinney > Rel:ms6      SCCA:26918  Delta:6.3
//   
//     09/12/97 14:11 <        Paul Kinney > Rel:ms6      SCCA:22677  Delta:6.2

fmuEdge::fmuEdge
(
    nodClass* pCFirstNode,
    nodClass* pCSecondNode,
    fmuMaster *pCMaster
) :
    pCStartNode ( pCFirstNode ),
    pCEndNode   ( pCSecondNode ),
    pCMidNode   ( NULL ),

    iNumberElems( 0 ),
    mId( 0 ),
    iMarked( 0 ),
    mIsFrozen( false )
{
    apElems[0] = NULL;
    apElems[1] = NULL;

    if(!pCStartNode->add_edge( this ) && pCMaster)
    {
        pCMaster->abort(true);
    }

    if(!pCEndNode->add_edge( this ) && pCMaster)
    {
        pCMaster->abort(true);
    }

    if( pCMaster )
    {
        const int id = pCMaster->global_edges()->length();
        pCMaster->global_edges()->append( this );
        this->id( id );
    }
}


// SDRC FE edge Member Function Descriptor Block -----
//
// fmuEdge::fmuEdge
// class - fmuEdge
// date 04-22-97  SDRC/P. Kinney
//
// Description:
//
// constructor with nodes and a elem supplied
//
// Access:
//
// fmuEdge( pCFirstNode, pCSecondNode, elem )
//
// Input Parameters:
//
// pCFirstNode                    nodClass*         first node of the edge
// pCSecondNode                   nodClass*         first node of the edge
// elem                           elmBase*          associated elem
//
// Output Parameters:
//
// none
//
// Return Code:
//
// address of new entity
//
// Error conditions:
//
// insufficient memory
//
// Additional comments:
//
// Edges are given labels that match their index position in the edge
// master global array. This is done for speed and is possible as edges
// are not output.
//
// Global Functions Referenced:
//
// new
//
// Global Data Referenced:
//
// none
//
// History:
//
// Created 4-22-97 as part of mesh valence analyzer and mesh cleaning
//
// ------------------------------------------
///

fmuEdge::fmuEdge
(
    nodClass* pCFirstNode,
    nodClass* pCSecondNode, 
    elmBase* elem,
    fmuMaster *pCMaster
) :
    pCStartNode ( pCFirstNode ),
    pCEndNode   ( pCSecondNode ),
    pCMidNode   ( NULL ),

    iNumberElems( 1 ),
    mId( 0 ),
    iMarked( 0 ),
    mIsFrozen( false )
{
    apElems[0] = elem;
    apElems[1] = NULL;

    if(!pCStartNode->add_edge( this ) && pCMaster)
    {
        pCMaster->abort(true);
    }

    if(!pCEndNode->add_edge( this ) && pCMaster)
    {
        pCMaster->abort(true);
    }

    if( pCMaster )
    {
        const int id = pCMaster->global_edges()->length();
        pCMaster->global_edges()->append( this );
        this->id( id );
    }
}


// SDRC FE edge Member Function Descriptor Block -----
//
// fmuEdge::other_node
// class - fmuEdge
// date 04-22-97  SDRC/P. Kinney
//
// Description:
//
// Get the node at the other end of this edge.
//
// Access:
//
// other_node( the_node )
//
// Input Parameters:
//
// the_node                    nodClass*         the node at one end
//
// Output Parameters:
//
// none
//
// Return Code:
//
// The node at the other end, NULL if "the_node" is not in the edge
//
// Error conditions:
//
// none
//
// Additional comments:
//
// none
//
// Global Functions Referenced:
//
// none
//
// Global Data Referenced:
//
// none
//
// History:
//
// Created 4-22-97 as part of mesh valence analyzer and mesh cleaning
//
// ------------------------------------------
///

nodClass* fmuEdge::other_node( nodClass *pCNode ) const
{
    if( pCStartNode == pCNode )
        return pCEndNode;
    if( pCEndNode == pCNode )
        return pCStartNode;
   return NULL;
}


// SDRC FE edge Member Function Descriptor Block -----
//
// fmuEdge::add_elem
// class - fmuEdge
// date 04-22-97  SDRC/P. Kinney
//
// Description:
//
// Add a referenced elem to an edge
//
// Access:
//
// add_elem( new_elem )
//
// Input Parameters:
//
// new_elem                    elmBase*          new elem to be added
//
// Output Parameters:
//
// none
//
// Return Code:
//
// none
//
// Error conditions:
//
// none
//
// Additional comments:
//
// If the elem is already listed, it is ignored.
//
// Global Functions Referenced:
//
// none
//
// Global Data Referenced:
//
// none
//
// History:
//
// Created 4-22-97 as part of mesh valence analyzer and mesh cleaning
//
// ------------------------------------------
///

void fmuEdge::add_elem( elmBase* pCNewElem )
{
    if(iNumberElems == 2)
        return;
    if( apElems[0] == pCNewElem )
        return;
    if( apElems[1] == pCNewElem )
        return;

    apElems[ iNumberElems ] = pCNewElem;
    ++iNumberElems;
}


// SDRC FE edge Member Function Descriptor Block -----
//
// fmuEdge::delete_elem
// class - fmuEdge
// date 04-22-97  SDRC/P. Kinney
//
// Description:
//
// remove the reference to a elem
//
// Access:
//
// delete_elem( pCOldElem )
//
// Input Parameters:
//
// pCOldElem                    elmBase*         elem to remove
//
// Output Parameters:
//
// none
//
// Return Code:
//
// none
//
// Error conditions:
//
// none
//
// Additional comments:
//
// If the elem is not already associated with the elem, it is ignored
//
// Global Functions Referenced:
//
// none
//
// Global Data Referenced:
//
// none
//
// History:
//
// Created 4-22-97 as part of mesh valence analyzer and mesh cleaning
//
// ------------------------------------------
///

void fmuEdge::delete_elem( elmBase* pCOldElem )
{
   if( apElems[0] == pCOldElem )
   {
      apElems[0] = apElems[1];
      apElems[1] = NULL;
      --iNumberElems;
   }
   if( apElems[1] == pCOldElem )
   {
      apElems[1] = NULL;
      --iNumberElems;
   }
}


// SDRC FE edge Member Function Descriptor Block -----
//
// fmuEdge::elems
// class - fmuEdge
// date 04-22-97  SDRC/P. Kinney
//
// Description:
//
// return pointers to the two referenced elems of the ege
//
// Access:
//
// elems( first_elem, second_elem )
//
// Input Parameters:
//
// none
//
// Output Parameters:
//
// first_elem                    elmBase**        the first elem
// second_elem                   elmBase**        the second elem
//
// Return Code:
//
// none
//
// Error conditions:
//
// none
//
// Additional comments:
//
// none
//
// Global Functions Referenced:
//
// none
//
// Global Data Referenced:
//
// none
//
// History:
//
// Created 4-22-97 as part of mesh valence analyzer and mesh cleaning
//
// ------------------------------------------
///

void fmuEdge::elems( elmBase** ppCFirstElem, elmBase** ppCSecondElem ) const
{
   *ppCFirstElem = apElems[0];
   *ppCSecondElem = apElems[1];
   return;
}


void fmuEdge::quads( elmQClass** ppCFirstQuad, elmQClass** ppCSecondQuad ) const
{
   *ppCFirstQuad = (elmQClass*) apElems[0];
   *ppCSecondQuad = (elmQClass*) apElems[1];
   return;
}


// SDRC FE edge Member Function Descriptor Block -----
//
// fmuEdge::other_elem
// class - fmuEdge
// date 04-22-97  SDRC/P. Kinney
//
// Description:
//
// return the elem on the other side of this edge
//
// Access:
//
// other_elem( known_elem )
//
// Input Parameters:
//
// known_elem                elmBase*        the known elem
//
// Output Parameters:
//
// none
//
// Return Code:
//
// other_elem               elmBase*         the elem across the edge
//
// Error conditions:
//
// none
//
// Additional comments:
//
// If the input elem is not adjacent to the edge, NULL is returned.
// If there is no elem on the other side of the edge, NULL is returned.
//
// Global Functions Referenced:
//
// none
//
// Global Data Referenced:
//
// None
//
// History:
//
// Created 4-22-97 as part of mesh valence analyzer and mesh cleaning
//
// ------------------------------------------
///

elmBase* fmuEdge::other_elem( elmBase* pCKnownElem ) const
{
    if( pCKnownElem == apElems[0] )
        return apElems[1];
    if( pCKnownElem == apElems[1] )
        return apElems[0];
    return NULL;
}


elmQClass* fmuEdge::other_quad( elmQClass* pCKnownQuad ) const
{
    elmBase *pCOtherQuad = NULL;

    if( pCKnownQuad == apElems[0] )
        pCOtherQuad = apElems[1];
    else if( pCKnownQuad == apElems[1] )
        pCOtherQuad = apElems[0];

    if( !pCOtherQuad )
        return NULL;

    if( pCOtherQuad->type() == SAM_QUAD )
        return (elmQClass*) pCOtherQuad;

    return NULL;
}

elmBase* fmuEdge::shared_elem( fmuEdge* pCEdge ) const
{
    if( apElems[0] == pCEdge->apElems[0] )
        return apElems[0];
    if( apElems[0] == pCEdge->apElems[1] )
        return apElems[0];
    if( apElems[1] == pCEdge->apElems[0] )
        return apElems[1];
    if( apElems[1] == pCEdge->apElems[1] )
        return apElems[1];
    return NULL;
}

elmQClass* fmuEdge::shared_quad( fmuEdge* pCEdge ) const
{
    elmBase *pCOtherElem = NULL;
    if( apElems[0] == pCEdge->apElems[0] )
        pCOtherElem = apElems[0];
    else if( apElems[0] == pCEdge->apElems[1] )
        pCOtherElem = apElems[0];
    else if( apElems[1] == pCEdge->apElems[0] )
        pCOtherElem = apElems[1];
    else if( apElems[1] == pCEdge->apElems[1] )
        pCOtherElem = apElems[1];

    if( !pCOtherElem )
        return NULL;

    if( pCOtherElem->type() == SAM_QUAD )
        return (elmQClass*) pCOtherElem;

    return NULL;
}


// SDRC FE edge Member Function Descriptor Block -----
//
// fmuEdge::draw_me
// class - fmuEdge
// date 04-22-97  SDRC/P. Kinney
//
// Description:
//
// draw this edge on the screen
//
// Access:
//
// draw_me()
//
// Input Parameters:
//
// none
//
// Output Parameters:
//
// none
//
// Return Code:
//
// none
//
// Error conditions:
//
// none
//
// Additional comments:
//
// none
//
// Global Functions Referenced:
//
// none
//
// Global Data Referenced:
//
// none
//
// History:
//
// Created 4-22-97 as part of mesh valence analyzer and mesh cleaning
//
// ------------------------------------------
///

void fmuEdge::draw_me( int iColor ) const
{
    pCStartNode->move_to_me();
    pCEndNode->draw_to_me( iColor );
}


void fmuEdge::highlight_me() const
{
   //Draw a circle at the middle of the edge
      //int iConst_1 = 1;
      //double dSX, dSY, dSZ, dEX, dEY, dEZ;
   //double dCoords[3];
      //dSX = pCStartNode->node_x();
      //dSY = pCStartNode->node_y();
      //dSZ = pCStartNode->node_z();
      //dEX = pCEndNode->node_x();
      //dEY = pCEndNode->node_y();
      //dEZ = pCEndNode->node_z();
      //dCoords[0] = ( dSX + dEX ) / 2.0;
      //dCoords[1] = ( dSY + dEY ) / 2.0;
      //dCoords[2] = ( dSZ + dEZ ) / 2.0;

   //rerota ( &iConst_1, &dCoords[0], &dCoords[1], &dCoords[2], TRANSF_1.atr,

   //gg3igmk( gg_kVP_MASK_ALL, NULL, NULL, dp_kCLR_WHITE, gg_eMRKR_CIRCLE, 1,
}

void fmuEdge::dump( FILE* ps ) const
{
    static const bool printPointerValue = false;

    if( ps == NULL )
        return;

    fprintf( ps, "edge (%8d,%2d)\n", mId, iNumberElems );
    fprintf( ps, "elements...\n" );
    for( size_t ie=0; ie<iNumberElems; ++ie )
    {
        if( printPointerValue )
            fprintf( ps, "%p, ", apElems[ie] );
        fprintf( ps, "%8d\n", apElems[ie]->id() );
    }
    fprintf( ps, "nodes: start,end,mid...\n" );
    if( printPointerValue )
        fprintf( ps, "%p, ", pCStartNode );
    fprintf( ps, "%8d\n", pCStartNode->id() );
    if( printPointerValue )
        fprintf( ps, "%p, ", pCEndNode );
    fprintf( ps, "%8d\n", pCEndNode->id() );
    if( printPointerValue )
        fprintf( ps, "%p, ", pCMidNode );
    fprintf( ps, "%8d\n", (pCMidNode?pCMidNode->id():-1) );
}

// ----------------------------------------------------------------------------
namespace edgeDynArray
{
void dump( FILE* ps, const fmuEdgeDynArray& dynEdges )
{
    size_t nE = dynEdges.length();
    for( size_t in=0; in<nE; ++in )
    {
        fmuEdge* pE = dynEdges[in];
        pE->dump( ps );
    }
}

void dump( int context, const fmuEdgeDynArray& dynEdges )
{
    static const size_t len = 128;
    static char buf[len];
    SNPRINTF( buf, len, "DEBUG_edgeDynArray_%u.txt", context );
    samFileNameUtility sfn( buf, CSamEnvVar::getSamReportPath() );
    FILE* ps = fopen( sfn.PathName(), "wbc" );
    dump( ps, dynEdges );
}
} // end namespace edgeDynArray
