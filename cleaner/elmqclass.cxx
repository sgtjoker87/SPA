/*

*/

#include "cleaner/elmqclass.hxx"
#include "cleaner/elmqarray.hxx"
#include "cleaner/nodclass.hxx"
#include "cleaner/fmuedge.hxx"
#include "cleaner/fmumaster.hxx"

#include "wksp/CSamEnvVar.hxx"

#include "MeshingDebugger/MeshingDebugger.h"
#include "MeshingDebugger/CMeshDebugOptions.hxx"

#include "FileImportExport/samFileNameUtility.h"
#include "FileImportExport/snprintf_macro.h"

#include <cstdio>
#include <cstddef>

static int kfCleanerDrawElementOn=1;
static int kfCleanerDrawElementOff=0;

// Do you want to draw the elements yes or no ?
static int kfCleanerDrawElementDfltStatus=kfCleanerDrawElementOff;



static int fDrawElement = kfCleanerDrawElementDfltStatus;

// SDRC FE quad Member Function Descriptor Block -----
//
// elmQClass::elmQClass
// class - elmQClass
// date 04-22-97  SDRC/P. Kinney
//
// Description:
//
// constructor
//
// Access:
//
// elmQClass()
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
// pointer to new entity
//
// Error conditions:
//
// insufficient memory
//
// Additional comments:
//
// This is the default constructor, it doesn't do you much good.
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
//     09/08/00 10:25 <          Mingwu Lai> Rel:ms9      SCCA:78384  Delta:9.3
//     
//     03/31/00 17:55 <               st qa> Rel:ms9      SCCA:71457  Delta:9.2
//     
//     11/06/97 11:55 <        Paul Kinney > Rel:ms6      SCCA:26918  Delta:6.3
//   
//     09/12/97 14:11 <        Paul Kinney > Rel:ms6      SCCA:22677  Delta:6.2





// SDRC FE quad Member Function Descriptor Block -----
//
// elmQClass::elmQClass
// class - elmQClass
// date 04-22-97  SDRC/P. Kinney
//
// Description:
//
// useful constructor
//
// Access:
//
// elmQClass( pCFirstNode, pCSecondNode, pCThirdNode, pCFourthNode, pCMaster )
//
// Input Parameters:
//
// pCFirstNode              nodClass*         the first node of the quad
// pCSecondNode             nodClass*         the second node of the quad
// pCThirdNode              nodClass*         the third node of the quad
// pCFourthNode             nodClass*         the fourth node of the quad
// pCMaster                 fmuMaster*        global list of entities,
//                                            defaults to NULL
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
// If no edge currently exists between a pair of nodes, it is created.
// The nodes input are assumed to be in order counterclockwise around
// the element normal.
//
// Quads are given labels that match their index position in the quad
// master global array. This is done for speed and is possible as quads
// are relabeled on output.
//
// Global Functions Referenced:
//
// new
//
// Global Data Referenced:
//
// none
//
// ------------------------------------------
///

elmQClass::elmQClass
(
    nodClass* pCFirstNode,
    nodClass* pCSecondNode,
    nodClass* pCThirdNode,
    nodClass* pCFourthNode,
    fmuMaster *pCMaster
) :
    elmBase()
{
    apCNodes[0] = pCFirstNode;
    apCNodes[1] = pCSecondNode;
    apCNodes[2] = pCThirdNode;
    apCNodes[3] = pCFourthNode;

    for( int ii=0, jj=0; ii < 4; ++ii )
    {
        jj = ii + 1;
        if( jj == 4 )
            jj = 0;

        apCEdges[ii] = apCNodes[ii]->shared_edge( apCNodes[jj] );
        if( apCEdges[ii] == NULL )
            apCEdges[ii] = new fmuEdge( apCNodes[ii], apCNodes[jj], this, pCMaster );
        else
            apCEdges[ii]->add_elem( this );
    }

    if( pCMaster )
    {
        const int id = pCMaster->global_quads()->length();
        pCMaster->global_quads()->append( this );
        elmBase::id( id );
    }
}


// SDRC FE quad Member Function Descriptor Block -----
//
// elmQClass::~elmQClass
// class - elmQClass
// date 04-22-97  SDRC/P. Kinney
//
// Description:
//
// destructor
//
// Access:
//
// delete elmqclass
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
// references to edges are deleted
//
// Global Functions Referenced:
//
// none
//
// Global Data Referenced:
//
// none
//
// ------------------------------------------
///

elmQClass::~elmQClass()
{
    apCEdges[0]->delete_elem( this );
    apCEdges[1]->delete_elem( this );
    apCEdges[2]->delete_elem( this );
    apCEdges[3]->delete_elem( this );
}


// SDRC FE Member Function Descriptor Block -----
//
// elmQClass::nodes
// class - elmQClass
// date 04-22-97  SDRC/P. Kinney
//
// Description:
//
// Retrieve an array of nodes for the quad. The nodes are in order
// counterclockwise around the element normal.
//
// Access:
//
// nodes( apCNodeList )
//
// Input Parameters:
//
// apCNodeList          nodclass*[4]          Array of pointers to nodes
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
// ------------------------------------------
///

void elmQClass::nodes( nodClass *apCNodeList[4] ) const
{
    apCNodeList[0]  = apCNodes[0];
    apCNodeList[1]  = apCNodes[1];
    apCNodeList[2]  = apCNodes[2];
    apCNodeList[3]  = apCNodes[3];
}


// SDRC FE quad Member Function Descriptor Block -----
//
// elmQClass::prev_node
// class - elmQClass
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
// mindexCurrentNode
//
// ------------------------------------------
///

nodClass* elmQClass::prev_node( nodClass* pCThisNode ) const
{
    if( apCNodes[ mindexCurrentNode ] != pCThisNode )
        if( !find_index( pCThisNode ) )
            return NULL;

    if( mindexCurrentNode == 0 )
        return apCNodes[3];

    return apCNodes[ mindexCurrentNode-1 ];
}


// SDRC FE quad Member Function Descriptor Block -----
//
// elmQClass::next_node
// class - elmQClass
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
// mindexCurrentNode
//
// ------------------------------------------
///

nodClass* elmQClass::next_node( nodClass* pCThisNode ) const
{
    if( apCNodes[ mindexCurrentNode ] != pCThisNode )
        if( !find_index( pCThisNode ) )
            return NULL;

    if( mindexCurrentNode == 3 )
        return apCNodes[0];

    return apCNodes[ mindexCurrentNode+1 ];
}


// SDRC FE quad Member Function Descriptor Block -----
//
// elmQClass::opposite_node
// class - elmQClass
// date 04-22-97  SDRC/P. Kinney
//
// Description:
//
// get the node that is diagonally across the quad from the one passed in.
//
// Access:
//
// opposite_node( pCThisNode )
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
// output_node                  nodClass*         the opposite node
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
// mindexCurrentNode 
//
// ------------------------------------------
///

nodClass* elmQClass::opposite_node( nodClass* pCThisNode ) const
{
    if( apCNodes[ mindexCurrentNode ] != pCThisNode )
        if( !find_index( pCThisNode ) )
            return NULL;

    const int iOppositeIndex = mindexCurrentNode + 2;
    if( iOppositeIndex > 3 )
        return apCNodes[ iOppositeIndex-4 ];

    return apCNodes[ iOppositeIndex ];
}


// SDRC FE quad Member Function Descriptor Block -----
//
// elmQClass::draw_me
// class - elmQClass
// date 04-22-97  SDRC/P. Kinney
//
// Description:
//
// display this quad on the screen
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
// ------------------------------------------
///

void elmQClass::draw_me( int iColor ) const
{
    /*apCEdges[0]->draw_me( iColor );
    apCEdges[1]->draw_me( iColor );
    apCEdges[2]->draw_me( iColor );
    apCEdges[3]->draw_me( iColor );*/
    int i=0, nNodes = 4, iElemId=0, fResponse=0;
    int aiNodeId[4] = {0};
    double dx = 0.0, dy = 0.0, adXY[8] = {0.0};
    char scPrompt[400];
    nodClass *apCNodes[4]={0};

    bool debugCleaner = CMeshDebugOptions::debugCleaner();

    if (!debugCleaner)
        return;

    iElemId = this->mId;

    this->nodes(apCNodes);

    for(i=0; i < 4; i++)
    {
        aiNodeId[i] = apCNodes[i]->id();
        dx = apCNodes[i]->node_u();
        dy = apCNodes[i]->node_v();
        adXY[2*i]   = apCNodes[i]->node_x();
        adXY[2*i+1] = apCNodes[i]->node_y();
    }


    if (fDrawElement == 0 || fDrawElement == 2)
    {
        sprintf(scPrompt, "Display cleaner element generation?");
        fResponse = MeshingDebugger::prompt(scPrompt,1);
            if (fResponse == 6) // Yes, draw them and prompt each time
                fDrawElement = 2;
            else if (fResponse == 7) // No, never draw them
            {
                    fDrawElement = -1;
                    goto Exit;
            }
            else // Continue, to draw without prompting
            {
                    fDrawElement = 1;
            }
    }
    else if (fDrawElement == -1)
    {
                goto Exit;
    }

    if (fDrawElement > 0)
    {
        // iEntType is set to Color Choice for cleaner
        // Usually -7 = Yellow -8 = Blue. Only 2 choices
        MeshingDebugger::display2dMesh( iColor,0,nNodes, aiNodeId, adXY, 
                                        1, &iElemId, &nNodes, aiNodeId);
    }

Exit:;
}


// SDRC FE quad Member Function Descriptor Block -----
//
// elmQClass::find_index
// class - elmQClass
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
// true if the pCTestNode is in the quad
// false if it is not
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
// mindexCurrentNode
//
// ------------------------------------------
///

bool elmQClass::find_index( nodClass* pCTestNode ) const
{
    mindexCurrentNode = -1;
    for( int ii = 0; ii < 4 && mindexCurrentNode == -1; ++ii )
        if( pCTestNode == apCNodes[ii] )
            mindexCurrentNode = ii;

    if( mindexCurrentNode > -1 )
        return true;

    mindexCurrentNode = 0;
    return false;
}

fmuEdge* elmQClass::opposite_edge(fmuEdge *pCEdge) const
{
    int iEdgeIndex = -1;

    for( int i=0; i<4; i++)
    {
        if(apCEdges[i] == pCEdge)
        {
            iEdgeIndex = i;
            break;
        }
    }

    if(iEdgeIndex == -1)
        return NULL;

    if(iEdgeIndex > 1)
        return apCEdges[iEdgeIndex-2];

    return apCEdges[iEdgeIndex+2];
}

void elmQClass::edges( fmuEdge *apCEdgeList[4] ) const
{
    apCEdgeList[0] = apCEdges[0];
    apCEdgeList[1] = apCEdges[1];
    apCEdgeList[2] = apCEdges[2];
    apCEdgeList[3] = apCEdges[3];
}

SAM_ElementShape elmQClass::type() const
{
    return SAM_QUAD;
}

// ----------------------------------------------------------------------------
namespace elemDynArray
{
void dump( FILE* ps, const elmQClassDynArray& dynQuads )
{
    size_t nQ = dynQuads.length();
    for( size_t ie=0; ie<nQ; ++ie )
    {
        elmBase* pQ = dynQuads[ie];
        if( pQ != NULL )
            pQ->dump( ps );
    }
}

void dump( int context, const elmQClassDynArray& dynQuads )
{
    static const size_t len = 128;
    static char buf[len];
    SNPRINTF( buf, len, "DEBUG_quadDynArray_%u.txt", context );
    samFileNameUtility sfn( buf, CSamEnvVar::getSamReportPath() );
    FILE* ps = fopen( sfn.PathName(), "wbc" );
    dump( ps, dynQuads );
}
} // end namespace elemDynArray
