/*

*/

#include "cleaner/nodclass.hxx"
#include "cleaner/elmqclass.hxx"
#include "cleaner/fmuedge.hxx"
#include "cleaner/fmuvector.hxx"
#include "cleaner/fmumaster.hxx"

#include "wksp/CSamEnvVar.hxx"

#include "FileImportExport/samFileNameUtility.h"
#include "FileImportExport/snprintf_macro.h"

#include <cstdlib>
#include <cstdio>
#include <cmath>

//=======================================================================================
// SDRC FE node Member Function Descriptor Block -----
//
// nodClass::nodClass default
// class - nodClass
// date 04-22-97  SDRC/P. Kinney
//
// Description:
//
// default constructor
//
// Access:
//
// nodClass()
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
// pointer to new class
//
// Error conditions:
//
// insufficient memory
//
// Additional comments:
//
// If used, the node created is at the origin
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
//     02/23/01 14:01 <            Hui Xiao> Rel:ms9      SCCA:88085  Delta:9.3
//     
//     09/08/00 10:25 <          Mingwu Lai> Rel:ms9      SCCA:78384  Delta:9.2
//     
//     09/01/98 13:20 <               st qa> Rel:ms7      SCCA:37372  Delta:7.3
//     
//     03/26/98 16:40 <         Paul Kinney> Rel:ms6a     SCCA:32091  Delta:61.2
//     
//     12/02/97 16:32 <        Paul Kinney > Rel:ms6      SCCA:28618  Delta:6.5
//     
//     10/27/97 17:19 <         Ravi Boppe > Rel:ms6      SCCA:25909  Delta:6.4
//     
//     10/07/97 17:41 <         Ravi Boppe > Rel:ms6      SCCA:24455  Delta:6.3
//   
//     09/12/97 14:11 <        Paul Kinney > Rel:ms6      SCCA:22677  Delta:6.2

nodClass::nodClass() :
    edgeList(),
    pCProjectNode( NULL ),
    pCSourceNode( NULL ),
    dGradient( 0.0 ),
    dNodeX( 0.0 ),
    dNodeY( 0.0 ),
    dNodeZ( 0.0 ),
    dNodeU( -1.0 ),
    dNodeV( -1.0 ),
    mId( 0 ),
    miType( 0 ),
    lPromoted( false ),
    lDegenerateLoopPoint( false ),
    lCleanList( false ),
    lSmoothList( false ),
    lSmoothHP( false ),
    lHardpoint( false ),
    lBoundary( false )
{
}

// SDRC FE node Member Function Descriptor Block -----
//
// nodClass::nodClass standard
// class - nodClass
// date 04-22-97  SDRC/P. Kinney
//
// Description:
//
// constructor that sets initial coordinates
//
// Access:
//
// nodClass( newX,
//           newY,
//           newZ )
//
// Input Parameters:
//
// newX          double          initial value for X coordinate
// newY          double          initial value for Y coordinate
// newZ          double          initial value for Z coordinate
//
// Output Parameters:
//
// none
//
// Return Code:
//
// pointer to new class
//
// Error conditions:
//
// insufficient memory
//
// Additional comments:
//
// Only the coordinates are set
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

nodClass::nodClass( double dNewX,
                    double dNewY,
                    double dNewZ,
                    fmuMaster *pCMaster,
                    double dNewU,
                    double dNewV) :
    edgeList(),
    pCProjectNode( NULL ),
    pCSourceNode( NULL ),
    dGradient( 0.0 ),
    dNodeX( dNewX ),
    dNodeY( dNewY ),
    dNodeZ( dNewZ ),
    dNodeU( dNewU ),
    dNodeV( dNewV ),
    mId( 0 ),
    miType( 0 ),
    lPromoted( false ),
    lDegenerateLoopPoint( false ),
    lCleanList( false ),
    lSmoothList( false ),
    lSmoothHP( false ),
    lHardpoint( false ),
    lBoundary( false ),
    lLockNode (false)
{
    if( pCMaster )
    {
        pCMaster->global_interior()->append( this );
        const int id = pCMaster->highest_node_id() + 1;
        this->id( id );
        pCMaster->highest_node_id( id );
    }
}

// SDRC FE node Member Function Descriptor Block -----
//
// nodClass::nodClass constructor
// class - nodClass
// date 04-22-97  SDRC/P. Kinney
//
// Description:
//
// constructor from a vector
//
// Access:
//
// nodClass( vector )
//
// Input Parameters:
//
// vector                    fmuVector          vector to use as coordinates
//
// Output Parameters:
//
// none
//
// Return Code:
//
// pointer to new class
//
// Error conditions:
//
// insufficient memory
//
// Additional comments:
//
// The vector is used as the coordinates of the node.
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

nodClass::nodClass( fmuVector& CVector,
                    fmuMaster *pCMaster,
                    double dNewU,
                    double dNewV) :
    edgeList(),
    pCProjectNode( NULL ),
    pCSourceNode( NULL ),
    dGradient( 0.0 ),
    dNodeX( CVector.x() ),
    dNodeY( CVector.y() ),
    dNodeZ( CVector.z() ),
    dNodeU( dNewU ),
    dNodeV( dNewV ),
    mId( 0 ),
    miType( 0 ),
    lPromoted( false ),
    lDegenerateLoopPoint( false ),
    lCleanList( false ),
    lSmoothList( false ),
    lSmoothHP( false ),
    lHardpoint( false ),
    lBoundary( false ),
    lLockNode (false)
{
    if( pCMaster )
    {
        pCMaster->global_interior()->append( this );
        const int id = pCMaster->highest_node_id() + 1;
        this->id( id );
        pCMaster->highest_node_id( id );
    }
}

// SDRC FE node Member Function Descriptor Block -----
//
// nodClass::~nodClass
// class - nodClass
// date 04-22-97  SDRC/P. Kinney
//
// Description:
//
// destructor
//
// Access:
//
// delete nodeclass
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
// node still references edges. Message is printed
//
// Additional comments:
//
// none
//
// Global Functions Referenced:
//
// delete
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

nodClass::~nodClass()
{
    if( edgeList.length() > 0 )
    {
        printf("Destructor called for referenced node\n");
    }
    pCProjectNode = NULL;
    pCSourceNode = NULL;
}

// SDRC FE node Member Function Descriptor Block -----
//
// nodClass::set
// class - nodClass
// date 04-22-97  SDRC/P. Kinney
//
// Description:
//
// Change the coordinates of a node from a vector
//
// Access:
//
// set( vector )
//
// Input Parameters:
//
// vector                     fmuVector          new coordinates
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


void nodClass::set( fmuVector& CVector )
{
    dNodeX = CVector.x();
    dNodeY = CVector.y();
    dNodeZ = CVector.z();
}

// SDRC FE node Member Function Descriptor Block -----
//
// nodClass::add_edge
// class - nodClass
// date 04-22-97  SDRC/P. Kinney
//
// Description:
//
// Add an edge to this node
//
// Access:
//
// add_edge( new_edge )
//
// Input Parameters:
//
// new_edge                    fmuEdge*          new edge to add
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
// If this edge is already in the list of edges, it is ignored.
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

bool nodClass::add_edge( fmuEdge* pCNewEdge )
{
    int iNumberEdges;

    iNumberEdges = edgeList.length();   
    for( int ii = 0; ii < iNumberEdges; ++ii)
    {
        if( edgeList[ii] == pCNewEdge)
            return true;
    }

    edgeList.append(pCNewEdge);

    return true;
}

// SDRC FE node Member Function Descriptor Block -----
//
// nodClass::delete_edge
// class - nodClass
// date 04-22-97  SDRC/P. Kinney
//
// Description:
//
// Removed a referenced edge
//
// Access:
//
// delete_edge( old_edge )
//
// Input Parameters:
//
// old_edge                    fmuEdge*          edge to remove
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
// If edge is not found, it is ignored
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

void nodClass::delete_edge( fmuEdge* pCOldEdge )
{
    const int iNumberEdges = edgeList.length();

    int ii;
    for( ii = 0; ii < iNumberEdges && edgeList[ii] != pCOldEdge; ++ii);

    if( edgeList[ii] == pCOldEdge )
    {
        edgeList.remove(ii);
    }
}

int nodClass::number_edges() const
{
    return edgeList.length(); 
}

fmuEdge* nodClass::first_edge() const
{
    if( edgeList.length() == 0 )
        return NULL;
    return edgeList[0];
}

void nodClass::edge_list( fmuEdgeDynArray& CEdges ) const
{
    const int iNumberEdges = edgeList.length();
    for( int ii = 0; ii < iNumberEdges; ++ii )
        CEdges[ii] = edgeList[ii];
}


// SDRC FE node Member Function Descriptor Block -----
//
// nodClass::shared_edge
// class - nodClass
// date 04-22-97  SDRC/P. Kinney
//
// Description:
//
// Find the edge that exists between two nodes
//
// Access:
//
// shared_edge( other )
//
// Input Parameters:
//
// other                    nodClass*          the other node
//
// Output Parameters:
//
// none
//
// Return Code:
//
// The edge between the two nodes, NULL if there is none
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

fmuEdge* nodClass::shared_edge( nodClass* pCNode ) const
{
    fmuEdge *pCEdge;
    const int iNumberEdges = edgeList.length();   

    for( int ii = 0; ii < iNumberEdges; ++ii )
    {
        pCEdge = edgeList[ii];
        if( pCEdge->start_node() == pCNode )
            return pCEdge;
        if( pCEdge->end_node() == pCNode )
            return pCEdge;
    }

    return NULL;
}

// SDRC FE node Member Function Descriptor Block -----
//
// nodClass::draw_me
// class - nodClass
// date 04-22-97  SDRC/P. Kinney
//
// Description:
//
// Display a marker on the screen at the node coordinates
//
// Access:
//
// draw_me
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

void nodClass::draw_me( int iColor ) const
{
   //int iConst_1 = 1;
   //double dCoords[3];
   //dCoords[0] = dNodeX;
   //dCoords[1] = dNodeY;
   //dCoords[2] = dNodeZ;
   //rerota ( &iConst_1, &dCoords[0], &dCoords[1], &dCoords[2], TRANSF_1.atr,
   //                                                           TRANSF_1.ash);

   //gg3igmk( gg_kVP_MASK_ALL, NULL, NULL, iColor, gg_eMRKR_CIRCLE, 1,
   //                                                dCoords, NULL );
}

// SDRC FE node Member Function Descriptor Block -----
//
// nodClass::smooth_me
// class - nodClass
// date 04-22-97  SDRC/P. Kinney
//
// Description:
//
// Apply laplacian smoothing to change the position of the node
//
// Access:
//
// smooth_me( pCNormal )
//
// Input Parameters:
//
// pCNormal          fmuVector*          normal vector of infinite plane
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
// If the node is connected to hardpoints by 3 or more edges, don't
// bother smoothing.
// If the node is connected to hardpoints by two edges, has four or
// less total edges, and the distance between it and one of the
// hardpoints is considerably less than it and one of the other points,
// use only those hardpoints to compute the new position and then mark this
// node as a hardpoint after smoothing as its position will not improve.
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

bool nodClass::smooth_me( fmuVector *pCNormal )
{
    int ii;
    int nHardpt = 0;
    double dX = 0.0, dY = 0.0, dZ = 0.0;
    double dXhp = 0.0, dYhp = 0.0, dZhp = 0.0;
    double dDivisor, dHardLen, dSoftLen;
    const double kSMOOTH_TOL = 0.0002;
    bool lChanged = false;
    bool lMakeHard = false;
    nodClass *pCNode;

    fmuEdge *pCSoftEdge = NULL;

    int iNumberEdges = edgeList.length();   

    if( lHardpoint )
        return false;

    if (lLockNode)
        return false;

    // JC - 12/03/2013 - use a vector instead of a hardcoded size of 10
    VEC_fmuEdge pCHardEdge(iNumberEdges);

    nHardpt = 0;

    dDivisor = (double) iNumberEdges;
    for( ii = 0; ii < iNumberEdges; ++ii)
    {
        pCNode = edgeList[ii]->other_node( this );
        if( pCNode->boundary() )
        {
            dXhp += pCNode->dNodeX;
            dYhp += pCNode->dNodeY;
            dZhp += pCNode->dNodeZ;
            pCHardEdge[nHardpt] = edgeList[ii];
            ++nHardpt;
        }
        else
        {
            pCSoftEdge = edgeList[ii];
        }
        dX += pCNode->dNodeX;
        dY += pCNode->dNodeY;
        dZ += pCNode->dNodeZ;
    }

    if( nHardpt == 2 &&
        iNumberEdges <= 4 &&
        pCSoftEdge &&
        !pCHardEdge[0]->shared_elem( pCHardEdge[1] ) )
    {
        fmuVector CHardVec( pCHardEdge[0]->start_node(), pCHardEdge[0]->end_node() );
        fmuVector CSoftVec( pCSoftEdge->start_node(), pCSoftEdge->end_node() );
        dHardLen = CHardVec.length();
        dSoftLen = CSoftVec.length();
        if( dHardLen * 3.0 < dSoftLen )
            lMakeHard = true;
    }

    //We do *not* want the hardpoint edges to share an element when
    //doing smoothing of nodes down the center of a flange.
    if( nHardpt == 2 &&
        iNumberEdges >= 4 &&
        !lMakeHard &&
        !pCHardEdge[0]->shared_elem( pCHardEdge[1] ) )
    {
        //Make the hardpoint edges count double by averaging in the separate
        //hardpoint average
        dX /= dDivisor;
        dY /= dDivisor;
        dZ /= dDivisor;
        dDivisor = (double) nHardpt;
        dXhp /= dDivisor;
        dYhp /= dDivisor;
        dZhp /= dDivisor;
        dX = ( dX + dXhp ) / 2.0;
        dY = ( dY + dYhp ) / 2.0;
        dZ = ( dZ + dZhp ) / 2.0;
        if( fabs( dNodeX - dX ) > kSMOOTH_TOL )
            lChanged = true;
        if( fabs( dNodeY - dY ) > kSMOOTH_TOL )
            lChanged = true;
        if( fabs( dNodeZ - dZ ) > kSMOOTH_TOL )
            lChanged = true;

        dNodeX = dX;
        dNodeY = dY;
        dNodeZ = dZ;
    }

    else if( lMakeHard )
    {
        dDivisor = (double) nHardpt;
        dXhp /= dDivisor;
        dYhp /= dDivisor;
        dZhp /= dDivisor;

        if( fabs( dNodeX - dXhp ) > kSMOOTH_TOL )
            lChanged = true;
        if( fabs( dNodeY - dYhp ) > kSMOOTH_TOL )
            lChanged = true;
        if( fabs( dNodeZ - dZhp ) > kSMOOTH_TOL )
            lChanged = true;

        dNodeX = dXhp;
        dNodeY = dYhp;
        dNodeZ = dZhp;

        //Check lengths again before actually freezing the node as smoothing
        //may have moved it enough to not be necessary.
        fmuVector CHardVec( pCHardEdge[0]->start_node(), pCHardEdge[0]->end_node() );
        fmuVector CSoftVec( pCSoftEdge->start_node(), pCSoftEdge->end_node() );
        dHardLen = CHardVec.length();
        dSoftLen = CSoftVec.length();
        if( dHardLen * 3.0 < dSoftLen )
        {
            lHardpoint = true;
            lSmoothHP = true;
        }
    }
    else
    {
        fmuVector CAnglePoint;

        dX /= dDivisor;
        dY /= dDivisor;
        dZ /= dDivisor;

        if( nHardpt == 1 )
        {
            if( equal_angles( pCHardEdge[0], pCNormal, &CAnglePoint ) )
            {
                dX = (dX + CAnglePoint.x() ) / 2.0;
                dY = (dY + CAnglePoint.y() ) / 2.0;
                dZ = (dZ + CAnglePoint.z() ) / 2.0;
            }
        }

        if( fabs( dNodeX - dX ) > kSMOOTH_TOL )
            lChanged = true;
        if( fabs( dNodeY - dY ) > kSMOOTH_TOL )
            lChanged = true;
        if( fabs( dNodeZ - dZ ) > kSMOOTH_TOL )
            lChanged = true;

        dNodeX = dX;
        dNodeY = dY;
        dNodeZ = dZ;
    }

    return lChanged;
}


// SDRC FE node Member Function Descriptor Block -----
//
// nodClass::equal_angles
// class - nodClass
// date 04-22-97  SDRC/P. Kinney
//
// Description:
//
// return a smoothing correction based on equal angles.
//
// Access:
//
// equal_angles( *pCEdge, *pCNormal, *pCAnglePoint )
//
// Input Parameters:
//
// pCEdge               fmuEdge*               Edge between this and a
//                                             boundary node
// pCNormal             fmuVector*             Normal of infinite plane
//
// Output Parameters:
//
// pCAnglePoint         fmuVector*             Recommended position of node
//
// Return Code:
//
// bool flag is returned to indicate if the correction can be used.
//
// Error conditions:
//
// none
//
// Additional comments:
//
// This routine should only be called if there is only one edge connecting
// "this" node to the mesh boundary. The goal of the routine is to compute
// a recommended position of "this" node as a correction to the standard
// smoothed value. This position is chosen to create equal angles on either
// side of the edge with the mesh boundary:
//
//     this       angle R-N-this = angle this-N-L
//      |
//      |         length N-this = average( length NR & length NL)
//  L---N---R
//
//  This correction is only flagged as useful if:
//  The elements in these angles are both quads.
//  Both L and R are on the mesh boundary
//  The angle R-N-L is greater than ... (otherwise the correction isn't needed)
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
// Created 11-07-97 as part of mesh valence analyzer and mesh cleaning
//
// ------------------------------------------
///

bool nodClass::equal_angles( fmuEdge *pCEdge, fmuVector *pCNormal,
                                                 fmuVector *pCAnglePoint )
{
    const double dANGLE_TOL = 3.4906585;  //200 degrees

    nodClass *pCNode = pCEdge->start_node();
    if( pCNode == this )
        pCNode = pCEdge->end_node();

    elmBase *pCFirstElem=0, *pCSecondElem=0;
    elmQClass *pCQuad1=0, *pCQuad2=0, *pCQuadL=0, *pCQuadR=0;

    pCEdge->elems( &pCFirstElem, &pCSecondElem );

    if (!pCFirstElem || !pCSecondElem)
        return false;

    if( pCFirstElem->type() == SAM_TRIA )
        return false;
    if( pCSecondElem->type() == SAM_TRIA )
        return false;

    pCQuad1 = (elmQClass*) pCFirstElem;
    pCQuad2 = (elmQClass*) pCSecondElem;

    if (!pCQuad1 || !pCQuad2)
        return false;

    if( pCNode == pCQuad1->next_node( this ) )
    {
        pCQuadR = pCQuad1;
        pCQuadL = pCQuad2;
    }
    else
    {
        pCQuadR = pCQuad2;
        pCQuadL = pCQuad1;
    }

    nodClass *pCNodeR = pCQuadR->opposite_node( this );
    if( !pCNodeR->boundary() )
        return false;
    nodClass *pCNodeL = pCQuadL->opposite_node( this );
    if( !pCNodeL->boundary() )
        return false;

    fmuVector CVecR( pCNode, pCNodeR );
    fmuVector CVecL( pCNode, pCNodeL );

    double dAngle = pCNormal->angle( CVecR, CVecL );
    if( dAngle < dANGLE_TOL )
        return false;

    dAngle /= 2.0;

    double dLength = (CVecR.length() + CVecL.length() ) / 2.0;

    fmuVector CVecNew = pCNormal->rotate( CVecR, dAngle );
    CVecNew.unitize();
    CVecNew *= dLength;
    fmuVector CVecNode( pCNode );

    *pCAnglePoint = CVecNode + CVecNew;

    return true;
}


// SDRC FE node Member Function Descriptor Block -----
//
// nodClass::size
// class - nodClass
// date 04-22-97  SDRC/P. Kinney
//
// Description:
//
// return the recommended edge length for edges near the node
//
// Access:
//
// size()
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
// size_value               double          the recommended size
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

double nodClass::size()
{
    return 1.0;
}


void nodClass::move_to_me() const
{
   //int iConst_1 = 1;
   //double dX, dY, dZ;
   //float rCoords[3];
   //dX = dNodeX;
   //dY = dNodeY;
   //dZ = dNodeZ;
   //rerota ( &iConst_1, &dX, &dY, &dZ, TRANSF_1.atr, TRANSF_1.ash);

   //rCoords[0] = dX;
   //rCoords[1] = dY;
   //rCoords[2] = dZ;

   //CALL_FORTRAN (spgmov) (  rCoords );
}


void nodClass::draw_to_me( int color ) const
{
   //int iConst_1 = 1;
   //double dX, dY, dZ;
   //float rCoords[3];
   //dX = dNodeX;
   //dY = dNodeY;
   //dZ = dNodeZ;
   //rerota ( &iConst_1, &dX, &dY, &dZ, TRANSF_1.atr, TRANSF_1.ash);

   //rCoords[0] = dX;
   //rCoords[1] = dY;
   //rCoords[2] = dZ;

   //CALL_FORTRAN (spgpln) ( &color, &iConst_1, &iConst_1, rCoords );
}

void nodClass::surround_nodes(nodClassDynArray &nodeList)
{
    nodClass *pcNode;
    int iNumberEdges = edgeList.length();   

    if( iNumberEdges > 0 )
    {
        for( int i=0; i<iNumberEdges; ++i )
        {
            if(edgeList[i] != NULL)
            {
                pcNode = edgeList[i]->other_node(this);
                nodeList.append(pcNode);
            }
        }
    }
}

void nodClass::dump( FILE* ps ) const
{
    static const bool printPointerValue = false;

    if( ps == NULL )
        return;

    size_t nE = edgeList.length();
    fprintf( ps, "node (%8d,%2d,%2d)\n", mId, miType, nE );
    fprintf( ps, "edges...\n" );
    for( size_t ie=0; ie<nE; ++ie )
    {
        fmuEdge* pE = edgeList[ie];
        if( printPointerValue )
            fprintf( ps, "%p, ", pE );
        fprintf( ps, "%8d\n", pE->id() );
    }
    fprintf( ps, "%12.4f %12.4f %12.4f %12.4f %12.4f\n",
            dNodeX, dNodeY, dNodeZ, dNodeU, dNodeV );
}

// ----------------------------------------------------------------------------
namespace nodeDynArray
{
void dump( FILE* ps, const nodClassDynArray& dynNodes )
{
    size_t nN = dynNodes.length();
    for( size_t in=0; in<nN; ++in )
    {
        nodClass* pN = dynNodes[in];
        if( pN != NULL )
            pN->dump( ps );
    }
}

void dump( int context, const nodClassDynArray& dynNodes )
{
    static const size_t len = 128;
    static char buf[len];
    SNPRINTF( buf, len, "DEBUG_nodDynArray_%u.txt", context );
    samFileNameUtility sfn( buf, CSamEnvVar::getSamReportPath() );
    FILE* ps = fopen( sfn.PathName(), "wbc" );
    dump( ps, dynNodes );
}
} // end namespace nodeDynArray
