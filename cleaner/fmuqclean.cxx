/*

*/

#define NM_TEST_QUAD_QUALITY 0



#define DEBUG_DumpNodeArray 0

#include "cleaner/fmuqclean.hxx"
#include "cleaner/nodclass.hxx"
#include "cleaner/fmuedge.hxx"
#include "cleaner/elmqclass.hxx"
#include "cleaner/fmumaster.hxx"
#include "cleaner/fmuvector.hxx"
#include "mesh2d/sam_gradpolyn.x"
#include "quality_check/sam_elemquality.hxx"
#include "MeshingDebugger/MeshingDebugger.h"
#include "MeshingDebugger/CMeshDebugOptions.hxx"
#include <cstdio>
#include <cmath>

const double fmu_kTWO_VALENT_ANGLE = 2.0943951;  // 120 degrees
const double fmu_kONE_VALENT_ANGLE = 4.1887902;  // 240 degrees
const double fmu_kZERO_VALENT_ANGLE = 5.7595865; // 320 degrees
const double fmu_kRIGHT_LOW_TOL = 1.3962634;     //  80 degrees
const double fmu_kRIGHT_HIGH_TOL = 1.7453293;    // 100 degrees
const double fmu_kSHAPE_ANGLE_TOL = 2.7925268;   // 160 degrees
const double fmu_kEXTREME_ANGLE_TOL = 5.4977871; // 315 degrees
const double fmu_kBNDY_ANGLE_TOL = 2.5307274;    // 145 degrees
const double fmu_kBNDY_SKEW_TOL = 0.6981317;     //  40 degrees
const double fmu_kBIG_SQUEEZE_ANGLE = 6.1086524; // 350 degrees
const double fmu_kSMALL_SQUEEZE_ANGLE = 0.43633231; // 25 degrees
const double fmu_kSMALL_BNDY_ANGLE_TOL = 0.034906585; // 2 degrees
const double fmu_kBIG_BNDY_ANGLE_TOL = 6.2482787; // 358 degrees
const double fmu_kSMALL_SHAPE_TOL = 0.52359878;  //  30 degrees
const double fmu_kHIGH_PI_TOL = 3.3161256;       // 190 degrees
const int ARRAY_SIZE = 100; //maximum 100 nodes connected to one node

const int dp_kCLR_YELLOW = -7;
const int dp_kCLR_BLUE = -8;



QuadCleanTool::QuadCleanTool()
{
    lDrawCleaning = true;
    lDrawHunting = false;
    lDrawSize = false;
}


//  quad mesh cleaning Member Function Descriptor Block -----
//
// QuadCleanTool::destructor
// class - QuadCleanTool
// 
//
// Description:
//
// destructor for the class.
//
// Access:
//
// ~QuadMeshClean()
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
// Not a whole lot to do in this destructor
//
// ------------------------------------------
QuadCleanTool::~QuadCleanTool()
{
}

//  quad mesh cleaning Member Function Descriptor Block -----
//
// QuadCleanTool::CleanMesh
// class - QuadCleanTool
// 
//
// Description:
//
// Clean a quad mesh. This is where the actual work gets done.
//
// Access:
//
// CleanMesh( pCMeshMaster, CPlaneNormal )
//
// Input Parameters:
//
// pCMeshMaster        fmuMaster*     Lists of boundary nodes, interior nodes,
//                                    edges, and elements of the mesh. The
//                                    contents of these lists may be different
//                                    on return.
// CPlaneNormal        fmuVector&     Normal vector to the infinite plane used
//                                    when meshing.
//
// Output Parameters:
//
// none
//
// Return Code:
//
// iRet      int        values:
//                          fmu_kCLEAN_NO_ACTION    = the mesh is not changed
//                          fmu_kCLEAN_MESH_CHANGED = the mesh is changed
//
// Error conditions:
//
// none
//
// Additional comments:
//
// This process cleans only quad elements. If there are tria elements in the
// mesh, it is assumed that the nodes around them have been included in a
// boundary loop so they can be treated as a hole in the mesh.
//
// Nodes that are on the boundary should be tagged as both hardpt and
// boundary. Nodes that are not to be moved or deleted should be tagged
// as hardpt. Edges that have a boundary node at both ends should be
// tagged as frozen.
//
// There are distinct parts to the mesh cleaning process:
//     valence_clean     perform cleaning operations based on the pattern
//                       of node valences around an interior node.
//     boundary_clean    perform cleaning of bad element angles along the
//                       boundary and of valence patterns around a boundary
//                       node.
//     shape_clean       perform cleaning of bad element angles.
//
// Each of these major parts calls several action routines. Various routines
// may be used in more than one major part.
//
// A major part of this routine is a loop over all the parts of the cleaning
// process. If no action is taken in one loop, the loop terminates. The loop
// also terminates after 3 iterations to avoid the risk of infinite loops
// due to cycling actions. All of the cleaning actions that will be done are
// done within the three iterations.
//
// Several groups of routines have matching call sequences, descriptions,
// and additional comments. To save file space, the comments for one group
// are not repeated in each routine of the group. Individual routines may
// still have comments appropriate to the routine.
// For valence_x_clean routines, see valence_3_clean.
// For action routines see clean_3434.
//
// ------------------------------------------
int QuadCleanTool::CleanMesh(   fmuMaster* pCMeshMaster,
                                fmuVector& CPlaneNormal )
{
    int iStatus = fmu_kCLEAN_NO_ACTION;
    int iCallStatus;
    int fSize = 1, fValency = 1, fBoundary = 1, fShape = 1;

    bool lMeshChanged;

    try
    {
        pCMaster = pCMeshMaster;

        CNormal = CPlaneNormal;

        load_size();

#if JC_DEBUG
        int nNodesOnThisBoundary=0;
        int nTotalBoundaryNodes=0 ;
        int nTotalInteriorNodes=pCMaster->global_interior()->length();

        //Boundary analysis
        nodClassDynArrayDynArray *boundaries = pCMaster->global_boundary();
        int nBoundaries = boundaries->length();
        for( int ii = 0; ii < nBoundaries; ++ii )
        {
            nodClassDynArray *boundary = (*boundaries)[ii];
            nNodesOnThisBoundary=boundary->length();

            nTotalBoundaryNodes += nNodesOnThisBoundary;
        }

        if ( (0 == nTotalInteriorNodes) && (nTotalBoundaryNodes <= 5) )
        {
            // Nothing much to do here!
            // We should not be coming to quad cleaning with triangles 
            // All triangles should have been isolated prior to quad cleaning.
            iStatus = fmu_kCLEAN_NO_ACTION;
            return iStatus;
        }
#endif

        smooth_mesh( pCMaster->global_interior() );

pCMaster->dump( 31 );

        if( lDrawCleaning )
        {
            draw_mesh();
            draw_close();
        }

        //Master loop. Execute this no more than 3 times.
        lMeshChanged = true;
#if 1
        for( int iCaseCount = 0; iCaseCount < 3 && lMeshChanged; ++iCaseCount )
#endif
        {
            lMeshChanged = false;
            //EraseElements();
            if (fSize == 1)
            {
                // Do size cleaning first - Nilanjan M
                iCallStatus = size_clean();

pCMaster->dump( 321 );

                if( iCallStatus == fmu_kCLEAN_ABORT )
                    return fmu_kCLEAN_ABORT;
                if( iCallStatus == fmu_kCLEAN_NO_MEMORY )
                    return fmu_kCLEAN_NO_MEMORY;
                if( iCallStatus == fmu_kCLEAN_MESH_CHANGED )
                    lMeshChanged = true;

                smooth_mesh( pCMaster->global_interior() );

pCMaster->dump( 322 );

                if( lDrawCleaning )
                {
                    draw_mesh();
                    draw_close();
                }
            }

            if (fValency == 1)
            {
pCMaster->dump( 331 );

                // Try valence cleaning next - Nilanjan M
                iCallStatus = valence_clean();

pCMaster->dump( 332 );

                if( iCallStatus == fmu_kCLEAN_ABORT )
                    return fmu_kCLEAN_ABORT;
                if( iCallStatus == fmu_kCLEAN_NO_MEMORY )
                    return fmu_kCLEAN_NO_MEMORY;
                if( iCallStatus == fmu_kCLEAN_MESH_CHANGED )
                    lMeshChanged = true;

                smooth_mesh( pCMaster->global_interior() );

pCMaster->dump( 333 );

                if( lDrawCleaning )
                {
                    draw_mesh();
                    draw_close();
                }
            }

            if (fBoundary == 1)
            {
                // Boundary cleaning third - Nilanjan M
                iCallStatus = boundary_clean();

pCMaster->dump( 341 );

                if( iCallStatus == fmu_kCLEAN_ABORT )
                    return fmu_kCLEAN_ABORT;
                if( iCallStatus == fmu_kCLEAN_NO_MEMORY )
                    return fmu_kCLEAN_NO_MEMORY;
                if( iCallStatus == fmu_kCLEAN_MESH_CHANGED )
                    lMeshChanged = true;

                smooth_mesh( pCMaster->global_interior() );

pCMaster->dump( 342 );

                if( lDrawCleaning )
                {
                    draw_mesh();
                    draw_close();
                }
            }

            if (fShape == 1)
            {
                // Finally shape cleaning - Nilanjan M
                iCallStatus = shape_clean();

pCMaster->dump( 351 );

                if( iCallStatus == fmu_kCLEAN_ABORT )
                    return fmu_kCLEAN_ABORT;
                if( iCallStatus == fmu_kCLEAN_NO_MEMORY )
                    return fmu_kCLEAN_NO_MEMORY;
                if( iCallStatus == fmu_kCLEAN_MESH_CHANGED )
                    lMeshChanged = true;

                smooth_mesh( pCMaster->global_interior() );

pCMaster->dump( 352 );

                if( lDrawCleaning )
                {
                    draw_mesh();
                    draw_close();
                }
            }

            if( lMeshChanged )
                iStatus = fmu_kCLEAN_MESH_CHANGED;
        }

        if( lDrawCleaning )
        {
            //Erase all the temporary graphics
            //gg3pger();
            MeshingDebugger::remove2dDisplay();
        }

        unload_size();
    }
    catch(...)
    {
        // If anything went wrong exit softly ... 
        iStatus = fmu_kCLEAN_NO_ACTION;
    }

    return iStatus;
}


int QuadCleanTool::size_clean()
{

    lDrawSize = 1;

    int ii;
    int iStatus = fmu_kCLEAN_NO_ACTION;
    int iEdgeStatus = fmu_kCLEAN_NO_ACTION;
    double dRatio;

    fmuEdge *pCEdge;
    nodClass *pCNodeA = NULL;
    nodClass *pCNodeB = NULL;

    const double kLOW_RATIO_TOL = 0.30;
    const double kHIGH_RATIO_TOL = 2.50;

    fmuEdgeDynArray *pCEdgeList = pCMaster->global_edges();
    int iInitialLength = pCEdgeList->length();
    int iMaxLength = (int) ( (double)(pCEdgeList->length()) * 1.5 );

    for( ii = 0; ii < pCEdgeList->length() &&
                 ( pCEdgeList->length() < iMaxLength ||
                 ii <= iInitialLength ); ++ii )
    {
        pCEdge = (*pCEdgeList)[ii];
        if( !pCEdge )
            continue;
        if( pCEdge->frozen() )
            continue;
        pCNodeA = pCNodeB = NULL;
        pCNodeA = pCEdge->start_node();
        pCNodeB = pCEdge->end_node();

        // Skip this edge if either of the end node
        // objects are missing - Nilanjan M
        if(!pCNodeA)
            continue;
        if(!pCNodeB)
            continue;

        if( lDrawHunting )
        {
            //draw_mesh();
            //draw_close();
            //pCEdge->highlight_me();
        }

        dRatio = edge_ratio( pCEdge );
        if( dRatio < kLOW_RATIO_TOL )
        {
            iEdgeStatus = size_small_clean( pCNodeA, pCNodeB );
        }
        else if ( dRatio > kHIGH_RATIO_TOL )
        {
            iEdgeStatus = size_big_clean( pCEdge, pCNodeA, pCNodeB );
        }

        if( iEdgeStatus == fmu_kCLEAN_ABORT )
            return fmu_kCLEAN_ABORT;
        if( iEdgeStatus == fmu_kCLEAN_NO_MEMORY )
            return fmu_kCLEAN_NO_MEMORY;
        if( iEdgeStatus == fmu_kCLEAN_MESH_CHANGED )
            iStatus = iEdgeStatus;
    }

    lDrawSize = 0;

    return iStatus;

}


int QuadCleanTool::size_small_clean( nodClass *pCNodeA, nodClass *pCNodeB )
{

    int iEdgeStatus = fmu_kCLEAN_NO_ACTION;
    nodClassDynArray CNodes(ARRAY_SIZE), CNodes1(ARRAY_SIZE);
    fmuEdgeDynArray CEdges(ARRAY_SIZE);
    elmQClassDynArray CQuads(ARRAY_SIZE);

    if( !pCNodeA->hardpt() )
    {
        int iValence = neighbors( pCNodeA, CNodes, CEdges, CQuads );
        if ( iValence == 3 )
        {	
            // Check if CNodes(20) forms a real 6 sided polygon 

            for( int jj = 0; jj< CNodes.length(); jj++)
            {
                if(CNodes1.find(CNodes[jj])>0)
                {
                    /* Do not error out, just be quiet and carry on */
                    break;
                }
                else
                {
                    CNodes1.append(CNodes[jj]);
                }
            }
            for( int jj = 0; jj < 3; ++jj )
            {
                iEdgeStatus = size_small_action( pCNodeA, CNodes, CEdges, CQuads, 3 ); 
                if( iEdgeStatus != fmu_kCLEAN_NO_ACTION )
                    return iEdgeStatus;
                CNodes.step(2);
                CEdges.step();
                CQuads.step();
            }
            for( int jj = 0; jj < 3; ++jj )
            {
                iEdgeStatus = size_small_action( pCNodeA, CNodes, CEdges, CQuads, 4 ); 
                if( iEdgeStatus != fmu_kCLEAN_NO_ACTION )
                    return iEdgeStatus;
                CNodes.step(2);
                CEdges.step();
                CQuads.step();
            }
        }
    }

    CNodes.clear();
    CEdges.clear();
    CQuads.clear();
    if( !pCNodeB->hardpt() )
    {
        int iValence = neighbors( pCNodeB, CNodes, CEdges, CQuads );
        if( iValence == 3 )
        {
            for( int jj = 0; jj < 3; ++jj )
            {
                iEdgeStatus = size_small_action( pCNodeB, CNodes, CEdges, CQuads, 3 ); 
                if( iEdgeStatus != fmu_kCLEAN_NO_ACTION )
                    return iEdgeStatus;
                CNodes.step(2);
                CEdges.step();
                CQuads.step();
            }
            for( int jj = 0; jj < 3; ++jj )
            {
                iEdgeStatus = size_small_action( pCNodeB, CNodes, CEdges, CQuads, 4 ); 
                if( iEdgeStatus != fmu_kCLEAN_NO_ACTION )
                    return iEdgeStatus;
                CNodes.step(2);
                CEdges.step();
                CQuads.step();
            }
        }
    }

    return fmu_kCLEAN_NO_ACTION;
}


int QuadCleanTool::size_small_action( nodClass *pCCenterNode,
                                          nodClassDynArray &CNodes,
                                          fmuEdgeDynArray &CEdges,
                                          elmQClassDynArray &CQuads,
                                          int iValence )
{
   
   int	iCallStatus;
   nodClass	*pCNode0 = 0, *pCNode1 = 0, *pCNode2 = 0, *pCNode3 = 0, 
			*pCNode4 = 0,
			*pCNode5 = 0;

   pCNode1 = CNodes.next();

   if( pCNode1->hardpt() ) return fmu_kCLEAN_NO_ACTION;

   if( get_valence( pCNode1 ) != iValence ) return fmu_kCLEAN_NO_ACTION;

   pCNode0 = CNodes.get();
   pCNode2 = CNodes.next(2);
   pCNode3 = CNodes.next(3);
   pCNode4 = CNodes.next(4);
   pCNode5 = CNodes.next(5);

   // Check for NULL pointers - Nilanjan M
	if (!pCNode0)
		return fmu_kCLEAN_ABORT;
	if (!pCNode1)
		return fmu_kCLEAN_ABORT;
	if (!pCNode2)
		return fmu_kCLEAN_ABORT;
	if (!pCNode3)
		return fmu_kCLEAN_ABORT;
	if (!pCNode4)
		return fmu_kCLEAN_ABORT;
	if (!pCNode5)
		return fmu_kCLEAN_ABORT;
    if (pCNode1 == pCNode3 || pCNode1 == pCNode5)    
        // This case can lead to collapsed or arrowhead quads or ones with with repeated nodes
        return fmu_kCLEAN_NO_ACTION;


   if( get_valence( pCNode0 ) == 3 &&
       pCNode4->hardpt() &&
       pCNode5->hardpt() ) return fmu_kCLEAN_NO_ACTION;

   if( get_valence( pCNode2 ) == 3 &&
       pCNode3->hardpt() &&
       pCNode4->hardpt() ) return fmu_kCLEAN_NO_ACTION;

   iCallStatus = across_remove_face( pCCenterNode, CNodes, CEdges, CQuads );

   if( iCallStatus != fmu_kCLEAN_MESH_CHANGED ) return iCallStatus;

   //Return the value from valence_2_check if an error.
   //Otherwise return mesh changed, as it already has.
   iCallStatus = valence_2_check( CNodes );
   if( iCallStatus == fmu_kCLEAN_ABORT ) return fmu_kCLEAN_ABORT;
   if( iCallStatus == fmu_kCLEAN_NO_MEMORY ) return fmu_kCLEAN_NO_MEMORY;

   return fmu_kCLEAN_MESH_CHANGED;
  
}


int QuadCleanTool::size_big_clean( fmuEdge *pCEdge, nodClass *pCNodeA,
                                                        nodClass *pCNodeB )
{
  
   int	 i = 0;
   elmQClass *pCQuadA = NULL, *pCQuadB = NULL;
   nodClassDynArray CRingNodes(ARRAY_SIZE);
   nodClassDynArray CDelNodes(ARRAY_SIZE);
   fmuEdgeDynArray CDelEdges(ARRAY_SIZE);
   elmQClassDynArray CDelQuads(ARRAY_SIZE);

   //Idiot proofing - Nilanjan M
   if (!pCNodeA)
	   return(fmu_kCLEAN_ABORT);
   if (!pCNodeB)
		return(fmu_kCLEAN_ABORT);

   pCEdge->quads( &pCQuadA, &pCQuadB );

   // If the 2 quads on this edge were not found 
   // abort cleaning, input mesh is defective.
   // - Nilanjan M
	if (!pCQuadA)
	 return(fmu_kCLEAN_ABORT);

	if (!pCQuadB)
	 return(fmu_kCLEAN_ABORT);

   CDelQuads[0] = pCQuadA;
   CDelQuads[1] = pCQuadB;
   CDelEdges[0] = pCEdge;

   // Initialize
   for(i=0; i < 6; i++)
	CRingNodes[i] = NULL;
  
   if( pCQuadA->prev_node( pCNodeA ) == pCNodeB )
   {
      CRingNodes[0] = pCNodeA;
      CRingNodes[1] = pCQuadA->next_node( pCNodeA );
      CRingNodes[2] = pCQuadA->opposite_node( pCNodeA );
      CRingNodes[3] = pCNodeB;
      CRingNodes[4] = pCQuadB->opposite_node( pCNodeA );
      CRingNodes[5] = pCQuadB->prev_node( pCNodeA );
   }
   else
   {
      CRingNodes[0] = pCNodeA;
      CRingNodes[1] = pCQuadB->next_node( pCNodeA );
      CRingNodes[2] = pCQuadB->opposite_node( pCNodeA );
      CRingNodes[3] = pCNodeB;
      CRingNodes[4] = pCQuadA->opposite_node( pCNodeA );
      CRingNodes[5] = pCQuadA->prev_node( pCNodeA );
   }

   // Check for invalid nodes - Nilanjan M
   for(i=0; i < 6; i++)
   {
	   if (CRingNodes[i] == NULL)
		   return(fmu_kCLEAN_ABORT);
   }

   return fill_choice( CDelNodes, CDelEdges, CDelQuads, CRingNodes );
  
}



//  quad mesh cleaning Member Function Descriptor Block -----
//
// QuadCleanTool::valence_clean
// class - QuadCleanTool
// 
//
// Description:
//
// Check for valence patterns in the mesh. Action routines are called
// to improve the mesh.
//
// Access:
//
// valence_clean()
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
// fmu_kCLEAN_NO_ACTION    = the mesh is as it was when the process started.
// fmu_kCLEAN_MESH_CHANGED = the mesh has been changed
//
// Error conditions:
//
// none
//
// Additional comments:
//
// Make a copy of the list of interior nodes. The original is not used as
// some nodes may need to be cleaned more than once. Some actions will put
// nodes back in the list. This, of course, does not need to be done for
// nodes that have not yet been processed. Since the list is not sorted, it
// is faster to query the node to determine if it is in the list rather
// than querying the list to determine if it contains the node. For this
// purpose, the node has a flag that is set when it is put in the list and
// unset it is when removed from the list. This is not ideal in terms
// of C++ classes, but is done only for speed.
//
// Various actions will require nodes to be cleaned again. These are put
// at the end of the cleaning list. Positions in the array are not reused.
//
// Lists are made of nodes, edges, and quads surrounding the node to be
// cleaned. These lists are passed to the pattern checking routines.
//
// ------------------------------------------
int QuadCleanTool::valence_clean()
{
    int iStatus = fmu_kCLEAN_NO_ACTION;
    int iCallStatus = 0;
    int iValence;
    nodClass *pCNode = 0;
    nodClassDynArray CNodes(ARRAY_SIZE);
    fmuEdgeDynArray CEdges(ARRAY_SIZE);
    elmQClassDynArray CQuads(ARRAY_SIZE);

    CCleaningList = (*pCMaster->global_interior());

#if DEBUG_DumpNodeArray
    nodeDynArray::dump( pCMaster->FaceId(), CCleaningList );
#endif

    //Mark all the nodes as being members of the list
    for( int ii = 0; ii < CCleaningList.length(); ii++)
        if( CCleaningList[ii] )
            CCleaningList[ii]->in_clean_list(true);

    //An occasional cycle between cases can happen. Set a limit to avoid
    //infinite loops, but make sure all nodes are checked at least once.
    int iOrigLength = CCleaningList.length();
    int iMaxLength = (int) ( (double)(CCleaningList.length()) * 1.5 );

    for( int ii = 0;
        ii < CCleaningList.length() && ( CCleaningList.length() < iMaxLength || ii <= iOrigLength );
        ++ ii )
    {
        CNodes.clear();
        CEdges.clear();
        CQuads.clear();

        pCNode = CCleaningList[ii];
        if( !pCNode )
            continue;

        int id = pCNode->id(); // debugging

        pCNode->in_clean_list(false);
        CCleaningList[ii] = NULL;

        if( lDrawHunting )
            draw_hunting( pCNode, ii );

        iValence = neighbors( pCNode, CNodes, CEdges, CQuads );

        if( iValence == 2 )
            iCallStatus = valence_2_clean( pCNode );

        else if( iValence == 3 )
            iCallStatus = valence_3_clean( pCNode, CNodes, CEdges, CQuads );

        else if (iValence == 4 )
            iCallStatus = valence_4_clean( pCNode, CNodes, CEdges, CQuads );

        else if (iValence == 5 )
            iCallStatus = valence_5_clean( pCNode, CNodes, CEdges, CQuads );

        else if (iValence > 5 )
            iCallStatus = valence_6_clean( pCNode, CNodes, CEdges, CQuads );

        if( iCallStatus == fmu_kCLEAN_ABORT )
            return fmu_kCLEAN_ABORT;
        if( iCallStatus == fmu_kCLEAN_NO_MEMORY )
            return fmu_kCLEAN_NO_MEMORY;
        if( iCallStatus == fmu_kCLEAN_MESH_CHANGED )
        {
            iStatus = iCallStatus;
            if( lDrawCleaning )
            {
                draw_mesh();
                draw_close();
            }
        }
    }

    CCleaningList.clear();

    return iStatus;
}


//  quad mesh cleaning Member Function Descriptor Block -----
//
// QuadCleanTool::valence_3_clean
// class - QuadCleanTool
// 
//
// Description:
//
// Check around 3 valent nodes for known valence patterns.
//
// Access:
//
// valence_3_clean( pCNode, CNodes, CEdges, CQuads )
//
// Input Parameters:
//
// pCNode          nodClass*          the center node of the pattern
// CNodes          nodClassDynArray&  the nodes, around the center node
// CEdges          fmuEdgeDynArray&   the edges, around the center node
// CQuads          elmQClassDynArray& the quads, around the center node
//
// Output Parameters:
//
// none
//
// Return Code:
//
// fmu_kCLEAN_NO_ACTION    = the mesh is as it was when the process started.
// fmu_kCLEAN_MESH_CHANGED = the mesh has been changed
//
// Error conditions:
//
// none
//
// Additional comments:
//
// The input lists of nodes, edges, and quads follow this pattern:
//
// n5----n4---n3    Where nc is the center node, the one being examined.
//  |    |    |     n0 to n7 are nodes 0 to 7, such that n0 and nc share
//  |    |    |     an edge. For edges, e0 connects nodes nc to n0, edge
// n6----nc---n2    e1 connectes nc to n2, edge e2 connects nc to n4,
//  |    |    |     edge e3 connects nc to n6. Quad q0 is made of nodes
//  |    |    |     nc, n0, n1, n2. Quad q1 has nodes nc, n2, n3, n4, and
// n7----n0---n1    so forth. There are twice as many nodes as edges or quads.
//
// There is a cooresponding nomenclature for nodes with valence other
// than 4.
//
// Arrays are constructed of the valence of each neighboring node and of
// whether the node is a hardpoint.
//
// Patterns are checked with whatever node is at the start of the lists.
// If no pattern is found, the lists are rotated one quad and edge and
// two nodes. The valence and hardpoint lists are rotated to match. The
// patterns are checked again with this new starting node. If all starting
// positions are checked with no pattern matched, the mesh is not changed.
//
// Valence pattern case nomenclature:
// The first node in the pattern is the valence of the center node. This
// is followed by a dash. Starting with node 0, the valence of each node
// is given. The values have these meaning:
// 3  = valence of 3
// 4- = valence of 4 or less
// 4  = valence of 4
// 4+ = valence of 4 or more
// 5  = valence of 5 or more
// 0  = valence is ignored, the node won't be affected by the mesh change.
// A blank or missing valence means the same as "0".
//
// Many valence patterns have a mirror image case. Many times the pattern
// for this new case starts such that the first node involved in the change
// is at position n0 and going counterclockwise rather than appearing to
// go clockwise through the "back" of the standard pattern. The mirror
// to pattern 4-03453500 may be designated 4-00535430 or 4-53543000.
// ------------------------------------------
int QuadCleanTool::valence_3_clean( nodClass *pCNode,
                                        nodClassDynArray &CNodes,
                                        fmuEdgeDynArray &CEdges,
                                        elmQClassDynArray &CQuads )
{
   
   int iStatus = fmu_kCLEAN_NO_ACTION;
   int iCallStatus;
   int ii;
   const int kNODE_COUNT = 6;
   int aiValence[kNODE_COUNT], aiHardpt[kNODE_COUNT];

   for( ii = 0; ii < kNODE_COUNT; ++ii )
   {
      aiValence[ii] = get_valence( CNodes[ii] );
      aiHardpt[ii] = CNodes[ii]->hardpt();
   }

   for( ii = 0; ii < 3; ++ii )
   {
      //This case checks for number of edges as valence doesn't distinguish
      //between a virtual edge on a flat side and in inside corner
      if( CNodes.get()->number_edges() > 3 &&
          aiValence[1] == 3 &&
          CNodes.next(2)->number_edges() > 3 &&
          !pCNode->hardpt() )
      {
         // case 3-4+34+
         iCallStatus = clean_3434( pCNode, CNodes, CEdges, CQuads );
         if( iCallStatus != fmu_kCLEAN_NO_ACTION )  return iCallStatus;
      }

      if( aiValence[0] == 3 &&
          aiValence[1] >= 5 &&
          aiValence[2] == 4 &&
          aiValence[3] >= 4 &&
          aiValence[4] == 4 &&
          aiValence[5] == 4 &&
          !aiHardpt[0]      &&
          !aiHardpt[4]      &&
          !pCNode->hardpt() )
      {
         // case 3-3544+44
         iCallStatus = double_three( pCNode, CNodes, CEdges, CQuads, 1 );
         if( iCallStatus != fmu_kCLEAN_NO_ACTION )  return iCallStatus;
      }

      if( aiValence[0] == 3 &&
          aiValence[1] == 4 &&
          aiValence[2] == 4 &&
          aiValence[3] >= 4 &&
          aiValence[4] == 4 &&
          aiValence[5] >= 5 &&
          !aiHardpt[0]      &&
          !aiHardpt[2]      &&
          !pCNode->hardpt() )
      {
         // case 3-3444+45
         iCallStatus = double_three( pCNode, CNodes, CEdges, CQuads, 2 );
         if( iCallStatus != fmu_kCLEAN_NO_ACTION )  return iCallStatus;
      }

      if( aiValence[0] == 3 &&
          aiValence[1] >= 4 &&
          aiValence[2] >= 4 &&
          aiValence[3] >= 4 &&
          aiValence[4] >= 4 &&
          aiValence[5] >= 4 &&
          !pCNode->hardpt()   )
      {
         // case 3-34+4+4+4+4+
         iCallStatus = double_three( pCNode, CNodes, CEdges, CQuads, 3 );
         if( iCallStatus != fmu_kCLEAN_NO_ACTION )  return iCallStatus;
      }

      if( aiValence[0] == 3 &&
          aiValence[1] == 3 &&
          !aiHardpt[0]      &&
          !aiHardpt[1]      &&
          !pCNode->hardpt() )
      {
         // case 3-33
         iCallStatus = triple_three( pCNode, CNodes, CEdges, CQuads );
         if( iCallStatus != fmu_kCLEAN_NO_ACTION )  return iCallStatus;
      }

      if( aiValence[0] == 4 &&
          aiValence[1] == 4 &&
          aiValence[2] == 4 &&
          aiValence[3] == 4 &&
          aiValence[4] == 4 &&
          aiValence[5] == 4 &&
          !aiHardpt[0] &&
          !pCNode->hardpt()   )
      {
         // case 3-444444
         iCallStatus = squeeze_double_row( pCNode, CNodes, CEdges, CQuads );
         if( iCallStatus != fmu_kCLEAN_NO_ACTION )  return iCallStatus;
      }

      rotate_list( aiValence, kNODE_COUNT, 2 );
      rotate_list( aiHardpt, kNODE_COUNT, 2 );
      CNodes.step(2);
      CEdges.step();
      CQuads.step();
   }

   return iStatus;
   
}


int QuadCleanTool::valence_4_clean( nodClass *pCNode,
                                        nodClassDynArray &CNodes,
                                        fmuEdgeDynArray &CEdges,
                                        elmQClassDynArray &CQuads )
{
   
   int iStatus = fmu_kCLEAN_NO_ACTION;
   int iCallStatus;
   int ii;
   const int kNODE_COUNT = 8;
   int aiValence[kNODE_COUNT], aiHardpt[kNODE_COUNT];

   for( ii = 0; ii < kNODE_COUNT; ++ii )
   {
	  // Check for NULL pointers - Nilanjan M
	  if (CNodes[ii] == NULL)
		  return(fmu_kCLEAN_NO_ACTION);
      aiValence[ii] = get_valence( CNodes[ii] );
      aiHardpt[ii] = CNodes[ii]->hardpt();
   }

   for( ii = 0; ii < 4; ++ii )
   {
      if( aiValence[0] == 3 &&
          aiValence[1] == 3 &&
          aiValence[2] == 3 &&
          aiValence[3] >= 4 &&
          aiValence[7] >= 4 &&
          !aiHardpt[0] &&
          !aiHardpt[1] &&
          !aiHardpt[2] )
      {
         // case 4-3334+0004+
         iCallStatus = clean_4333( pCNode, CNodes, CEdges, CQuads );
         if( iCallStatus != fmu_kCLEAN_NO_ACTION )  return iCallStatus;
      }


      if( aiValence[0] == 3 &&
          aiValence[1] >= 5 &&
          aiValence[2] == 4 &&
          aiValence[3] == 3 &&
          aiValence[4] >= 4 && 
          !aiHardpt[0]      &&
          !aiHardpt[3]      &&
          !pCNode->hardpt() )
      {
         // case 4-35435
         iCallStatus = clean_435435( pCNode, CNodes, CEdges, CQuads );
         if( iCallStatus != fmu_kCLEAN_NO_ACTION )  return iCallStatus;
      }

      if( aiValence[0] >= 4 &&
          aiValence[1] == 3 &&
          aiValence[2] == 4 &&
          aiValence[3] >= 5 &&
          aiValence[4] == 3 &&
          !aiHardpt[1]      &&
          !aiHardpt[4]      &&
          !pCNode->hardpt()  )
      {
         // case 4-53453
         iCallStatus = clean_453453( pCNode, CNodes, CEdges, CQuads );
         if( iCallStatus != fmu_kCLEAN_NO_ACTION )  return iCallStatus;
      }

      if( aiValence[0] == 3 &&
          aiValence[1] == 4 &&
          aiValence[2] == 4 &&
          aiValence[3] == 3 &&
          aiValence[4] >= 5 &&
          aiValence[5] == 4 &&
          aiValence[6] == 4 &&
          aiValence[7] >= 5 &&
          !aiHardpt[0] &&
          !aiHardpt[2] &&
          !aiHardpt[3] &&
          !pCNode->hardpt() )
      {
         // case 4-34435445
         iCallStatus = clean_434435445( pCNode, CNodes, CEdges, CQuads );
         if( iCallStatus != fmu_kCLEAN_NO_ACTION )  return iCallStatus;
      }

      if( aiValence[0] >= 5 &&
          aiValence[1] == 3 &&
          aiValence[2] == 4 &&
          aiValence[3] == 4 &&
          aiValence[4] == 3 &&
          aiValence[5] >= 5 &&
          aiValence[6] == 4 &&
          aiValence[7] == 4 &&
          !aiHardpt[1] &&
          !aiHardpt[2] &&
          !aiHardpt[4] &&
          !pCNode->hardpt() )
      {
         // case 4-53443544
         iCallStatus = clean_453443544( pCNode, CNodes, CEdges, CQuads );
         if( iCallStatus != fmu_kCLEAN_NO_ACTION )  return iCallStatus;
      }

      if( aiValence[0] == 3 &&
          aiValence[1] == 4 &&
          aiValence[2] == 4 &&
          aiValence[3] == 3 &&
          aiValence[4] == 4 &&
          !aiHardpt[0] &&
          !aiHardpt[1] &&
          !aiHardpt[2] &&
          !aiHardpt[3] &&
          !aiHardpt[4] &&
          !pCNode->hardpt() )
      {
         // case 4-34434
         iCallStatus = clean_434434( pCNode, CNodes, CEdges, CQuads,
                                                     aiValence );
         if( iCallStatus != fmu_kCLEAN_NO_ACTION )  return iCallStatus;
      }

      if( aiValence[0] == 4 &&
          aiValence[1] == 3 &&
          aiValence[2] == 4 &&
          aiValence[3] == 4 &&
          aiValence[4] == 3 &&
          !aiHardpt[1] &&
          !aiHardpt[2] &&
          !aiHardpt[3] &&
          !aiHardpt[4] &&
          !pCNode->hardpt() )
      {
         // case 4-43443
         iCallStatus = clean_443443( pCNode, CNodes, CEdges, CQuads,
                                                     aiValence );
         if( iCallStatus != fmu_kCLEAN_NO_ACTION )  return iCallStatus;
      }

      if( aiValence[0] == 3 &&
          aiValence[1] >= 5 &&
          aiValence[2] == 4 &&
          aiValence[3] >= 5 &&
          aiValence[4] == 3 &&
          aiValence[6] >= 4 &&
          !aiHardpt[0] &&
          !aiHardpt[4] &&
          !aiHardpt[6] &&
          !pCNode->hardpt() )
      {
         // case 4-3545304+0
         iCallStatus = clean_435453( pCNode, CNodes, CEdges, CQuads );
         if( iCallStatus != fmu_kCLEAN_NO_ACTION )  return iCallStatus;
      }

      if( aiValence[0] == 4 &&
          aiValence[1] == 3 &&
          aiValence[3] == 4 &&
          aiValence[4] == 4 &&
          aiValence[5] == 3 &&
          aiValence[7] == 4 &&
          ( ( aiValence[2] >= 5 && aiValence[6] >= 4 ) ||
            ( aiValence[2] >= 4 && aiValence[6] >= 5 ) ) &&
          !pCNode->hardpt() )
      {
         // case 4-434+44354 or 4-4354434+4
         iCallStatus = clean_443544354( pCNode, CNodes, CEdges, CQuads, 1 );
         if( iCallStatus != fmu_kCLEAN_NO_ACTION )  return iCallStatus;
      }

      if( aiValence[1] == 3 &&
          aiValence[2] == 4 &&
          aiValence[3] == 4 &&
          aiValence[5] == 3 &&
          aiValence[6] == 4 &&
          aiValence[7] == 4 &&
          ( ( aiValence[0] >= 5 && aiValence[4] >= 4 ) ||
            ( aiValence[0] >= 4 && aiValence[4] >= 5 ) ) &&
          !pCNode->hardpt() )
      {
         // case 4-4+3445344 or 4-53444+344
         iCallStatus = clean_443544354( pCNode, CNodes, CEdges, CQuads, 2 );
         if( iCallStatus != fmu_kCLEAN_NO_ACTION )  return iCallStatus;
      }

      if( aiValence[0] == 3 &&
          aiValence[2] == 4 &&
          aiValence[3] == 4 &&
          aiValence[4] == 3 &&
          aiValence[6] == 4 &&
          aiValence[7] == 4 &&
          ( ( aiValence[1] >= 5 && aiValence[5] >= 4 ) ||
            ( aiValence[1] >= 4 && aiValence[5] >= 5 ) ) &&
          !aiHardpt[0] &&
          !aiHardpt[4] &&
          !pCNode->hardpt() )
      {
         // case 4-34+443544 or 4-354434+44
         iCallStatus = clean_435443544( pCNode, CNodes, CEdges, CQuads, 1 );
         if( iCallStatus != fmu_kCLEAN_NO_ACTION )  return iCallStatus;
      }

      if( aiValence[0] == 3 &&
          aiValence[1] == 4 &&
          aiValence[2] == 4 &&
          aiValence[4] == 3 &&
          aiValence[5] == 4 &&
          aiValence[6] == 4 &&
          ( ( aiValence[3] >= 5 && aiValence[7] >= 4 ) ||
            ( aiValence[3] >= 4 && aiValence[7] >= 5 ) ) &&
          !aiHardpt[0] &&
          !aiHardpt[4] &&
          !pCNode->hardpt() )
      {
         // case 4-3444+3445 or 4-34453444+
         iCallStatus = clean_435443544( pCNode, CNodes, CEdges, CQuads, 2 );
         if( iCallStatus != fmu_kCLEAN_NO_ACTION )  return iCallStatus;
      }

      if( aiValence[0] == 3 &&
          aiValence[1] == 4 &&
          aiValence[3] == 4 &&
          aiValence[4] == 3 &&
          aiValence[5] == 4 &&
          aiValence[7] == 4 &&
          ( ( aiValence[2] >= 5 && aiValence[6] >= 4 ) ||
            ( aiValence[2] >= 4 && aiValence[6] >= 5 ) ) &&
          !aiHardpt[0]      &&
          !aiHardpt[4]      &&
          !pCNode->hardpt()  )
      {
         // case 4-344+43454 or 4-3454344+4
         iCallStatus = double_triangle( pCNode, CNodes, CEdges, CQuads,
                                                     aiValence, aiHardpt );
         if( iCallStatus != fmu_kCLEAN_NO_ACTION )  return iCallStatus;
      }

      if( aiValence[0] == 3 &&
          aiValence[1] == 4 &&
          aiValence[2] >= 5 &&
          aiValence[3] == 3 &&
          !aiHardpt[2] )
      {
         // case 4-3453
         iCallStatus = rotate_edge( pCNode, CNodes, CEdges, CQuads, 2 );
         if( iCallStatus != fmu_kCLEAN_NO_ACTION )  return iCallStatus;
      }

      if( aiValence[1] == 3 &&
          aiValence[2] >= 5 &&
          aiValence[3] == 4 &&
          aiValence[4] == 3 &&
          !aiHardpt[2] )
      {
         // case 4-03543
         iCallStatus = rotate_edge( pCNode, CNodes, CEdges, CQuads, 1 );
         if( iCallStatus != fmu_kCLEAN_NO_ACTION )  return iCallStatus;
      }

      if( aiValence[1] == 3 &&
          aiValence[2] >= 5 &&
          aiValence[3] == 3 &&
          !aiHardpt[2] )
      {
         // case 4-0353
         iCallStatus = reverse_open_face( pCNode, CNodes, CEdges, CQuads );
         if( iCallStatus != fmu_kCLEAN_NO_ACTION )  return iCallStatus;
      }

      if( aiValence[0] == 3 &&
          aiValence[1] == 4 &&
          aiValence[2] >= 5 &&
          aiValence[3] == 4 &&
          aiValence[4] == 3 &&
          aiValence[5] >= 5 &&
          aiValence[6] == 4 &&
          aiValence[7] == 4 &&
          !pCNode->hardpt() &&
          !aiHardpt[4] )
      {
         // case 4-34543544
         iCallStatus = clean_434543544( pCNode, CNodes, CEdges, CQuads );
         if( iCallStatus != fmu_kCLEAN_NO_ACTION )  return iCallStatus;
      }

      if( aiValence[0] == 3 &&
          aiValence[1] == 4 &&
          aiValence[2] >= 4 &&
          aiValence[3] == 5 &&
          aiValence[4] == 3 &&
          aiValence[5] >= 4 &&
          aiValence[6] == 5 &&
          aiValence[7] == 4 &&
          !pCNode->hardpt() &&
          !aiHardpt[4] )
      {
         // case 4-34453454 mirror of 4-34543544 
         iCallStatus = clean_434453454( pCNode, CNodes, CEdges, CQuads );
         if( iCallStatus != fmu_kCLEAN_NO_ACTION )  return iCallStatus;
      }
 
      rotate_list( aiValence, kNODE_COUNT, 2 );
      rotate_list( aiHardpt, kNODE_COUNT, 2 );
      CNodes.step(2);
      CEdges.step();
      CQuads.step();
   }

   return iStatus;

}



int QuadCleanTool::valence_5_clean( nodClass *pCNode,
                                        nodClassDynArray &CNodes,
                                        fmuEdgeDynArray &CEdges,
                                        elmQClassDynArray &CQuads )
{

   int iStatus = fmu_kCLEAN_NO_ACTION;
   int iCallStatus;
   int ii;
   const int kNODE_COUNT = 10;
   int aiValence[kNODE_COUNT], aiHardpt[kNODE_COUNT];

	// Check for NULL pointers - Nilanjan M
	if (pCNode == NULL)
		return(fmu_kCLEAN_NO_ACTION);

   for( ii = 0; ii < kNODE_COUNT; ++ii )
   {	  
	  // Check for NULL pointers - Nilanjan M
	  if (CNodes[ii] == NULL)
		  return(fmu_kCLEAN_NO_ACTION);
      aiValence[ii] = get_valence( CNodes[ii] );
      aiHardpt[ii] = CNodes[ii]->hardpt();
   }

   for( ii = 0; ii < 5; ++ii )
   {
      if( aiValence[0] == 3 &&
          aiValence[5] == 3 &&
          !pCNode->hardpt() )
      {
         // case 5-300003
         iCallStatus = clean_5300003( pCNode, CNodes, CEdges, CQuads );
         if( iCallStatus != fmu_kCLEAN_NO_ACTION )  return iCallStatus;
      }

      if( aiValence[0] == 3 &&
          aiValence[1] >= 4 &&
          aiValence[2] == 4 &&
          aiValence[3] >= 4 &&
          aiValence[4] == 3  )
      {
         // case 5-34+44+3
         iCallStatus = across_open_face( pCNode, CNodes, CEdges, CQuads );
         if( iCallStatus != fmu_kCLEAN_NO_ACTION )  return iCallStatus;
      }

      if( aiValence[0] >= 5 &&
          aiValence[1] == 4 &&
          aiValence[2] == 4 &&
          aiValence[3] == 3   )
      {
         // case 5-5443
         iCallStatus = clean_55443( pCNode, CNodes, CEdges, CQuads );
         if( iCallStatus != fmu_kCLEAN_NO_ACTION )  return iCallStatus;
      }

      if( aiValence[1] == 3 &&
          aiValence[2] == 4 &&
          aiValence[3] == 4 &&
          aiValence[4] >= 5   )
      {
         // case 5-03445
         iCallStatus = clean_503445( pCNode, CNodes, CEdges, CQuads );
         if( iCallStatus != fmu_kCLEAN_NO_ACTION )  return iCallStatus;
      }

      if( aiValence[0] == 4 &&
          aiValence[1] == 3 &&
          aiValence[2] == 4 &&
          aiValence[3] == 4 &&
          aiValence[4] == 4 &&
          aiValence[8] == 4 &&
          aiValence[9] == 4 &&
          !aiHardpt[1]      &&
          !aiHardpt[2] )
      {
         // case 5-4344400044
         iCallStatus = zip_unzip( pCNode, CNodes, CEdges, CQuads, aiValence,
                                                              aiHardpt );
         if( iCallStatus != fmu_kCLEAN_NO_ACTION )  return iCallStatus;
      }

      if( aiValence[0] == 3 &&
          aiValence[1] >= 4 &&
          aiValence[2] >= 4 &&
          aiValence[3] == 3 &&
          !aiHardpt[2] )
      {
         // case 5-34+4+3
         iCallStatus = rotate_edge( pCNode, CNodes, CEdges, CQuads, 2 );
         if( iCallStatus != fmu_kCLEAN_NO_ACTION )  return iCallStatus;
      }

      if( aiValence[1] == 3 &&
          aiValence[2] >= 4 &&
          aiValence[3] >= 4 &&
          aiValence[4] == 3 &&
          !aiHardpt[2] )
      {
         // case 5-034+4+3
         iCallStatus = rotate_edge( pCNode, CNodes, CEdges, CQuads, 1 );
         if( iCallStatus != fmu_kCLEAN_NO_ACTION )  return iCallStatus;
      }

      if( aiValence[0] >= 5 &&
          aiValence[1] == 3 &&
          aiValence[2] == 3 &&
          aiValence[8] <= 4 &&
          aiValence[9] <= 4 &&
          aiHardpt[3] &&
          aiHardpt[4] )
      {
         // case 5-533000004-4-
         //Note that nodes 3 and 4 and center or 0 are hardpoints,
         //otherwise this should be handled by double three
         CNodes.step(-2);
         CEdges.step(-1);
         CQuads.step(-1);
         if( aiHardpt[0] && !pCNode->hardpt() )
         {
            iCallStatus = rotate_edge( pCNode, CNodes, CEdges, CQuads, 2 );
            if( iCallStatus != fmu_kCLEAN_NO_ACTION )  return iCallStatus;
         }
         else if( !aiHardpt[0] && pCNode->hardpt() )
         {
            iCallStatus = rotate_edge( pCNode, CNodes, CEdges, CQuads, 1 );
            if( iCallStatus != fmu_kCLEAN_NO_ACTION )  return iCallStatus;
         }
         //Reset lists since action was not done.
         CNodes.step(2);
         CEdges.step(1);
         CQuads.step(1);
      }


      if( aiValence[0] == 3 &&
          aiValence[1] == 3 &&
          aiValence[2] >= 5 &&
          aiValence[3] <= 4 &&
          aiValence[4] <= 4 &&
          aiHardpt[8] &&
          aiHardpt[9] )
      {
         // case 5-3354-4-
         //Note that nodes 8 and 9 and center or 2 are hardpoints,
         //otherwise this should be handled by double three
         if( aiHardpt[2] && !pCNode->hardpt() )
         {
            iCallStatus = rotate_edge( pCNode, CNodes, CEdges, CQuads, 1 );
            if( iCallStatus != fmu_kCLEAN_NO_ACTION )  return iCallStatus;
         }
         else if( !aiHardpt[2] && pCNode->hardpt() )
         {
            iCallStatus = rotate_edge( pCNode, CNodes, CEdges, CQuads, 2 );
            if( iCallStatus != fmu_kCLEAN_NO_ACTION )  return iCallStatus;
         }
      }

      if( aiValence[0] >= 5 &&
          aiValence[1] == 4 &&
          aiValence[2] == 4 &&
          aiValence[3] == 4 &&
          aiValence[4] == 3 &&
          ( !aiHardpt[0] || !pCNode->hardpt() ) )
      {
         iCallStatus = clean_554443( pCNode, CNodes, CEdges, CQuads );
         if( iCallStatus != fmu_kCLEAN_NO_ACTION )  return iCallStatus;
      }

      if( aiValence[0] == 3 &&
          aiValence[1] == 4 &&
          aiValence[2] == 4 &&
          aiValence[3] == 4 &&
          aiValence[4] >= 5 &&
          ( !aiHardpt[4] || !pCNode->hardpt() ) )
      {
         iCallStatus = clean_534445( pCNode, CNodes, CEdges, CQuads );
         if( iCallStatus != fmu_kCLEAN_NO_ACTION )  return iCallStatus;
      }

      rotate_list( aiValence, kNODE_COUNT, 2 );
      rotate_list( aiHardpt, kNODE_COUNT, 2 );
      CNodes.step(2);
      CEdges.step();
      CQuads.step();
   }

   return iStatus;
 
}


int QuadCleanTool::valence_6_clean( nodClass *pCNode,
                                        nodClassDynArray &CNodes,
                                        fmuEdgeDynArray &CEdges,
                                        elmQClassDynArray &CQuads )
{
  
   int iStatus = fmu_kCLEAN_NO_ACTION;
   int iCallStatus;
   int ii;
   //Array size and node count are separate as this routine
   //is used for valence of 6 and greater.
   //  const int ARRAY_SIZE = 20;
   int node_count = CNodes.length();
   int aiValence[ARRAY_SIZE], aiHardpt[ARRAY_SIZE];

	// Check for NULL pointers - Nilanjan M
	if (pCNode == NULL)
		return(fmu_kCLEAN_NO_ACTION);

   for( ii = 0; ii < node_count; ++ii )
   {
	  // Check for NULL pointers - Nilanjan M
	  if (CNodes[ii] == NULL)
		  return(fmu_kCLEAN_NO_ACTION);
      aiValence[ii] = get_valence( CNodes[ii] );
      aiHardpt[ii] = CNodes[ii]->hardpt();
   }

   for( ii = 0; ii < node_count; ++ii )
   {

      if( aiValence[0] == 3 &&
          aiValence[2] >= 4 &&
          aiValence[3] == 3 &&
          !aiHardpt[2] )
      {
         // case 6-304+3
         iCallStatus = rotate_edge( pCNode, CNodes, CEdges, CQuads, 2 );
         if( iCallStatus != fmu_kCLEAN_NO_ACTION )  return iCallStatus;
      }

      if( aiValence[1] == 3 &&
          aiValence[2] >= 4 &&
          aiValence[4] == 3 &&
          !aiHardpt[2] )
      {
         // case 6-034+03
         iCallStatus = rotate_edge( pCNode, CNodes, CEdges, CQuads, 1 );
         if( iCallStatus != fmu_kCLEAN_NO_ACTION )  return iCallStatus;
      }

      if( aiValence[0] == 3 &&
          aiValence[1] >= 4 &&
          aiValence[2] == 4 &&
          aiValence[3] >= 4 &&
          aiValence[4] == 3  )
      {
         // case 6-34+44+3
         iCallStatus = across_open_face( pCNode, CNodes, CEdges, CQuads );
         if( iCallStatus != fmu_kCLEAN_NO_ACTION )  return iCallStatus;
      }

      if( aiValence[0] == 3 &&
          aiValence[1] == 4 &&
          aiValence[2] == 4 &&
          aiValence[3] == 4 &&
          aiValence[4] == 4 &&
          aiValence[5] == 4 &&
          aiValence[6] == 3  )
      {
         // case 6-3444443
         iCallStatus = crack_open_middle( pCNode, CNodes, CEdges, CQuads );
         if( iCallStatus != fmu_kCLEAN_NO_ACTION )  return iCallStatus;
      }

      rotate_list( aiValence, node_count, 2 );
      rotate_list( aiHardpt, node_count, 2 );
      CNodes.step(2);
      CEdges.step();
      CQuads.step();
   }

   //Last chance efforts. Check all possible orientations of each before
   //going on to the next.
   for( ii = 0; ii < node_count; ++ii )
   {

      if( aiValence[0] == 3 &&
          !aiHardpt[6] )
      {
         // case 6-3 without opposite hardpoint
         iCallStatus = crack_open_middle( pCNode, CNodes, CEdges, CQuads );
         if( iCallStatus != fmu_kCLEAN_NO_ACTION )  return iCallStatus;
      }

      rotate_list( aiValence, node_count, 2 );
      rotate_list( aiHardpt, node_count, 2 );
      CNodes.step(2);
      CEdges.step();
      CQuads.step();
   }

   for( ii = 0; ii < node_count; ++ii )
   {

      if( aiValence[0] == 4 &&
          aiValence[6] == 4 &&
          !aiHardpt[0] &&
          !aiHardpt[6] )
      {
         // case 6-4
         iCallStatus = crack_open_middle( pCNode, CNodes, CEdges, CQuads );
         if( iCallStatus != fmu_kCLEAN_NO_ACTION )  return iCallStatus;
      }

      rotate_list( aiValence, node_count, 2 );
      rotate_list( aiHardpt, node_count, 2 );
      CNodes.step(2);
      CEdges.step();
      CQuads.step();
   }

   for( ii = 0; ii < node_count; ++ii )
   {

      if( aiValence[0] == 3 &&
          aiValence[6] < 5)
      {
         // case 6-3 with opposite hardpoint
         iCallStatus = crack_open_middle( pCNode, CNodes, CEdges, CQuads );
         if( iCallStatus != fmu_kCLEAN_NO_ACTION )  return iCallStatus;
      }

      rotate_list( aiValence, node_count, 2 );
      rotate_list( aiHardpt, node_count, 2 );
      CNodes.step(2);
      CEdges.step();
      CQuads.step();
   }

   return iStatus;
 
}

//  quad mesh cleaning Member Function Descriptor Block -----
//
// QuadCleanTool::clean_3434
// class - QuadCleanTool
// 
//
// Description:
//
// Handle valence pattern 3-434. Much of this description is true for
// all of the action rotines.
//
// Access:
//
// clean_3434( pCCenterNode, CNodes, CEdges, CQuads )
//
// Input Parameters:
//
// pCCenterNode    nodClass*          the center node of the pattern
// CNodes          nodClassDynArray&  the nodes, around the center node
// CEdges          fmuEdgeDynArray&   the edges, around the center node
// CQuads          elmQClassDynArray& the quads, around the center node
//
// Output Parameters:
//
// none
//
// Return Code:
//
// fmu_kCLEAN_NO_ACTION    = the mesh is as it was when the process started.
// fmu_kCLEAN_MESH_CHANGED = the mesh has been changed
//
// Error conditions:
//
// none
//
// Additional comments:
//
// A standard variable naming convention is used in the action routines.
// pCNode0 will refer to the node in "CNodes.get()". This is "n0" in the
// valence pattern descriptions. pCNode1 is in "CNodes.next()" and pCNode2
// is in "CNodes.next(2)", etc. In the same manner, edges are named pCEdge0,
// pCEdge1, etc. and quads are pCQuad0, etc. Edges around the perimiter of
// the pattern are designated by the nodes of the edge. Edge pCEdge01 is
// the edge from pCNode0 to pCNode1.
//
// Several action routines use nodes, edges, and faces outside of the
// original circle. These nodes and quads will frequently be given
// designations of pCNodeA and pCQuadB. Their position is dependent on
// the case. pCEdge3B is the edge from pCNode3 to pCNodeB.
//
// Action routines may be full routines that delete part of the mesh and
// create a replacement mesh. Others may simply do more checking before
// calling a full routine. A full action routine will contain the following
// parts:
//    Pull the required nodes, edges, and quads out of the input lists,
//       checking for various factors that would prevent the action from
//       being completed. This must include checks for whether the edges
//       to be deleted are frozen. It is assumed that there are no quads
//       on the other side of a frozen edge.
//    Draw the quads to be deleted, within the block "if( lDrawCleaning )"
//       See examples on how to open and close the graphics. Try to select
//       a different color combination for each action so they can be
//       distinguished (to the practiced eye) by their colors.
//    Actually perform the action. Delete the affected quads, edges, and
//       nodes, and create new ones. The deletes must be in the order of
//       quads, edges, nodes to avoid deleting a node that has an edge
//       referring to it. The routines "delete_this_node", "delete_this_edge"
//       and "delete_this_quad" must be used to ensure that various lists
//       are kept up to date.
//    Smooth the nodes that were involved. Input to smoothing does not
//       have to specify all nodes as large moves will smooth neighbors.
//       If the action is executed only within valence cleaning, smoothing
//       should be done within the "after" draw to save time.
//    Draw the new quads, similar to the "before" block.
//    If the action is a part of valence clean, call "reclean" for the
//       nodes that had their valence changed. Again, not all nodes have
//       to be listed as all the nodes around each one in the call are also
//       added to the list to be recleaned.
//
// ------------------------------------------
int QuadCleanTool::clean_3434( nodClass *pCCenterNode,
                                   nodClassDynArray &CNodes,
                                   fmuEdgeDynArray &CEdges,
                                   elmQClassDynArray &CQuads )
{
   
   nodClass *pCNode0 = 0;
   nodClass *pCNode1 = 0;
   nodClass *pCNode2 = 0;
   nodClass *pCNode3 = 0;
   nodClass *pCNode4 = 0;
   nodClass *pCNode5 = 0;

   // Idiot-proofing - Nilanjan M
   if (!pCCenterNode)
	   return(fmu_kCLEAN_NO_ACTION);

   if( pCCenterNode->hardpt() ) return fmu_kCLEAN_NO_ACTION;

   pCNode0 = CNodes.get();
   pCNode1 = CNodes.next();
   pCNode2 = CNodes.next(2);
   pCNode3 = CNodes.next(3);
   pCNode4 = CNodes.next(4);
   pCNode5 = CNodes.next(5);

   // Check for NULL pointers - Nilanjan M
   if (	pCNode0 == NULL ||
		pCNode1 == NULL	||
		pCNode2 == NULL	||
		pCNode3 == NULL ||
		pCNode4 == NULL ||
		pCNode5 == NULL)
	return(fmu_kCLEAN_NO_ACTION);

   if( pCNode5->hardpt() &&
       pCNode0->hardpt() &&
       pCNode1->hardpt() )
   {
      fmuVector CVec_01( pCNode0, pCNode1 );
      fmuVector CVec_05( pCNode0, pCNode5 );
      if( CNormal.angle( CVec_01, CVec_05 ) >= fmu_kSHAPE_ANGLE_TOL )
        return fmu_kCLEAN_NO_ACTION;
   }

   if( pCNode1->hardpt() &&
       pCNode2->hardpt() &&
       pCNode3->hardpt() )
   {
      fmuVector CVec_21( pCNode2, pCNode1 );
      fmuVector CVec_23( pCNode2, pCNode3 );
      if( CNormal.angle( CVec_23, CVec_21 ) >= fmu_kSHAPE_ANGLE_TOL )
        return fmu_kCLEAN_NO_ACTION;
   }

   if( pCNode0->hardpt() )
   {
      fmuVector CVec_41( pCNode4, pCNode1 );
      fmuVector CVec_50( pCNode5, pCNode0 );
      if( CVec_41.length() > 1.3 * CVec_50.length() )
          return fmu_kCLEAN_NO_ACTION;
   }

   return across_remove_face( pCCenterNode, CNodes, CEdges, CQuads );
 
}


int QuadCleanTool::across_remove_face( nodClass *pCCenterNode,
                                           nodClassDynArray &CNodes,
                                           fmuEdgeDynArray &CEdges,
                                           elmQClassDynArray &CQuads )
{
    int ii=0;
    nodClass *PCMergeNode = 0;

    // Idiot-proofing - Nilanjan M
    if (!pCCenterNode)
        return(fmu_kCLEAN_NO_ACTION);

    //All edges and all quads must not be frozen. Test edges now. Test
    //quads if any quads are set to frozen.

    for( ii = 0; ii < CEdges.length(); ++ii )
    {
        if( CEdges[ii]->frozen() )
            return fmu_kCLEAN_NO_ACTION;
    }

    PCMergeNode = CNodes.next();
    if (!PCMergeNode)
        return(fmu_kCLEAN_NO_ACTION);

    //Delete all the quads and around the center node and rebuild around
    //the merge node.
    int nQuads = CQuads.length();
    if( lDrawCleaning )
    {
        //draw_open();
        // Refresh debug graphics
        //EraseElements();
        for( int ii = 0; ii < nQuads; ++ii )
        {
            CQuads[ii]->draw_me( dp_kCLR_YELLOW ); 
        }
        //draw_close();
    }
    
    // Check the quality of the existing elements
    nodClass *apCNodes[4]={0};
    double dQS=0.0, dQE=0.0;
    for( int ii = 0; ii < nQuads; ++ii )
    {
        CQuads[ii]->nodes(apCNodes);
        // Compute Jac ratio 
        dQE = test_quad_quality( apCNodes[0]->node_x(), apCNodes[0]->node_y(),
                                 apCNodes[1]->node_x(), apCNodes[1]->node_y(),
                                 apCNodes[2]->node_x(), apCNodes[2]->node_y(),
                                 apCNodes[3]->node_x(), apCNodes[3]->node_y());

        dQS += dQE;
    }

    dQS /= double(nQuads);

    // Check quality of future elements 
    double dQSum=0.0, dQ=0.0;
    CNodes.step(2);
    for( int ii = 0; ii < (nQuads-1); ++ii )
    {
        // Compute Jac ratio of the would-be elements
        dQ = test_quad_quality( PCMergeNode->node_x(), PCMergeNode->node_y(),
                                       (CNodes.get())->node_x(), (CNodes.get())->node_y(),
                                       (CNodes.next())->node_x(), (CNodes.next())->node_y(),
                                       (CNodes.next(2))->node_x(), (CNodes.next(2))->node_y());
        dQSum += dQ;

        CNodes.step(2);
    }

    dQSum /= double(nQuads-1);


    for( int ii = 0; ii < nQuads; ++ii )
    {
        delete_this_quad( CQuads[ii] );
        CQuads[ii] = NULL;
    }
    for( int ii = 0; ii < CEdges.length(); ++ii)
    {
        if( delete_this_edge( CEdges[ii] ) == fmu_kCLEAN_ABORT )
            return fmu_kCLEAN_ABORT;
    }

    if( delete_this_node( pCCenterNode ) == fmu_kCLEAN_ABORT )
        return fmu_kCLEAN_ABORT;

    //Make one fewer quads than what we started with.
    --nQuads;

    CNodes.step(2);

    for( int ii = 0; ii < nQuads; ++ii )
    {
        CQuads[ii] = new elmQClass( PCMergeNode, CNodes.get(), CNodes.next(),
                                    CNodes.next(2), pCMaster );
        if( !CQuads[ii] )
            return fmu_kCLEAN_NO_MEMORY;

        CNodes.step(2);
    }

    smooth_mesh( &CNodes );

    dQSum = 0.0; 
    CNodes.step(2);
    for( int ii = 0; ii < nQuads; ++ii )
    {
        // Compute Jac ratio of the would-be elements
        dQ = test_quad_quality( PCMergeNode->node_x(), PCMergeNode->node_y(),
                                       (CNodes.get())->node_x(), (CNodes.get())->node_y(),
                                       (CNodes.next())->node_x(), (CNodes.next())->node_y(),
                                       (CNodes.next(2))->node_x(), (CNodes.next(2))->node_y());
        dQSum += dQ;

        CNodes.step(2);
    }

    dQSum /= double(nQuads-1);
    if( lDrawCleaning )
    {
        //draw_mesh();
        // Refresh debug graphics
        //EraseElements();
        for( ii = 0; ii < nQuads; ++ii )
        {
            CQuads[ii]->draw_me( dp_kCLR_BLUE );
        }
        //draw_close();
    }

    for( int ii = 0; ii < CNodes.length(); ++ii )
        reclean( CNodes[ii] );

    return fmu_kCLEAN_MESH_CHANGED;
}


int QuadCleanTool::double_three( nodClass *pCCenterNode,
                                     nodClassDynArray &CNodes,
                                     fmuEdgeDynArray &CEdges,
                                     elmQClassDynArray &CQuads,
                                     int iSide )
{
   
   int ii, iCallStatus = fmu_kCLEAN_NO_ACTION;
   nodClassDynArray CDeleteNodes(ARRAY_SIZE);
   nodClass *pCNode0 = 0;
   nodClass *pCNode1 = 0;
   nodClass *pCNode2 = 0;
   nodClass *pCNode3 = 0;
   nodClass *pCNode4 = 0;
   nodClass *pCNode5 = 0;

   //All edges must not be frozen
   for( ii = 0; ii < CEdges.length(); ++ii )
   {
	  if (CEdges[ii] == NULL)
		return fmu_kCLEAN_NO_ACTION;
      if( CEdges[ii]->frozen() ) 
		return fmu_kCLEAN_NO_ACTION;
   }

   pCNode0 = CNodes.get();
   pCNode1 = CNodes.next();
   pCNode2 = CNodes.next(2);
   pCNode3 = CNodes.next(3);
   pCNode4 = CNodes.next(4);
   pCNode5 = CNodes.prev();

   // Check for NULL pointers - Nilanjan M
   if (	pCNode0 == NULL ||
		pCNode1 == NULL	||
		pCNode2 == NULL	||
		pCNode3 == NULL ||
		pCNode4 == NULL ||
		pCNode5 == NULL)
	return(fmu_kCLEAN_NO_ACTION);

   if(pCNode0->hardpt()) return fmu_kCLEAN_NO_ACTION;

   //Get the quad on the other side of the second 3 valent node
   fmuEdge  *pCEdge01 = 0;
   pCEdge01 = pCNode0->shared_edge( pCNode1 );
   if( pCEdge01 == NULL ||
	   pCEdge01->frozen() ) return fmu_kCLEAN_NO_ACTION;
   fmuEdge  *pCEdge05 = 0;
   pCEdge05 = pCNode0->shared_edge( pCNode5 );
   if( pCEdge05 == NULL ||
	   pCEdge05->frozen() ) return fmu_kCLEAN_NO_ACTION;
   elmQClass *pCOuterQuad = 0;
   pCOuterQuad = pCEdge01->shared_quad( pCEdge05 );
   if (pCOuterQuad == NULL)
		return fmu_kCLEAN_NO_ACTION;
   nodClass *pCOuterNode = 0;
   pCOuterNode = pCOuterQuad->opposite_node( pCNode0 );
   if (pCOuterNode == NULL)
	   return fmu_kCLEAN_NO_ACTION;

   if( iSide == 1 &&
       pCOuterNode->hardpt() &&
       pCNode1->hardpt() &&
       pCNode2->hardpt() )
   {
      fmuVector CVec_12( pCNode1, pCNode2 );
      fmuVector CVec_1O( pCNode1, pCOuterNode );
      if( CNormal.angle( CVec_12, CVec_1O ) >= fmu_kSHAPE_ANGLE_TOL )
        return fmu_kCLEAN_NO_ACTION;
   }

   if( iSide == 1 &&
       pCNode3->hardpt() &&
       pCNode4->hardpt() &&
       pCNode5->hardpt() )
   {
      fmuVector CVec_45( pCNode4, pCNode5 );
      fmuVector CVec_43( pCNode4, pCNode3 );
      if( CNormal.angle( CVec_45, CVec_43 ) >= fmu_kSHAPE_ANGLE_TOL )
        return fmu_kCLEAN_NO_ACTION;
   }

   if( iSide == 2 &&
       pCNode1->hardpt() &&
       pCNode2->hardpt() &&
       pCNode3->hardpt() )
   {
      fmuVector CVec_23( pCNode2, pCNode3 );
      fmuVector CVec_21( pCNode2, pCNode1 );
      if( CNormal.angle( CVec_23, CVec_21 ) >= fmu_kSHAPE_ANGLE_TOL )
        return fmu_kCLEAN_NO_ACTION;
   }

   if( iSide == 2 &&
       pCNode4->hardpt() &&
       pCNode5->hardpt() &&
       pCOuterNode->hardpt() )
   {
      fmuVector CVec_5O( pCNode5, pCOuterNode );
      fmuVector CVec_54( pCNode5, pCNode4 );
      if( CNormal.angle( CVec_5O, CVec_54 ) >= fmu_kSHAPE_ANGLE_TOL )
        return fmu_kCLEAN_NO_ACTION;
   }

   //CNodes is now the outer loop of nodes
   CNodes.replace( pCOuterNode );
   //These are the interior nodes that will disappear
   CDeleteNodes[0] = pCCenterNode;
   CDeleteNodes[1] = pCNode0;
   //edges and quads will disappear
   CEdges.append( pCEdge01 );
   CEdges.append( pCEdge05 );
   CQuads.append( pCOuterQuad );

   if( iSide == 1 )
   {
      CNodes.step(-1);
      iCallStatus = fill_two( CDeleteNodes, CEdges, CQuads, CNodes );
   }
   else if( iSide == 2 )
   {
      CNodes.step(1);
      iCallStatus = fill_two( CDeleteNodes, CEdges, CQuads, CNodes );
   }
   else if( iSide == 3 )
   {
      iCallStatus = fill_choice( CDeleteNodes, CEdges, CQuads, CNodes );
   }

   return iCallStatus;
   
}


int QuadCleanTool::triple_three( nodClass *pCCenterNode,
                                     nodClassDynArray &CNodes,
                                     fmuEdgeDynArray &CEdges,
                                     elmQClassDynArray &CQuads )
{
    // Check for a valid center node
    if (pCCenterNode == NULL)
        return fmu_kCLEAN_NO_ACTION;

    nodClass* pCNode0 = CNodes.get();
    nodClass* pCNode1 = CNodes.next();
    nodClass* pCNode2 = CNodes.next(2);
    nodClass* pCNode3 = CNodes.next(3);
    nodClass* pCNode4 = CNodes.next(4);
    nodClass* pCNode5 = CNodes.next(5);

    if (pCNode0 == NULL ||
        pCNode1 == NULL ||
        pCNode2 == NULL ||
        pCNode3 == NULL ||
        pCNode4 == NULL ||
        pCNode5 == NULL)
        return fmu_kCLEAN_NO_ACTION;

    fmuEdge* pCEdge0 = CEdges.get();
    fmuEdge* pCEdge1 = CEdges.next();
    fmuEdge* pCEdge2 = CEdges.next(2);

    if( pCEdge0 == NULL ||
        pCEdge0->frozen() )
        return fmu_kCLEAN_NO_ACTION;
    if( pCEdge1 == NULL ||
        pCEdge1->frozen() )
        return fmu_kCLEAN_NO_ACTION;
    if( pCEdge1 == NULL ||
        pCEdge2->frozen() )
        return fmu_kCLEAN_NO_ACTION;

    elmQClass* pCQuad0 = CQuads.get();
    elmQClass* pCQuad1 = CQuads.next();
    elmQClass* pCQuad2 = CQuads.next(2);

    if( pCQuad0 == NULL )
        return fmu_kCLEAN_NO_ACTION;
    if( pCQuad1 == NULL )
        return fmu_kCLEAN_NO_ACTION;
    if( pCQuad2 == NULL )
        return fmu_kCLEAN_NO_ACTION;

    fmuEdge *pCEdge01 = pCNode0->shared_edge ( pCNode1 );
    if( pCEdge01 == NULL)
        return fmu_kCLEAN_NO_ACTION;
    if( pCEdge01->frozen() )
        return fmu_kCLEAN_NO_ACTION;
    elmQClass *pCQuadA = pCEdge01->other_quad( pCQuad0 );
    if( pCQuadA == NULL)
        return fmu_kCLEAN_NO_ACTION;
    nodClass *pCNodeA = pCQuadA->opposite_node ( pCNode0 );
    if( pCNodeA == NULL)
        return fmu_kCLEAN_NO_ACTION;

    fmuEdge *pCEdge1A = pCNode1->shared_edge ( pCNodeA );
    if( pCEdge1A == NULL)
        return fmu_kCLEAN_NO_ACTION;
    if( pCEdge1A->frozen() )
        return fmu_kCLEAN_NO_ACTION;
    elmQClass *pCQuadB = pCEdge1A->other_quad( pCQuadA );
    if( pCQuadB == NULL)
        return fmu_kCLEAN_NO_ACTION;
    nodClass *pCNodeB = pCQuadB->opposite_node ( pCNode1 );
    if( pCNodeB == NULL)
        return fmu_kCLEAN_NO_ACTION;

    fmuEdge *pCEdge12 = pCNode1->shared_edge ( pCNode2 );
    if( pCEdge12 == NULL)
        return fmu_kCLEAN_NO_ACTION;
    fmuEdge *pCEdge05 = pCNode0->shared_edge ( pCNode5 );
    if( pCEdge05 == NULL)
        return fmu_kCLEAN_NO_ACTION;
    if( pCEdge05->frozen() )
        return fmu_kCLEAN_NO_ACTION;
    if( pCEdge12->frozen() )
        return fmu_kCLEAN_NO_ACTION;

    int iPattern = 0;
    fmuVector CVec45( pCNode4, pCNode5 );
    fmuVector CVec43( pCNode4, pCNode3 );
    if( CNormal.angle( CVec45, CVec43 ) > fmu_kSHAPE_ANGLE_TOL )
        iPattern = 1;
    else
    {
        fmuVector CVecAB( pCNodeA, pCNodeB );
        fmuVector CVecA5( pCNodeA, pCNode5 );
        if( CNormal.angle( CVecAB, CVecA5 ) > fmu_kSHAPE_ANGLE_TOL )
            iPattern = 1;
    }
    if( iPattern == 1 )
    {
        fmuVector CVec5A( pCNode5, pCNodeA );
        if( CNormal.angle( CVec5A, -CVec45 ) > fmu_kSHAPE_ANGLE_TOL )
            return fmu_kCLEAN_NO_ACTION;
    }

    if( lDrawCleaning )
    {
        //draw_open();
        // Refresh debug graphics
        //EraseElements();
        pCQuad0->draw_me( dp_kCLR_YELLOW );
        pCQuad1->draw_me( dp_kCLR_YELLOW );
        pCQuad2->draw_me( dp_kCLR_YELLOW );
        pCQuadA->draw_me( dp_kCLR_YELLOW );
        pCQuadB->draw_me( dp_kCLR_YELLOW );
        //draw_close();
    }

    delete_this_quad( pCQuad0 ); pCQuad0 = NULL;
    delete_this_quad( pCQuad1 ); pCQuad1 = NULL;
    delete_this_quad( pCQuad2 ); pCQuad2 = NULL;
    delete_this_quad( pCQuadA ); pCQuadA = NULL;
    delete_this_quad( pCQuadB ); pCQuadB = NULL;

    if( delete_this_edge( pCEdge0 ) == fmu_kCLEAN_ABORT )
        return fmu_kCLEAN_ABORT;
    if( delete_this_edge( pCEdge1 ) == fmu_kCLEAN_ABORT )
        return fmu_kCLEAN_ABORT;
    if( delete_this_edge( pCEdge2 ) == fmu_kCLEAN_ABORT )
        return fmu_kCLEAN_ABORT;
    if( delete_this_edge( pCEdge01 ) == fmu_kCLEAN_ABORT )
        return fmu_kCLEAN_ABORT;
    if( delete_this_edge( pCEdge05 ) == fmu_kCLEAN_ABORT )
        return fmu_kCLEAN_ABORT;
    if( delete_this_edge( pCEdge12 ) == fmu_kCLEAN_ABORT )
        return fmu_kCLEAN_ABORT;
    if( delete_this_edge( pCEdge1A ) == fmu_kCLEAN_ABORT )
        return fmu_kCLEAN_ABORT;

    if( delete_this_node( pCCenterNode ) == fmu_kCLEAN_ABORT )
        return fmu_kCLEAN_ABORT;
    if( delete_this_node( pCNode1 ) == fmu_kCLEAN_ABORT )
        return fmu_kCLEAN_ABORT;
    if( iPattern == 0 )
    {
        if( delete_this_node( pCNode0 ) == fmu_kCLEAN_ABORT )
            return fmu_kCLEAN_ABORT;
    }

    //Build 2 or 3 new quads
    elmQClass *pCQuadC = 0;
    if( iPattern == 0 )
    {
        pCQuadA = new elmQClass( pCNode2, pCNode3, pCNode4, pCNode5, pCMaster );
        if( !pCQuadA )
            return fmu_kCLEAN_NO_MEMORY;
        pCQuadB = new elmQClass( pCNode5, pCNodeA, pCNodeB, pCNode2, pCMaster );
        if( !pCQuadB )
            return fmu_kCLEAN_NO_MEMORY;
    }
    else
    {
        pCQuadA = new elmQClass( pCNode2, pCNode3, pCNode4, pCNode0, pCMaster );
        if( !pCQuadA )
            return fmu_kCLEAN_NO_MEMORY;
        pCQuadB = new elmQClass( pCNode0, pCNodeA, pCNodeB, pCNode2, pCMaster );
        if( !pCQuadB )
            return fmu_kCLEAN_NO_MEMORY;
        pCQuadC = new elmQClass( pCNode4, pCNode5, pCNodeA, pCNode0, pCMaster );
        if( !pCQuadC )
            return fmu_kCLEAN_NO_MEMORY;
    }

    if( lDrawCleaning )
    {
        // Refresh debug graphics
        //EraseElements();
        nodClassDynArray CSmoothList;
        CSmoothList.append( pCNode2 );
        CSmoothList.append( pCNode3 );
        CSmoothList.append( pCNode4 );
        CSmoothList.append( pCNode5 );
        CSmoothList.append( pCNodeA );
        CSmoothList.append( pCNodeB );
        smooth_mesh( &CSmoothList );
        //draw_mesh();
        pCQuadA->draw_me( dp_kCLR_BLUE );
        pCQuadB->draw_me( dp_kCLR_BLUE );
        if( iPattern == 1 ) pCQuadB->draw_me( dp_kCLR_BLUE );
        //draw_close();
    }

    reclean( pCNode2 );
    reclean( pCNode3 );
    reclean( pCNode4 );
    reclean( pCNode5 );
    reclean( pCNodeA );
    reclean( pCNodeB );
    if( iPattern == 1 )
        reclean( pCNode0 );

    return fmu_kCLEAN_MESH_CHANGED;
}


int QuadCleanTool::squeeze_double_row(nodClass *pCCenterNode,
                                          nodClassDynArray &CNodes,
                                          fmuEdgeDynArray &CEdges,
                                          elmQClassDynArray &CQuads )
{
    //This is a mesh size adjustment hiding as a valence case. This case
    //is 3-444444. Check around the central node for some elements with
    //bad aspect ratios. This might be done in either direction as the
    //elements may be alongside a boundary. If the elements would have
    //better aspect ratios if 3 rows were reduced to one row, start squeezing
    //out elements. Continue until a non-4-valent node interferes or the
    //aspect ratio needs no more correcting.

    int iDebugFlag = 0;
    bool lProceeding = true;
    bool lMeshChanged = false;
    double dRatio, dAngle;
    double dRATIO_TOL = 1.60;
    nodClassDynArray CSmoothList;
    nodClass *pCNode0 = 0, *pCNode1 = 0, *pCNode2 = 0, *pCNode3 = 0,
        *pCNode4 = 0, *pCNode5 = 0;
    fmuEdge  *pCEdge0 = 0, *pCEdge1 = 0, *pCEdge2 = 0;
    elmQClass *pCQuad0 = 0, *pCQuad1 = 0, *pCQuad2 = 0;

    //If the edges are frozen, no action regardless of orientation
    if( CEdges.get()->frozen() ) 
        return fmu_kCLEAN_NO_ACTION;
    if( CEdges.next()->frozen() ) 
        return fmu_kCLEAN_NO_ACTION;
    if( CEdges.next(2)->frozen() ) 
        return fmu_kCLEAN_NO_ACTION;

    for( int ii = 0; ii < 3; ++ii )
    {
        //The original orientation is irrelevant. This will check them all anyway
        CNodes.step(2);
        CEdges.step();
        CQuads.step();
        pCNode0 = CNodes.get();
        pCNode1 = CNodes.next();
        pCNode2 = CNodes.next(2);
        pCNode3 = CNodes.next(3);
        pCNode4 = CNodes.next(4);
        pCNode5 = CNodes.next(5);
        pCEdge0 = CEdges.get();
        pCEdge1 = CEdges.next();
        pCEdge2 = CEdges.next(2);
        pCQuad0 = CQuads.get();
        pCQuad1 = CQuads.next();
        pCQuad2 = CQuads.next(2);

        if (	pCNode0 == NULL ||
            pCNode1 == NULL ||
            pCNode2 == NULL ||
            pCNode3 == NULL ||
            pCNode4 == NULL ||
            pCNode5 == NULL ||
            pCEdge0 == NULL ||
            pCEdge1 == NULL ||
            pCEdge2 == NULL ||
            pCQuad0 == NULL ||
            pCQuad1 == NULL ||
            pCQuad2 == NULL)
            return fmu_kCLEAN_NO_ACTION;

        while( lProceeding )
        {
            lProceeding = false;
            fmuEdge *pCEdge12 = pCNode1->shared_edge( pCNode2 );
            if (pCEdge12 == NULL) 
                return fmu_kCLEAN_NO_ACTION;

            fmuEdge *pCEdge23 = pCNode2->shared_edge( pCNode3 );
            if (pCEdge23 == NULL) 
                return fmu_kCLEAN_NO_ACTION;

            if( !pCEdge12->frozen() && !pCEdge23->frozen() &&
                !pCNode2->hardpt() && !pCNode1->hardpt() && !pCNode0->hardpt() )
            {
                elmQClass *pCQuadA = pCEdge12->other_quad( pCQuad0 );
                if (pCQuadA == NULL) 
                    return fmu_kCLEAN_NO_ACTION;

                nodClass *pCNodeA = pCQuadA->next_node( pCNode1 );
                if (pCNodeA == NULL) 
                    return fmu_kCLEAN_NO_ACTION;

                fmuVector CVec5A( pCNode5, pCNodeA );
                fmuVector CVec54( pCNode5, pCNode4 );
                fmuVector CVec12( pCNode1, pCNode2 );
                //Use the average of vector 54 & vector 12.
                dRatio = ( CVec5A.length() * 2.0 ) / ( CVec54.length() + CVec12.length() );
                if( dRatio < dRATIO_TOL )
                {
                    elmQClass *pCQuadB = pCEdge23->other_quad( pCQuad1 );
                    if (pCQuadB == NULL) 
                        return fmu_kCLEAN_NO_ACTION;

                    nodClass *pCNodeB = pCQuadB->next_node( pCNode2 );
                    if (pCNodeB == NULL)
                        return fmu_kCLEAN_NO_ACTION;

                    nodClass *pCNodeC = pCQuadB->opposite_node( pCNode2 );
                    if (pCNodeC == NULL) 
                        return fmu_kCLEAN_NO_ACTION;

                    fmuEdge *pCEdge2B = pCNode2->shared_edge( pCNodeB );
                    if (pCEdge2B == NULL) 
                        return fmu_kCLEAN_NO_ACTION;

                    dAngle = 0.0;
                    if( pCNode4->hardpt() && pCNode3->hardpt() &&
                        pCNodeC->hardpt() )
                    {
                        fmuVector CVec34( pCNode3, pCNode4 );
                        fmuVector CVec3C( pCNode3, pCNodeC );
                        dAngle = CNormal.angle( CVec34, CVec3C );
                    }
                    if( !pCEdge2B->frozen() && dAngle < fmu_kSHAPE_ANGLE_TOL )
                    {
                        if( lDrawCleaning )
                        {
                            //draw_open();
                            // Refresh debug graphics
                            //EraseElements();
                            pCQuad0->draw_me( dp_kCLR_YELLOW );
                            pCQuad1->draw_me( dp_kCLR_YELLOW );
                            pCQuad2->draw_me( dp_kCLR_YELLOW );
                            pCQuadA->draw_me( dp_kCLR_YELLOW );
                            pCQuadB->draw_me( dp_kCLR_YELLOW );
                            //draw_close();
                        }
                        delete_this_quad( pCQuad0 ); pCQuad0 = NULL;
                        delete_this_quad( pCQuad1 ); pCQuad1 = NULL;
                        delete_this_quad( pCQuad2 ); pCQuad2 = NULL;
                        delete_this_quad( pCQuadA ); pCQuadA = NULL;
                        delete_this_quad( pCQuadB ); pCQuadB = NULL;

                        if( delete_this_edge( pCEdge0 ) == fmu_kCLEAN_ABORT )
                            return fmu_kCLEAN_ABORT;
                        if( delete_this_edge( pCEdge1 ) == fmu_kCLEAN_ABORT )
                            return fmu_kCLEAN_ABORT;
                        if( delete_this_edge( pCEdge2 ) == fmu_kCLEAN_ABORT )
                            return fmu_kCLEAN_ABORT;
                        if( delete_this_edge( pCEdge12 ) == fmu_kCLEAN_ABORT )
                            return fmu_kCLEAN_ABORT;
                        if( delete_this_edge( pCEdge23 ) == fmu_kCLEAN_ABORT )
                            return fmu_kCLEAN_ABORT;
                        if( delete_this_edge( pCEdge2B ) == fmu_kCLEAN_ABORT )
                            return fmu_kCLEAN_ABORT;

                        if( delete_this_node( pCCenterNode ) == fmu_kCLEAN_ABORT )
                            return fmu_kCLEAN_ABORT;
                        if( delete_this_node( pCNode2 ) == fmu_kCLEAN_ABORT )
                            return fmu_kCLEAN_ABORT;

                        //Three new quads, assigned to variables to be ready for
                        //the next round of squeezing
                        pCQuad1 = new elmQClass( pCNode0, pCNode1, pCNode4, pCNode5, pCMaster );
                        if( !pCQuad1 )
                            return fmu_kCLEAN_NO_MEMORY;
                        pCQuad0 = new elmQClass( pCNode1, pCNodeA, pCNodeB, pCNode4, pCMaster );
                        if( !pCQuad0 )
                            return fmu_kCLEAN_NO_MEMORY;
                        pCQuad2 = new elmQClass( pCNode4, pCNodeB, pCNodeC, pCNode3, pCMaster );
                        if( !pCQuad2 )
                            return fmu_kCLEAN_NO_MEMORY;

                        CSmoothList.clear();
                        CSmoothList.append( pCNode0 );
                        CSmoothList.append( pCNode1 );
                        CSmoothList.append( pCNode3 );
                        CSmoothList.append( pCNode4 );
                        CSmoothList.append( pCNode5 );
                        CSmoothList.append( pCNodeA );
                        CSmoothList.append( pCNodeB );
                        CSmoothList.append( pCNodeC );
                        smooth_mesh( &CSmoothList );
                        if( lDrawCleaning )
                        {
                            //draw_mesh();
                            // Refresh debug graphics
                            //EraseElements();
                            pCQuad0->draw_me( dp_kCLR_BLUE);
                            pCQuad1->draw_me( dp_kCLR_BLUE );
                            pCQuad2->draw_me( dp_kCLR_BLUE );
                            //draw_close();
                        }
                        reclean( pCNodeC );
                        reclean( pCNode4 );
                        lMeshChanged = true;

                        //Rearrange nodes, edges, quads so the squeezing might
                        //be propagated down the row

                        pCNode3 = pCNode4;
                        pCNode4 = pCNode5;
                        pCCenterNode = pCNode0;
                        pCNode2 = pCNode1;
                        pCEdge1 = pCCenterNode->shared_edge( pCNode2 );
                        pCEdge2 = pCCenterNode->shared_edge( pCNode4 );

                        // Check for NULL pointers - NM
                        if (pCEdge1 == NULL || pCEdge2 == NULL)
                            return fmu_kCLEAN_ABORT;

                        if( pCEdge1->frozen() || pCEdge2->frozen() )
                            return fmu_kCLEAN_MESH_CHANGED;

                        pCQuad0 = pCEdge1->other_quad( pCQuad1 );
                        if (pCQuad0 == NULL)
                            return fmu_kCLEAN_ABORT;

                        pCNode0 = pCQuad0->next_node( pCCenterNode );
                        if (pCNode0 == NULL)
                            return fmu_kCLEAN_ABORT;

                        pCEdge0 = pCCenterNode->shared_edge( pCNode0 );
                        if (pCEdge0 == NULL)
                            return fmu_kCLEAN_ABORT;

                        pCNode1 = pCQuad0->opposite_node( pCCenterNode );
                        if (pCNode1 == NULL)
                            return fmu_kCLEAN_ABORT;

                        pCQuad2 = pCEdge2->other_quad( pCQuad1 );
                        if (pCQuad2 == NULL)
                            return fmu_kCLEAN_ABORT;

                        pCNode5 = pCQuad2->opposite_node( pCCenterNode );
                        if (pCNode5 == NULL)
                            return fmu_kCLEAN_ABORT;

                        if( get_valence( pCNode0 ) != 4 ||
                            get_valence( pCNode1 ) != 4 ||
                            get_valence( pCNode5 ) != 4 )
                            return fmu_kCLEAN_MESH_CHANGED;
                        lProceeding = true;
                        if( iDebugFlag )
                        {
                            //draw_open();
                            // Refresh debug graphics
                            //EraseElements();
                            pCQuad0->draw_me( dp_kCLR_BLUE );
                            pCQuad1->draw_me( dp_kCLR_BLUE );
                            pCQuad2->draw_me( dp_kCLR_BLUE );
                            //draw_close();
                        }
                    }
                }
            }
        }
        if( lMeshChanged )
            return fmu_kCLEAN_MESH_CHANGED;

        lProceeding = true;
        while( lProceeding )
        {
            lProceeding = false;
            fmuEdge *pCEdge45 = pCNode4->shared_edge( pCNode5 );
            fmuEdge *pCEdge34 = pCNode3->shared_edge( pCNode4 );
            // Check if the objects are valid - NM
            if (pCEdge45 == NULL || pCEdge34 == NULL)
                return fmu_kCLEAN_ABORT;

            if( !pCEdge45->frozen() && !pCEdge34->frozen()
                && !pCNode4->hardpt() && !pCNode5->hardpt() && !pCNode0->hardpt() )
            {
                elmQClass *pCQuadA = pCEdge45->other_quad( pCQuad2 );
                if (pCQuadA == NULL)
                    return fmu_kCLEAN_ABORT;

                nodClass *pCNodeA = pCQuadA->prev_node( pCNode5 );
                if (pCNodeA == NULL)
                    return fmu_kCLEAN_ABORT;

                fmuVector CVec1A( pCNode1, pCNodeA );
                fmuVector CVec12( pCNode1, pCNode2 );
                fmuVector CVec54( pCNode5, pCNode4 );
                //Use the average of vector 54 & vector 12.
                dRatio = ( CVec1A.length() * 2.0 ) / ( CVec54.length() + CVec12.length() );
                if( dRatio < dRATIO_TOL )
                {
                    elmQClass *pCQuadB = pCEdge34->other_quad( pCQuad1 );
                    if (pCQuadB == NULL)
                        return fmu_kCLEAN_ABORT;

                    nodClass *pCNodeB = pCQuadB->prev_node( pCNode4 );
                    if (pCNodeB == NULL)
                        return fmu_kCLEAN_ABORT;

                    nodClass *pCNodeC = pCQuadB->opposite_node( pCNode4 );
                    if (pCNodeC == NULL)
                        return fmu_kCLEAN_ABORT;

                    fmuEdge *pCEdge4B = pCNode4->shared_edge( pCNodeB );
                    if (pCEdge4B == NULL)
                        return fmu_kCLEAN_ABORT;

                    dAngle = 0.0;
                    if( pCNodeC->hardpt() && 
                        pCNode3->hardpt() &&
                        pCNode2->hardpt() )
                    {
                        fmuVector CVec3C( pCNode3, pCNodeC );
                        fmuVector CVec32( pCNode3, pCNode2 );
                        dAngle = CNormal.angle( CVec3C, CVec32 );
                    }
                    if( !pCEdge4B->frozen() && dAngle < fmu_kSHAPE_ANGLE_TOL )
                    {
                        if( lDrawCleaning )
                        {
                            //draw_open();
                            // Refresh debug graphics
                            //EraseElements();
                            pCQuad0->draw_me( dp_kCLR_YELLOW );
                            pCQuad1->draw_me( dp_kCLR_YELLOW );
                            pCQuad2->draw_me( dp_kCLR_YELLOW );
                            pCQuadA->draw_me( dp_kCLR_YELLOW );
                            pCQuadB->draw_me( dp_kCLR_YELLOW );
                            //draw_close();
                        }
                        delete_this_quad( pCQuad0 ); pCQuad0 = NULL;
                        delete_this_quad( pCQuad1 ); pCQuad1 = NULL;
                        delete_this_quad( pCQuad2 ); pCQuad2 = NULL;
                        delete_this_quad( pCQuadA ); pCQuadA = NULL;
                        delete_this_quad( pCQuadB ); pCQuadB = NULL;

                        if( delete_this_edge( pCEdge0 ) == fmu_kCLEAN_ABORT )
                            return fmu_kCLEAN_ABORT;
                        if( delete_this_edge( pCEdge1 ) == fmu_kCLEAN_ABORT )
                            return fmu_kCLEAN_ABORT;
                        if( delete_this_edge( pCEdge2 ) == fmu_kCLEAN_ABORT )
                            return fmu_kCLEAN_ABORT;
                        if( delete_this_edge( pCEdge34 ) == fmu_kCLEAN_ABORT )
                            return fmu_kCLEAN_ABORT;
                        if( delete_this_edge( pCEdge45 ) == fmu_kCLEAN_ABORT )
                            return fmu_kCLEAN_ABORT;
                        if( delete_this_edge( pCEdge4B ) == fmu_kCLEAN_ABORT )
                            return fmu_kCLEAN_ABORT;

                        if( delete_this_node( pCCenterNode ) == fmu_kCLEAN_ABORT )
                            return fmu_kCLEAN_ABORT;
                        if( delete_this_node( pCNode4 ) == fmu_kCLEAN_ABORT )
                            return fmu_kCLEAN_ABORT;

                        pCQuad1 = new elmQClass( pCNode0, pCNode1, pCNode2, pCNode5, pCMaster );
                        if( !pCQuad1 )
                            return fmu_kCLEAN_NO_MEMORY;
                        pCQuad0 = new elmQClass( pCNodeA, pCNode5, pCNode2, pCNodeB, pCMaster );
                        if( !pCQuad0 )
                            return fmu_kCLEAN_NO_MEMORY;
                        pCQuad2 = new elmQClass( pCNodeB, pCNode2, pCNode3, pCNodeC, pCMaster );
                        if( !pCQuad2 )
                            return fmu_kCLEAN_NO_MEMORY;

                        CSmoothList.clear();
                        CSmoothList.append( pCNode0 );
                        CSmoothList.append( pCNode1 );
                        CSmoothList.append( pCNode2 );
                        CSmoothList.append( pCNode3 );
                        CSmoothList.append( pCNode5 );
                        CSmoothList.append( pCNodeA );
                        CSmoothList.append( pCNodeB );
                        CSmoothList.append( pCNodeC );
                        smooth_mesh( &CSmoothList );
                        if( lDrawCleaning )
                        {
                            //draw_mesh();
                            // Refresh debug graphics
                            //EraseElements();
                            //pCQuad0->draw_me( dp_kCLR_PINK );
                            //pCQuad1->draw_me( dp_kCLR_GOLD );
                            //pCQuad2->draw_me( dp_kCLR_MAGNTA );
                            pCQuad0->draw_me( dp_kCLR_BLUE );
                            pCQuad1->draw_me( dp_kCLR_BLUE );
                            pCQuad2->draw_me( dp_kCLR_BLUE );
                           
                            //draw_close();
                        }
                        reclean( pCNodeC );
                        reclean( pCNode2 );
                        lMeshChanged = true;

                        //Rearrange nodes, edges, quads so the squeezing might
                        //be propagated down the row

                        pCNode3 = pCNode2;
                        pCNode2 = pCNode1;
                        pCCenterNode = pCNode0;
                        pCNode4 = pCNode5;
                        pCEdge1 = pCCenterNode->shared_edge( pCNode2 );
                        pCEdge2 = pCCenterNode->shared_edge( pCNode4 );
                        // Check for NULL pointers - NM
                        if (pCEdge1 == NULL || pCEdge2 == NULL) 
                            return fmu_kCLEAN_ABORT;

                        if( pCEdge1->frozen() || pCEdge2->frozen() )
                            return fmu_kCLEAN_MESH_CHANGED;
                        pCQuad0 = pCEdge1->other_quad( pCQuad1 );
                        if (pCQuad0 == NULL)
                            return fmu_kCLEAN_ABORT;
                        pCNode0 = pCQuad0->next_node( pCCenterNode );
                        if (pCNode0 == NULL)
                            return fmu_kCLEAN_ABORT;
                        pCEdge0 = pCCenterNode->shared_edge( pCNode0 );
                        if (pCEdge0 == NULL)
                            return fmu_kCLEAN_ABORT;
                        pCNode1 = pCQuad0->opposite_node( pCCenterNode );
                        if (pCNode1 == NULL)
                            return fmu_kCLEAN_ABORT;
                        pCQuad2 = pCEdge2->other_quad( pCQuad1 );
                        if (pCQuad2 == NULL)
                            return fmu_kCLEAN_ABORT;

                        pCNode5 = pCQuad2->opposite_node( pCCenterNode );
                        if (pCNode5 == NULL)
                            return fmu_kCLEAN_ABORT;

                        if( get_valence( pCNode0 ) != 4 ||
                            get_valence( pCNode1 ) != 4 ||
                            get_valence( pCNode5 ) != 4 )
                            return fmu_kCLEAN_MESH_CHANGED;
                        lProceeding = true;
                        if( iDebugFlag )
                        {
                            //draw_open();
                            // Refresh debug graphics
                            //EraseElements();
                            //pCQuad0->draw_me( dp_kCLR_BLUE );
                            //pCQuad1->draw_me( dp_kCLR_RED );
                            //pCQuad2->draw_me( dp_kCLR_BLUE );
                            pCQuad0->draw_me( dp_kCLR_BLUE );
                            pCQuad1->draw_me( dp_kCLR_BLUE );
                            pCQuad2->draw_me( dp_kCLR_BLUE );
                          
                            //draw_close();
                        }
                    }
                }
            }
        }
        if( lMeshChanged )
            return fmu_kCLEAN_MESH_CHANGED;
    }

    return fmu_kCLEAN_NO_ACTION;
}


int QuadCleanTool::clean_4333( nodClass *pCCenterNode,
    nodClassDynArray &CNodes,
    fmuEdgeDynArray &CEdges,
    elmQClassDynArray &CQuads )
{

    int iCallStatus;
    nodClassDynArray CDeleteNodes(10);
    nodClassDynArray CRingNodes(10);
    fmuEdgeDynArray CDeleteEdges(10);
    elmQClassDynArray CDeleteQuads(10);
    nodClass *pCNode0 = 0;
    nodClass *pCNode1 = 0;
    nodClass *pCNode2 = 0;
    nodClass *pCNode3 = 0;
    nodClass *pCNode4 = 0;
    nodClass *pCNode6 = 0;
    nodClass *pCNode7 = 0;
    fmuEdge  *pCEdge0 = 0;
    fmuEdge  *pCEdge1 = 0;
    elmQClass *pCQuad0 = 0;
    elmQClass *pCQuad1 = 0;
    elmQClass *pCQuad3 = 0;

    pCNode0 = CNodes.get();
    pCNode1 = CNodes.next();
    pCNode2 = CNodes.next(2);
    pCNode3 = CNodes.next(3);
    pCNode4 = CNodes.next(4);
    pCNode6 = CNodes.next(6);
    pCNode7 = CNodes.next(7);
    pCEdge0 = CEdges.get();
    pCEdge1 = CEdges.next();
    pCQuad0 = CQuads.get();
    pCQuad1 = CQuads.next();
    pCQuad3 = CQuads.next(3);

    if (	pCNode0 == NULL ||
        pCNode1 == NULL ||
        pCNode2 == NULL ||
        pCNode3 == NULL ||
        pCNode4 == NULL ||
        pCNode6 == NULL ||
        pCNode7 == NULL ||
        pCEdge0 == NULL ||
        pCEdge1 == NULL ||
        pCQuad0 == NULL ||
        pCQuad1 == NULL ||
        pCQuad3 == NULL)
        return fmu_kCLEAN_NO_ACTION;

    if( pCEdge0->frozen() )
        return fmu_kCLEAN_NO_ACTION;
    if( pCEdge1->frozen() ) 
        return fmu_kCLEAN_NO_ACTION;

    fmuEdge *pCEdge01 = pCNode0->shared_edge( pCNode1 );
    if( pCEdge01== NULL)
        return fmu_kCLEAN_NO_ACTION;

    if( pCEdge01->frozen() ) 
        return fmu_kCLEAN_NO_ACTION;

    fmuEdge *pCEdge12 = pCNode1->shared_edge( pCNode2 );
    if( pCEdge12== NULL) 
        return fmu_kCLEAN_NO_ACTION;

    if( pCEdge12->frozen() ) 
        return fmu_kCLEAN_NO_ACTION;

    fmuEdge *pCEdge23 = pCNode2->shared_edge( pCNode3 );
    if( pCEdge23== NULL) 
        return fmu_kCLEAN_NO_ACTION;
    if( pCEdge23->frozen() ) 
        return fmu_kCLEAN_NO_ACTION;

    fmuEdge *pCEdge07 = pCNode0->shared_edge( pCNode7 );
    if( pCEdge07== NULL) 
        return fmu_kCLEAN_NO_ACTION;
    if( pCEdge07->frozen() ) 
        return fmu_kCLEAN_NO_ACTION;

    elmQClass *pCQuadA = pCEdge01->other_quad( pCQuad0 );
    if( pCQuadA== NULL) 
        return fmu_kCLEAN_NO_ACTION;

    nodClass *pCNodeA = pCQuadA->opposite_node( pCNode0 );
    if( pCNodeA== NULL) 
        return fmu_kCLEAN_NO_ACTION;

    fmuEdge *pCEdge1A = pCNode1->shared_edge( pCNodeA );
    if( pCEdge1A == NULL) 
        return fmu_kCLEAN_NO_ACTION;
    if( pCEdge1A->frozen() ) 
        return fmu_kCLEAN_NO_ACTION;

    elmQClass *pCQuadB = pCEdge12->other_quad( pCQuad0 );
    if( pCQuadB== NULL) 
        return fmu_kCLEAN_NO_ACTION;

    //Quads to delete
    CDeleteQuads[0] = pCQuad0;
    CDeleteQuads[1] = pCQuad1;
    CDeleteQuads[2] = pCQuad3;
    CDeleteQuads[3] = pCQuadA;
    CDeleteQuads[4] = pCQuadB;
    //Edges to delete
    CDeleteEdges[0] = pCEdge0;
    CDeleteEdges[1] = pCEdge1;
    CDeleteEdges[2] = pCEdge01;
    CDeleteEdges[3] = pCEdge12;
    CDeleteEdges[4] = pCEdge23;
    CDeleteEdges[5] = pCEdge07;
    CDeleteEdges[6] = pCEdge1A;
    //Nodes to delete
    CDeleteNodes[0] = pCNode0;
    CDeleteNodes[1] = pCNode1;
    CDeleteNodes[2] = pCNode2;
    //Ring nodes
    CRingNodes[0] = pCCenterNode;
    CRingNodes[1] = pCNode6;
    CRingNodes[2] = pCNode7;
    CRingNodes[3] = pCNodeA;
    CRingNodes[4] = pCNode3;
    CRingNodes[5] = pCNode4;

    iCallStatus = fill_two( CDeleteNodes, CDeleteEdges, CDeleteQuads, CRingNodes );

    return iCallStatus;
}


int QuadCleanTool::clean_435435( nodClass *pCCenterNode,
                                     nodClassDynArray &CNodes,
                                     fmuEdgeDynArray &CEdges,
                                     elmQClassDynArray &CQuads ) 
{
   
   int ii, nCount;
   nodClass *pCNode0 = 0;
   nodClass *pCNode1 = 0;
   nodClass *pCNode2 = 0;
   nodClass *pCNode3 = 0;
   nodClass *pCNode4 = 0;
   nodClass *pCNode5 = 0;
   nodClass *pCNode6 = 0;
   nodClass *pCNode7 = 0;

   for( ii = 0; ii < CEdges.length(); ++ii )
   {
	  if (CEdges[ii] == NULL) return fmu_kCLEAN_NO_ACTION;
      if( CEdges[ii]->frozen() ) return fmu_kCLEAN_NO_ACTION;
   }

   pCNode0 = CNodes.get();
   pCNode1 = CNodes.next();
   pCNode2 = CNodes.next(2);
   pCNode3 = CNodes.next(3);
   pCNode4 = CNodes.next(4);
   pCNode5 = CNodes.next(5);
   pCNode6 = CNodes.next(6);
   pCNode7 = CNodes.next(7);

   // Check to see if all nodes exist - NM
	  if (	pCNode0 == NULL ||
			pCNode1 == NULL ||
			pCNode2 == NULL ||
			pCNode3 == NULL ||
			pCNode4 == NULL ||
			pCNode5 == NULL ||
			pCNode6 == NULL ||
			pCNode7 == NULL)
			return fmu_kCLEAN_NO_ACTION;


   if( pCNode4->hardpt() && pCNode6->hardpt() ) return fmu_kCLEAN_NO_ACTION;

   //Find the quad on the other side of the 3 valent node
   fmuEdge  *pCEdge01 = 0;
   pCEdge01 = pCNode0->shared_edge( pCNode1 );
   if( pCEdge01 == NULL ) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge01->frozen() ) return fmu_kCLEAN_NO_ACTION;

   fmuEdge  *pCEdge07 = 0;
   pCEdge07 = pCNode0->shared_edge( pCNode7 );
	if( pCEdge07 == NULL ) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge07->frozen() ) return fmu_kCLEAN_NO_ACTION;

   elmQClass *pCOuterQuad = 0;
   pCOuterQuad = pCEdge01->shared_quad( pCEdge07 );
	if( pCOuterQuad == NULL ) return fmu_kCLEAN_NO_ACTION;

   nodClass *pCOuterNode = 0;
   pCOuterNode = pCOuterQuad->opposite_node( pCNode0 );
	if( pCOuterNode == NULL ) return fmu_kCLEAN_NO_ACTION;

   if( pCOuterNode->hardpt() &&
       pCNode1->hardpt() &&
       pCNode2->hardpt() )
   {
      fmuVector CVec_12( pCNode1, pCNode2 );
      fmuVector CVec_1O( pCNode1, pCOuterNode );
      if( CNormal.angle( CVec_12, CVec_1O ) >= fmu_kSHAPE_ANGLE_TOL )
        return fmu_kCLEAN_NO_ACTION;
   }

   if( pCNode3->hardpt() &&
       pCNode4->hardpt() &&
       pCNode5->hardpt() )
   {
      fmuVector CVec_45( pCNode4, pCNode5 );
      fmuVector CVec_43( pCNode4, pCNode3 );
      if( CNormal.angle( CVec_45, CVec_43 ) >= fmu_kSHAPE_ANGLE_TOL )
        return fmu_kCLEAN_NO_ACTION;
   }

   if( lDrawCleaning )
   {
   //draw_open();
       // Refresh debug graphics
        //EraseElements();
        CQuads.get()->draw_me( dp_kCLR_YELLOW );
        CQuads.next()->draw_me( dp_kCLR_YELLOW );
        CQuads.next(2)->draw_me( dp_kCLR_YELLOW );
        CQuads.next(3)->draw_me( dp_kCLR_YELLOW );
        pCOuterQuad->draw_me( dp_kCLR_YELLOW );
   //draw_close();
   }

   //Delete all the quads and edges around the center node along
   //with the outer quad found above.
   nCount = CQuads.length();
   for( ii = 0; ii < nCount; ++ii )
   {
      delete_this_quad( CQuads[ii] );
   }
   delete_this_quad( pCOuterQuad );
   nCount = CEdges.length();
   for( ii = 0; ii < nCount; ++ii )
   {
      if( delete_this_edge( CEdges[ii] ) == fmu_kCLEAN_ABORT )
         return fmu_kCLEAN_ABORT;
   }
   if( delete_this_edge( pCEdge01 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge07 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;

   if( delete_this_node( pCCenterNode ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_node( pCNode0 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;

   //Rebuild: 3 new quads
   elmQClass *pCQuad1 = new elmQClass( pCOuterNode, pCNode1, pCNode2,
                                                    pCNode7, pCMaster );
   if( !pCQuad1 ) return fmu_kCLEAN_NO_MEMORY;
   elmQClass *pCQuad2 = new elmQClass( pCNode7, pCNode2, pCNode3, pCNode6,
                                                                pCMaster );
   if( !pCQuad2 ) return fmu_kCLEAN_NO_MEMORY;
   elmQClass *pCQuad3 = new elmQClass( pCNode6, pCNode3, pCNode4, pCNode5,
                                                                pCMaster );
   if( !pCQuad3 ) return fmu_kCLEAN_NO_MEMORY;

   if( lDrawCleaning )
   {
      nodClassDynArray CSmoothList;
      CSmoothList.append( pCNode1 );
      CSmoothList.append( pCNode2 );
      CSmoothList.append( pCNode3 );
      CSmoothList.append( pCNode4 );
      CSmoothList.append( pCNode5 );
      CSmoothList.append( pCNode6 );
      CSmoothList.append( pCNode7 );
      CSmoothList.append( pCOuterNode );
      smooth_mesh( &CSmoothList );
      //draw_mesh();
      //pCQuad1->draw_me( dp_kCLR_BLUE );
      //pCQuad2->draw_me( dp_kCLR_YELLOW );
      //pCQuad3->draw_me( dp_kCLR_CYAN );
           // Refresh debug graphics
            //EraseElements();
            pCQuad1->draw_me( dp_kCLR_BLUE );
            pCQuad2->draw_me( dp_kCLR_BLUE );
            pCQuad3->draw_me( dp_kCLR_BLUE );
                           
      //draw_close();
   }

   //The node in the current loop position has been deleted. All the
   //rest need to be cleaned.
   for( ii = 1; ii < 8; ++ii ) reclean( CNodes.next(ii) );
   reclean( pCOuterNode );

   return fmu_kCLEAN_MESH_CHANGED;   
   
}


int QuadCleanTool::clean_453453( nodClass *pCCenterNode,
                                     nodClassDynArray &CNodes,
                                     fmuEdgeDynArray &CEdges,
                                     elmQClassDynArray &CQuads )
{
   
   int ii, nCount;
   nodClass *pCNode0 = 0;
   nodClass *pCNode1 = 0;
   nodClass *pCNode2 = 0;
   nodClass *pCNode3 = 0;
   nodClass *pCNode4 = 0;
   nodClass *pCNode5 = 0;
   nodClass *pCNode6 = 0;
   nodClass *pCNode7 = 0;

   for( ii = 0; ii < CEdges.length(); ++ii )
   {
	  if( CEdges[ii] == NULL ) return fmu_kCLEAN_NO_ACTION;
      if( CEdges[ii]->frozen() ) return fmu_kCLEAN_NO_ACTION;
   }

   pCNode0 = CNodes.get();
   pCNode1 = CNodes.next();
   pCNode2 = CNodes.next(2);
   pCNode3 = CNodes.next(3);
   pCNode4 = CNodes.next(4);
   pCNode5 = CNodes.next(5);
   pCNode6 = CNodes.next(6);
   pCNode7 = CNodes.next(7);

   // Check to see if all nodes exist - NM
	  if (	pCNode0 == NULL ||
			pCNode1 == NULL ||
			pCNode2 == NULL ||
			pCNode3 == NULL ||
			pCNode4 == NULL ||
			pCNode5 == NULL ||
			pCNode6 == NULL ||
			pCNode7 == NULL)
			return fmu_kCLEAN_NO_ACTION;


   if( pCNode0->hardpt() && pCNode6->hardpt() ) return fmu_kCLEAN_NO_ACTION;

   //Find the quad on the other side of the 3 valent node
   fmuEdge  *pCEdge43 = 0;
   pCEdge43 = pCNode4->shared_edge( pCNode3 );
   if( pCEdge43 == NULL ) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge43->frozen() ) return fmu_kCLEAN_NO_ACTION;

   fmuEdge  *pCEdge45 = 0;
   pCEdge45 = pCNode4->shared_edge( pCNode5 );
   if( pCEdge45 == NULL ) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge45->frozen() ) return fmu_kCLEAN_NO_ACTION;

   elmQClass *pCOuterQuad = 0;
   pCOuterQuad = pCEdge43->shared_quad( pCEdge45 );
   if( pCOuterQuad == NULL ) return fmu_kCLEAN_NO_ACTION;

   nodClass *pCOuterNode = 0;
   pCOuterNode = pCOuterQuad->opposite_node( pCNode4 );
   if( pCOuterNode == NULL ) return fmu_kCLEAN_NO_ACTION;

   if( pCNode2->hardpt() &&
       pCNode3->hardpt() &&
       pCOuterNode->hardpt() )
   {
      fmuVector CVec_3O( pCNode3, pCOuterNode );
      fmuVector CVec_32( pCNode3, pCNode2 );
      if( CNormal.angle( CVec_3O, CVec_32 ) >= fmu_kSHAPE_ANGLE_TOL )
        return fmu_kCLEAN_NO_ACTION;
   }

   if( pCNode7->hardpt() &&
       pCNode0->hardpt() &&
       pCNode1->hardpt() )
   {
      fmuVector CVec_01( pCNode0, pCNode1 );
      fmuVector CVec_07( pCNode0, pCNode7 );
      if( CNormal.angle( CVec_01, CVec_07 ) >= fmu_kSHAPE_ANGLE_TOL )
        return fmu_kCLEAN_NO_ACTION;
   }

   if( lDrawCleaning )
   {
      //draw_open();
        CQuads.get()->draw_me( dp_kCLR_YELLOW );
        CQuads.next()->draw_me( dp_kCLR_YELLOW );
        CQuads.next(2)->draw_me( dp_kCLR_YELLOW );
        CQuads.next(3)->draw_me( dp_kCLR_YELLOW );
        pCOuterQuad->draw_me( dp_kCLR_YELLOW );
   //draw_close();
   }

   //Delete all the CQuads and edges around the center node along
   //with the outer quad found above.
   nCount = CQuads.length();
   for( ii = 0; ii < nCount; ++ii )
   {
      delete_this_quad( CQuads[ii] );
   }
   delete_this_quad( pCOuterQuad );
   nCount = CEdges.length();
   for( ii = 0; ii < nCount; ++ii )
   {
      if( delete_this_edge( CEdges[ii] ) == fmu_kCLEAN_ABORT )
         return fmu_kCLEAN_ABORT;
   }
   if( delete_this_edge( pCEdge43 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge45 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_node( pCCenterNode ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_node( pCNode4 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;

   //Rebuild: 3 new CQuads
   elmQClass *pCQuad1 = new elmQClass( pCNode1, pCNode6, pCNode7, pCNode0,
                                                                 pCMaster );
   if( !pCQuad1 ) return fmu_kCLEAN_NO_MEMORY;
   elmQClass *pCQuad2 = new elmQClass( pCNode2, pCNode5, pCNode6, pCNode1,
                                                                 pCMaster );
   if( !pCQuad2 ) return fmu_kCLEAN_NO_MEMORY;
   elmQClass *pCQuad3 = new elmQClass( pCNode3, pCOuterNode, pCNode5, pCNode2,
                                                                 pCMaster );
   if( !pCQuad3 ) return fmu_kCLEAN_NO_MEMORY;

   if( lDrawCleaning )
   {
      nodClassDynArray CSmoothList;
      CSmoothList.append( pCNode0 );
      CSmoothList.append( pCNode1 );
      CSmoothList.append( pCNode2 );
      CSmoothList.append( pCNode3 );
      CSmoothList.append( pCNode5 );
      CSmoothList.append( pCNode6 );
      CSmoothList.append( pCNode7 );
      CSmoothList.append( pCOuterNode );
      smooth_mesh( &CSmoothList );
      //draw_mesh();
      //pCQuad1->draw_me( dp_kCLR_BLUE );
      //pCQuad2->draw_me( dp_kCLR_GOLD );
      //pCQuad3->draw_me( dp_kCLR_LTBLUE );
            // Refresh debug graphics
            //EraseElements();
            pCQuad1->draw_me( dp_kCLR_BLUE );
            pCQuad2->draw_me( dp_kCLR_BLUE );
            pCQuad3->draw_me( dp_kCLR_BLUE );
                           
      //draw_close();
   }

   //The node in the next(4) position of the loop has been deleted. All the
   //rest need to be cleaned.
   CNodes.step(4);
   for( ii = 1; ii < 8; ++ii ) reclean( CNodes.next(ii) );
   reclean( pCOuterNode );

   return fmu_kCLEAN_MESH_CHANGED;   
   
}


int QuadCleanTool::double_triangle( nodClass *pCCenterNode,
                                        nodClassDynArray &CNodes,
                                        fmuEdgeDynArray &CEdges,
                                        elmQClassDynArray &CQuads,
                                        int aiValence[8],
                                        int aiHardpt[8] )
{
  
   int ii, nCount;
   nodClass *middle_node = NULL;
   elmQClass *pCQuad1 = 0, *pCQuad2 = 0;

   for( ii = 0; ii < CEdges.length(); ++ii )
   {
	  if( CEdges[ii] == NULL ) return fmu_kCLEAN_NO_ACTION;
      if( CEdges[ii]->frozen() ) return fmu_kCLEAN_NO_ACTION;
   }

   nodClass *pCNode0 = 0;
   pCNode0 = CNodes.get();
   nodClass *pCNode1 = 0;
   pCNode1 = CNodes.next();
   nodClass *pCNode2 = 0;
   pCNode2 = CNodes.next(2);
   nodClass *pCNode3 = 0;
   pCNode3 = CNodes.next(3);
   nodClass *pCNode4 = 0;
   pCNode4 = CNodes.next(4);
   nodClass *pCNode5 = 0;
   pCNode5 = CNodes.next(5);
   nodClass *pCNode6 = 0;
   pCNode6 = CNodes.next(6);
   nodClass *pCNode7 = 0;
   pCNode7 = CNodes.next(7);

   
   // Check to see if all nodes exist - NM
	  if (	pCNode0 == NULL ||
			pCNode1 == NULL ||
			pCNode2 == NULL ||
			pCNode3 == NULL ||
			pCNode4 == NULL ||
			pCNode5 == NULL ||
			pCNode6 == NULL ||
			pCNode7 == NULL)
			return fmu_kCLEAN_NO_ACTION;

   //Check for special combinations where this action won't work.
   //Prevent triangular shaped quad on the border
   if( aiValence[2] == 4 && aiHardpt[2] ) return fmu_kCLEAN_NO_ACTION;
   if( aiValence[6] == 4 && aiHardpt[6] ) return fmu_kCLEAN_NO_ACTION;
   //Prevent stretched CQuads
   if( aiHardpt[0] && ( aiHardpt[3] || aiHardpt[5] ) )
      return fmu_kCLEAN_NO_ACTION;
   if( aiHardpt[4] && ( aiHardpt[1] || aiHardpt[7] ) )
      return fmu_kCLEAN_NO_ACTION;
   //Check angles to avoid flattened angles
   fmuVector CVec23( pCNode2, pCNode3 );
   fmuVector CVec21( pCNode2, pCNode1 );
   if( CNormal.angle( CVec23, CVec21 ) > fmu_kSHAPE_ANGLE_TOL )
      return fmu_kCLEAN_NO_ACTION;
   fmuVector CVec67( pCNode6, pCNode7 );
   fmuVector CVec65( pCNode6, pCNode5 );
   if( CNormal.angle( CVec67, CVec65 ) > fmu_kSHAPE_ANGLE_TOL )
      return fmu_kCLEAN_NO_ACTION;

   //Only one interior node can be a hardpoint, if none are,
   // keep pCCenterNode.
   if( aiHardpt[0] )  middle_node = pCNode0;
   if( aiHardpt[4] )
   {
      if( middle_node ) return fmu_kCLEAN_NO_ACTION;
      else              middle_node = pCNode4;
   }
   if( pCCenterNode->hardpt() && middle_node ) return fmu_kCLEAN_NO_ACTION;
   if( !middle_node ) middle_node = pCCenterNode;

   //Find the CQuads on the other side of the 3 valent node
   fmuEdge  *pCEdge01 = pCNode0->shared_edge( pCNode1 );
   if( pCEdge01->frozen() ) return fmu_kCLEAN_NO_ACTION;
   fmuEdge  *pCEdge07 = pCNode0->shared_edge( pCNode7 );
   if( pCEdge07->frozen() ) return fmu_kCLEAN_NO_ACTION;
   elmQClass *pCQuadA = pCEdge01->shared_quad( pCEdge07 );
   nodClass *pCNodeA = pCQuadA->opposite_node( pCNode0 );

   fmuEdge  *pCEdge43 = pCNode4->shared_edge( pCNode3 );
   if( pCEdge43->frozen() ) return fmu_kCLEAN_NO_ACTION;
   fmuEdge  *pCEdge45 = pCNode4->shared_edge( pCNode5 );
   if( pCEdge45->frozen() ) return fmu_kCLEAN_NO_ACTION;
   elmQClass *pCQuadB = pCEdge43->shared_quad( pCEdge45 );
   nodClass *pCNodeB = pCQuadB->opposite_node( pCNode4 );

   if( lDrawCleaning )
   {
   //draw_open();
        // Refresh debug graphics
        //EraseElements();
        CQuads.get()->draw_me( dp_kCLR_YELLOW );
        CQuads.next()->draw_me( dp_kCLR_YELLOW );
        CQuads.next(2)->draw_me( dp_kCLR_YELLOW );
        CQuads.next(3)->draw_me( dp_kCLR_YELLOW );
        pCQuadA->draw_me( dp_kCLR_YELLOW );
        pCQuadB->draw_me( dp_kCLR_YELLOW );
   //draw_close();
   }

   //Delete all the CQuads and edges around the center node along
   //with the outer quad found above.
   nCount = CQuads.length();
   for( ii = 0; ii < nCount; ++ii )
   {
      delete_this_quad( CQuads[ii] );
   }
   delete_this_quad( pCQuadA );
   delete_this_quad( pCQuadB );

   nCount = CEdges.length();
   for( ii = 0; ii < nCount; ++ii )
   {
      if( delete_this_edge( CEdges[ii] ) == fmu_kCLEAN_ABORT )
         return fmu_kCLEAN_ABORT;
   }
   if( delete_this_edge( pCEdge01 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge07 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge43 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge45 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;

   if( middle_node != pCCenterNode )
      if( delete_this_node( pCCenterNode ) == fmu_kCLEAN_ABORT )
         return fmu_kCLEAN_ABORT;
   if( middle_node != pCNode0 )
      if( delete_this_node( pCNode0 ) == fmu_kCLEAN_ABORT )
         return fmu_kCLEAN_ABORT;
   if( middle_node != pCNode4 )
      if( delete_this_node( pCNode4 ) == fmu_kCLEAN_ABORT )
         return fmu_kCLEAN_ABORT;

   //Build 4 new CQuads
   pCQuadA = new elmQClass( middle_node, pCNode7, pCNodeA, pCNode1,
                                                                pCMaster );
   if( !pCQuadA ) return fmu_kCLEAN_NO_MEMORY;
   pCQuadB = new elmQClass( middle_node, pCNode1, pCNode2, pCNode3,
                                                                pCMaster );
   if( !pCQuadB ) return fmu_kCLEAN_NO_MEMORY;
   pCQuad1 = new elmQClass( middle_node, pCNode3, pCNodeB, pCNode5,
                                                                pCMaster );
   if( !pCQuad1 ) return fmu_kCLEAN_NO_MEMORY;
   pCQuad2 = new elmQClass( middle_node, pCNode5, pCNode6, pCNode7,
                                                                pCMaster );
   if( !pCQuad2 ) return fmu_kCLEAN_NO_MEMORY;

   if( lDrawCleaning )
   {
      nodClassDynArray CSmoothList;
      CSmoothList.append( pCNodeA );
      CSmoothList.append( pCNodeB );
      CSmoothList.append( pCNode1 );
      CSmoothList.append( pCNode2 );
      CSmoothList.append( pCNode3 );
      CSmoothList.append( pCNode5 );
      CSmoothList.append( pCNode6 );
      CSmoothList.append( pCNode7 );
      CSmoothList.append( middle_node );
      smooth_mesh( &CSmoothList );
      //draw_mesh();
      //pCQuad1->draw_me( dp_kCLR_BLUE );
      //pCQuad2->draw_me( dp_kCLR_BLUE );
      //pCQuadA->draw_me( dp_kCLR_GOLD );
      //pCQuadB->draw_me( dp_kCLR_RED );
          // Refresh debug graphics
          //EraseElements();
          pCQuad1->draw_me( dp_kCLR_BLUE );
          pCQuad2->draw_me( dp_kCLR_BLUE );
          pCQuadA->draw_me( dp_kCLR_BLUE );
          pCQuadB->draw_me( dp_kCLR_BLUE );
      //draw_close();
   }

   reclean( pCNode1 );
   reclean( pCNode2 );
   reclean( pCNode3 );
   reclean( pCNode5 );
   reclean( pCNode6 );
   reclean( pCNode7 );

   return fmu_kCLEAN_MESH_CHANGED;   
}


int QuadCleanTool::clean_434435445( nodClass *pCCenterNode,
                                        nodClassDynArray &CNodes,
                                        fmuEdgeDynArray &CEdges,
                                        elmQClassDynArray &CQuads )
{
  
   nodClass *pCNode0 = 0;
   pCNode0 = CNodes.get();
   nodClass *pCNode1 = 0;
   pCNode1 = CNodes.next();
   nodClass *pCNode2 = 0;
   pCNode2 = CNodes.next(2);
   nodClass *pCNode3 = 0;
   pCNode3 = CNodes.next(3);
   nodClass *pCNode4 = 0;
   pCNode4 = CNodes.next(4);
   nodClass *pCNode5 = 0;
   pCNode5 = CNodes.next(5);
   nodClass *pCNode6 = 0;
   pCNode6 = CNodes.next(6);
   nodClass *pCNode7 = 0;
   pCNode7 = CNodes.next(7);

   
   // Check to see if all nodes exist - NM
	  if (	pCNode0 == NULL ||
			pCNode1 == NULL ||
			pCNode2 == NULL ||
			pCNode3 == NULL ||
			pCNode4 == NULL ||
			pCNode5 == NULL ||
			pCNode6 == NULL ||
			pCNode7 == NULL)
			return fmu_kCLEAN_NO_ACTION;

//   nodClass *pCNode0 = CNodes.get();
//   nodClass *pCNode1 = CNodes.next();
//   nodClass *pCNode2 = CNodes.next(2);
//   nodClass *pCNode3 = CNodes.next(3);
//   nodClass *pCNode4 = CNodes.next(4);
//   nodClass *pCNode5 = CNodes.next(5);
//   nodClass *pCNode6 = CNodes.next(6);
//   nodClass *pCNode7 = CNodes.next(7);

   fmuEdge  *pCEdge0 = 0;
   pCEdge0 = CEdges.get();
   fmuEdge  *pCEdge1 = 0;
   pCEdge1 = CEdges.next();
   fmuEdge  *pCEdge2 = 0;
   pCEdge2 = CEdges.next(2);
   fmuEdge  *pCEdge3 = 0;
   pCEdge3 = CEdges.next(3);

   // Null pointer check - NM
   if( !pCEdge0 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCEdge1 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCEdge2 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCEdge3 ) return fmu_kCLEAN_NO_ACTION;


   elmQClass *pCQuad0 = 0;
   pCQuad0 = CQuads.get();
   elmQClass *pCQuad1 = 0;
   pCQuad1 = CQuads.next();
   elmQClass *pCQuad2 = 0;
   pCQuad2 = CQuads.next(2);
   elmQClass *pCQuad3 = 0;
   pCQuad3 = CQuads.next(3);

   // Null pointer check - NM
   if( !pCQuad0 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCQuad1 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCQuad2 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCQuad3 ) return fmu_kCLEAN_NO_ACTION;

   if( pCEdge0->frozen() ) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge1->frozen() ) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge2->frozen() ) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge3->frozen() ) return fmu_kCLEAN_NO_ACTION;

   fmuEdge *pCEdge01 = pCNode0->shared_edge( pCNode1 );
   if( pCEdge01->frozen() ) return fmu_kCLEAN_NO_ACTION;
   fmuEdge *pCEdge07 = pCNode0->shared_edge( pCNode7 );
   if( pCEdge07->frozen() ) return fmu_kCLEAN_NO_ACTION;

   elmQClass *pCQuadA = 0;
   pCQuadA = pCEdge01->other_quad( pCQuad0 );
   nodClass *pCNodeA = 0;
   pCNodeA = pCQuadA->opposite_node( pCNode0 );

   if( lDrawCleaning )
   {
   //draw_open();
        //pCQuadA->draw_me( dp_kCLR_BLUE );
   //pCQuad3->draw_me( dp_kCLR_RED );
   //pCQuad2->draw_me( dp_kCLR_YELLOW );
   //pCQuad1->draw_me( dp_kCLR_MAGNTA );
   //pCQuad0->draw_me( dp_kCLR_BLUE );
        // Refresh debug graphics
        //EraseElements();
        pCQuadA->draw_me( dp_kCLR_YELLOW );
        pCQuad3->draw_me( dp_kCLR_YELLOW );
        pCQuad2->draw_me( dp_kCLR_YELLOW );
        pCQuad1->draw_me( dp_kCLR_YELLOW );
        pCQuad0->draw_me( dp_kCLR_YELLOW );
   //draw_close();
   }

   //Have 5 quads, 6 edges, 2 nodes to delete
   delete_this_quad( pCQuad0 );
   delete_this_quad( pCQuad1 );
   delete_this_quad( pCQuad2 );
   delete_this_quad( pCQuad3 );
   delete_this_quad( pCQuadA );

   if( delete_this_edge( pCEdge0 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge1 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge2 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge3 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge01 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge07 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;

   if( delete_this_node( pCCenterNode ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_node( pCNode0 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;

   //Now make 3 new quads
   pCQuad0 = new elmQClass( pCNode6, pCNode7, pCNodeA, pCNode1, pCMaster );
   if( !pCQuad0 ) return fmu_kCLEAN_NO_MEMORY;
   pCQuad1 = new elmQClass( pCNode6, pCNode1, pCNode2, pCNode3, pCMaster );
   if( !pCQuad1 ) return fmu_kCLEAN_NO_MEMORY;
   pCQuad2 = new elmQClass( pCNode5, pCNode6, pCNode3, pCNode4, pCMaster );
   if( !pCQuad2 ) return fmu_kCLEAN_NO_MEMORY;

   if( lDrawCleaning )
   {
      nodClassDynArray CSmoothList;
      CSmoothList.append( pCNodeA );
      CSmoothList.append( pCNode5 );
      CSmoothList.append( pCNode6 );
      CSmoothList.append( pCNode7 );
      CSmoothList.append( pCNode1 );
      CSmoothList.append( pCNode4 );
      CSmoothList.append( pCNode3 );
      CSmoothList.append( pCNode2 );
      smooth_mesh( &CSmoothList );
      //draw_mesh();
      //pCQuad2->draw_me( dp_kCLR_YELLOW );
      //pCQuad1->draw_me( dp_kCLR_BLUE );
      //pCQuad0->draw_me( dp_kCLR_BLUE );
      // Refresh debug graphics
        //EraseElements();
        pCQuad2->draw_me( dp_kCLR_BLUE );
        pCQuad1->draw_me( dp_kCLR_BLUE );
        pCQuad0->draw_me( dp_kCLR_BLUE );
      //draw_close();
   }

   reclean( pCNode1 );
   reclean( pCNode4 );
   reclean( pCNode5 );
   reclean( pCNode6 );
   reclean( pCNode7 );
   reclean( pCNode2 );
   reclean( pCNode3 );

   return fmu_kCLEAN_MESH_CHANGED;   
}


int QuadCleanTool::clean_453443544( nodClass *pCCenterNode,
                                        nodClassDynArray &CNodes,
                                        fmuEdgeDynArray &CEdges,
                                        elmQClassDynArray &CQuads )
{
  
   nodClass *pCNode0 = CNodes.get();
   nodClass *pCNode1 = CNodes.next();
   nodClass *pCNode2 = CNodes.next(2);
   nodClass *pCNode3 = CNodes.next(3);
   nodClass *pCNode4 = CNodes.next(4);
   nodClass *pCNode5 = CNodes.next(5);
   nodClass *pCNode6 = CNodes.next(6);
   nodClass *pCNode7 = CNodes.next(7);
   fmuEdge  *pCEdge0 = CEdges.get();
   fmuEdge  *pCEdge1 = CEdges.next();
   fmuEdge  *pCEdge2 = CEdges.next(2);
   fmuEdge  *pCEdge3 = CEdges.next(3);
   elmQClass *pCQuad0 = CQuads.get();
   elmQClass *pCQuad1 = CQuads.next();
   elmQClass *pCQuad2 = CQuads.next(2);
   elmQClass *pCQuad3 = CQuads.next(3);

   // Check to see if all nodes exist - NM
	  if (	pCNode0 == NULL ||
			pCNode1 == NULL ||
			pCNode2 == NULL ||
			pCNode3 == NULL ||
			pCNode4 == NULL ||
			pCNode5 == NULL ||
			pCNode6 == NULL ||
			pCNode7 == NULL)
			return fmu_kCLEAN_NO_ACTION;

   // Check to see if all nodes exist  - NM
   if( !pCEdge0 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCEdge1 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCEdge2 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCEdge3 ) return fmu_kCLEAN_NO_ACTION;

   // Check to see if all nodes exist - NM
   if( !pCQuad0 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCQuad1 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCQuad2 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCQuad3 ) return fmu_kCLEAN_NO_ACTION;

   if( pCEdge0->frozen() ) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge1->frozen() ) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge2->frozen() ) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge3->frozen() ) return fmu_kCLEAN_NO_ACTION;


   fmuEdge *pCEdge34 = 0;
   pCEdge34 = pCNode3->shared_edge( pCNode4 );
   if (!pCEdge34) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge34->frozen() ) return fmu_kCLEAN_NO_ACTION;

   fmuEdge *pCEdge45 = 0;
   pCEdge45 = pCNode4->shared_edge( pCNode5 );
   if (!pCEdge45) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge45->frozen() ) return fmu_kCLEAN_NO_ACTION;

   elmQClass *pCQuadA = 0;
   pCQuadA = pCEdge34->other_quad( pCQuad1 );
   nodClass *pCNodeA = 0;
   pCNodeA = pCQuadA->opposite_node( pCNode4 );

   if( lDrawCleaning )
   {
   //draw_open();
        //pCQuadA->draw_me( dp_kCLR_BLUE );
   //pCQuad3->draw_me( dp_kCLR_RED );
   //pCQuad2->draw_me( dp_kCLR_YELLOW );
   //pCQuad1->draw_me( dp_kCLR_LTMGNT );
   //pCQuad0->draw_me( dp_kCLR_BLUE );
       // Refresh debug graphics
        //EraseElements();
        pCQuadA->draw_me( dp_kCLR_YELLOW );
        pCQuad3->draw_me( dp_kCLR_YELLOW );
        pCQuad2->draw_me( dp_kCLR_YELLOW );
        pCQuad1->draw_me( dp_kCLR_YELLOW );
        pCQuad0->draw_me( dp_kCLR_YELLOW );
   //draw_close();
   }

   //Have 5 quads, 6 edges, 2 nodes to delete
   delete_this_quad( pCQuad0 );
   delete_this_quad( pCQuad1 );
   delete_this_quad( pCQuad2 );
   delete_this_quad( pCQuad3 );
   delete_this_quad( pCQuadA );

   if( delete_this_edge( pCEdge0 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge1 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge2 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge3 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge34 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge45 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;

   if( delete_this_node( pCCenterNode ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_node( pCNode4 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;

   //Now make 3 new quads
   pCQuad0 = new elmQClass( pCNode6, pCNode7, pCNode0, pCNode1, pCMaster );
   if( !pCQuad0 ) return fmu_kCLEAN_NO_MEMORY;
   pCQuad1 = new elmQClass( pCNode6, pCNode1, pCNode2, pCNode3, pCMaster );
   if( !pCQuad1 ) return fmu_kCLEAN_NO_MEMORY;
   pCQuad2 = new elmQClass( pCNode5, pCNode6, pCNode3, pCNodeA, pCMaster );
   if( !pCQuad2 ) return fmu_kCLEAN_NO_MEMORY;

   if( lDrawCleaning )
   {
      nodClassDynArray CSmoothList;
      CSmoothList.append( pCNodeA );
      CSmoothList.append( pCNode5 );
      CSmoothList.append( pCNode6 );
      CSmoothList.append( pCNode7 );
      CSmoothList.append( pCNode0 );
      CSmoothList.append( pCNode1 );
      CSmoothList.append( pCNode3 );
      CSmoothList.append( pCNode2 );
      smooth_mesh( &CSmoothList );
      //draw_mesh();
      //pCQuad2->draw_me( dp_kCLR_YELLOW );
      //pCQuad1->draw_me( dp_kCLR_BLUE );
      //pCQuad0->draw_me( dp_kCLR_BLUE );
      // Refresh debug graphics
        //EraseElements();
        pCQuad2->draw_me( dp_kCLR_BLUE );
        pCQuad1->draw_me( dp_kCLR_BLUE );
        pCQuad0->draw_me( dp_kCLR_BLUE );
      //draw_close();
   }

   reclean( pCNode0 );
   reclean( pCNode3 );
   reclean( pCNode5 );
   reclean( pCNode6 );
   reclean( pCNode7 );
   reclean( pCNode1 );
   reclean( pCNode2 );
   reclean( pCNodeA );

   return fmu_kCLEAN_MESH_CHANGED;   
}


int QuadCleanTool::clean_434434( nodClass *pCCenterNode,
                                     nodClassDynArray &CNodes,
                                     fmuEdgeDynArray &CEdges,
                                     elmQClassDynArray &CQuads,
                                     int aiValence[8] )
{

   int iPattern;
   nodClass *pCNode0 = CNodes.get();
   nodClass *pCNode1 = CNodes.next();
   nodClass *pCNode2 = CNodes.next(2);
   nodClass *pCNode3 = CNodes.next(3);
   nodClass *pCNode4 = CNodes.next(4);
   nodClass *pCNode5 = CNodes.next(5);
   nodClass *pCNode6 = CNodes.next(6);
   nodClass *pCNode7 = CNodes.next(7);
   fmuEdge  *pCEdge0 = CEdges.get();
   fmuEdge  *pCEdge1 = CEdges.next();
   fmuEdge  *pCEdge2 = CEdges.next(2);
   fmuEdge  *pCEdge3 = CEdges.next(3);
   elmQClass *pCQuad0 = CQuads.get();
   elmQClass *pCQuad1 = CQuads.next();
   elmQClass *pCQuad2 = CQuads.next(2);
   elmQClass *pCQuad3 = CQuads.next(3);

   // Check to see if all nodes exist - NM
	  if (	pCNode0 == NULL ||
			pCNode1 == NULL ||
			pCNode2 == NULL ||
			pCNode3 == NULL ||
			pCNode4 == NULL ||
			pCNode5 == NULL ||
			pCNode6 == NULL ||
			pCNode7 == NULL)
			return fmu_kCLEAN_NO_ACTION;

   // Check to see if all nodes exist  - NM
   if( !pCEdge0 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCEdge1 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCEdge2 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCEdge3 ) return fmu_kCLEAN_NO_ACTION;

   // Check to see if all nodes exist - NM
   if( !pCQuad0 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCQuad1 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCQuad2 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCQuad3 ) return fmu_kCLEAN_NO_ACTION;


   if( pCEdge0->frozen() ) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge1->frozen() ) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge2->frozen() ) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge3->frozen() ) return fmu_kCLEAN_NO_ACTION;

   //The pattern used to replace what is deleted depends where
   //the 5 valent nodes are. Whether the fix proceeds depends
   //on the valence of various other nodes.
   if(      aiValence[5] == 4 &&
            aiValence[6] == 5 &&
            aiValence[7] == 4 ) iPattern = 1;

   else if( aiValence[5] == 4 &&
            aiValence[6] == 4 &&
            aiValence[7] == 5 ) iPattern = 2;

   else if( aiValence[5] == 5 &&
            aiValence[6] == 4 &&
            aiValence[7] == 5 ) iPattern = 3;

   else if( aiValence[5] == 5 &&
            aiValence[6] == 4 &&
            aiValence[7] == 4 ) iPattern = 4;

   else if( aiValence[5] == 4 &&
            aiValence[6] == 4 &&
            aiValence[7] == 4 )
   {
      return fix_double_trans( pCCenterNode, CNodes, CEdges, CQuads,
                                                             aiValence[0] );
   }

   else return fmu_kCLEAN_NO_ACTION;

   fmuEdge *pCEdge01 = 0;
   pCEdge01 = pCNode0->shared_edge( pCNode1 );
   if (!pCEdge01) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge01->frozen() ) return fmu_kCLEAN_NO_ACTION;

   fmuEdge *pCEdge12 = 0;
   pCEdge12 = pCNode1->shared_edge( pCNode2 );
	if (!pCEdge12) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge12->frozen() ) return fmu_kCLEAN_NO_ACTION;

   fmuEdge *pCEdge23 = 0;
   pCEdge23 = pCNode2->shared_edge( pCNode3 );
	if (!pCEdge23) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge23->frozen() ) return fmu_kCLEAN_NO_ACTION;

   fmuEdge *pCEdge34 = 0;
   pCEdge34 = pCNode3->shared_edge( pCNode4 );
   if (!pCEdge34) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge34->frozen() ) return fmu_kCLEAN_NO_ACTION;

   fmuEdge *pCEdge07 = 0;
   pCEdge07 = pCNode0->shared_edge( pCNode7 );
	if (!pCEdge07) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge07->frozen() ) return fmu_kCLEAN_NO_ACTION;

   elmQClass *pCQuadA = 0;
   pCQuadA = pCEdge01->other_quad( pCQuad0 );
   nodClass *pCNodeA = 0;
   pCNodeA = pCQuadA->opposite_node( pCNode0 );

   elmQClass *pCQuadB = 0;
   pCQuadB = pCEdge12->other_quad( pCQuad0 );
   nodClass *pCNodeB = 0;
   if(!pCQuadB) return fmu_kCLEAN_NO_ACTION;
   pCNodeB = pCQuadB->opposite_node( pCNode2 );
   nodClass *pCNodeC = 0;
   pCNodeC = pCQuadB->prev_node( pCNode2 );
   fmuEdge *pCEdge2C = 0;
   pCEdge2C = pCNode2->shared_edge( pCNodeC );
   if( pCEdge2C->frozen() ) return fmu_kCLEAN_NO_ACTION;

   elmQClass *pCQuadC = pCEdge23->other_quad( pCQuad1 );
   nodClass *pCNodeD = pCQuadC->opposite_node( pCNode2 );
   fmuEdge *pCEdge3D = pCNode3->shared_edge( pCNodeD );
   if( pCEdge3D->frozen() ) return fmu_kCLEAN_NO_ACTION;

   elmQClass *pCQuadD = pCEdge34->other_quad( pCQuad1 );
   nodClass *pCNodeE = pCQuadD->opposite_node( pCNode3 );

   fmuEdge *pCEdge1A = 0, *pCEdge1B = 0, *pCEdge45 = 0, *pCEdge4E = 0;
   elmQClass *pCQuadE = 0, *pCQuadF = 0;
   nodClass *pCNodeF = 0, *pCNodeG = 0;
   if( iPattern == 1 )
   {
      //Node C must be 4 or 5 valent
      if( get_valence( pCNodeC ) < 4 ) return fmu_kCLEAN_NO_ACTION;
      //Check for possible angle problems
      if( pCNode5->hardpt() &&
          pCNode6->hardpt() &&
          pCNode7->hardpt() )
      {
         fmuVector CVec_67( pCNode6, pCNode7 );
         fmuVector CVec_65( pCNode6, pCNode5 );
         if( CNormal.angle( CVec_67, CVec_65 ) >= fmu_kSHAPE_ANGLE_TOL )
            return fmu_kCLEAN_NO_ACTION;
      }
      if( pCNodeB->hardpt() &&
          pCNodeC->hardpt() &&
          pCNodeD->hardpt() )
      {
         fmuVector CVec_CD( pCNodeC, pCNodeD );
         fmuVector CVec_CB( pCNodeC, pCNodeB );
         if( CNormal.angle( CVec_CD, CVec_CB ) >= fmu_kSHAPE_ANGLE_TOL )
            return fmu_kCLEAN_NO_ACTION;
      }
   }
   else if( iPattern == 2 )
   {
      //Node D must be 4 or 5 valent, node E must be 3 or 4 valent
      if( get_valence( pCNodeD ) < 4 ) return fmu_kCLEAN_NO_ACTION;
      if( get_valence( pCNodeE ) > 4 ) return fmu_kCLEAN_NO_ACTION;
      //Check for possible angle problems
      if( pCNode6->hardpt() &&
          pCNode7->hardpt() &&
          pCNodeA->hardpt() )
      {
         fmuVector CVec_7A( pCNode7, pCNodeA );
         fmuVector CVec_76( pCNode7, pCNode6 );
         if( CNormal.angle( CVec_7A, CVec_76 ) >= fmu_kSHAPE_ANGLE_TOL )
            return fmu_kCLEAN_NO_ACTION;
      }
      if( pCNodeC->hardpt() &&
          pCNodeD->hardpt() &&
          pCNodeE->hardpt() )
      {
         fmuVector CVec_DE( pCNodeD, pCNodeE );
         fmuVector CVec_DC( pCNodeD, pCNodeC );
         if( CNormal.angle( CVec_DE, CVec_DC ) >= fmu_kSHAPE_ANGLE_TOL )
            return fmu_kCLEAN_NO_ACTION;
      }
   }
   else if( iPattern == 3 )
   {
      //Node B must be 4 or 5 valent, use an alternate quad
      if( get_valence( pCNodeB ) < 4 ) return fmu_kCLEAN_NO_ACTION;
      pCEdge1A = pCNode1->shared_edge( pCNodeA );
		if(!pCEdge1A) return fmu_kCLEAN_NO_ACTION;
      if( pCEdge1A->frozen() ) return fmu_kCLEAN_NO_ACTION;
		
      pCEdge1B = pCNode1->shared_edge( pCNodeB );
		if(!pCEdge1B) return fmu_kCLEAN_NO_ACTION;
      if( pCEdge1B->frozen() ) return fmu_kCLEAN_NO_ACTION;

      pCQuadE = pCEdge1A->other_quad( pCQuadA );

		if(!pCQuadE) return fmu_kCLEAN_NO_ACTION;
      pCNodeF = pCQuadE->opposite_node( pCNode1 );
   }
   else if( iPattern == 4 )
   {
      //Node A must be 5 valent, node D must be 4 or 5 valent
      //Two additional quads in this pattern
      if( get_valence( pCNodeA ) != 5 ) return fmu_kCLEAN_NO_ACTION;
      if( get_valence( pCNodeD ) > 4 ) return fmu_kCLEAN_NO_ACTION;
      pCEdge1A = pCNode1->shared_edge( pCNodeA );

		if(!pCEdge1A) return fmu_kCLEAN_NO_ACTION;
      if( pCEdge1A->frozen() ) return fmu_kCLEAN_NO_ACTION;

      pCEdge1B = pCNode1->shared_edge( pCNodeB );
		if(!pCEdge1B) return fmu_kCLEAN_NO_ACTION;
      if( pCEdge1B->frozen() ) return fmu_kCLEAN_NO_ACTION;

      pCQuadE = pCEdge1A->other_quad( pCQuadA );

		if(!pCQuadE) return fmu_kCLEAN_NO_ACTION;
      pCNodeF = pCQuadE->opposite_node( pCNode1 );
      pCEdge45 = pCNode4->shared_edge( pCNode5 );
		if(!pCEdge45) return fmu_kCLEAN_NO_ACTION;
      if( pCEdge45->frozen() ) return fmu_kCLEAN_NO_ACTION;

      pCEdge4E = pCNode4->shared_edge( pCNodeE );
		if(!pCEdge4E) return fmu_kCLEAN_NO_ACTION;
      if( pCEdge4E->frozen() ) return fmu_kCLEAN_NO_ACTION;

		if(!pCEdge45) return fmu_kCLEAN_NO_ACTION;
      pCQuadF = pCEdge45->other_quad( pCQuad2 );
      pCNodeG = pCQuadF->opposite_node( pCNode4 );
   }

   if( lDrawCleaning )
   {
   //draw_open();
   //pCQuadD->draw_me( dp_kCLR_PINK );
   //pCQuadC->draw_me( dp_kCLR_GOLD );
   //pCQuadB->draw_me( dp_kCLR_LTBLUE );
        //pCQuadA->draw_me( dp_kCLR_BLUE );
   //pCQuad3->draw_me( dp_kCLR_RED );
   //pCQuad2->draw_me( dp_kCLR_YELLOW );
   //pCQuad1->draw_me( dp_kCLR_MAGNTA );
   //pCQuad0->draw_me( dp_kCLR_BLUE );
       // Refresh debug graphics
        //EraseElements(); 
        pCQuadD->draw_me( dp_kCLR_YELLOW );
        pCQuadC->draw_me( dp_kCLR_YELLOW );
        pCQuadB->draw_me( dp_kCLR_YELLOW );
        pCQuadA->draw_me( dp_kCLR_YELLOW );
        pCQuad3->draw_me( dp_kCLR_YELLOW );
        pCQuad2->draw_me( dp_kCLR_YELLOW );
        pCQuad1->draw_me( dp_kCLR_YELLOW );
        pCQuad0->draw_me( dp_kCLR_YELLOW );
   //if( iPattern == 3 ) pCQuadE->draw_me( dp_kCLR_LTMGNT );
        if( iPattern == 3 ) pCQuadE->draw_me( dp_kCLR_YELLOW );
        else if( iPattern == 4 )
        {
   //pCQuadE->draw_me( dp_kCLR_LTMGNT );
   //pCQuadF->draw_me( dp_kCLR_DKGRN );
            pCQuadE->draw_me( dp_kCLR_YELLOW );
            pCQuadF->draw_me( dp_kCLR_YELLOW );
        }
   //draw_close();
   }

   delete_this_quad( pCQuad0 );
   delete_this_quad( pCQuad1 );
   delete_this_quad( pCQuad2 );
   delete_this_quad( pCQuad3 );
   delete_this_quad( pCQuadA );
   delete_this_quad( pCQuadB );
   delete_this_quad( pCQuadC );
   if( iPattern < 3 )       delete_this_quad( pCQuadD );
   else if( iPattern == 3 ) delete_this_quad( pCQuadE );
   else if( iPattern == 4 )
   {
      delete_this_quad( pCQuadD );
      delete_this_quad( pCQuadE );
      delete_this_quad( pCQuadF );
   }

   if( delete_this_edge( pCEdge0 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge1 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge2 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge3 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge01 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge12 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge23 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge07 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge2C ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( iPattern < 3 )
   {
      if( delete_this_edge( pCEdge34 ) == fmu_kCLEAN_ABORT )
         return fmu_kCLEAN_ABORT;
      if( delete_this_edge( pCEdge3D ) == fmu_kCLEAN_ABORT )
         return fmu_kCLEAN_ABORT;
   }
   else if( iPattern == 3 )
   {
      if( delete_this_edge( pCEdge1A ) == fmu_kCLEAN_ABORT )
         return fmu_kCLEAN_ABORT;
      if( delete_this_edge( pCEdge1B ) == fmu_kCLEAN_ABORT )
         return fmu_kCLEAN_ABORT;
   }
   else if( iPattern == 4 )
   {
      if( delete_this_edge( pCEdge1A ) == fmu_kCLEAN_ABORT )
         return fmu_kCLEAN_ABORT;
      if( delete_this_edge( pCEdge1B ) == fmu_kCLEAN_ABORT )
         return fmu_kCLEAN_ABORT;
      if( delete_this_edge( pCEdge34 ) == fmu_kCLEAN_ABORT )
         return fmu_kCLEAN_ABORT;
      if( delete_this_edge( pCEdge3D ) == fmu_kCLEAN_ABORT )
         return fmu_kCLEAN_ABORT;
      if( delete_this_edge( pCEdge45 ) == fmu_kCLEAN_ABORT )
         return fmu_kCLEAN_ABORT;
      if( delete_this_edge( pCEdge4E ) == fmu_kCLEAN_ABORT )
         return fmu_kCLEAN_ABORT;
   }

   if( delete_this_node( pCCenterNode ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_node( pCNode0 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_node( pCNode2 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( iPattern < 3 )
   {
      if( delete_this_node( pCNode3 ) == fmu_kCLEAN_ABORT )
         return fmu_kCLEAN_ABORT;
   }
   else if( iPattern == 3 )
   {
      if( delete_this_node( pCNode1 ) == fmu_kCLEAN_ABORT )
         return fmu_kCLEAN_ABORT;
   }
   else if( iPattern == 4 )
   {
      if( delete_this_node( pCNode1 ) == fmu_kCLEAN_ABORT )
         return fmu_kCLEAN_ABORT;
      if( delete_this_node( pCNode3 ) == fmu_kCLEAN_ABORT )
         return fmu_kCLEAN_ABORT;
      if( delete_this_node( pCNode4 ) == fmu_kCLEAN_ABORT )
         return fmu_kCLEAN_ABORT;
   }

   //Now make 4 new quads, depending on the pattern
   if( iPattern == 1 )
   {
      pCQuad0 = new elmQClass( pCNode6, pCNode7, pCNode4, pCNode5, pCMaster );
      if( !pCQuad0 ) return fmu_kCLEAN_NO_MEMORY;
      pCQuad1 = new elmQClass( pCNode4, pCNode7, pCNodeA, pCNode1, pCMaster );
      if( !pCQuad1 ) return fmu_kCLEAN_NO_MEMORY;
      pCQuad2 = new elmQClass( pCNodeE, pCNode4, pCNode1, pCNodeD, pCMaster );
      if( !pCQuad2 ) return fmu_kCLEAN_NO_MEMORY;
      pCQuad3 = new elmQClass( pCNodeD, pCNode1, pCNodeB, pCNodeC, pCMaster );
      if( !pCQuad3 ) return fmu_kCLEAN_NO_MEMORY;
   }
   else if( iPattern == 2 )
   {
      pCQuad0 = new elmQClass( pCNode6, pCNode7, pCNodeA, pCNode1, pCMaster );
      if( !pCQuad0 ) return fmu_kCLEAN_NO_MEMORY;
      pCQuad1 = new elmQClass( pCNode4, pCNode5, pCNode6, pCNode1, pCMaster );
      if( !pCQuad1 ) return fmu_kCLEAN_NO_MEMORY;
      pCQuad2 = new elmQClass( pCNode4, pCNode1, pCNodeB, pCNodeC, pCMaster );
      if( !pCQuad2 ) return fmu_kCLEAN_NO_MEMORY;
      pCQuad3 = new elmQClass( pCNodeE, pCNode4, pCNodeC, pCNodeD, pCMaster );
      if( !pCQuad3 ) return fmu_kCLEAN_NO_MEMORY;
   }
   else if( iPattern == 3 )
   {
      pCQuad0 = new elmQClass( pCNode6, pCNode7, pCNodeA, pCNode3, pCMaster );
      if( !pCQuad0 ) return fmu_kCLEAN_NO_MEMORY;
      pCQuad1 = new elmQClass( pCNode4, pCNode5, pCNode6, pCNode3, pCMaster );
      if( !pCQuad1 ) return fmu_kCLEAN_NO_MEMORY;
      pCQuad2 = new elmQClass( pCNodeA, pCNodeF, pCNodeB, pCNodeC, pCMaster );
      if( !pCQuad2 ) return fmu_kCLEAN_NO_MEMORY;
      pCQuad3 = new elmQClass( pCNodeA, pCNodeC, pCNodeD, pCNode3, pCMaster );
      if( !pCQuad3 ) return fmu_kCLEAN_NO_MEMORY;
      pCNode1 = pCNode3;      //Make things easier below
   }
   else if( iPattern == 4 )
   {
      pCQuad0 = new elmQClass( pCNode7, pCNodeA, pCNodeF, pCNodeB, pCMaster );
      if( !pCQuad0 ) return fmu_kCLEAN_NO_MEMORY;
      pCQuad1 = new elmQClass( pCNode6, pCNode7, pCNodeB, pCNodeC, pCMaster );
      if( !pCQuad1 ) return fmu_kCLEAN_NO_MEMORY;
      pCQuad2 = new elmQClass( pCNode6, pCNodeC, pCNodeD, pCNodeE, pCMaster );
      if( !pCQuad2 ) return fmu_kCLEAN_NO_MEMORY;
      pCQuad3 = new elmQClass( pCNode5, pCNode6, pCNodeE, pCNodeG, pCMaster );
      if( !pCQuad3 ) return fmu_kCLEAN_NO_MEMORY;
      pCNode1 = pCNodeF;      //Make things easier below
      pCNode4 = pCNodeG;
   }


   if( lDrawCleaning )
   {
      nodClassDynArray CSmoothList;
      CSmoothList.append( pCNodeA );
      CSmoothList.append( pCNodeB );
      CSmoothList.append( pCNodeD );
      CSmoothList.append( pCNodeE );
      CSmoothList.append( pCNodeC );
      CSmoothList.append( pCNode5 );
      CSmoothList.append( pCNode6 );
      CSmoothList.append( pCNode7 );
      CSmoothList.append( pCNode1 );
      CSmoothList.append( pCNode4 );
      smooth_mesh( &CSmoothList );
      //draw_mesh();
      //pCQuad3->draw_me( dp_kCLR_RED );
      //pCQuad2->draw_me( dp_kCLR_YELLOW );
      //pCQuad1->draw_me( dp_kCLR_BLUE );
      //pCQuad0->draw_me( dp_kCLR_BLUE );
        // Refresh debug graphics
        //EraseElements();
        pCQuad3->draw_me( dp_kCLR_BLUE );
        pCQuad2->draw_me( dp_kCLR_BLUE );
        pCQuad1->draw_me( dp_kCLR_BLUE );
        pCQuad0->draw_me( dp_kCLR_BLUE );
      //draw_close();
   }

   reclean( pCNode1 );
   reclean( pCNode4 );
   reclean( pCNode5 );
   reclean( pCNode6 );
   reclean( pCNode7 );
   reclean( pCNodeB );
   reclean( pCNodeC );
   reclean( pCNodeD );

   return fmu_kCLEAN_MESH_CHANGED;   

}


int QuadCleanTool::clean_443443( nodClass *pCCenterNode,
                                     nodClassDynArray &CNodes,
                                     fmuEdgeDynArray &CEdges,
                                     elmQClassDynArray &CQuads,
                                     int aiValence[8] )
{

   int iPattern;
   nodClass *pCNode0 = CNodes.get();
   nodClass *pCNode1 = CNodes.next();
   nodClass *pCNode2 = CNodes.next(2);
   nodClass *pCNode3 = CNodes.next(3);
   nodClass *pCNode4 = CNodes.next(4);
   nodClass *pCNode5 = CNodes.next(5);
   nodClass *pCNode6 = CNodes.next(6);
   nodClass *pCNode7 = CNodes.next(7);
   fmuEdge  *pCEdge0 = CEdges.get();
   fmuEdge  *pCEdge1 = CEdges.next();
   fmuEdge  *pCEdge2 = CEdges.next(2);
   fmuEdge  *pCEdge3 = CEdges.next(3);
   elmQClass *pCQuad0 = CQuads.get();
   elmQClass *pCQuad1 = CQuads.next();
   elmQClass *pCQuad2 = CQuads.next(2);
   elmQClass *pCQuad3 = CQuads.next(3);

   // Check to see if all nodes exist - NM
	  if (	pCNode0 == NULL ||
			pCNode1 == NULL ||
			pCNode2 == NULL ||
			pCNode3 == NULL ||
			pCNode4 == NULL ||
			pCNode5 == NULL ||
			pCNode6 == NULL ||
			pCNode7 == NULL)
			return fmu_kCLEAN_NO_ACTION;

   // Check to see if all nodes exist  - NM
   if( !pCEdge0 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCEdge1 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCEdge2 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCEdge3 ) return fmu_kCLEAN_NO_ACTION;

   // Check to see if all nodes exist - NM
   if( !pCQuad0 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCQuad1 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCQuad2 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCQuad3 ) return fmu_kCLEAN_NO_ACTION;

   if( pCEdge0->frozen() ) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge1->frozen() ) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge2->frozen() ) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge3->frozen() ) return fmu_kCLEAN_NO_ACTION;

   //The pattern used to replace what is deleted depends where
   //the 5 valent nodes are. Whether the fix proceeds depends
   //on the valence of various other nodes.
   if(      aiValence[5] == 4 &&
            aiValence[6] == 5 &&
            aiValence[7] == 4 ) iPattern = 1;

   else if( aiValence[5] == 5 &&
            aiValence[6] == 4 &&
            aiValence[7] == 4 ) iPattern = 2;

   else if( aiValence[5] == 5 &&
            aiValence[6] == 4 &&
            aiValence[7] == 5 ) iPattern = 3;

   else if( aiValence[5] == 4 &&
            aiValence[6] == 4 &&
            aiValence[7] == 5 ) iPattern = 4;

   else if( aiValence[5] == 4 &&
            aiValence[6] == 4 &&
            aiValence[7] == 4 )
   {
      return fix_double_trans( pCCenterNode, CNodes, CEdges, CQuads,
                                                             aiValence[0] );
   }

   else return fmu_kCLEAN_NO_ACTION;

   fmuEdge *pCEdge01 = 0;
   pCEdge01 = pCNode0->shared_edge( pCNode1 );
	if(!pCEdge01) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge01->frozen() ) return fmu_kCLEAN_NO_ACTION;

   fmuEdge *pCEdge12 = 0;
   pCEdge12 = pCNode1->shared_edge( pCNode2 );
	if(!pCEdge12) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge12->frozen() ) return fmu_kCLEAN_NO_ACTION;

   fmuEdge *pCEdge23 = 0;
   pCEdge23 = pCNode2->shared_edge( pCNode3 );
	if(!pCEdge23) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge23->frozen() ) return fmu_kCLEAN_NO_ACTION;

   fmuEdge *pCEdge34 = 0;
   pCEdge34 = pCNode3->shared_edge( pCNode4 );
	if(!pCEdge34) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge34->frozen() ) return fmu_kCLEAN_NO_ACTION;

   fmuEdge *pCEdge45 = 0;
   pCEdge45 = pCNode4->shared_edge( pCNode5 );
	if(!pCEdge45) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge45->frozen() ) return fmu_kCLEAN_NO_ACTION;

   elmQClass *pCQuadA = pCEdge01->other_quad( pCQuad0 );
	if(!pCQuadA) return fmu_kCLEAN_NO_ACTION;
   nodClass *pCNodeA = pCQuadA->opposite_node( pCNode1 );
   nodClass *pCNodeB = pCQuadA->prev_node( pCNode1 );

   fmuEdge *pCEdge1B = pCNode1->shared_edge( pCNodeB );
   if( pCEdge1B->frozen() ) return fmu_kCLEAN_NO_ACTION;

   elmQClass *pCQuadB = pCEdge12->other_quad( pCQuad0 );
	if(!pCQuadB) return fmu_kCLEAN_NO_ACTION;
   nodClass *pCNodeC = pCQuadB->opposite_node( pCNode1 );
   fmuEdge *pCEdge2C = pCNode2->shared_edge( pCNodeC );
	if(!pCEdge2C) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge2C->frozen() ) return fmu_kCLEAN_NO_ACTION;


   elmQClass *pCQuadC = pCEdge23->other_quad( pCQuad1 );

	if(!pCQuadC) return fmu_kCLEAN_NO_ACTION;
   nodClass *pCNodeD = pCQuadC->opposite_node( pCNode2 );

   elmQClass *pCQuadD = pCEdge34->other_quad( pCQuad1 );
	if(!pCQuadD) return fmu_kCLEAN_NO_ACTION;
   nodClass *pCNodeE = pCQuadD->opposite_node( pCNode4 );

   fmuEdge *pCEdge3D = 0, *pCEdge3E = 0, *pCEdge07 = 0, *pCEdge0A = 0;
   elmQClass *pCQuadE = 0, *pCQuadF = 0;
   nodClass *pCNodeF = 0, *pCNodeG = 0;
   if( iPattern == 1 )
   {
      //Node C must be 4 or 5 valent
      if( get_valence( pCNodeC ) < 4 ) return fmu_kCLEAN_NO_ACTION;
      //Check for possible angle problems
      if( pCNode5->hardpt() &&
          pCNode6->hardpt() &&
          pCNode7->hardpt() )
      {
         fmuVector CVec_67( pCNode6, pCNode7 );
         fmuVector CVec_65( pCNode6, pCNode5 );
         if( CNormal.angle( CVec_67, CVec_65 ) >= fmu_kSHAPE_ANGLE_TOL )
            return fmu_kCLEAN_NO_ACTION;
      }
      if( pCNodeB->hardpt() &&
          pCNodeC->hardpt() &&
          pCNodeD->hardpt() )
      {
         fmuVector CVec_CD( pCNodeC, pCNodeD );
         fmuVector CVec_CB( pCNodeC, pCNodeB );
         if( CNormal.angle( CVec_CD, CVec_CB ) >= fmu_kSHAPE_ANGLE_TOL )
            return fmu_kCLEAN_NO_ACTION;
      }
   }
   else if( iPattern == 2 )
   {
      //Node B must be 4 or 5 valent, node A must be 3 or 4 valent
      if( get_valence( pCNodeB ) < 4 ) return fmu_kCLEAN_NO_ACTION;
      if( get_valence( pCNodeA ) > 4 ) return fmu_kCLEAN_NO_ACTION;
      //Check for possible angle problems
      if( pCNodeE->hardpt() &&
          pCNode5->hardpt() &&
          pCNode6->hardpt() )
      {
         fmuVector CVec_56( pCNode5, pCNode6 );
         fmuVector CVec_5E( pCNode5, pCNodeE );
         if( CNormal.angle( CVec_56, CVec_5E ) >= fmu_kSHAPE_ANGLE_TOL )
            return fmu_kCLEAN_NO_ACTION;
      }
      if( pCNodeA->hardpt() &&
          pCNodeB->hardpt() &&
          pCNodeC->hardpt() )
      {
         fmuVector CVec_BC( pCNodeB, pCNodeC );
         fmuVector CVec_BA( pCNodeB, pCNodeA );
         if( CNormal.angle( CVec_BC, CVec_BA ) >= fmu_kSHAPE_ANGLE_TOL )
            return fmu_kCLEAN_NO_ACTION;
      }
   }
   else if( iPattern == 3 )
   {
      //Node D must be 4 or 5 valent, use an alternate quad
      if( get_valence( pCNodeD ) < 4 ) return fmu_kCLEAN_NO_ACTION;
		
      pCEdge3D = pCNode3->shared_edge( pCNodeD );

		if(!pCEdge3D) return fmu_kCLEAN_NO_ACTION;
      if( pCEdge3D->frozen() ) return fmu_kCLEAN_NO_ACTION;

      pCEdge3E = pCNode3->shared_edge( pCNodeE );
		if(!pCEdge3E) return fmu_kCLEAN_NO_ACTION;
      if( pCEdge3E->frozen() ) return fmu_kCLEAN_NO_ACTION;
      pCQuadE = pCEdge3D->other_quad( pCQuadC );

		if(!pCQuadE) return fmu_kCLEAN_NO_ACTION;
      pCNodeF = pCQuadE->opposite_node( pCNode3 );
   }
   else if( iPattern == 4 )
   {
      //Node E must be 5 valent, node B must be 4 or 5 valent
      //Two additional quads in this pattern
      if( get_valence( pCNodeE ) != 5 ) return fmu_kCLEAN_NO_ACTION;
      if( get_valence( pCNodeB ) > 4 ) return fmu_kCLEAN_NO_ACTION;
      pCEdge3D = pCNode3->shared_edge( pCNodeD );
		if(!pCEdge3D) return fmu_kCLEAN_NO_ACTION;
      if( pCEdge3D->frozen() ) return fmu_kCLEAN_NO_ACTION;
      pCEdge3E = pCNode3->shared_edge( pCNodeE );
		if(!pCEdge3E) return fmu_kCLEAN_NO_ACTION;
      if( pCEdge3E->frozen() ) return fmu_kCLEAN_NO_ACTION;
      pCQuadE = pCEdge3D->other_quad( pCQuadC );
		if(!pCQuadE) return fmu_kCLEAN_NO_ACTION;
      pCNodeF = pCQuadE->opposite_node( pCNode3 );
      pCEdge07 = pCNode0->shared_edge( pCNode7 );
      if( pCEdge07->frozen() ) return fmu_kCLEAN_NO_ACTION;
      pCEdge0A = pCNode0->shared_edge( pCNodeA );
		if(!pCEdge0A) return fmu_kCLEAN_NO_ACTION;
      if( pCEdge0A->frozen() ) return fmu_kCLEAN_NO_ACTION;
      pCQuadF = pCEdge07->other_quad( pCQuad3 );
		if(!pCQuadF) return fmu_kCLEAN_NO_ACTION;
      pCNodeG = pCQuadF->opposite_node( pCNode0 );
   }

   if( lDrawCleaning )
   {
   //draw_open();
   //pCQuadD->draw_me( dp_kCLR_PINK );
   //pCQuadC->draw_me( dp_kCLR_ORANGE );
   //pCQuadB->draw_me( dp_kCLR_LTBLUE );
   //pCQuadA->draw_me( dp_kCLR_BLUE );
   //pCQuad3->draw_me( dp_kCLR_RED );
   //pCQuad2->draw_me( dp_kCLR_YELLOW );
   //pCQuad1->draw_me( dp_kCLR_LTMGNT );
   //pCQuad0->draw_me( dp_kCLR_BLUE );
       // Refresh debug graphics
        //EraseElements();
        pCQuadD->draw_me( dp_kCLR_YELLOW );
        pCQuadC->draw_me( dp_kCLR_YELLOW );
        pCQuadB->draw_me( dp_kCLR_YELLOW );
        pCQuadA->draw_me( dp_kCLR_YELLOW );
        pCQuad3->draw_me( dp_kCLR_YELLOW );
        pCQuad2->draw_me( dp_kCLR_YELLOW );
        pCQuad1->draw_me( dp_kCLR_YELLOW );
        pCQuad0->draw_me( dp_kCLR_YELLOW );
   //if( iPattern == 3 ) pCQuadE->draw_me( dp_kCLR_CYAN );
        if( iPattern == 3 ) pCQuadE->draw_me( dp_kCLR_YELLOW );
        else if( iPattern == 4 )
        {
   //pCQuadE->draw_me( dp_kCLR_DKGRN );
   //pCQuadF->draw_me( dp_kCLR_LTMGNT );
            pCQuadE->draw_me( dp_kCLR_YELLOW );
            pCQuadF->draw_me( dp_kCLR_YELLOW );
        }
   //draw_close();
   }

   delete_this_quad( pCQuad0 );
   delete_this_quad( pCQuad1 );
   delete_this_quad( pCQuad2 );
   delete_this_quad( pCQuad3 );
   delete_this_quad( pCQuadB );
   delete_this_quad( pCQuadC );
   delete_this_quad( pCQuadD );
   if( iPattern < 3 )       delete_this_quad( pCQuadA );
   else if( iPattern == 3 ) delete_this_quad( pCQuadE );
   else if( iPattern == 4 )
   {
      delete_this_quad( pCQuadA );
      delete_this_quad( pCQuadE );
      delete_this_quad( pCQuadF );
   }

   if( delete_this_edge( pCEdge0 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge1 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge2 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge3 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge12 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge23 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge34 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge45 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge2C ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( iPattern < 3 )
   {
      if( delete_this_edge( pCEdge01 ) == fmu_kCLEAN_ABORT )
         return fmu_kCLEAN_ABORT;
      if( delete_this_edge( pCEdge1B ) == fmu_kCLEAN_ABORT )
         return fmu_kCLEAN_ABORT;
   }
   else if( iPattern == 3 )
   {
      if( delete_this_edge( pCEdge3D ) == fmu_kCLEAN_ABORT )
         return fmu_kCLEAN_ABORT;
      if( delete_this_edge( pCEdge3E ) == fmu_kCLEAN_ABORT )
         return fmu_kCLEAN_ABORT;
   }
   else if( iPattern == 4 )
   {
      if( delete_this_edge( pCEdge01 ) == fmu_kCLEAN_ABORT )
         return fmu_kCLEAN_ABORT;
      if( delete_this_edge( pCEdge07 ) == fmu_kCLEAN_ABORT )
         return fmu_kCLEAN_ABORT;
      if( delete_this_edge( pCEdge0A ) == fmu_kCLEAN_ABORT )
         return fmu_kCLEAN_ABORT;
      if( delete_this_edge( pCEdge1B ) == fmu_kCLEAN_ABORT )
         return fmu_kCLEAN_ABORT;
      if( delete_this_edge( pCEdge3D ) == fmu_kCLEAN_ABORT )
         return fmu_kCLEAN_ABORT;
      if( delete_this_edge( pCEdge3E ) == fmu_kCLEAN_ABORT )
         return fmu_kCLEAN_ABORT;
   }

   if( delete_this_node( pCCenterNode ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_node( pCNode2 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_node( pCNode4 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( iPattern < 3 )
   {
      if( delete_this_node( pCNode1 ) == fmu_kCLEAN_ABORT )
	return fmu_kCLEAN_ABORT;
   }
   else if( iPattern == 3 )
   {
      if( delete_this_node( pCNode3 ) == fmu_kCLEAN_ABORT )
         return fmu_kCLEAN_ABORT;
   }
   else if( iPattern == 4 )
   {
      if( delete_this_node( pCNode0 ) == fmu_kCLEAN_ABORT )
         return fmu_kCLEAN_ABORT;
      if( delete_this_node( pCNode1 ) == fmu_kCLEAN_ABORT )
         return fmu_kCLEAN_ABORT;
      if( delete_this_node( pCNode3 ) == fmu_kCLEAN_ABORT )
         return fmu_kCLEAN_ABORT;
   }

   //Now make 4 new quads, depending on the pattern
   if( iPattern == 1 )
   {
      pCQuad0 = new elmQClass( pCNode6, pCNode7, pCNode0, pCNode5, pCMaster );
      if( !pCQuad0 ) return fmu_kCLEAN_NO_MEMORY;
      pCQuad1 = new elmQClass( pCNode5, pCNode0, pCNode3, pCNodeE, pCMaster );
      if( !pCQuad1 ) return fmu_kCLEAN_NO_MEMORY;
      pCQuad2 = new elmQClass( pCNode0, pCNodeA, pCNodeB, pCNode3, pCMaster );
      if( !pCQuad2 ) return fmu_kCLEAN_NO_MEMORY;
      pCQuad3 = new elmQClass( pCNode3, pCNodeB, pCNodeC, pCNodeD, pCMaster );
      if( !pCQuad3 ) return fmu_kCLEAN_NO_MEMORY;
   }
   else if( iPattern == 2 )
   {
      pCQuad0 = new elmQClass( pCNode5, pCNode6, pCNode3, pCNodeE, pCMaster );
      if( !pCQuad0 ) return fmu_kCLEAN_NO_MEMORY;
      pCQuad1 = new elmQClass( pCNode6, pCNode7, pCNode0, pCNode3, pCMaster );
      if( !pCQuad1 ) return fmu_kCLEAN_NO_MEMORY;
      pCQuad2 = new elmQClass( pCNode3, pCNode0, pCNodeC, pCNodeD, pCMaster );
      if( !pCQuad2 ) return fmu_kCLEAN_NO_MEMORY;
      pCQuad3 = new elmQClass( pCNode0, pCNodeA, pCNodeB, pCNodeC, pCMaster );
      if( !pCQuad3 ) return fmu_kCLEAN_NO_MEMORY;
   }
   else if( iPattern == 3 )
   {
      pCQuad0 = new elmQClass( pCNode6, pCNode7, pCNode0, pCNode1, pCMaster );
      if( !pCQuad0 ) return fmu_kCLEAN_NO_MEMORY;
      pCQuad1 = new elmQClass( pCNodeE, pCNode5, pCNode6, pCNode1, pCMaster );
      if( !pCQuad1 ) return fmu_kCLEAN_NO_MEMORY;
      pCQuad2 = new elmQClass( pCNodeE, pCNode1, pCNodeB, pCNodeC, pCMaster );
      if( !pCQuad2 ) return fmu_kCLEAN_NO_MEMORY;
      pCQuad3 = new elmQClass( pCNodeE, pCNodeC, pCNodeD, pCNodeF, pCMaster );
      if( !pCQuad3 ) return fmu_kCLEAN_NO_MEMORY;
      pCNode3 = pCNode1;      //Make things easier below
   }
   else if( iPattern == 4 )
   {
      pCQuad0 = new elmQClass( pCNode6, pCNode7, pCNodeG, pCNodeA, pCMaster );
      if( !pCQuad0 ) return fmu_kCLEAN_NO_MEMORY;
      pCQuad1 = new elmQClass( pCNode6, pCNodeA, pCNodeB, pCNodeC, pCMaster );
      if( !pCQuad1 ) return fmu_kCLEAN_NO_MEMORY;
      pCQuad2 = new elmQClass( pCNode5, pCNode6, pCNodeC, pCNodeD, pCMaster );
      if( !pCQuad2 ) return fmu_kCLEAN_NO_MEMORY;
      pCQuad3 = new elmQClass( pCNodeE, pCNode5, pCNodeD, pCNodeF, pCMaster );
      if( !pCQuad3 ) return fmu_kCLEAN_NO_MEMORY;
      pCNode0 = pCNodeF;      //Make things easier below
      pCNode3 = pCNodeG;
   }


   if( lDrawCleaning )
   {
      nodClassDynArray CSmoothList;
      CSmoothList.append( pCNodeA );
      CSmoothList.append( pCNodeB );
      CSmoothList.append( pCNodeD );
      CSmoothList.append( pCNodeE );
      CSmoothList.append( pCNodeC );
      CSmoothList.append( pCNode5 );
      CSmoothList.append( pCNode6 );
      CSmoothList.append( pCNode7 );
      CSmoothList.append( pCNode0 );
      CSmoothList.append( pCNode3 );
      smooth_mesh( &CSmoothList );
      //draw_mesh();
      // Refresh debug graphics
        //EraseElements();
      //pCQuad3->draw_me( dp_kCLR_RED );
      //pCQuad2->draw_me( dp_kCLR_YELLOW );
      //pCQuad1->draw_me( dp_kCLR_BLUE );
      //pCQuad0->draw_me( dp_kCLR_BLUE );
        pCQuad3->draw_me( dp_kCLR_BLUE );
        pCQuad2->draw_me( dp_kCLR_BLUE );
        pCQuad1->draw_me( dp_kCLR_BLUE );
        pCQuad0->draw_me( dp_kCLR_BLUE );
      //draw_close();
   }

   reclean( pCNode0 );
   reclean( pCNode3 );
   reclean( pCNode5 );
   reclean( pCNode6 );
   reclean( pCNode7 );
   reclean( pCNodeB );
   reclean( pCNodeC );
   reclean( pCNodeD );

   return fmu_kCLEAN_MESH_CHANGED;
 
}


int QuadCleanTool::fix_double_trans( nodClass *pCCenterNode,
                                         nodClassDynArray &CNodes,
                                         fmuEdgeDynArray &CEdges,
                                         elmQClassDynArray &CQuads,
                                         int iValence0 )
{

   nodClass *pCNode0 = CNodes.get();
   nodClass *pCNode1 = CNodes.next();
   nodClass *pCNode2 = CNodes.next(2);
   nodClass *pCNode3 = CNodes.next(3);
   nodClass *pCNode4 = CNodes.next(4);
   nodClass *pCNode5 = CNodes.next(5);
   nodClass *pCNode6 = CNodes.next(6);
   nodClass *pCNode7 = CNodes.next(7);
   fmuEdge  *pCEdge0 = CEdges.get();
   fmuEdge  *pCEdge1 = CEdges.next();
   fmuEdge  *pCEdge2 = CEdges.next(2);
   fmuEdge  *pCEdge3 = CEdges.next(3);
   elmQClass *pCQuad0 = CQuads.get();
   elmQClass *pCQuad1 = CQuads.next();
   elmQClass *pCQuad2 = CQuads.next(2);
   elmQClass *pCQuad3 = CQuads.next(3);

   // Check to see if all nodes exist - NM
	  if (	pCNode0 == NULL ||
			pCNode1 == NULL ||
			pCNode2 == NULL ||
			pCNode3 == NULL ||
			pCNode4 == NULL ||
			pCNode5 == NULL ||
			pCNode6 == NULL ||
			pCNode7 == NULL)
			return fmu_kCLEAN_NO_ACTION;

   // Check to see if all nodes exist  - NM
   if( !pCEdge0 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCEdge1 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCEdge2 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCEdge3 ) return fmu_kCLEAN_NO_ACTION;

   // Check to see if all nodes exist - NM
   if( !pCQuad0 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCQuad1 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCQuad2 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCQuad3 ) return fmu_kCLEAN_NO_ACTION;

   if( pCEdge0->frozen() ) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge1->frozen() ) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge2->frozen() ) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge3->frozen() ) return fmu_kCLEAN_NO_ACTION;

   fmuEdge *pCEdge01 = pCNode0->shared_edge( pCNode1 );
	if(!pCEdge01) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge01->frozen() ) return fmu_kCLEAN_NO_ACTION;

   fmuEdge *pCEdge12 = pCNode1->shared_edge( pCNode2 );
	if(!pCEdge12) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge12->frozen() ) return fmu_kCLEAN_NO_ACTION;

   fmuEdge *pCEdge23 = pCNode2->shared_edge( pCNode3 );
	if(!pCEdge23) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge23->frozen() ) return fmu_kCLEAN_NO_ACTION;

   fmuEdge *pCEdge34 = pCNode3->shared_edge( pCNode4 );
	if(!pCEdge34) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge34->frozen() ) return fmu_kCLEAN_NO_ACTION;

   fmuEdge *pCEdge45 = pCNode4->shared_edge( pCNode5 );
	if(!pCEdge45) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge45->frozen() ) return fmu_kCLEAN_NO_ACTION;

   fmuEdge *pCEdge07 = pCNode0->shared_edge( pCNode7 );
	if(!pCEdge07) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge07->frozen() ) return fmu_kCLEAN_NO_ACTION;


   elmQClass *pCQuadA = pCEdge07->other_quad( pCQuad3 );
	if(!pCQuadA) return fmu_kCLEAN_NO_ACTION;

   nodClass *pCNodeA = pCQuadA->opposite_node( pCNode0 );
   elmQClass *pCQuadC = pCEdge12->other_quad( pCQuad0 );

	if(!pCQuadC) return fmu_kCLEAN_NO_ACTION;
   nodClass *pCNodeC = pCQuadC->next_node( pCNode1 );
   nodClass *pCNodeD = pCQuadC->opposite_node( pCNode1 );
   fmuEdge *pCEdge1C = pCNode1->shared_edge( pCNodeC );
	if(!pCEdge1C) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge1C->frozen() ) return fmu_kCLEAN_NO_ACTION;

   elmQClass *pCQuadB = pCEdge1C->other_quad( pCQuadC );
   nodClass *pCNodeB = pCQuadB->opposite_node( pCNode1 );
   elmQClass *pCQuadD = pCEdge23->other_quad( pCQuad1 );
	if(!pCQuadD) return fmu_kCLEAN_NO_ACTION;
   nodClass *pCNodeE = pCQuadD->opposite_node( pCNode2 );

   fmuEdge *pCEdge2D = pCNode2->shared_edge( pCNodeD );
	if(!pCEdge2D) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge2D->frozen() ) return fmu_kCLEAN_NO_ACTION;

   fmuEdge *pCEdge3E = pCNode3->shared_edge( pCNodeE );
	if(!pCEdge3E) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge3E->frozen() ) return fmu_kCLEAN_NO_ACTION;

   elmQClass *pCQuadE = pCEdge3E->other_quad( pCQuadD );
	if(!pCQuadE) return fmu_kCLEAN_NO_ACTION;
   nodClass *pCNodeF = pCQuadE->opposite_node( pCNode3 );
   elmQClass *pCQuadF = pCEdge45->other_quad( pCQuad2 );
	if(!pCQuadF) return fmu_kCLEAN_NO_ACTION;

   nodClass *pCNodeG = pCQuadF->opposite_node( pCNode4 );
   fmuEdge *pCEdgeQAQB, *pCEdgeQEQF;
   int iValenceA = get_valence( pCNodeA );
   int iValenceB = get_valence( pCNodeB );
   int iValenceF = get_valence( pCNodeF );
   int iValenceG = get_valence( pCNodeG );

   //This is the valence of node0. There are two mirror patterns to handle
   if( iValence0 == 3 )
   {
      pCEdgeQAQB = pCNode1->shared_edge( pCNodeA );
      if( pCEdgeQAQB->frozen() ) return fmu_kCLEAN_NO_ACTION;
      pCEdgeQEQF = pCNode4->shared_edge( pCNodeF );
      if( pCEdgeQEQF->frozen() ) return fmu_kCLEAN_NO_ACTION;
      if( ( iValenceA != 5 && iValenceF != 5 ) ||
          iValenceB != 4 || iValenceG != 4 ) return fmu_kCLEAN_NO_ACTION;
      if( iValenceA == 4 && pCNodeA->hardpt() ) return fmu_kCLEAN_NO_ACTION;
      if( iValenceF == 4 && pCNodeF->hardpt() ) return fmu_kCLEAN_NO_ACTION;
   }
   else
   {
      pCEdgeQAQB = pCNode0->shared_edge( pCNodeB );
      if( pCEdgeQAQB->frozen() ) return fmu_kCLEAN_NO_ACTION;
      pCEdgeQEQF = pCNode3->shared_edge( pCNodeG );
      if( pCEdgeQEQF->frozen() ) return fmu_kCLEAN_NO_ACTION;
      if( ( iValenceB != 5 && iValenceG != 5 ) ||
          iValenceA != 4 || iValenceF != 4 ) return fmu_kCLEAN_NO_ACTION;
      if( iValenceB == 4 && pCNodeB->hardpt() ) return fmu_kCLEAN_NO_ACTION;
      if( iValenceG == 4 && pCNodeG->hardpt() ) return fmu_kCLEAN_NO_ACTION;
   }

   if( lDrawCleaning )
   {
   //draw_open();
   //pCQuadF->draw_me( dp_kCLR_MAGNTA );
   //pCQuadE->draw_me( dp_kCLR_LTBLUE );
   //pCQuadD->draw_me( dp_kCLR_LTMGNT );
   //pCQuadC->draw_me( dp_kCLR_BLUE );
   //pCQuadB->draw_me( dp_kCLR_PINK );
   //pCQuadA->draw_me( dp_kCLR_GOLD );
   //pCQuad3->draw_me( dp_kCLR_RED );
   //pCQuad2->draw_me( dp_kCLR_YELLOW );
   //pCQuad1->draw_me( dp_kCLR_DKGRN );
   //pCQuad0->draw_me( dp_kCLR_BLUE );
       // Refresh debug graphics
        //EraseElements();
        pCQuadF->draw_me( dp_kCLR_YELLOW );
        pCQuadE->draw_me( dp_kCLR_YELLOW );
        pCQuadD->draw_me( dp_kCLR_YELLOW );
        pCQuadC->draw_me( dp_kCLR_YELLOW );
        pCQuadB->draw_me( dp_kCLR_YELLOW );
        pCQuadA->draw_me( dp_kCLR_YELLOW );
        pCQuad3->draw_me( dp_kCLR_YELLOW );
        pCQuad2->draw_me( dp_kCLR_YELLOW );
        pCQuad1->draw_me( dp_kCLR_YELLOW );
        pCQuad0->draw_me( dp_kCLR_YELLOW );
   //draw_close();
   }

   delete_this_quad( pCQuad0 );
   delete_this_quad( pCQuad1 );
   delete_this_quad( pCQuad2 );
   delete_this_quad( pCQuad3 );
   delete_this_quad( pCQuadA );
   delete_this_quad( pCQuadB );
   delete_this_quad( pCQuadC );
   delete_this_quad( pCQuadD );
   delete_this_quad( pCQuadE );
   delete_this_quad( pCQuadF );

   if( delete_this_edge( pCEdge0 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge1 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge2 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge3 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge01 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge12 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge23 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge34 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge45 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge07 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge1C ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge2D ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge3E ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdgeQAQB ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdgeQEQF ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;

   if( delete_this_node( pCCenterNode ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_node( pCNode0 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_node( pCNode1 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_node( pCNode2 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_node( pCNode3 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_node( pCNode4 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;

   pCQuad0 = new elmQClass( pCNodeA, pCNodeB, pCNodeC, pCNode7, pCMaster );
   if( !pCQuad0 ) return fmu_kCLEAN_NO_MEMORY;
   pCQuad1 = new elmQClass( pCNode7, pCNodeC, pCNodeD, pCNode6, pCMaster );
   if( !pCQuad1 ) return fmu_kCLEAN_NO_MEMORY;
   pCQuad2 = new elmQClass( pCNode6, pCNodeD, pCNodeE, pCNode5, pCMaster );
   if( !pCQuad2 ) return fmu_kCLEAN_NO_MEMORY;
   pCQuad3 = new elmQClass( pCNode5, pCNodeE, pCNodeF, pCNodeG, pCMaster );
   if( !pCQuad3 ) return fmu_kCLEAN_NO_MEMORY;

   if( lDrawCleaning )
   {
      nodClassDynArray CSmoothList;
      CSmoothList.append( pCNodeA );
      CSmoothList.append( pCNodeC );
      CSmoothList.append( pCNode6 );
      CSmoothList.append( pCNodeE );
      CSmoothList.append( pCNodeG );
      smooth_mesh( &CSmoothList );
      //draw_mesh();
      //pCQuad3->draw_me( dp_kCLR_CYAN );
      //pCQuad2->draw_me( dp_kCLR_YELLOW );
      //pCQuad1->draw_me( dp_kCLR_DKGRN );
      //pCQuad0->draw_me( dp_kCLR_GOLD );
        // Refresh debug graphics
        //EraseElements();
        pCQuad3->draw_me( dp_kCLR_BLUE );
        pCQuad2->draw_me( dp_kCLR_BLUE );
        pCQuad1->draw_me( dp_kCLR_BLUE );
        pCQuad0->draw_me( dp_kCLR_BLUE );
       
      //draw_close();
   }

   reclean( pCNodeB );
   reclean( pCNodeD );
   reclean( pCNodeF );
   reclean( pCNode5 );
   reclean( pCNode7 );

   return fmu_kCLEAN_MESH_CHANGED;
}


int QuadCleanTool::clean_435453( nodClass *pCCenterNode,
                                     nodClassDynArray &CNodes,
                                     fmuEdgeDynArray &CEdges,
                                     elmQClassDynArray &CQuads )
{

   nodClass *pCNode0 = CNodes.get();
   nodClass *pCNode1 = CNodes.next();
   nodClass *pCNode2 = CNodes.next(2);
   nodClass *pCNode3 = CNodes.next(3);
   nodClass *pCNode4 = CNodes.next(4);
   nodClass *pCNode5 = CNodes.next(5);
   nodClass *pCNode6 = CNodes.next(6);
   nodClass *pCNode7 = CNodes.next(7);
   fmuEdge  *pCEdge0 = CEdges.get();
   fmuEdge  *pCEdge1 = CEdges.next();
   fmuEdge  *pCEdge2 = CEdges.next(2);
   fmuEdge  *pCEdge3 = CEdges.next(3);
   elmQClass *pCQuad0 = CQuads.get();
   elmQClass *pCQuad1 = CQuads.next();
   elmQClass *pCQuad2 = CQuads.next(2);
   elmQClass *pCQuad3 = CQuads.next(3);

   // Check to see if all nodes exist - NM
	  if (	pCNode0 == NULL ||
			pCNode1 == NULL ||
			pCNode2 == NULL ||
			pCNode3 == NULL ||
			pCNode4 == NULL ||
			pCNode5 == NULL ||
			pCNode6 == NULL ||
			pCNode7 == NULL)
			return fmu_kCLEAN_NO_ACTION;

   // Check to see if all nodes exist  - NM
   if( !pCEdge0 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCEdge1 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCEdge2 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCEdge3 ) return fmu_kCLEAN_NO_ACTION;

   // Check to see if all nodes exist - NM
   if( !pCQuad0 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCQuad1 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCQuad2 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCQuad3 ) return fmu_kCLEAN_NO_ACTION;

   if( pCEdge0->frozen() ) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge1->frozen() ) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge2->frozen() ) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge3->frozen() ) return fmu_kCLEAN_NO_ACTION;

   fmuEdge *pCEdge01 = pCNode0->shared_edge( pCNode1 );
   if( pCEdge01->frozen() ) return fmu_kCLEAN_NO_ACTION;
   fmuEdge *pCEdge07 = pCNode0->shared_edge( pCNode7 );
   if( pCEdge07->frozen() ) return fmu_kCLEAN_NO_ACTION;
   fmuEdge *pCEdge34 = pCNode3->shared_edge( pCNode4 );
   if( pCEdge34->frozen() ) return fmu_kCLEAN_NO_ACTION;
   fmuEdge *pCEdge45 = pCNode4->shared_edge( pCNode5 );
   if( pCEdge45->frozen() ) return fmu_kCLEAN_NO_ACTION;

   elmQClass *pCQuadA = pCEdge01->other_quad( pCQuad0 );
   nodClass *pCNodeA = pCQuadA->opposite_node( pCNode0 );
   elmQClass *pCQuadB = pCEdge34->other_quad( pCQuad1 );
   nodClass *pCNodeB = pCQuadB->opposite_node( pCNode4 );

   if( lDrawCleaning )
   {
   //draw_open();
   //pCQuadB->draw_me( dp_kCLR_DKGRN );
   //pCQuadA->draw_me( dp_kCLR_CYAN );
   //pCQuad3->draw_me( dp_kCLR_YELLOW );
   //pCQuad2->draw_me( dp_kCLR_LTBLUE );
   //pCQuad1->draw_me( dp_kCLR_RED );
   //pCQuad0->draw_me( dp_kCLR_BLUE );
        // Refresh debug graphics
        //EraseElements();
        pCQuadB->draw_me( dp_kCLR_YELLOW );
        pCQuadA->draw_me( dp_kCLR_YELLOW );
        pCQuad3->draw_me( dp_kCLR_YELLOW );
        pCQuad2->draw_me( dp_kCLR_YELLOW );
        pCQuad1->draw_me( dp_kCLR_YELLOW );
        pCQuad0->draw_me( dp_kCLR_YELLOW );
        
   //draw_close();
   }

   delete_this_quad( pCQuad0 );
   delete_this_quad( pCQuad1 );
   delete_this_quad( pCQuad2 );
   delete_this_quad( pCQuad3 );
   delete_this_quad( pCQuadA );
   delete_this_quad( pCQuadB );

   if( delete_this_edge( pCEdge0 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge1 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge2 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge3 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge01 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge07 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge34 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge45 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;

   if( delete_this_node( pCCenterNode ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_node( pCNode0 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_node( pCNode4 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;

   //Three new quads
   pCQuad0 = new elmQClass( pCNode2, pCNode7, pCNodeA, pCNode1, pCMaster );
   if( !pCQuad0 ) return fmu_kCLEAN_NO_MEMORY;
   pCQuad1 = new elmQClass( pCNode2, pCNode5, pCNode6, pCNode7, pCMaster );
   if( !pCQuad1 ) return fmu_kCLEAN_NO_MEMORY;
   pCQuad2 = new elmQClass( pCNode3, pCNodeB, pCNode5, pCNode2, pCMaster );
   if( !pCQuad2 ) return fmu_kCLEAN_NO_MEMORY;

   if( lDrawCleaning )
   {
      nodClassDynArray CSmoothList;
      CSmoothList.append( pCNode1 );
      CSmoothList.append( pCNode2 );
      CSmoothList.append( pCNode3 );
      CSmoothList.append( pCNode5 );
      CSmoothList.append( pCNode6 );
      CSmoothList.append( pCNode7 );
      CSmoothList.append( pCNodeA );
      CSmoothList.append( pCNodeB );
      smooth_mesh( &CSmoothList );
      //draw_mesh();
      //pCQuad0->draw_me( dp_kCLR_CYAN );
      //pCQuad1->draw_me( dp_kCLR_RED );
      //pCQuad2->draw_me( dp_kCLR_DKGRN );
        // Refresh debug graphics
        //EraseElements();
        pCQuad2->draw_me( dp_kCLR_BLUE );
        pCQuad1->draw_me( dp_kCLR_BLUE );
        pCQuad0->draw_me( dp_kCLR_BLUE );
      
      //draw_close();
   }

   reclean( pCNode1 );
   reclean( pCNode2 );
   reclean( pCNode3 );
   reclean( pCNode5 );
   reclean( pCNode6 );
   reclean( pCNode7 );

   return fmu_kCLEAN_MESH_CHANGED;
   
}


int QuadCleanTool::clean_443544354( nodClass *pCCenterNode,
                                        nodClassDynArray &CNodes,
                                        fmuEdgeDynArray &CEdges,
                                        elmQClassDynArray &CQuads,
                                        int iSide )
{
   
   nodClass *pCNode0 = CNodes.get();
   nodClass *pCNode1 = CNodes.next();
   nodClass *pCNode2 = CNodes.next(2);
   nodClass *pCNode3 = CNodes.next(3);
   nodClass *pCNode4 = CNodes.next(4);
   nodClass *pCNode5 = CNodes.next(5);
   nodClass *pCNode6 = CNodes.next(6);
   nodClass *pCNode7 = CNodes.next(7);
   fmuEdge  *pCEdge0 = CEdges.get();
   fmuEdge  *pCEdge1 = CEdges.next();
   fmuEdge  *pCEdge2 = CEdges.next(2);
   fmuEdge  *pCEdge3 = CEdges.next(3);
   elmQClass *pCQuad0 = CQuads.get();
   elmQClass *pCQuad1 = CQuads.next();
   elmQClass *pCQuad2 = CQuads.next(2);
   elmQClass *pCQuad3 = CQuads.next(3);

   // Check to see if all nodes exist - NM
	  if (	pCNode0 == NULL ||
			pCNode1 == NULL ||
			pCNode2 == NULL ||
			pCNode3 == NULL ||
			pCNode4 == NULL ||
			pCNode5 == NULL ||
			pCNode6 == NULL ||
			pCNode7 == NULL)
			return fmu_kCLEAN_NO_ACTION;

   // Check to see if all nodes exist  - NM
   if( !pCEdge0 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCEdge1 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCEdge2 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCEdge3 ) return fmu_kCLEAN_NO_ACTION;

   // Check to see if all nodes exist - NM
   if( !pCQuad0 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCQuad1 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCQuad2 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCQuad3 ) return fmu_kCLEAN_NO_ACTION;

   if( pCEdge0->frozen() ) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge1->frozen() ) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge2->frozen() ) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge3->frozen() ) return fmu_kCLEAN_NO_ACTION;

   if( iSide == 1 &&
       pCNode1->hardpt() &&
       pCNode2->hardpt() &&
       pCNode3->hardpt() )
   {
      fmuVector CVec_23( pCNode2, pCNode3 );
      fmuVector CVec_21( pCNode2, pCNode1 );
      if( CNormal.angle( CVec_23, CVec_21 ) >= fmu_kSHAPE_ANGLE_TOL )
        return fmu_kCLEAN_NO_ACTION;
   }

   if( iSide == 1 &&
       pCNode5->hardpt() &&
       pCNode6->hardpt() &&
       pCNode7->hardpt() )
   {
      fmuVector CVec_67( pCNode6, pCNode7 );
      fmuVector CVec_65( pCNode6, pCNode5 );
      if( CNormal.angle( CVec_67, CVec_65 ) >= fmu_kSHAPE_ANGLE_TOL )
        return fmu_kCLEAN_NO_ACTION;
   }

   if( iSide == 2 &&
       pCNode1->hardpt() &&
       pCNode0->hardpt() &&
       pCNode7->hardpt() )
   {
      fmuVector CVec_01( pCNode0, pCNode1 );
      fmuVector CVec_07( pCNode0, pCNode7 );
      if( CNormal.angle( CVec_01, CVec_07 ) >= fmu_kSHAPE_ANGLE_TOL )
        return fmu_kCLEAN_NO_ACTION;
   }

   if( iSide == 2 &&
       pCNode3->hardpt() &&
       pCNode4->hardpt() &&
       pCNode5->hardpt() )
   {
      fmuVector CVec_45( pCNode4, pCNode5 );
      fmuVector CVec_43( pCNode4, pCNode3 );
      if( CNormal.angle( CVec_45, CVec_43 ) >= fmu_kSHAPE_ANGLE_TOL )
        return fmu_kCLEAN_NO_ACTION;
   }

   if( lDrawCleaning )
   {
   //draw_open();
   //pCQuad0->draw_me( dp_kCLR_RED );
   //pCQuad1->draw_me( dp_kCLR_MAGNTA);
   //pCQuad2->draw_me( dp_kCLR_LTMGNT );
   //pCQuad3->draw_me( dp_kCLR_PINK );
        // Refresh debug graphics
        //EraseElements();
        pCQuad0->draw_me( dp_kCLR_YELLOW );
        pCQuad1->draw_me( dp_kCLR_YELLOW );
        pCQuad2->draw_me( dp_kCLR_YELLOW );
        pCQuad3->draw_me( dp_kCLR_YELLOW );
       
   //draw_close();
   }

   delete_this_quad( pCQuad0 );
   delete_this_quad( pCQuad1 );
   delete_this_quad( pCQuad2 );
   delete_this_quad( pCQuad3 );

   if( delete_this_edge( pCEdge0 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge1 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge2 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge3 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;

   if( delete_this_node( pCCenterNode ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;

   //Build 3 new quads
   if( iSide == 1 )
   {
      pCQuad0 = new elmQClass( pCNode6, pCNode7, pCNode0, pCNode5, pCMaster );
      if( !pCQuad0 ) return fmu_kCLEAN_NO_MEMORY;
      pCQuad1 = new elmQClass( pCNode5, pCNode0, pCNode1, pCNode4, pCMaster );
      if( !pCQuad1 ) return fmu_kCLEAN_NO_MEMORY;
      pCQuad2 = new elmQClass( pCNode4, pCNode1, pCNode2, pCNode3, pCMaster );
      if( !pCQuad2 ) return fmu_kCLEAN_NO_MEMORY;
   }
   else if( iSide == 2 )
   {
      pCQuad0 = new elmQClass( pCNode6, pCNode7, pCNode0, pCNode1, pCMaster );
      if( !pCQuad0 ) return fmu_kCLEAN_NO_MEMORY;
      pCQuad1 = new elmQClass( pCNode6, pCNode1, pCNode2, pCNode5, pCMaster );
      if( !pCQuad1 ) return fmu_kCLEAN_NO_MEMORY;
      pCQuad2 = new elmQClass( pCNode2, pCNode3, pCNode4, pCNode5, pCMaster );
      if( !pCQuad2 ) return fmu_kCLEAN_NO_MEMORY;
   }

   if( lDrawCleaning )
   {
      nodClassDynArray CSmoothList;
      CSmoothList.append( pCNode0 );
      CSmoothList.append( pCNode1 );
      CSmoothList.append( pCNode2 );
      CSmoothList.append( pCNode3 );
      CSmoothList.append( pCNode4 );
      CSmoothList.append( pCNode5 );
      CSmoothList.append( pCNode6 );
      CSmoothList.append( pCNode7 );
      smooth_mesh( &CSmoothList );
      //draw_mesh();
      //pCQuad0->draw_me( dp_kCLR_RED );
      //pCQuad1->draw_me( dp_kCLR_PINK );
      //pCQuad2->draw_me( dp_kCLR_MAGNTA );
        // Refresh debug graphics
        //EraseElements();
        pCQuad0->draw_me( dp_kCLR_BLUE );
        pCQuad1->draw_me( dp_kCLR_BLUE );
        pCQuad2->draw_me( dp_kCLR_BLUE );

      //draw_close();
   }

   reclean( pCNode0 );
   reclean( pCNode1 );
   reclean( pCNode2 );
   reclean( pCNode4 );
   reclean( pCNode5 );
   reclean( pCNode6 );

   return fmu_kCLEAN_MESH_CHANGED;

}


int QuadCleanTool::clean_435443544( nodClass *pCCenterNode,
                                        nodClassDynArray &CNodes,
                                        fmuEdgeDynArray &CEdges,
                                        elmQClassDynArray &CQuads,
                                        int iSide )
{
   
   nodClass *pCNode0 = CNodes.get();
   nodClass *pCNode1 = CNodes.next();
   nodClass *pCNode2 = CNodes.next(2);
   nodClass *pCNode3 = CNodes.next(3);
   nodClass *pCNode4 = CNodes.next(4);
   nodClass *pCNode5 = CNodes.next(5);
   nodClass *pCNode6 = CNodes.next(6);
   nodClass *pCNode7 = CNodes.next(7);
   fmuEdge  *pCEdge0 = CEdges.get();
   fmuEdge  *pCEdge1 = CEdges.next();
   fmuEdge  *pCEdge2 = CEdges.next(2);
   fmuEdge  *pCEdge3 = CEdges.next(3);
   elmQClass *pCQuad0 = CQuads.get();
   elmQClass *pCQuad1 = CQuads.next();
   elmQClass *pCQuad2 = CQuads.next(2);
   elmQClass *pCQuad3 = CQuads.next(3);
   // Check to see if all nodes exist - NM
	  if (	pCNode0 == NULL ||
			pCNode1 == NULL ||
			pCNode2 == NULL ||
			pCNode3 == NULL ||
			pCNode4 == NULL ||
			pCNode5 == NULL ||
			pCNode6 == NULL ||
			pCNode7 == NULL)
			return fmu_kCLEAN_NO_ACTION;

   // Check to see if all nodes exist  - NM
   if( !pCEdge0 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCEdge1 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCEdge2 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCEdge3 ) return fmu_kCLEAN_NO_ACTION;

   // Check to see if all nodes exist - NM
   if( !pCQuad0 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCQuad1 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCQuad2 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCQuad3 ) return fmu_kCLEAN_NO_ACTION;


   if( pCEdge0->frozen() ) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge1->frozen() ) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge2->frozen() ) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge3->frozen() ) return fmu_kCLEAN_NO_ACTION;

   fmuEdge *pCEdge01 = pCNode0->shared_edge( pCNode1 );
   if( pCEdge01->frozen() ) return fmu_kCLEAN_NO_ACTION;
   fmuEdge *pCEdge07 = pCNode0->shared_edge( pCNode7 );
   if( pCEdge07->frozen() ) return fmu_kCLEAN_NO_ACTION;
   fmuEdge *pCEdge34 = pCNode3->shared_edge( pCNode4 );
   if( pCEdge34->frozen() ) return fmu_kCLEAN_NO_ACTION;
   fmuEdge *pCEdge45 = pCNode4->shared_edge( pCNode5 );
   if( pCEdge45->frozen() ) return fmu_kCLEAN_NO_ACTION;

   elmQClass *pCQuadA = pCEdge01->other_quad( pCQuad0 );
	if(!pCQuadA) return fmu_kCLEAN_NO_ACTION;
   nodClass *pCNodeA = pCQuadA->opposite_node( pCNode0 );
   elmQClass *pCQuadB = pCEdge34->other_quad( pCQuad1 );
	if(!pCQuadB) return fmu_kCLEAN_NO_ACTION;
   nodClass *pCNodeB = pCQuadB->opposite_node( pCNode4 );

   if( iSide == 1 &&
       pCNodeA->hardpt() &&
       pCNode1->hardpt() &&
       pCNode2->hardpt() )
   {
      fmuVector CVec_12( pCNode1, pCNode2 );
      fmuVector CVec_1A( pCNode1, pCNodeA );
      if( CNormal.angle( CVec_12, CVec_1A ) >= fmu_kSHAPE_ANGLE_TOL )
        return fmu_kCLEAN_NO_ACTION;
   }

   if( iSide == 1 &&
       pCNodeB->hardpt() &&
       pCNode5->hardpt() &&
       pCNode6->hardpt() )
   {
      fmuVector CVec_56( pCNode5, pCNode6 );
      fmuVector CVec_5B( pCNode5, pCNodeB );
      if( CNormal.angle( CVec_56, CVec_5B ) >= fmu_kSHAPE_ANGLE_TOL )
        return fmu_kCLEAN_NO_ACTION;
   }

   if( iSide != 1 &&
       pCNode2->hardpt() &&
       pCNode3->hardpt() &&
       pCNodeB->hardpt() )
   {
      fmuVector CVec_3B( pCNode3, pCNodeB );
      fmuVector CVec_32( pCNode3, pCNode2 );
      if( CNormal.angle( CVec_3B, CVec_32 ) >= fmu_kSHAPE_ANGLE_TOL )
        return fmu_kCLEAN_NO_ACTION;
   }

   if( iSide != 1 &&
       pCNode6->hardpt() &&
       pCNode7->hardpt() &&
       pCNodeA->hardpt() )
   {
      fmuVector CVec_7A( pCNode7, pCNodeA );
      fmuVector CVec_76( pCNode7, pCNode6 );
      if( CNormal.angle( CVec_7A, CVec_76 ) >= fmu_kSHAPE_ANGLE_TOL )
        return fmu_kCLEAN_NO_ACTION;
   }

   //Check for bad angles already in an element. This is not the best
   //way to clean it. Shape cleaning does better without changing the
   //mesh as there is a 3 valent node to work wit.
   if( pCNode7->hardpt() &&
       pCNodeA->hardpt() &&
       pCNode1->hardpt() )
   {
      fmuVector CVec_A1( pCNodeA, pCNode1 );
      fmuVector CVec_A7( pCNodeA, pCNode7 );
      if( CNormal.angle( CVec_A1, CVec_A7 ) >= fmu_kSHAPE_ANGLE_TOL )
        return fmu_kCLEAN_NO_ACTION;
   }

   if( pCNode3->hardpt() &&
       pCNodeB->hardpt() &&
       pCNode5->hardpt() )
   {
      fmuVector CVec_B5( pCNodeB, pCNode5 );
      fmuVector CVec_B3( pCNodeB, pCNode3 );
      if( CNormal.angle( CVec_B5, CVec_B3 ) >= fmu_kSHAPE_ANGLE_TOL )
        return fmu_kCLEAN_NO_ACTION;
   }

   //Check if the pattern will wrap around a corner, causing a distorted
   //element.
   if( iSide == 1 &&
       pCNodeA->hardpt() &&
       pCNode7->hardpt() &&
       pCNode2->hardpt() )
   {
      fmuVector CVec_7A( pCNode7, pCNodeA );
      fmuVector CVec_72( pCNode7, pCNode2 );
      if( CNormal.angle( CVec_7A, CVec_72 ) >= fmu_kSHAPE_ANGLE_TOL )
        return fmu_kCLEAN_NO_ACTION;
   }

   if( iSide == 1 &&
       pCNodeB->hardpt() &&
       pCNode3->hardpt() &&
       pCNode6->hardpt() )
   {
      fmuVector CVec_3B( pCNode3, pCNodeB );
      fmuVector CVec_36( pCNode3, pCNode6 );
      if( CNormal.angle( CVec_3B, CVec_36 ) >= fmu_kSHAPE_ANGLE_TOL )
        return fmu_kCLEAN_NO_ACTION;
   }

   if( iSide != 1 &&
       pCNode2->hardpt() &&
       pCNode5->hardpt() &&
       pCNodeB->hardpt() )
   {
      fmuVector CVec_52( pCNode5, pCNode2 );
      fmuVector CVec_5B( pCNode5, pCNodeB );
      if( CNormal.angle( CVec_52, CVec_5B ) >= fmu_kSHAPE_ANGLE_TOL )
        return fmu_kCLEAN_NO_ACTION;
   }

   if( iSide != 1 &&
       pCNode6->hardpt() &&
       pCNode1->hardpt() &&
       pCNodeA->hardpt() )
   {
      fmuVector CVec_16( pCNode1, pCNode6 );
      fmuVector CVec_1A( pCNode1, pCNodeA );
      if( CNormal.angle( CVec_16, CVec_1A ) >= fmu_kSHAPE_ANGLE_TOL )
        return fmu_kCLEAN_NO_ACTION;
   }

   if( lDrawCleaning )
   {
   //draw_open();
   //pCQuadB->draw_me( dp_kCLR_DKGRN );
   //pCQuadA->draw_me( dp_kCLR_BLUE );
   //pCQuad0->draw_me( dp_kCLR_BLUE );
   //pCQuad1->draw_me( dp_kCLR_GOLD);
   //pCQuad2->draw_me( dp_kCLR_LTBLUE );
   //pCQuad3->draw_me( dp_kCLR_YELLOW );

        // Refresh debug graphics
        //EraseElements();
        pCQuadB->draw_me( dp_kCLR_YELLOW );
        pCQuadA->draw_me( dp_kCLR_YELLOW );
        pCQuad0->draw_me( dp_kCLR_YELLOW );
        pCQuad1->draw_me( dp_kCLR_YELLOW );
        pCQuad2->draw_me( dp_kCLR_YELLOW );
        pCQuad3->draw_me( dp_kCLR_YELLOW );
       
   }

   delete_this_quad( pCQuad0 );
   delete_this_quad( pCQuad1 );
   delete_this_quad( pCQuad2 );
   delete_this_quad( pCQuad3 );
   delete_this_quad( pCQuadA );
   delete_this_quad( pCQuadB );

   if( delete_this_edge( pCEdge0 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge1 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge2 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge3 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge01 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge07 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge34 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge45 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;

   if( delete_this_node( pCCenterNode ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_node( pCNode0 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_node( pCNode4 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;

   //Three new quads
   if( iSide == 1 )
   {
      pCQuad0 = new elmQClass( pCNode2, pCNode7, pCNodeA, pCNode1, pCMaster );
      if( !pCQuad0 ) return fmu_kCLEAN_NO_MEMORY;
      pCQuad1 = new elmQClass( pCNode3, pCNode6, pCNode7, pCNode2, pCMaster );
      if( !pCQuad1 ) return fmu_kCLEAN_NO_MEMORY;
      pCQuad2 = new elmQClass( pCNode5, pCNode6, pCNode3, pCNodeB, pCMaster );
      if( !pCQuad2 ) return fmu_kCLEAN_NO_MEMORY;
   }
   else
   {
      pCQuad0 = new elmQClass( pCNode6, pCNode7, pCNodeA, pCNode1, pCMaster );
      if( !pCQuad0 ) return fmu_kCLEAN_NO_MEMORY;
      pCQuad1 = new elmQClass( pCNode5, pCNode6, pCNode1, pCNode2, pCMaster );
      if( !pCQuad1 ) return fmu_kCLEAN_NO_MEMORY;
      pCQuad2 = new elmQClass( pCNodeB, pCNode5, pCNode2, pCNode3, pCMaster );
      if( !pCQuad2 ) return fmu_kCLEAN_NO_MEMORY;
   }

   if( lDrawCleaning )
   {
      nodClassDynArray CSmoothList;
      CSmoothList.append( pCNode1 );
      CSmoothList.append( pCNode2 );
      CSmoothList.append( pCNode3 );
      CSmoothList.append( pCNode5 );
      CSmoothList.append( pCNode7 );
      CSmoothList.append( pCNode6 );
      CSmoothList.append( pCNodeA );
      CSmoothList.append( pCNodeB );
      smooth_mesh( &CSmoothList );
      //draw_mesh();
      //pCQuad0->draw_me( dp_kCLR_BLUE );
      //pCQuad1->draw_me( dp_kCLR_GOLD );
      //pCQuad2->draw_me( dp_kCLR_DKGRN );
        // Refresh debug graphics
        //EraseElements();
        pCQuad0->draw_me( dp_kCLR_BLUE );
        pCQuad1->draw_me( dp_kCLR_BLUE );
        pCQuad2->draw_me( dp_kCLR_BLUE );
       
      //draw_close();
   }

   reclean( pCNode1 );
   reclean( pCNode2 );
   reclean( pCNode3 );
   reclean( pCNode5 );
   reclean( pCNode6 );
   reclean( pCNode7 );

   return fmu_kCLEAN_MESH_CHANGED;
}


int QuadCleanTool::clean_434543544( nodClass *pCCenterNode,
                                        nodClassDynArray &CNodes,
                                        fmuEdgeDynArray &CEdges,
                                        elmQClassDynArray &CQuads)
{
  
   nodClass *pCNode0 = CNodes.get();
   nodClass *pCNode1 = CNodes.next();
   nodClass *pCNode2 = CNodes.next(2);
   nodClass *pCNode3 = CNodes.next(3);
   nodClass *pCNode4 = CNodes.next(4);
   nodClass *pCNode5 = CNodes.next(5);
   nodClass *pCNode6 = CNodes.next(6);
   nodClass *pCNode7 = CNodes.next(7);
   fmuEdge  *pCEdge0 = CEdges.get();
   fmuEdge  *pCEdge1 = CEdges.next();
   fmuEdge  *pCEdge2 = CEdges.next(2);
   fmuEdge  *pCEdge3 = CEdges.next(3);
   elmQClass *pCQuad0 = CQuads.get();
   elmQClass *pCQuad1 = CQuads.next();
   elmQClass *pCQuad2 = CQuads.next(2);
   elmQClass *pCQuad3 = CQuads.next(3);

   // Check to see if all nodes exist - NM
	  if (	pCNode0 == NULL ||
			pCNode1 == NULL ||
			pCNode2 == NULL ||
			pCNode3 == NULL ||
			pCNode4 == NULL ||
			pCNode5 == NULL ||
			pCNode6 == NULL ||
			pCNode7 == NULL)
			return fmu_kCLEAN_NO_ACTION;

   // Check to see if all nodes exist  - NM
   if( !pCEdge0 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCEdge1 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCEdge2 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCEdge3 ) return fmu_kCLEAN_NO_ACTION;

   // Check to see if all nodes exist - NM
   if( !pCQuad0 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCQuad1 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCQuad2 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCQuad3 ) return fmu_kCLEAN_NO_ACTION;


   if( pCEdge0->frozen() ) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge1->frozen() ) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge2->frozen() ) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge3->frozen() ) return fmu_kCLEAN_NO_ACTION;

   fmuEdge *pCEdge34 = pCNode3->shared_edge( pCNode4 );
	if(!pCEdge34) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge34->frozen() ) return fmu_kCLEAN_NO_ACTION;

   fmuEdge *pCEdge45 = pCNode4->shared_edge( pCNode5 );
	if(!pCEdge45) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge45->frozen() ) return fmu_kCLEAN_NO_ACTION;

   elmQClass *pCQuadA = pCEdge34->other_quad( pCQuad1 );
	if(!pCQuadA) return fmu_kCLEAN_NO_ACTION;
   nodClass *pCNodeA = pCQuadA->opposite_node( pCNode4 );
	if(!pCNodeA) return fmu_kCLEAN_NO_ACTION;

   if( pCNode5->hardpt() && pCNode6->hardpt() && pCNodeA->hardpt() )
   {
      fmuVector CVec56( pCNode5, pCNode6 );
      fmuVector CVec5A( pCNode5, pCNodeA );
      if( CNormal.angle( CVec56, CVec5A ) > fmu_kSHAPE_ANGLE_TOL )
         return fmu_kCLEAN_NO_ACTION;
   }

   if( lDrawCleaning )
   {
   //draw_open();
   //pCQuadA->draw_me( dp_kCLR_BLUE );
   //pCQuad0->draw_me( dp_kCLR_BLUE );
   //pCQuad1->draw_me( dp_kCLR_GOLD);
   //pCQuad2->draw_me( dp_kCLR_LTBLUE );
   //pCQuad3->draw_me( dp_kCLR_YELLOW );
        // Refresh debug graphics
        //EraseElements();
        pCQuadA->draw_me( dp_kCLR_YELLOW );
        pCQuad0->draw_me( dp_kCLR_YELLOW );
        pCQuad1->draw_me( dp_kCLR_YELLOW );
        pCQuad2->draw_me( dp_kCLR_YELLOW );
        pCQuad3->draw_me( dp_kCLR_YELLOW );
      
   //draw_close();
   }

   delete_this_quad( pCQuad0 );
   delete_this_quad( pCQuad1 );
   delete_this_quad( pCQuad2 );
   delete_this_quad( pCQuad3 );
   delete_this_quad( pCQuadA );

   if( delete_this_edge( pCEdge0 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge1 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge2 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge3 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge34 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge45 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;

   if( delete_this_node( pCCenterNode ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_node( pCNode4 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;

   //Three new quads
   pCQuad0 = new elmQClass( pCNode0, pCNode1, pCNode2, pCNode3, pCMaster );
   if( !pCQuad0 ) return fmu_kCLEAN_NO_MEMORY;
   pCQuad1 = new elmQClass( pCNode0, pCNode3, pCNode6, pCNode7, pCMaster );
   if( !pCQuad1 ) return fmu_kCLEAN_NO_MEMORY;
   pCQuad2 = new elmQClass( pCNode3, pCNodeA, pCNode5, pCNode6, pCMaster );
   if( !pCQuad2 ) return fmu_kCLEAN_NO_MEMORY;

   if( lDrawCleaning )
   {
      nodClassDynArray CSmoothList;
      CSmoothList.append( pCNode0 );
      CSmoothList.append( pCNode1 );
      CSmoothList.append( pCNode2 );
      CSmoothList.append( pCNode3 );
      CSmoothList.append( pCNodeA );
      CSmoothList.append( pCNode5 );
      CSmoothList.append( pCNode6 );
      CSmoothList.append( pCNode7 );
      smooth_mesh( &CSmoothList );
      //draw_mesh();
      //pCQuad0->draw_me( dp_kCLR_BLUE );
      //pCQuad1->draw_me( dp_kCLR_GOLD );
      //pCQuad2->draw_me( dp_kCLR_DKGRN );
        // Refresh debug graphics
        //EraseElements();
        pCQuad0->draw_me( dp_kCLR_BLUE );
        pCQuad1->draw_me( dp_kCLR_BLUE );
        pCQuad2->draw_me( dp_kCLR_BLUE );
        
      //draw_close();
   }

   reclean( pCNode0 );
   reclean( pCNode1 );
   reclean( pCNode2 );
   reclean( pCNode3 );
   reclean( pCNodeA );
   reclean( pCNode5 );
   reclean( pCNode6 );
   reclean( pCNode7 );

   return fmu_kCLEAN_MESH_CHANGED;
  
}


int QuadCleanTool::clean_434453454( nodClass *pCCenterNode,
                                        nodClassDynArray &CNodes,
                                        fmuEdgeDynArray &CEdges,
                                        elmQClassDynArray &CQuads)
{
   
   nodClass *pCNode0 = CNodes.get();
   nodClass *pCNode1 = CNodes.next();
   nodClass *pCNode2 = CNodes.next(2);
   nodClass *pCNode3 = CNodes.next(3);
   nodClass *pCNode4 = CNodes.next(4);
   nodClass *pCNode5 = CNodes.next(5);
   nodClass *pCNode6 = CNodes.next(6);
   nodClass *pCNode7 = CNodes.next(7);
   fmuEdge  *pCEdge0 = CEdges.get();
   fmuEdge  *pCEdge1 = CEdges.next();
   fmuEdge  *pCEdge2 = CEdges.next(2);
   fmuEdge  *pCEdge3 = CEdges.next(3);
   elmQClass *pCQuad0 = CQuads.get();
   elmQClass *pCQuad1 = CQuads.next();
   elmQClass *pCQuad2 = CQuads.next(2);
   elmQClass *pCQuad3 = CQuads.next(3);

   // Check to see if all nodes exist - NM
	  if (	pCNode0 == NULL ||
			pCNode1 == NULL ||
			pCNode2 == NULL ||
			pCNode3 == NULL ||
			pCNode4 == NULL ||
			pCNode5 == NULL ||
			pCNode6 == NULL ||
			pCNode7 == NULL)
			return fmu_kCLEAN_NO_ACTION;

   // Check to see if all nodes exist  - NM
   if( !pCEdge0 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCEdge1 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCEdge2 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCEdge3 ) return fmu_kCLEAN_NO_ACTION;

   // Check to see if all nodes exist - NM
   if( !pCQuad0 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCQuad1 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCQuad2 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCQuad3 ) return fmu_kCLEAN_NO_ACTION;

   if( pCEdge0->frozen() ) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge1->frozen() ) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge2->frozen() ) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge3->frozen() ) return fmu_kCLEAN_NO_ACTION;

   fmuEdge *pCEdge34 = pCNode3->shared_edge( pCNode4 );
	if(!pCEdge34) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge34->frozen() ) return fmu_kCLEAN_NO_ACTION;
   fmuEdge *pCEdge45 = pCNode4->shared_edge( pCNode5 );
	if(!pCEdge45) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge45->frozen() ) return fmu_kCLEAN_NO_ACTION;

   elmQClass *pCQuadA = pCEdge34->other_quad( pCQuad1 );
	if(!pCQuadA) return fmu_kCLEAN_NO_ACTION;
   nodClass *pCNodeA = pCQuadA->opposite_node( pCNode4 );
	if(!pCNodeA) return fmu_kCLEAN_NO_ACTION;

   if( pCNode2->hardpt() && pCNode3->hardpt() && pCNodeA->hardpt() )
   {
      fmuVector CVec32( pCNode3, pCNode2 );
      fmuVector CVec3A( pCNode3, pCNodeA );
      if( CNormal.angle( CVec3A, CVec32 ) > fmu_kSHAPE_ANGLE_TOL )
         return fmu_kCLEAN_NO_ACTION;
   }

   if( lDrawCleaning )
   {
   //draw_open();
   //pCQuadA->draw_me( dp_kCLR_BLUE );
   //pCQuad0->draw_me( dp_kCLR_BLUE );
   //pCQuad1->draw_me( dp_kCLR_GOLD);
   //pCQuad2->draw_me( dp_kCLR_LTBLUE );
   //pCQuad3->draw_me( dp_kCLR_YELLOW );
        
        // Refresh debug graphics
        //EraseElements();
        pCQuadA->draw_me( dp_kCLR_YELLOW );
        pCQuad0->draw_me( dp_kCLR_YELLOW );
        pCQuad1->draw_me( dp_kCLR_YELLOW );
        pCQuad2->draw_me( dp_kCLR_YELLOW );
        pCQuad3->draw_me( dp_kCLR_YELLOW );
        
   //draw_close();
   }

   delete_this_quad( pCQuad0 );
   delete_this_quad( pCQuad1 );
   delete_this_quad( pCQuad2 );
   delete_this_quad( pCQuad3 );
   delete_this_quad( pCQuadA );

   if( delete_this_edge( pCEdge0 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge1 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge2 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge3 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge34 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge45 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;

   if( delete_this_node( pCCenterNode ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_node( pCNode4 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;

   //Three new quads
   pCQuad0 = new elmQClass( pCNode0, pCNode1, pCNode2, pCNode5, pCMaster );
   if( !pCQuad0 ) return fmu_kCLEAN_NO_MEMORY;
   pCQuad1 = new elmQClass( pCNode0, pCNode5, pCNode6, pCNode7, pCMaster );
   if( !pCQuad1 ) return fmu_kCLEAN_NO_MEMORY;
   pCQuad2 = new elmQClass( pCNode3, pCNodeA, pCNode5, pCNode2, pCMaster );
   if( !pCQuad2 ) return fmu_kCLEAN_NO_MEMORY;

   if( lDrawCleaning )
   {
      nodClassDynArray CSmoothList;
      CSmoothList.append( pCNode0 );
      CSmoothList.append( pCNode1 );
      CSmoothList.append( pCNode2 );
      CSmoothList.append( pCNode3 );
      CSmoothList.append( pCNodeA );
      CSmoothList.append( pCNode5 );
      CSmoothList.append( pCNode6 );
      CSmoothList.append( pCNode7 );
      smooth_mesh( &CSmoothList );
      //draw_mesh();
      //pCQuad0->draw_me( dp_kCLR_BLUE );
      //pCQuad1->draw_me( dp_kCLR_GOLD );
      //pCQuad2->draw_me( dp_kCLR_DKGRN );
        // Refresh debug graphics
        //EraseElements();
        pCQuad0->draw_me( dp_kCLR_BLUE );
        pCQuad1->draw_me( dp_kCLR_BLUE );
        pCQuad2->draw_me( dp_kCLR_BLUE );
        
      //draw_close();
   }

   reclean( pCNode0 );
   reclean( pCNode1 );
   reclean( pCNode2 );
   reclean( pCNode3 );
   reclean( pCNodeA );
   reclean( pCNode5 );
   reclean( pCNode6 );
   reclean( pCNode7 );

   return fmu_kCLEAN_MESH_CHANGED;
   
}


int QuadCleanTool::reverse_open_face( nodClass *pCCenterNode,
                                          nodClassDynArray &CNodes,
                                          fmuEdgeDynArray &CEdges,
                                          elmQClassDynArray &CQuads )
{
  
   nodClassDynArray CDelNodes, CRingNodes;
   fmuEdgeDynArray CDelEdges;
   elmQClassDynArray CDelQuads;
   int iCallStatus;

   if(CEdges.next()->frozen() ) return fmu_kCLEAN_NO_ACTION;

   CDelQuads[0] = CQuads.get();
   CDelQuads[1] = CQuads.next();
   CDelEdges[0] = CEdges.next();
   CRingNodes[0] = pCCenterNode;
   CRingNodes[1] = CNodes.get();
   CRingNodes[2] = CNodes.next();
   CRingNodes[3] = CNodes.next(2);
   CRingNodes[4] = CNodes.next(3);
   CRingNodes[5] = CNodes.next(4);

   iCallStatus = fill_three( CDelNodes, CDelEdges, CDelQuads, CRingNodes );

   return iCallStatus;
   
}


int QuadCleanTool::clean_5300003( nodClass *pCCenterNode,
                                      nodClassDynArray &CNodes,
                                      fmuEdgeDynArray &CEdges,
                                      elmQClassDynArray &CQuads )
{
   
   nodClass *pCNodeA = 0, *pCNodeB = 0;
   elmQClass *pCQuadA = 0, *pCQuadB = 0;
   nodClass *pCNode0 = CNodes.get();
   nodClass *pCNode1 = CNodes.next();
   nodClass *pCNode2 = CNodes.next(2);
   nodClass *pCNode3 = CNodes.next(3);
   nodClass *pCNode4 = CNodes.next(4);
   nodClass *pCNode5 = CNodes.next(5);
   nodClass *pCNode6 = CNodes.next(6);
   nodClass *pCNode7 = CNodes.next(7);
   nodClass *pCNode8 = CNodes.next(8);
   nodClass *pCNode9 = CNodes.next(9);
   fmuEdge  *pCEdge1 = CEdges.next();
   fmuEdge  *pCEdge2 = CEdges.next(2);
   elmQClass *pCQuad0 = CQuads.get();
   elmQClass *pCQuad1 = CQuads.next();
   elmQClass *pCQuad2 = CQuads.next(2);

   // Check to see if all nodes exist - NM
	  if (	pCNode0 == NULL ||
			pCNode1 == NULL ||
			pCNode2 == NULL ||
			pCNode3 == NULL ||
			pCNode4 == NULL ||
			pCNode5 == NULL ||
			pCNode6 == NULL ||
			pCNode7 == NULL ||
			pCNode8 == NULL ||
			pCNode9 == NULL)
			return fmu_kCLEAN_NO_ACTION;

   // Check to see if all nodes exist  - NM
   if( !pCEdge1 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCEdge2 ) return fmu_kCLEAN_NO_ACTION;
 

   // Check to see if all nodes exist - NM
   if( !pCQuad0 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCQuad1 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCQuad2 ) return fmu_kCLEAN_NO_ACTION;


   if( pCEdge1->frozen() ) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge2->frozen() ) return fmu_kCLEAN_NO_ACTION;

   //Three quads are to be added. Make sure there is room so
   //squashed elements are avoided
   if( pCNode2->hardpt() && pCNode8->hardpt() )
   {
      double dLen2, dLen8, dLenSide;
      fmuVector CVecC2( pCCenterNode, pCNode2 );
      fmuVector CVecC8( pCCenterNode, pCNode8 );
      dLen2 = CVecC2.length();
      if( pCNode3->hardpt() )
      {
         fmuVector CVec23( pCNode2, pCNode3 );
         dLenSide = CVec23.length();
      }
      else
      {
         fmuVector CVec21( pCNode2, pCNode1 );
         dLenSide = CVec21.length();
      }
      if( dLen2 * 1.5 < dLenSide ) return fmu_kCLEAN_NO_ACTION;
      dLen8 = CVecC8.length();
      if( pCNode7->hardpt() )
      {
         fmuVector CVec87( pCNode8, pCNode7 );
         dLenSide = CVec87.length();
      }
      else
      {
         fmuVector CVec89( pCNode8, pCNode9 );
         dLenSide = CVec89.length();
      }
      if( dLen8 * 1.5 < dLenSide ) return fmu_kCLEAN_NO_ACTION;
   }

   if( lDrawCleaning )
   {
   //draw_open();
   //pCQuad0->draw_me( dp_kCLR_DKGRN );
   //pCQuad1->draw_me( dp_kCLR_CYAN );
   //pCQuad2->draw_me( dp_kCLR_BLUE );
        // Refresh debug graphics
        //EraseElements();
        pCQuad0->draw_me( dp_kCLR_YELLOW );
        pCQuad1->draw_me( dp_kCLR_YELLOW );
        pCQuad2->draw_me( dp_kCLR_YELLOW );
       
   //draw_close();
   }

   delete_this_quad( pCQuad0 );
   delete_this_quad( pCQuad1 );
   delete_this_quad( pCQuad2 );

   if( delete_this_edge( pCEdge1 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge2 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;

   //Need two new nodes
   two_middle_nodes( pCNode2, pCNode5, &pCNodeA, &pCNodeB );

   //Build 5 new quads
   pCQuad0 = new elmQClass( pCNode0, pCNode1, pCNode2, pCNodeA, pCMaster );
   if( !pCQuad0 ) return fmu_kCLEAN_NO_MEMORY;
   pCQuad1 = new elmQClass( pCNode2, pCNode3, pCNode4, pCNodeA, pCMaster );
   if( !pCQuad1 ) return fmu_kCLEAN_NO_MEMORY;
   pCQuad2 = new elmQClass( pCNodeA, pCNode4, pCNode5, pCNodeB, pCMaster );
   if( !pCQuad2 ) return fmu_kCLEAN_NO_MEMORY;
   pCQuadA = new elmQClass( pCNodeB, pCNode5, pCNode6, pCCenterNode,
                                                                pCMaster );
   if( !pCQuadA ) return fmu_kCLEAN_NO_MEMORY;
   pCQuadB = new elmQClass( pCCenterNode, pCNode0, pCNodeA, pCNodeB,
                                                                pCMaster );
   if( !pCQuadB ) return fmu_kCLEAN_NO_MEMORY;

   if( lDrawCleaning )
   {
      nodClassDynArray CSmoothList;
      CSmoothList.append( pCNode0 );
      CSmoothList.append( pCNode1 );
      CSmoothList.append( pCNode2 );
      CSmoothList.append( pCNode3 );
      CSmoothList.append( pCNode4 );
      CSmoothList.append( pCNode5 );
      CSmoothList.append( pCNode6 );
      CSmoothList.append( pCNodeA );
      CSmoothList.append( pCNodeB );
      CSmoothList.append( pCCenterNode );
      smooth_mesh( &CSmoothList );
      //draw_mesh();
      //pCQuad0->draw_me( dp_kCLR_DKGRN );
      //pCQuad1->draw_me( dp_kCLR_CYAN );
      //pCQuad2->draw_me( dp_kCLR_BLUE );
      //pCQuadA->draw_me( dp_kCLR_GOLD );
      //pCQuadB->draw_me( dp_kCLR_RED );
        // Refresh debug graphics
        //EraseElements();
        pCQuad0->draw_me( dp_kCLR_BLUE );
        pCQuad1->draw_me( dp_kCLR_BLUE );
        pCQuad2->draw_me( dp_kCLR_BLUE );
        pCQuadA->draw_me( dp_kCLR_BLUE );
        pCQuadB->draw_me( dp_kCLR_BLUE );
      
      //draw_close();
   }

   reclean( pCNode0 );
   reclean( pCNode5 );
   reclean( pCCenterNode );

   return fmu_kCLEAN_MESH_CHANGED;
}


int QuadCleanTool::clean_55443( nodClass *pCCenterNode,
                                    nodClassDynArray &CNodes,
                                    fmuEdgeDynArray &CEdges,
                                    elmQClassDynArray &CQuads )
{
   
   nodClass *pCNode0 = CNodes.get();
   nodClass *pCNode1 = CNodes.next();
   nodClass *pCNode2 = CNodes.next(2);
   nodClass *pCNode3 = CNodes.next(3);
   nodClass *pCNode4 = CNodes.next(4);
   fmuEdge  *pCEdge1 = CEdges.next();
   elmQClass *pCQuad0 = CQuads.get();
   elmQClass *pCQuad1 = CQuads.next();


   // Check to see if all nodes exist - NM
	  if (	pCNode0 == NULL ||
			pCNode1 == NULL ||
			pCNode2 == NULL ||
			pCNode3 == NULL ||
			pCNode4 == NULL )
			return fmu_kCLEAN_NO_ACTION;

   // Check to see if all nodes exist  - NM
   if( !pCEdge1 ) return fmu_kCLEAN_NO_ACTION;
  
   // Check to see if all nodes exist - NM
   if( !pCQuad0 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCQuad1 ) return fmu_kCLEAN_NO_ACTION;
 

   if( pCEdge1->frozen() ) return fmu_kCLEAN_NO_ACTION;

   fmuEdge *pCEdge01 = pCNode0->shared_edge( pCNode1 );
	if(!pCEdge01) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge01->frozen() ) return fmu_kCLEAN_NO_ACTION;
   
   elmQClass *pCQuadB = pCEdge01->other_quad( pCQuad0 );
	if(!pCQuadB) return fmu_kCLEAN_NO_ACTION;
   nodClass *pCNodeD = pCQuadB->opposite_node( pCNode0 );
   nodClass *pCNodeC = pCQuadB->next_node( pCNode0 );
    if(!pCNodeC) return fmu_kCLEAN_NO_ACTION;
   fmuEdge *pCEdge0C = pCNode0->shared_edge( pCNodeC );
	if(!pCEdge0C) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge0C->frozen() ) return fmu_kCLEAN_NO_ACTION;
   
   elmQClass *pCQuadA = pCEdge0C->other_quad( pCQuadB );
	if(!pCQuadA) return fmu_kCLEAN_NO_ACTION;
   nodClass *pCNodeB = pCQuadA->opposite_node( pCNode0 );
   nodClass *pCNodeA = pCQuadA->next_node( pCNode0 );
   if( get_valence( pCNodeA ) > 3 ) return fmu_kCLEAN_NO_ACTION;

   if( lDrawCleaning )
   {
   //draw_open();
   //pCQuad0->draw_me( dp_kCLR_BLUE );
   //pCQuad1->draw_me( dp_kCLR_CYAN );
   //pCQuadA->draw_me( dp_kCLR_GRBLUE );
   //pCQuadB->draw_me( dp_kCLR_LTBLUE );
       // Refresh debug graphics
        //EraseElements();
        pCQuad0->draw_me( dp_kCLR_YELLOW );
        pCQuad1->draw_me( dp_kCLR_YELLOW );
        pCQuadA->draw_me( dp_kCLR_YELLOW );
        pCQuadB->draw_me( dp_kCLR_YELLOW );
      
   //draw_close();
   }

   delete_this_quad( pCQuad0 );
   delete_this_quad( pCQuad1 );
   delete_this_quad( pCQuadA );
   delete_this_quad( pCQuadB );

   if( delete_this_edge( pCEdge1 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge01 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge0C ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;

   //Need two new nodes
   nodClass *pCNodeE, *pCNodeF;
   two_middle_nodes( pCNode0, pCNodeD, &pCNodeE, &pCNodeF );

   //Build 6 new quads
   elmQClass *pCQuadC, *pCQuadD;
   pCQuad0 = new elmQClass( pCNode0, pCNode3, pCNode4, pCCenterNode,
                                                                pCMaster );
   if( !pCQuad0 ) return fmu_kCLEAN_NO_MEMORY;
   pCQuad1 = new elmQClass( pCNode0, pCNodeE, pCNode2, pCNode3, pCMaster );
   if( !pCQuad1 ) return fmu_kCLEAN_NO_MEMORY;
   pCQuadA = new elmQClass( pCNode0, pCNodeA, pCNodeF, pCNodeE, pCMaster );
   if( !pCQuadA ) return fmu_kCLEAN_NO_MEMORY;
   pCQuadB = new elmQClass( pCNodeE, pCNodeF, pCNode1, pCNode2, pCMaster );
   if( !pCQuadB ) return fmu_kCLEAN_NO_MEMORY;
   pCQuadC = new elmQClass( pCNodeA, pCNodeB, pCNodeC, pCNodeF, pCMaster );
   if( !pCQuadC ) return fmu_kCLEAN_NO_MEMORY;
   pCQuadD = new elmQClass( pCNodeF, pCNodeC, pCNodeD, pCNode1, pCMaster );
   if( !pCQuadD ) return fmu_kCLEAN_NO_MEMORY;

   if( lDrawCleaning )
   {
      nodClassDynArray CSmoothList;
      CSmoothList.append( pCNode0 );
      CSmoothList.append( pCNode1 );
      CSmoothList.append( pCNode2 );
      CSmoothList.append( pCNode3 );
      CSmoothList.append( pCNode4 );
      CSmoothList.append( pCNodeA );
      CSmoothList.append( pCNodeB );
      CSmoothList.append( pCNodeC );
      CSmoothList.append( pCNodeD );
      CSmoothList.append( pCNodeE );
      CSmoothList.append( pCNodeF );
      CSmoothList.append( pCCenterNode );
      smooth_mesh( &CSmoothList );
      //draw_mesh();
      //pCQuad0->draw_me( dp_kCLR_BLUE );
      //pCQuad1->draw_me( dp_kCLR_CYAN );
      //pCQuadA->draw_me( dp_kCLR_GRBLUE );
      //pCQuadB->draw_me( dp_kCLR_LTBLUE );
      //pCQuadC->draw_me( dp_kCLR_GOLD );
      //pCQuadD->draw_me( dp_kCLR_YELLOW );
        // Refresh debug graphics
        //EraseElements();
        pCQuad0->draw_me( dp_kCLR_BLUE );
        pCQuad1->draw_me( dp_kCLR_BLUE );
        pCQuadA->draw_me( dp_kCLR_BLUE );
        pCQuadB->draw_me( dp_kCLR_BLUE );
        pCQuadC->draw_me( dp_kCLR_BLUE );
        pCQuadD->draw_me( dp_kCLR_BLUE );

      //draw_close();
   }

   reclean( pCNode0 );
   reclean( pCNode1 );
   reclean( pCNode2 );
   reclean( pCNode3 );
   reclean( pCNodeA );
   reclean( pCCenterNode );

   return fmu_kCLEAN_MESH_CHANGED;
  
}


int QuadCleanTool::clean_503445( nodClass *pCCenterNode,
                                     nodClassDynArray &CNodes,
                                     fmuEdgeDynArray &CEdges,
                                     elmQClassDynArray &CQuads )
{
   
   nodClass *pCNode0 = CNodes.get();
   nodClass *pCNode1 = CNodes.next();
   nodClass *pCNode2 = CNodes.next(2);
   nodClass *pCNode3 = CNodes.next(3);
   nodClass *pCNode4 = CNodes.next(4);
   fmuEdge  *pCEdge1 = CEdges.next();
   elmQClass *pCQuad0 = CQuads.get();
   elmQClass *pCQuad1 = CQuads.next();

   // Check to see if all nodes exist - NM
	  if (	pCNode0 == NULL ||
			pCNode1 == NULL ||
			pCNode2 == NULL ||
			pCNode3 == NULL ||
			pCNode4 == NULL)
			return fmu_kCLEAN_NO_ACTION;

   // Check to see if all nodes exist  - NM
   if( !pCEdge1 ) return fmu_kCLEAN_NO_ACTION;
 

   // Check to see if all nodes exist - NM
   if( !pCQuad0 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCQuad1 ) return fmu_kCLEAN_NO_ACTION;


   if( pCEdge1->frozen() ) return fmu_kCLEAN_NO_ACTION;

   fmuEdge *pCEdge34 = pCNode3->shared_edge( pCNode4 );
	if(!pCEdge34) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge34->frozen() ) return fmu_kCLEAN_NO_ACTION;

   elmQClass *pCQuadA = pCEdge34->other_quad( pCQuad1 );
	if(!pCQuadA) return fmu_kCLEAN_NO_ACTION;
   nodClass *pCNodeA = pCQuadA->opposite_node( pCNode4 );
   nodClass *pCNodeB = pCQuadA->prev_node( pCNode4 );
	if(!pCNodeB) return fmu_kCLEAN_NO_ACTION;
   fmuEdge *pCEdge4B = pCNode4->shared_edge( pCNodeB );
	if(!pCEdge4B) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge4B->frozen() ) return fmu_kCLEAN_NO_ACTION;

   elmQClass *pCQuadB = pCEdge4B->other_quad( pCQuadA );
	if(!pCQuadB) return fmu_kCLEAN_NO_ACTION;
   nodClass *pCNodeC = pCQuadB->opposite_node( pCNode4 );
   nodClass *pCNodeD = pCQuadB->prev_node( pCNode4 );
   if( get_valence( pCNodeD ) > 3 ) return fmu_kCLEAN_NO_ACTION;

   if( lDrawCleaning )
   {
   //draw_open();
   //pCQuad0->draw_me( dp_kCLR_BLUE );
   //pCQuad1->draw_me( dp_kCLR_CYAN );
   //pCQuadA->draw_me( dp_kCLR_GRBLUE );
   //pCQuadB->draw_me( dp_kCLR_LTBLUE );
        // Refresh debug graphics
        //EraseElements();
        pCQuad0->draw_me( dp_kCLR_YELLOW );
        pCQuad1->draw_me( dp_kCLR_YELLOW );
        pCQuadA->draw_me( dp_kCLR_YELLOW );
        pCQuadB->draw_me( dp_kCLR_YELLOW );
       
   //draw_close();
   }

   delete_this_quad( pCQuad0 );
   delete_this_quad( pCQuad1 );
   delete_this_quad( pCQuadA );
   delete_this_quad( pCQuadB );

   if( delete_this_edge( pCEdge1 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge34 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge4B ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;

   //Need two new nodes
   nodClass *pCNodeE, *pCNodeF;
   two_middle_nodes( pCNode4, pCNodeA, &pCNodeE, &pCNodeF );

   //Build 6 new quads
   elmQClass *pCQuadC, *pCQuadD;
   pCQuad0 = new elmQClass( pCCenterNode, pCNode0, pCNode1, pCNode4,
                                                                pCMaster );
   if( !pCQuad0 ) return fmu_kCLEAN_NO_MEMORY;
   pCQuad1 = new elmQClass( pCNode4, pCNode1, pCNode2, pCNodeE, pCMaster );
   if( !pCQuad1 ) return fmu_kCLEAN_NO_MEMORY;
   pCQuadA = new elmQClass( pCNode4, pCNodeE, pCNodeF, pCNodeD, pCMaster );
   if( !pCQuadA ) return fmu_kCLEAN_NO_MEMORY;
   pCQuadB = new elmQClass( pCNodeE, pCNode2, pCNode3, pCNodeF, pCMaster );
   if( !pCQuadB ) return fmu_kCLEAN_NO_MEMORY;
   pCQuadC = new elmQClass( pCNodeF, pCNode3, pCNodeA, pCNodeB, pCMaster );
   if( !pCQuadC ) return fmu_kCLEAN_NO_MEMORY;
   pCQuadD = new elmQClass( pCNodeD, pCNodeF, pCNodeB, pCNodeC, pCMaster );
   if( !pCQuadD ) return fmu_kCLEAN_NO_MEMORY;

   if( lDrawCleaning )
   {
      nodClassDynArray CSmoothList;
      CSmoothList.append( pCNode0 );
      CSmoothList.append( pCNode1 );
      CSmoothList.append( pCNode2 );
      CSmoothList.append( pCNode3 );
      CSmoothList.append( pCNode4 );
      CSmoothList.append( pCNodeA );
      CSmoothList.append( pCNodeB );
      CSmoothList.append( pCNodeC );
      CSmoothList.append( pCNodeD );
      CSmoothList.append( pCNodeE );
      CSmoothList.append( pCNodeF );
      CSmoothList.append( pCCenterNode );
      smooth_mesh( &CSmoothList );
      //draw_mesh();
      //pCQuad0->draw_me( dp_kCLR_BLUE );
      //pCQuad1->draw_me( dp_kCLR_CYAN );
      //pCQuadA->draw_me( dp_kCLR_GRBLUE );
      //pCQuadB->draw_me( dp_kCLR_LTBLUE );
      //pCQuadC->draw_me( dp_kCLR_GOLD );
      //pCQuadD->draw_me( dp_kCLR_YELLOW );
        // Refresh debug graphics
        //EraseElements();
        pCQuad0->draw_me( dp_kCLR_BLUE );
        pCQuad1->draw_me( dp_kCLR_BLUE );
        pCQuadA->draw_me( dp_kCLR_BLUE );
        pCQuadB->draw_me( dp_kCLR_BLUE );
        pCQuadC->draw_me( dp_kCLR_BLUE );
        pCQuadD->draw_me( dp_kCLR_BLUE );
       
      //draw_close();
   }

   reclean( pCNode0 );
   reclean( pCNode1 );
   reclean( pCNode2 );
   reclean( pCNode3 );
   reclean( pCNodeA );
   reclean( pCCenterNode );

   return fmu_kCLEAN_MESH_CHANGED;
   
}


int QuadCleanTool::zip_unzip( nodClass *pCCenterNode,
                                  nodClassDynArray &CNodes,
                                  fmuEdgeDynArray &CEdges,
                                  elmQClassDynArray &CQuads,
                                  int *aiValence, int *aiHardpt )
{
   
   int iDebugFlag = 0;
   nodClass *pCNode0 = CNodes.get();
   nodClass *pCNode1 = CNodes.next();
   nodClass *pCNode2 = CNodes.next(2);
   nodClass *pCNode3 = CNodes.next(3);
   nodClass *pCNode4 = CNodes.next(4);
   nodClass *pCNode5 = CNodes.next(5);
   nodClass *pCNode6 = CNodes.next(6);
   nodClass *pCNode7 = CNodes.next(7);
   nodClass *pCNode9 = CNodes.next(9);
   nodClass *pCNodeA, *pCNodeB, *pCNodeE;
   fmuEdge  *pCEdge0 = CEdges.get();
   fmuEdge  *pCEdge1 = CEdges.next();
   fmuEdge  *pCEdge2 = CEdges.next(2);
   fmuEdge  *pCEdge4 = CEdges.next(4);
   elmQClass *pCQuad0 = CQuads.get();
   elmQClass *pCQuad1 = CQuads.next();
   elmQClass *pCQuad2 = CQuads.next(2);
   elmQClass *pCQuadA = 0, *pCQuadB = 0, *pCQuadC = 0, *pCQuadD = 0;

   // Check to see if all entites exist - NM
	  if (	pCNode0 == NULL ||
			pCNode1 == NULL ||
			pCNode2 == NULL ||
			pCNode3 == NULL ||
			pCNode4 == NULL ||
			pCNode5 == NULL ||
			pCNode6 == NULL ||
			pCNode7 == NULL ||
			pCNode9 == NULL)
			return fmu_kCLEAN_NO_ACTION;

   // Check to see if all nodes exist  - NM
	if( !pCEdge0 ) return fmu_kCLEAN_NO_ACTION;
	if( !pCEdge1 ) return fmu_kCLEAN_NO_ACTION;
	if( !pCEdge2 ) return fmu_kCLEAN_NO_ACTION;
	if( !pCEdge4 ) return fmu_kCLEAN_NO_ACTION;
 

   // Check to see if all nodes exist - NM
   if( !pCQuad0 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCQuad1 ) return fmu_kCLEAN_NO_ACTION;
	if( !pCQuad2 ) return fmu_kCLEAN_NO_ACTION;


   if( pCEdge1->frozen() ) return fmu_kCLEAN_NO_ACTION;
   if( pCCenterNode->hardpt() ) return fmu_kCLEAN_NO_ACTION;

   if( iDebugFlag )
   {
   //draw_open();
   //pCQuad0->draw_me( dp_kCLR_BLUE );
   //pCQuad1->draw_me( dp_kCLR_CYAN );
   //pCQuad2->draw_me( dp_kCLR_LTBLUE );
        // Refresh debug graphics
        //EraseElements();
        pCQuad0->draw_me( dp_kCLR_YELLOW );
        pCQuad1->draw_me( dp_kCLR_YELLOW );
        pCQuad2->draw_me( dp_kCLR_YELLOW );
       
   //draw_close();
   }

   //Determine if zip or unzip is appropriate

   //Is unzip the thing to do? Are the edges opposite the 3 valent node
   //too long?
   double dRatio2 = edge_ratio( pCEdge2 );
   double dRatio4 = edge_ratio( pCEdge4 );
   if( dRatio2 > 1.6 && dRatio4 > 1.6 )
   {
      if( aiValence[5] != 4 || aiValence[6] != 4 ||  aiValence[7] != 4 )
                                                      return fmu_kCLEAN_NO_ACTION;
      if( aiHardpt[6] ) return fmu_kCLEAN_NO_ACTION;
      if( pCEdge2->frozen() ) return fmu_kCLEAN_NO_ACTION;

      if( lDrawCleaning )
      {
      //draw_open();
      //pCQuad0->draw_me( dp_kCLR_BLUE );
      //pCQuad1->draw_me( dp_kCLR_CYAN );
      //pCQuad2->draw_me( dp_kCLR_LTBLUE );
        // Refresh debug graphics
        //EraseElements();
        pCQuad0->draw_me( dp_kCLR_YELLOW );
        pCQuad1->draw_me( dp_kCLR_YELLOW );
        pCQuad2->draw_me( dp_kCLR_YELLOW );
    

      //draw_close();
      }

      delete_this_quad( pCQuad0 );
      delete_this_quad( pCQuad1 );
      delete_this_quad( pCQuad2 );

      if( delete_this_edge( pCEdge1 ) == fmu_kCLEAN_ABORT )
         return fmu_kCLEAN_ABORT;
      if( delete_this_edge( pCEdge2 ) == fmu_kCLEAN_ABORT )
         return fmu_kCLEAN_ABORT;

      //Five new quads
      two_middle_nodes( pCCenterNode, pCNode4, &pCNodeA, &pCNodeB );
      pCQuad0 = new elmQClass( pCCenterNode, pCNode0, pCNode1, pCNodeA,
                                                                   pCMaster );
      if( !pCQuad0 ) return fmu_kCLEAN_NO_MEMORY;
      pCQuad1 = new elmQClass( pCNodeA, pCNode1, pCNode2, pCNodeB, pCMaster );
      if( !pCQuad1 ) return fmu_kCLEAN_NO_MEMORY;
      pCQuad2 = new elmQClass( pCNodeB, pCNode2, pCNode3, pCNode4, pCMaster );
      if( !pCQuad2 ) return fmu_kCLEAN_NO_MEMORY;
      pCQuadA = new elmQClass( pCCenterNode, pCNodeA, pCNodeB, pCNode6,
                                                                   pCMaster );
      if( !pCQuadA ) return fmu_kCLEAN_NO_MEMORY;
      pCQuadB = new elmQClass( pCNodeB, pCNode4, pCNode5, pCNode6, pCMaster );
      if( !pCQuadB ) return fmu_kCLEAN_NO_MEMORY;

      if( lDrawCleaning )
      {
         nodClassDynArray CSmoothList;
         CSmoothList.append( pCNode0 );
         CSmoothList.append( pCNode1 );
         CSmoothList.append( pCNode2 );
         CSmoothList.append( pCNode3 );
         CSmoothList.append( pCNode4 );
         CSmoothList.append( pCNode5 );
         CSmoothList.append( pCNode6 );
         CSmoothList.append( pCCenterNode );
         CSmoothList.append( pCNodeA );
         CSmoothList.append( pCNodeB );
         smooth_mesh( &CSmoothList );
         //draw_mesh();
         //pCQuad0->draw_me( dp_kCLR_BLUE );
         //pCQuad1->draw_me( dp_kCLR_CYAN );
         //pCQuad2->draw_me( dp_kCLR_LTBLUE );
         //pCQuadA->draw_me( dp_kCLR_GOLD );
         //pCQuadB->draw_me( dp_kCLR_RED );
        // Refresh debug graphics
        //EraseElements();
        pCQuad0->draw_me( dp_kCLR_BLUE );
        pCQuad1->draw_me( dp_kCLR_BLUE );
        pCQuad2->draw_me( dp_kCLR_BLUE );
        pCQuadA->draw_me( dp_kCLR_BLUE );
        pCQuadB->draw_me( dp_kCLR_BLUE );
        
         //draw_close();
      }
      reclean( pCNode0 );
      reclean( pCNode1 );
      reclean( pCNode2 );
      reclean( pCNode3 );
      reclean( pCNode4 );
      reclean( pCNode5 );
      reclean( pCNode6 );
      reclean( pCCenterNode );

      return fmu_kCLEAN_MESH_CHANGED;
   }

   //Unzip is not appropriate, try zip instead. Reject as we go as appropriate
   if( pCEdge0->frozen() ) return fmu_kCLEAN_NO_ACTION;
   fmuEdge *pCEdge01 = pCNode0->shared_edge( pCNode1 );
	if( !pCEdge01 ) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge01->frozen() ) return fmu_kCLEAN_NO_ACTION;

   fmuEdge *pCEdge12 = pCNode1->shared_edge( pCNode2 );
	if( !pCEdge12 ) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge12->frozen() ) return fmu_kCLEAN_NO_ACTION;

   fmuEdge *pCEdge23 = pCNode2->shared_edge( pCNode3 );
	if( !pCEdge23 ) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge23->frozen() ) return fmu_kCLEAN_NO_ACTION;

   fmuEdge *pCEdge09 = pCNode0->shared_edge( pCNode9 );
	if( !pCEdge09 ) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge09->frozen() ) return fmu_kCLEAN_NO_ACTION;

   pCQuadA = pCEdge01->other_quad( pCQuad0 );
   pCQuadB = pCEdge12->other_quad( pCQuad0 );

	if( !pCQuadA ) return fmu_kCLEAN_NO_ACTION;
   pCNodeB = pCQuadA->opposite_node( pCNode1 );
	if( !pCNodeB ) return fmu_kCLEAN_NO_ACTION;
   if( pCNodeB->hardpt() ) return fmu_kCLEAN_NO_ACTION;
   if( get_valence( pCNodeB ) != 4 ) return fmu_kCLEAN_NO_ACTION;

   nodClass *pCNodeC = pCQuadA->prev_node( pCNode1 );
	if( !pCNodeC ) return fmu_kCLEAN_NO_ACTION;
   if( pCNodeC->hardpt() ) return fmu_kCLEAN_NO_ACTION;
   if( get_valence( pCNodeC ) != 4 ) return fmu_kCLEAN_NO_ACTION;

   nodClass *pCNodeD = pCQuadB->opposite_node( pCNode1 );
	if( !pCNodeD ) return fmu_kCLEAN_NO_ACTION;
   if( pCNodeD->hardpt() ) return fmu_kCLEAN_NO_ACTION;
   if( get_valence( pCNodeD ) != 4 ) return fmu_kCLEAN_NO_ACTION;

   fmuEdge *pCEdge0B = pCNode0->shared_edge( pCNodeB );
	if( !pCEdge0B ) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge0B->frozen() ) return fmu_kCLEAN_NO_ACTION;

   fmuEdge *pCEdge1C = pCNode1->shared_edge( pCNodeC );
	if( !pCEdge1C ) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge1C->frozen() ) return fmu_kCLEAN_NO_ACTION;

   fmuEdge *pCEdge2D = pCNode2->shared_edge( pCNodeD );
	if( !pCEdge2D ) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge2D->frozen() ) return fmu_kCLEAN_NO_ACTION;

   pCQuadC = pCEdge2D->other_quad( pCQuadB );
   pCQuadD = pCEdge0B->other_quad( pCQuadA );

	if( !pCQuadD ) return fmu_kCLEAN_NO_ACTION;
   pCNodeA = pCQuadD->prev_node( pCNodeB );
   if( get_valence( pCNodeA ) != 4 ) return fmu_kCLEAN_NO_ACTION;

	if( !pCQuadC ) return fmu_kCLEAN_NO_ACTION;
   pCNodeE = pCQuadC->next_node( pCNodeD );
   if( get_valence( pCNodeE ) != 4 ) return fmu_kCLEAN_NO_ACTION;

   if( iDebugFlag )
   {
   //draw_open();
   //pCQuadA->draw_me( dp_kCLR_BLUE );
   //pCQuadB->draw_me( dp_kCLR_GOLD );
   //pCQuadC->draw_me( dp_kCLR_RED );
   //pCQuadD->draw_me( dp_kCLR_YELLOW );
        // Refresh debug graphics
        //EraseElements();
        pCQuadA->draw_me( dp_kCLR_YELLOW );
        pCQuadB->draw_me( dp_kCLR_YELLOW );
        pCQuadC->draw_me( dp_kCLR_YELLOW );
        pCQuadD->draw_me( dp_kCLR_YELLOW );
   //draw_close();
   }

   double dRatio01 = edge_ratio( pCEdge01 );
   double dRatio12 = edge_ratio( pCEdge12 );
   double dRatio23 = edge_ratio( pCEdge23 );
   double dRatio09 = edge_ratio( pCEdge09 );

   if( dRatio01 > 0.55 && dRatio12 > 0.55 ) return fmu_kCLEAN_NO_ACTION;

   if( lDrawCleaning )
   {
   //draw_open();
   //pCQuad0->draw_me( dp_kCLR_BLUE );
   //pCQuad1->draw_me( dp_kCLR_CYAN );
   //pCQuadA->draw_me( dp_kCLR_LTBLUE );
   //pCQuadB->draw_me( dp_kCLR_GOLD );
   //pCQuadC->draw_me( dp_kCLR_RED );

        // Refresh debug graphics
        //EraseElements();
        pCQuad0->draw_me( dp_kCLR_YELLOW );
        pCQuad1->draw_me( dp_kCLR_YELLOW );
        pCQuadA->draw_me( dp_kCLR_YELLOW );
        pCQuadB->draw_me( dp_kCLR_YELLOW );
        pCQuadC->draw_me( dp_kCLR_YELLOW );
       
   //draw_close();
   }

   delete_this_quad( pCQuad0 );
   delete_this_quad( pCQuad1 );
   delete_this_quad( pCQuadA );
   delete_this_quad( pCQuadB );
   delete_this_quad( pCQuadC );

   if( delete_this_edge( pCEdge1 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge01 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge12 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge23 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge1C ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge2D ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;

   if( delete_this_node( pCNode1 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_node( pCNode2 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;

   pCQuad0 = new elmQClass( pCCenterNode, pCNode0, pCNode3, pCNode4,
                                                                pCMaster );
   if( !pCQuad0 ) return fmu_kCLEAN_NO_MEMORY;
   pCQuad1 = new elmQClass( pCNode0, pCNodeB, pCNodeC, pCNodeD, pCMaster );
   if( !pCQuad1 ) return fmu_kCLEAN_NO_MEMORY;
   pCQuadA = new elmQClass( pCNode0, pCNodeD, pCNodeE, pCNode3, pCMaster );
   if( !pCQuadA ) return fmu_kCLEAN_NO_MEMORY;

   if( lDrawCleaning )
   {
      nodClassDynArray CSmoothList;
      CSmoothList.append( pCNode0 );
      CSmoothList.append( pCNode3 );
      CSmoothList.append( pCNode4 );
      CSmoothList.append( pCNodeB );
      CSmoothList.append( pCNodeC );
      CSmoothList.append( pCNodeD );
      CSmoothList.append( pCCenterNode );
      smooth_mesh( &CSmoothList );
      //draw_mesh();
      //pCQuad0->draw_me( dp_kCLR_BLUE );
      //pCQuad1->draw_me( dp_kCLR_CYAN );
      //pCQuadA->draw_me( dp_kCLR_LTBLUE );
        // Refresh debug graphics
        //EraseElements();
        pCQuad0->draw_me( dp_kCLR_BLUE );
        pCQuad1->draw_me( dp_kCLR_BLUE );
        pCQuadA->draw_me( dp_kCLR_BLUE );
       
      //draw_close();
   }
   reclean( pCNode0 );
   reclean( pCNode3 );
   reclean( pCNodeB );
   reclean( pCNodeC );
   reclean( pCNodeD );

   return fmu_kCLEAN_MESH_CHANGED;
}


int QuadCleanTool::clean_554443( nodClass *pCCenterNode,
                                     nodClassDynArray &CNodes,
                                     fmuEdgeDynArray &CEdges,
                                     elmQClassDynArray &CQuads )
{
   
   int iSide = 0;
   nodClass *pCNode0 = CNodes.get();
   nodClass *pCNode1 = CNodes.next();
   nodClass *pCNode2 = CNodes.next(2);
   nodClass *pCNode3 = CNodes.next(3);
   nodClass *pCNode4 = CNodes.next(4);
   fmuEdge  *pCEdge1 = CEdges.next();
   elmQClass *pCQuad0 = CQuads.get();
   elmQClass *pCQuad1 = CQuads.next();


   // Check to see if all entites exist - NM
	  if (	pCNode0 == NULL ||
			pCNode1 == NULL ||
			pCNode2 == NULL ||
			pCNode3 == NULL ||
			pCNode4 == NULL)
			return fmu_kCLEAN_NO_ACTION;

   // Check to see if all edges exist  - NM
	if( !pCEdge1 ) return fmu_kCLEAN_NO_ACTION;
 

   // Check to see if all quads exist - NM
   if( !pCQuad0 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCQuad1 ) return fmu_kCLEAN_NO_ACTION;
	
   if( pCEdge1->frozen() ) return fmu_kCLEAN_NO_ACTION;
   fmuEdge *pCEdge01 = pCNode0->shared_edge( pCNode1 );
	if( !pCEdge01 ) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge01->frozen() ) return fmu_kCLEAN_NO_ACTION;

   elmQClass *pCQuadA = pCEdge01->other_quad( pCQuad0 );
	if( !pCQuadA ) return fmu_kCLEAN_NO_ACTION;
   nodClass *pCNodeA = pCQuadA->next_node( pCNode0 );
	if( !pCNodeA ) return fmu_kCLEAN_NO_ACTION;
   if( get_valence( pCNodeA ) > 4 ) return fmu_kCLEAN_NO_ACTION;
   nodClass *pCNodeB = pCQuadA->opposite_node( pCNode0 );
	if( !pCNodeB ) return fmu_kCLEAN_NO_ACTION;

   if( pCNode0->hardpt() ) iSide = 1;
   else if( pCCenterNode->hardpt() ) iSide = 2;
   else
   {
      if( pCNodeB->hardpt() && pCNode1->hardpt() && pCNode2->hardpt() )
      {
         fmuVector CVec12( pCNode1, pCNode2 );
         fmuVector CVec1B( pCNode1, pCNodeB );
         if( CNormal.angle( CVec12, CVec1B ) > fmu_kHIGH_PI_TOL ) iSide = 1;
      }
      if( pCNode1->hardpt() && pCNode2->hardpt() && pCNode3->hardpt() )
      {
         fmuVector CVec23( pCNode2, pCNode3 );
         fmuVector CVec21( pCNode2, pCNode1 );
         if( CNormal.angle( CVec23, CVec21 ) > fmu_kHIGH_PI_TOL ) iSide = 2;
      }
   }

   if( lDrawCleaning )
   {
   //draw_open();
   //pCQuad0->draw_me( dp_kCLR_PINK );
   //pCQuad1->draw_me( dp_kCLR_CYAN );
   //pCQuadA->draw_me( dp_kCLR_GRBLUE );
        // Refresh debug graphics
        //EraseElements();
        pCQuad0->draw_me( dp_kCLR_YELLOW );
        pCQuad1->draw_me( dp_kCLR_YELLOW );
        pCQuadA->draw_me( dp_kCLR_YELLOW );
      
   //draw_close();
   }

   delete_this_quad( pCQuad0 );
   delete_this_quad( pCQuad1 );
   delete_this_quad( pCQuadA );

   if( delete_this_edge( pCEdge1 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge01 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;

   //Need two new nodes
   nodClass *pCNodeC, *pCNodeD;
   two_middle_nodes( pCNodeA, pCNode4, &pCNodeC, &pCNodeD );

   //Build 5 new quads
   elmQClass *pCQuadB, *pCQuadC;
   pCQuad0 = new elmQClass( pCNodeA, pCNodeB, pCNode1, pCNodeC, pCMaster );
   if( !pCQuad0 ) return fmu_kCLEAN_NO_MEMORY;
   pCQuad1 = new elmQClass( pCNodeC, pCNode1, pCNode2, pCNodeD, pCMaster );
   if( !pCQuad1 ) return fmu_kCLEAN_NO_MEMORY;
   pCQuadA = new elmQClass( pCNodeD, pCNode2, pCNode3, pCNode4, pCMaster );
   if( !pCQuadA ) return fmu_kCLEAN_NO_MEMORY;

   if( iSide == 1 )
   {
      pCQuadB = new elmQClass( pCCenterNode, pCNode0, pCNodeA, pCNodeC,
                                                               pCMaster );
      if( !pCQuadB ) return fmu_kCLEAN_NO_MEMORY;
      pCQuadC = new elmQClass( pCCenterNode, pCNodeC, pCNodeD, pCNode4,
                                                               pCMaster );
      if( !pCQuadC ) return fmu_kCLEAN_NO_MEMORY;
   }
   else
   {
      pCQuadB = new elmQClass( pCNode0, pCNodeA, pCNodeC, pCNodeD, pCMaster );
      if( !pCQuadB ) return fmu_kCLEAN_NO_MEMORY;
      pCQuadC = new elmQClass( pCNode0, pCNodeD, pCNode4, pCCenterNode,
                                                                   pCMaster );
      if( !pCQuadC ) return fmu_kCLEAN_NO_MEMORY;
   }

   if( lDrawCleaning )
   {
      nodClassDynArray CSmoothList;
      CSmoothList.append( pCNodeC );
      CSmoothList.append( pCNodeD );
      CSmoothList.append( pCNode4 );
      CSmoothList.append( pCNode0 );
      CSmoothList.append( pCNodeA );
      CSmoothList.append( pCCenterNode );
      smooth_mesh( &CSmoothList );
      //draw_mesh();
      //pCQuad0->draw_me( dp_kCLR_PINK );
      //pCQuad1->draw_me( dp_kCLR_CYAN );
      //pCQuadA->draw_me( dp_kCLR_GRBLUE );
      //pCQuadB->draw_me( dp_kCLR_LTBLUE );
      //pCQuadC->draw_me( dp_kCLR_GOLD );
        // Refresh debug graphics
        //EraseElements();
        pCQuad0->draw_me( dp_kCLR_BLUE );
        pCQuad1->draw_me( dp_kCLR_BLUE );
        pCQuadA->draw_me( dp_kCLR_BLUE );
        pCQuadB->draw_me( dp_kCLR_BLUE );
        pCQuadC->draw_me( dp_kCLR_BLUE );
       
      //draw_close();
   }

   reclean( pCNodeC );
   reclean( pCNodeD );
   reclean( pCNode4 );
   reclean( pCNode0 );
   reclean( pCNodeA );
   reclean( pCCenterNode );

   return fmu_kCLEAN_MESH_CHANGED;
   
}


int QuadCleanTool::clean_534445( nodClass *pCCenterNode,
                                     nodClassDynArray &CNodes,
                                     fmuEdgeDynArray &CEdges,
                                     elmQClassDynArray &CQuads )
{
   
   int iSide = 0;
   nodClass *pCNode0 = CNodes.get();
   nodClass *pCNode1 = CNodes.next();
   nodClass *pCNode2 = CNodes.next(2);
   nodClass *pCNode3 = CNodes.next(3);
   nodClass *pCNode4 = CNodes.next(4);
   fmuEdge  *pCEdge1 = CEdges.next();
   elmQClass *pCQuad0 = CQuads.get();
   elmQClass *pCQuad1 = CQuads.next();

  // Check to see if all entites exist - NM
	  if (	pCNode0 == NULL ||
			pCNode1 == NULL ||
			pCNode2 == NULL ||
			pCNode3 == NULL ||
			pCNode4 == NULL)
			return fmu_kCLEAN_NO_ACTION;

   // Check to see if all edges exist  - NM
	if( !pCEdge1 ) return fmu_kCLEAN_NO_ACTION;
 

   // Check to see if all quads exist - NM
   if( !pCQuad0 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCQuad1 ) return fmu_kCLEAN_NO_ACTION;

   if( pCEdge1->frozen() ) return fmu_kCLEAN_NO_ACTION;
   fmuEdge *pCEdge34 = pCNode3->shared_edge( pCNode4 );
	if( !pCEdge34 ) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge34->frozen() ) return fmu_kCLEAN_NO_ACTION;

   elmQClass *pCQuadA = pCEdge34->other_quad( pCQuad1 );
	if( !pCQuadA ) return fmu_kCLEAN_NO_ACTION;
   nodClass *pCNodeA = pCQuadA->next_node( pCNode3 );
	if( !pCNodeA ) return fmu_kCLEAN_NO_ACTION;
   nodClass *pCNodeB = pCQuadA->opposite_node( pCNode3 );
	if( !pCNodeB ) return fmu_kCLEAN_NO_ACTION;
   if( get_valence( pCNodeB ) > 4 ) return fmu_kCLEAN_NO_ACTION;

   if( pCNode4->hardpt() ) iSide = 1;
   else if( pCCenterNode->hardpt() ) iSide = 2;
   else
   {
      if( pCNode2->hardpt() && pCNode3->hardpt() && pCNodeA->hardpt() )
      {
         fmuVector CVec3A( pCNode3, pCNodeA );
         fmuVector CVec32( pCNode3, pCNode2 );
         if( CNormal.angle( CVec3A, CVec32 ) > fmu_kHIGH_PI_TOL ) iSide = 1;
      }
      if( pCNode1->hardpt() && pCNode2->hardpt() && pCNode3->hardpt() )
      {
         fmuVector CVec23( pCNode2, pCNode3 );
         fmuVector CVec21( pCNode2, pCNode1 );
         if( CNormal.angle( CVec23, CVec21 ) > fmu_kHIGH_PI_TOL ) iSide = 2;
      }
   }

   if( lDrawCleaning )
   {
   //draw_open();
   //pCQuad0->draw_me( dp_kCLR_BLUE );
   //pCQuad1->draw_me( dp_kCLR_RED );
   //pCQuadA->draw_me( dp_kCLR_GRBLUE );
       // Refresh debug graphics
        //EraseElements();
        pCQuad0->draw_me( dp_kCLR_YELLOW );
        pCQuad1->draw_me( dp_kCLR_YELLOW );
        pCQuadA->draw_me( dp_kCLR_YELLOW );
      
   //draw_close();
   }

   delete_this_quad( pCQuad0 );
   delete_this_quad( pCQuad1 );
   delete_this_quad( pCQuadA );

   if( delete_this_edge( pCEdge1 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge34 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;

   //Need two new nodes
   nodClass *pCNodeC, *pCNodeD;
   two_middle_nodes( pCNode0, pCNodeB, &pCNodeC, &pCNodeD );

   //Build 5 new quads
   elmQClass *pCQuadB, *pCQuadC;
   pCQuad0 = new elmQClass( pCNode0, pCNode1, pCNode2, pCNodeC, pCMaster );
   if( !pCQuad0 ) return fmu_kCLEAN_NO_MEMORY;
   pCQuad1 = new elmQClass( pCNodeC, pCNode2, pCNode3, pCNodeD, pCMaster );
   if( !pCQuad1 ) return fmu_kCLEAN_NO_MEMORY;
   pCQuadA = new elmQClass( pCNodeD, pCNode3, pCNodeA, pCNodeB, pCMaster );
   if( !pCQuadA ) return fmu_kCLEAN_NO_MEMORY;

   if( iSide == 1 )
   {
      pCQuadB = new elmQClass( pCCenterNode, pCNode0, pCNodeC, pCNodeD,
                                                               pCMaster );
      if( !pCQuadB ) return fmu_kCLEAN_NO_MEMORY;
      pCQuadC = new elmQClass( pCCenterNode, pCNodeD, pCNodeB, pCNode4,
                                                               pCMaster );
      if( !pCQuadC ) return fmu_kCLEAN_NO_MEMORY;
   }
   else
   {
      pCQuadB = new elmQClass( pCNode4, pCNodeC, pCNodeD, pCNodeB, pCMaster );
      if( !pCQuadB ) return fmu_kCLEAN_NO_MEMORY;
      pCQuadC = new elmQClass( pCNode4, pCCenterNode, pCNode0, pCNodeC,
                                                                   pCMaster );
      if( !pCQuadC ) return fmu_kCLEAN_NO_MEMORY;
   }
   if( lDrawCleaning )
   {
      nodClassDynArray CSmoothList;
      CSmoothList.append( pCNodeC );
      CSmoothList.append( pCNodeD );
      CSmoothList.append( pCNode0 );
      CSmoothList.append( pCNodeB );
      CSmoothList.append( pCNode4 );
      CSmoothList.append( pCCenterNode );
      smooth_mesh( &CSmoothList );
      //draw_mesh();
      //pCQuad0->draw_me( dp_kCLR_BLUE );
      //pCQuad1->draw_me( dp_kCLR_RED );
      //pCQuadA->draw_me( dp_kCLR_GRBLUE );
      //pCQuadB->draw_me( dp_kCLR_LTBLUE );
      //pCQuadC->draw_me( dp_kCLR_GOLD );
        // Refresh debug graphics
        //EraseElements();
        pCQuad0->draw_me( dp_kCLR_BLUE );
        pCQuad1->draw_me( dp_kCLR_BLUE );
        pCQuadA->draw_me( dp_kCLR_BLUE );
        pCQuadB->draw_me( dp_kCLR_BLUE );
        pCQuadC->draw_me( dp_kCLR_BLUE );
        
      //draw_close();
   }

   reclean( pCNodeC );
   reclean( pCNodeD );
   reclean( pCNode4 );
   reclean( pCNode0 );
   reclean( pCNodeB );
   reclean( pCCenterNode );

   return fmu_kCLEAN_MESH_CHANGED;
  
}


int QuadCleanTool::rotate_edge( nodClass *pCCenterNode,
                                    nodClassDynArray &CNodes,
                                    fmuEdgeDynArray &CEdges,
                                    elmQClassDynArray &CQuads,
                                    int direction )
{
   
   int iCallStatus;
   nodClassDynArray CDeleteNodes(ARRAY_SIZE);
   nodClass *pCNode0 = CNodes.get();
   nodClass *pCNode1 = CNodes.next();
   nodClass *pCNode2 = CNodes.next(2);
   nodClass *pCNode3 = CNodes.next(3);
   nodClass *pCNode4 = CNodes.next(4);
   fmuEdge  *pCEdge1 = CEdges.next();
   elmQClass *pCQuad0 = CQuads.get();
   elmQClass *pCQuad1 = CQuads.next();

  // Check to see if all entites exist - NM
	  if (	pCNode0 == NULL ||
			pCNode1 == NULL ||
			pCNode2 == NULL ||
			pCNode3 == NULL ||
			pCNode4 == NULL)
			return fmu_kCLEAN_NO_ACTION;

   // Check to see if all edges exist  - NM
	if( !pCEdge1 ) return fmu_kCLEAN_NO_ACTION;
 

   // Check to see if all quads exist - NM
   if( !pCQuad0 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCQuad1 ) return fmu_kCLEAN_NO_ACTION;


   if( pCEdge1->frozen() ) return fmu_kCLEAN_NO_ACTION;

   //Only two quads and one edge are to be deleted, reset the lists
   CQuads.clear();
   CQuads[0] = pCQuad0;
   CQuads[1] = pCQuad1;
   CEdges.clear();
   CEdges[0] = pCEdge1;

   //No nodes to delete, there are only 6 nodes in the ring. Current
   //position depends on the direction of the new cut. Hint: nodes.get()
   //is not necessarily nodes[0].
   CNodes.clear();
   CNodes[0] = pCNode0;
   CNodes[1] = pCNode1;
   CNodes[2] = pCNode2;
   CNodes[3] = pCNode3;
   CNodes[4] = pCNode4;
   CNodes[5] = pCCenterNode;
   CNodes.set( 0 );
   if( direction == 1 ) CNodes.step();

   iCallStatus = fill_two( CDeleteNodes, CEdges, CQuads, CNodes );

   return iCallStatus;
  
}


int QuadCleanTool::across_open_face(    nodClass *pCCenterNode,
                                        nodClassDynArray &CNodes,
                                        fmuEdgeDynArray &CEdges,
                                        elmQClassDynArray &CQuads )
{
    nodClass *pCNode0 = CNodes.get();
    nodClass *pCNode1 = CNodes.next();
    nodClass *pCNode2 = CNodes.next(2);
    nodClass *pCNode3 = CNodes.next(3);
    nodClass *pCNode4 = CNodes.next(4);

    // Check to see if all entites exist - NM
    if( pCNode0 == NULL ||
        pCNode1 == NULL ||
        pCNode2 == NULL ||
        pCNode3 == NULL ||
        pCNode4 == NULL)
        return fmu_kCLEAN_NO_ACTION;

    // Check to see if all edges exist  - NM
    fmuEdge  *pCEdge1 = CEdges.next();

    if( !pCEdge1 ) 
        return fmu_kCLEAN_NO_ACTION;
    // or are frozen
    if( pCEdge1->frozen() )
        return fmu_kCLEAN_NO_ACTION;

    // Check to see if all quads exist - NM
    elmQClass *pCQuad0 = CQuads.get();
    elmQClass *pCQuad1 = CQuads.next();

    if( !pCQuad0 ) 
        return fmu_kCLEAN_NO_ACTION;
    if( !pCQuad1 ) 
        return fmu_kCLEAN_NO_ACTION;

    //Only two quads and one edge are to be deleted, reset the lists
    CQuads.clear();
    CQuads.resize( 2 );
    CQuads[0] = pCQuad0;
    CQuads[1] = pCQuad1;

    CEdges.clear();
    CEdges.resize( 1 );
    CEdges[0] = pCEdge1;

    //No nodes to delete, there are only 6 nodes in the ring. Current
    //position depends on the direction of the new cut. Hint: nodes.get()
    //is not necessarily nodes[0].
    CNodes.clear();
    CNodes.resize( 6 );
    CNodes[0] = pCNode0;
    CNodes[1] = pCNode1;
    CNodes[2] = pCNode2;
    CNodes[3] = pCNode3;
    CNodes[4] = pCNode4;
    CNodes[5] = pCCenterNode;
    CNodes.set( 0 );

    nodClassDynArray CDeleteNodes( 0 );
    int iCallStatus = fill_three( CDeleteNodes, CEdges, CQuads, CNodes );

    return iCallStatus;
}


int QuadCleanTool::crack_open_middle( nodClass *pCCenterNode,
                                          nodClassDynArray &CNodes,
                                          fmuEdgeDynArray &CEdges,
                                          elmQClassDynArray &CQuads )
{
   
   double dQ=0.0; 
   double dNewX=0.0, dNewY=0.0, dNewX1=0.0, dNewY1=0.0;
   bool bLockNode = false;

   nodClass *pCNode0 = CNodes.get();
   nodClass *pCNode1 = CNodes.next();
   nodClass *pCNode2 = CNodes.next(2);
   nodClass *pCNode3 = CNodes.next(3);
   nodClass *pCNode4 = CNodes.next(4);
   nodClass *pCNode5 = CNodes.next(5);
   nodClass *pCNode6 = CNodes.next(6);
   fmuEdge  *pCEdge1 = CEdges.next();
   fmuEdge  *pCEdge2 = CEdges.next(2);
   elmQClass *pCQuad0 = CQuads.get();
   elmQClass *pCQuad1 = CQuads.next();
   elmQClass *pCQuad2 = CQuads.next(2);

  // Check to see if all entites exist - NM
	  if (	pCNode0 == NULL ||
			pCNode1 == NULL ||
			pCNode2 == NULL ||
			pCNode3 == NULL ||
			pCNode4 == NULL ||
			pCNode5 == NULL ||
			pCNode6 == NULL)
			return fmu_kCLEAN_NO_ACTION;

   // Check to see if all edges exist  - NM
	if( !pCEdge1 ) return fmu_kCLEAN_NO_ACTION;
    if( !pCEdge2 ) return fmu_kCLEAN_NO_ACTION;

   // Check to see if all quads exist - NM
   if( !pCQuad0 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCQuad1 ) return fmu_kCLEAN_NO_ACTION;
   if( !pCQuad2 ) return fmu_kCLEAN_NO_ACTION;

   if( pCEdge1->frozen() ) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge2->frozen() ) return fmu_kCLEAN_NO_ACTION;

   if( lDrawCleaning )
   {
   //draw_open();
   //pCQuad0->draw_me( dp_kCLR_RED );
   //pCQuad1->draw_me( dp_kCLR_YELLOW );
   //pCQuad2->draw_me( dp_kCLR_BLUE );
    // Refresh debug graphics
    //EraseElements();
    pCQuad0->draw_me( dp_kCLR_YELLOW );
    pCQuad1->draw_me( dp_kCLR_YELLOW );
    pCQuad2->draw_me( dp_kCLR_YELLOW );
   //draw_close();
   }

   delete_this_quad( pCQuad0 );
   delete_this_quad( pCQuad1 );
   delete_this_quad( pCQuad2 );

   if( delete_this_edge( pCEdge1 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge2 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;

   //Need a new node in the middle
   nodClass *pCNewNode = one_middle_node( pCCenterNode, pCNode3 );
   dNewX = pCNewNode->node_x();
   dNewY = pCNewNode->node_y();

   //Four new quads
   pCQuad0 = new elmQClass( pCNewNode, pCNode0, pCNode1, pCNode2, pCMaster );
   if( !pCQuad0 ) return fmu_kCLEAN_NO_MEMORY;
   pCQuad1 = new elmQClass( pCNewNode, pCNode2, pCNode3, pCNode4, pCMaster );
   if( !pCQuad1 ) return fmu_kCLEAN_NO_MEMORY;
   pCQuad2 = new elmQClass( pCNewNode, pCNode4, pCNode5, pCNode6, pCMaster );
   if( !pCQuad2 ) return fmu_kCLEAN_NO_MEMORY;
   elmQClass *pCQuadA = new elmQClass( pCNewNode, pCNode6, pCCenterNode, pCNode0, pCMaster );
   if( !pCQuadA ) return fmu_kCLEAN_NO_MEMORY;

   if( lDrawCleaning )
   {
      nodClassDynArray CSmoothList;
      CSmoothList.append( pCCenterNode );
      CSmoothList.append( pCNewNode );
      smooth_mesh( &CSmoothList );
      //draw_mesh();
      //pCQuad0->draw_me( dp_kCLR_RED );
      //pCQuad1->draw_me( dp_kCLR_YELLOW );
      //pCQuad2->draw_me( dp_kCLR_BLUE );
      //pCQuadA->draw_me( dp_kCLR_LTBLUE );
        // Refresh debug graphics
        //EraseElements();
        pCQuad0->draw_me( dp_kCLR_BLUE );
        pCQuad1->draw_me( dp_kCLR_BLUE );
        pCQuad2->draw_me( dp_kCLR_BLUE );
        pCQuadA->draw_me( dp_kCLR_BLUE );
      //draw_close();
   }

   // Check Quad quality
   dQ = test_quad_quality(pCQuad0);
   if (dQ < 0.1)
       bLockNode = true;

   dQ = test_quad_quality(pCQuad1);
   if (dQ < 0.1)
       bLockNode = true;

   dQ = test_quad_quality(pCQuad2);
   if (dQ < 0.1)
       bLockNode = true;

   dQ = test_quad_quality(pCQuadA);
   if (dQ < 0.1)
       bLockNode = true;

   if (bLockNode)
   {
       pCNewNode->set(dNewX, dNewY, 0.0);
       pCNewNode->lock_node();
   }
   reclean( pCCenterNode );
   reclean( pCNewNode );

   return fmu_kCLEAN_MESH_CHANGED;
  
}

int QuadCleanTool::boundary_clean()
{
    int iDebugFlag = 0;
    int iStatus = fmu_kCLEAN_NO_ACTION;
    int iCallStatus;
    int ii, jj;
    int iBoundaryValence, iSurfaceValence, iFullValence;
    double dAngle;

    nodClass *pCNode;
    nodClassDynArray CBackstepNodes;
    nodClassDynArray CNodes;
    fmuEdgeDynArray CEdges;
    elmQClassDynArray CQuads;

    nodClassDynArrayDynArray *boundaries = pCMaster->global_boundary();
    for( ii = 0; ii < boundaries->length(); ++ii )
    {
        nodClassDynArray *boundary = (*boundaries)[ii];

        //Hardpoint or "anchor" nodes were single node loops
        //Nodes in a weldline were changed to interior, so
        //the loop is empty. Skip these cases.
      if( boundary->empty() )
            continue;

        for( jj = 0; jj < boundary->length(); ++jj )
        {
            boundary->set( jj );
            pCNode = boundary->get();
            if( iDebugFlag )
            {
                //draw_mesh();
                //draw_close();
                //pCPrevNode = boundary->prev();
                //pCNextNode = boundary->next();
                //pCNode->draw_me( dp_kCLR_MAGNTA );
                //pCPrevNode->draw_me( dp_kCLR_YELLOW );
                //pCNextNode->draw_me( dp_kCLR_LTBLUE );
            }

            CNodes.clear();
            CEdges.clear();
            CQuads.clear();

            //iBoundaryValence = number of quads on this surface at this node
            iBoundaryValence = boundary_neighbors(  boundary,
                                                    CNodes, CEdges, CQuads );

            if( iBoundaryValence == 0) 
            {
                boundary->step();
                continue;
            }
            if( iBoundaryValence == -1)
                return fmu_kCLEAN_ABORT;

            //iSurfaceValence = number of edges from node excluding boundary
            iSurfaceValence = iBoundaryValence - 1;

            //iFullValence = virtual boundary valence of node
            iFullValence = get_valence( pCNode );

            dAngle = angle_on_boundary( boundary );
            if( dAngle > fmu_kBNDY_ANGLE_TOL && iSurfaceValence == 0 )
            {
                //Handle flat angles - an angle that needs an edge out of it
                iCallStatus = bndy_flat_clean( pCNode, CQuads[0], boundary );
                if( iCallStatus == fmu_kCLEAN_NO_MEMORY )
                    return fmu_kCLEAN_NO_MEMORY;
                if( iCallStatus == fmu_kCLEAN_ABORT )
                    return fmu_kCLEAN_ABORT;
                if( iCallStatus == fmu_kCLEAN_MESH_CHANGED )
                    iStatus = iCallStatus;
            }

            else
            {
                //Check for strange valence patterns and flange patterns
                iCallStatus = bndy_valence_clean(   pCNode,iFullValence, 
                                                    CNodes, CEdges, CQuads );
                if( iCallStatus == fmu_kCLEAN_NO_MEMORY )
                    return fmu_kCLEAN_NO_MEMORY;
                if( iCallStatus == fmu_kCLEAN_ABORT )
                    return fmu_kCLEAN_ABORT;

                //An operation needs to propagate against the direction of testing
                //but don't fall into infinite loops.
                if( iCallStatus == fmu_kCLEAN_BACKSTEP )
                {
                    if( CBackstepNodes.find( pCNode ) == -1 )
                    {
                        CBackstepNodes.append( pCNode );
                        if( jj > 4 ) 
                            jj -= 3;
                        else
                            jj = -1;
                    }
                    iStatus = fmu_kCLEAN_MESH_CHANGED;
                }
                if( iCallStatus == fmu_kCLEAN_MESH_CHANGED )
                    iStatus = iCallStatus;
            }

            boundary->step();
        }
    }

    return iStatus;
}


int QuadCleanTool::bndy_flat_clean( nodClass *pCCenterNode,
                                    elmQClass *pCQuad,
                                    nodClassDynArray *pCBoundary )
{
   //This routine assumes that there is only one quad at this node
   //on the surface.

    if( pCQuad == NULL )
        return fmu_kCLEAN_NO_ACTION;

   nodClass *pCNextNode = pCQuad->next_node( pCCenterNode );
   nodClass *pCPrevNode = pCQuad->prev_node( pCCenterNode );
   nodClass *pCOppNode  = pCQuad->opposite_node( pCCenterNode );

	if (!pCNextNode)
	   return fmu_kCLEAN_NO_ACTION;
	if (!pCPrevNode)
	   return fmu_kCLEAN_NO_ACTION;
	if (!pCOppNode)
	   return fmu_kCLEAN_NO_ACTION;

   int iCallStatus;
   nodClassDynArray CNodes;
   fmuEdgeDynArray CEdges;
   elmQClassDynArray CQuads;

   //Check first for extreme aspect ratio in this triagular quad as the fix
   //overrides fixes based on valence
   fmuVector CVecPO( pCPrevNode, pCOppNode );
   fmuVector CVecNO( pCNextNode, pCOppNode );
   if (CVecPO.length() <= kdfToler)
	    return fmu_kCLEAN_NO_ACTION;

   double dRatio = CVecNO.length() / CVecPO.length();
   if( dRatio > 2.0 )
   {
      iCallStatus = shape_rotate_next( pCQuad, pCCenterNode, pCNextNode,
                                       pCOppNode, pCPrevNode, false );
      if( iCallStatus != fmu_kCLEAN_NO_ACTION )
            return iCallStatus;
   }
   else if( dRatio < 0.5 )
   {
      iCallStatus = shape_rotate_prev( pCQuad, pCCenterNode, pCNextNode,
                                       pCOppNode, pCPrevNode, false );
      if( iCallStatus != fmu_kCLEAN_NO_ACTION )
            return iCallStatus;
   }

   int iOppValence  = get_valence( pCOppNode );
   int iNextValence = get_valence( pCNextNode );
   int iPrevValence = get_valence( pCPrevNode );

   if( iOppValence == 5 )
   {
      boundary_neighbors( pCBoundary, CNodes, CEdges, CQuads );
      iCallStatus = bndy_flat_5( pCCenterNode, CNodes, CEdges, CQuads );
      if( iCallStatus != fmu_kCLEAN_NO_ACTION )
            return iCallStatus;
   }

   if( iNextValence >= 5 && iPrevValence == 4 )
   {
      iCallStatus = shape_544( pCQuad, pCCenterNode, pCNextNode, pCOppNode,
                                                                 pCPrevNode );
      if( iCallStatus != fmu_kCLEAN_NO_ACTION )
            return iCallStatus;
   }

   if( iPrevValence >= 5 && iNextValence == 4)
   {
      iCallStatus = shape_445( pCQuad, pCCenterNode, pCNextNode, pCOppNode,
                                                                 pCPrevNode );
      if( iCallStatus != fmu_kCLEAN_NO_ACTION )
            return iCallStatus;
   }

   //The opposite node is of valence 4.
   if( iOppValence == 4 )
   {
      iCallStatus = bndy_flat_4( pCBoundary );
      if( iCallStatus != fmu_kCLEAN_NO_ACTION )
            return iCallStatus;
   }

   return fmu_kCLEAN_NO_ACTION;
}


int QuadCleanTool::bndy_flat_4( nodClassDynArray *pCBoundary )
{
   
   int ii = 0;
   int iCallStatus = 0;
   int iStartIndex = 0, iDirection = 0, iScan = 0, nSurfaceEdges = 0, nSteps = 0;
   int iPrevScan = 0, iNextScan = 0;
   bool lWorking = true;
   nodClass *pCNode;
   nodClassDynArray CNodes;
   fmuEdgeDynArray CEdges;
   elmQClassDynArray CQuads;

   //This routine rotates edges, propagating them along the boundary
   //until it finds a 5 valent node. The first step is to find the
   //best direction for propagation. The search is interrupted if
   //there is a node with its only edges on the boundary (another
   //3 valent or a corner).
   if( pCBoundary->length() < 8 ) return fmu_kCLEAN_NO_ACTION;

   iStartIndex = pCBoundary->step(0);
   for( iDirection = -1; iDirection < 2; iDirection+=2 )
   {
      lWorking = true;
      for( ii = 1; ii < 4 && lWorking; ++ii )
      {
         iScan = 0;
         pCBoundary->step( iDirection );
         CNodes.clear();
         CEdges.clear();
         CQuads.clear();
         pCNode = pCBoundary->get();
         boundary_neighbors( pCBoundary, CNodes, CEdges, CQuads );
         nSurfaceEdges = CQuads.length() - 1;
         if( iDirection == 1 )
         {
            nSteps = CQuads.length() - 2;
            if( nSteps > 0 )
            {
               CNodes.step(2 * nSteps );
               CEdges.step( nSteps );
               CQuads.step( nSteps );
            }
         }
         if( nSurfaceEdges == 0 )
         {
            iScan = -1;
            lWorking = false;
         }
         if( nSurfaceEdges >= 2 )
         {
            fmuVector CVec23( CNodes.next(2), CNodes.next(3) );
            fmuVector CVec21( CNodes.next(2), CNodes.next() );
            if( CNormal.angle( CVec23, CVec21 ) > fmu_kSHAPE_ANGLE_TOL )
            {
               iScan = -1;
               lWorking = false;
               continue;
            }
            fmuVector CVec12( CNodes.next(), CNodes.next(2) );
            fmuVector CVec10( CNodes.next(), CNodes.get() );
            if( CNormal.angle( CVec12, CVec10 ) > fmu_kSHAPE_ANGLE_TOL )
            {
               iScan = -1;
               lWorking = false;
               continue;
            }
            fmuVector CVecC0( pCNode, CNodes.get() );
            fmuVector CVecC4( pCNode, CNodes.next(4) );
            if( CNormal.angle( CVecC0, CVecC4 ) > fmu_kSHAPE_ANGLE_TOL )
            {
               iScan = -1;
               lWorking = false;
               continue;
            }
            fmuVector CVec4C( CNodes.next(4), pCNode );
            fmuVector CVec43( CNodes.next(4), CNodes.next(3) );
            if( CNormal.angle( CVec4C, CVec43 ) > fmu_kSHAPE_ANGLE_TOL )
            {
               iScan = -1;
               lWorking = false;
               continue;
            }
            iScan = ii;
            lWorking = false;
         }
      }
      pCBoundary->set( iStartIndex );
      if( iDirection == -1 ) iPrevScan = iScan;
      else                   iNextScan = iScan;
   }

   if( iPrevScan <= 0 && iNextScan <= 0 ) return fmu_kCLEAN_NO_ACTION;

   if( iPrevScan > 0 && ( iNextScan <= 0 || iNextScan > iPrevScan ) )
   {
      pCNode = pCBoundary->get();
      for( ii = 0; ii < iPrevScan; ++ii )
      {
         CNodes.clear();
         CEdges.clear();
         CQuads.clear();
         pCBoundary->step(-1);
         pCNode = pCBoundary->get();
         boundary_neighbors( pCBoundary, CNodes, CEdges, CQuads );
         iCallStatus = rotate_edge( pCNode, CNodes, CEdges, CQuads, 0 );
         if( iCallStatus == fmu_kCLEAN_NO_MEMORY ) return fmu_kCLEAN_NO_MEMORY;
         if( iCallStatus == fmu_kCLEAN_ABORT ) return fmu_kCLEAN_ABORT;
      }
   }
   else
   {
      pCNode = pCBoundary->get();
      for( ii = 0; ii < iNextScan; ++ii )
      {
         CNodes.clear();
         CEdges.clear();
         CQuads.clear();
         pCBoundary->step(1);
         pCNode = pCBoundary->get();
         boundary_neighbors( pCBoundary, CNodes, CEdges, CQuads );
         nSteps = CQuads.length() - 2;
         if( nSteps > 0 )
         {
            CNodes.step(2 * nSteps );
            CEdges.step( nSteps );
            CQuads.step( nSteps );
         }
         iCallStatus = rotate_edge( pCNode, CNodes, CEdges, CQuads, 1 );
         if( iCallStatus == fmu_kCLEAN_NO_MEMORY ) return fmu_kCLEAN_NO_MEMORY;
         if( iCallStatus == fmu_kCLEAN_ABORT ) return fmu_kCLEAN_ABORT;
      }
   }
   return fmu_kCLEAN_MESH_CHANGED;
  
}


int QuadCleanTool::bndy_flat_5( nodClass *pCCenterNode,
                                nodClassDynArray &CNodes,
                                fmuEdgeDynArray &CEdges,
                                elmQClassDynArray &CQuads )
{
    if( CNodes.length() <= 0 )
        return fmu_kCLEAN_NO_ACTION;

    if( CEdges.length() <= 0 )
        return fmu_kCLEAN_NO_ACTION;

    if( CQuads.length() <= 0 )
        return fmu_kCLEAN_NO_ACTION;

    nodClass *pCNode0 = CNodes.get();
    nodClass *pCNode1 = CNodes.next();
    nodClass *pCNode2 = CNodes.next(2);
    fmuEdge *pCEdge0  = CEdges.get();
    elmQClass *pCQuad0 = CQuads.get();

    // Check to see if all entites exist - NM
    if( pCNode0 == NULL ||
        pCNode1 == NULL ||
        pCNode2 == NULL )
        return fmu_kCLEAN_NO_ACTION;

    // Check to see if all edges exist  - NM
    if( !pCEdge0 ) 
        return fmu_kCLEAN_NO_ACTION;

    // Check to see if all quads exist - NM
    if( !pCQuad0 ) 
        return fmu_kCLEAN_NO_ACTION;

    fmuEdge *pCEdge01 = pCNode0->shared_edge ( pCNode1 );
    if( !pCEdge01 ) 
        return fmu_kCLEAN_NO_ACTION;
    if( pCEdge01->frozen() ) 
        return fmu_kCLEAN_NO_ACTION;

    elmQClass *pCQuadA = pCEdge01->other_quad( pCQuad0 );
    if( !pCQuadA ) 
        return fmu_kCLEAN_NO_ACTION;

    nodClass *pCNodeA = pCQuadA->opposite_node ( pCNode1 );
    if( !pCNodeA ) 
        return fmu_kCLEAN_NO_ACTION;

    //Use a different pattern if node 0 doesn't need an edge into it.
    fmuVector CVec0( pCNode0, pCCenterNode );
    fmuVector CVecA( pCNode0, pCNodeA );
    double dAngle = CNormal.angle( CVecA, CVec0 );

    //If the angle on one side is within tolerance, test the other side
    //to see if that one is outside of tolerance
    if( dAngle < fmu_kSHAPE_ANGLE_TOL )
    {
        fmuEdge *pCEdge12 = pCNode1->shared_edge ( pCNode2 );
        if( !pCEdge12 ) 
            return fmu_kCLEAN_NO_ACTION;
        if( pCEdge12->frozen() )
            return fmu_kCLEAN_NO_ACTION;

        elmQClass *pCQuadC = pCEdge12->other_quad( pCQuad0 );
        if( !pCQuadC ) 
            return fmu_kCLEAN_NO_ACTION;

        nodClass *pCNodeE = pCQuadC->opposite_node ( pCNode1 );
        if( !pCNodeE ) 
            return fmu_kCLEAN_NO_ACTION;

        fmuVector CVec2( pCNode2, pCCenterNode );
        fmuVector CVecE( pCNode2, pCNodeE );
        dAngle = CNormal.angle( CVec2, CVecE );
        if( dAngle < fmu_kSHAPE_ANGLE_TOL )
            return fmu_kCLEAN_NO_ACTION;
    }

    nodClass *pCNodeB = pCQuadA->prev_node ( pCNode1 );
    fmuEdge *pCEdge1B = pCNode1->shared_edge ( pCNodeB );
    if( !pCEdge1B ) 
        return fmu_kCLEAN_NO_ACTION;
    if( pCEdge1B->frozen() ) 
        return fmu_kCLEAN_NO_ACTION;
    elmQClass *pCQuadB = pCEdge1B->other_quad( pCQuadA );
    if( !pCQuadB ) 
        return fmu_kCLEAN_NO_ACTION;
    nodClass *pCNodeC = pCQuadB->opposite_node ( pCNode1 );
    nodClass *pCNodeD = pCQuadB->prev_node ( pCNode1 );

    if( lDrawCleaning )
    {
        //draw_open();
        //pCQuad0->draw_me( dp_kCLR_BLUE );
        //pCQuadA->draw_me( dp_kCLR_GRBLUE );
        //pCQuadB->draw_me( dp_kCLR_LTBLUE );
        // Refresh debug graphics
        //EraseElements();
        pCQuad0->draw_me( dp_kCLR_YELLOW );
        pCQuadA->draw_me( dp_kCLR_YELLOW );
        pCQuadB->draw_me( dp_kCLR_YELLOW );
        //draw_close();
    }


    // needs to be this high because of the 'goto's below...
    nodClassDynArray CSmoothList;
    int fStatus = fmu_kCLEAN_NO_ACTION;

    // Test the element quality of the 5 would-be quads
    nodClass *pCNodeX=NULL, *pCNodeY=NULL;
    double dQ=0.0;
    double adNodeX[2]={0.0}, adNodeY[2]={0.0};
    // Get locations of 2 would-be nodes
    this->get_two_middle_node_location(pCNode1, pCNodeB, adNodeX, adNodeY);

    // Compute Jac ratio of the would-be elements
    dQ = test_quad_quality( pCCenterNode->node_x(), pCCenterNode->node_y(),
                            adNodeX[0],             adNodeX[1],
                            pCNode1->node_x(),      pCNode1->node_y(),
                            pCNode2->node_x(),      pCNode2->node_y());
    if (dQ < 0.1)
        goto Error;
    dQ = test_quad_quality( pCCenterNode->node_x(), pCCenterNode->node_y(), 
                            pCNode0->node_x(),      pCNode0->node_y(),
                            adNodeY[0],             adNodeY[1],
                            adNodeX[0],             adNodeX[1]);
    if (dQ < 0.1)
        goto Error;
    dQ = test_quad_quality( pCNode0->node_x(), pCNode0->node_y(),
                            pCNodeA->node_x(), pCNodeA->node_y(),
                            pCNodeB->node_x(), pCNodeB->node_y(), 
                            adNodeY[0],        adNodeY[1]);
    if (dQ < 0.1)
        goto Error;
    dQ = test_quad_quality( adNodeY[0],        adNodeY[1],
                            pCNodeB->node_x(), pCNodeB->node_y(),
                            pCNodeC->node_x(), pCNodeC->node_y(), 
                            pCNodeD->node_x(), pCNodeD->node_y());
    if (dQ < 0.1)
        goto Error;
    dQ = test_quad_quality( adNodeX[0],        adNodeX[1],
                            adNodeY[0],        adNodeY[1], 
                            pCNodeD->node_x(), pCNodeD->node_y(), 
                            pCNode1->node_x(), pCNode1->node_y());
    if (dQ < 0.1)
        goto Error;

    //Need two new nodes
    two_middle_nodes( pCNode1, pCNodeB, &pCNodeX, &pCNodeY );

    delete_this_quad( pCQuad0 );
    delete_this_quad( pCQuadA );
    delete_this_quad( pCQuadB );

    if( delete_this_edge( pCEdge01 ) == fmu_kCLEAN_ABORT )
        return fmu_kCLEAN_ABORT;
    if( delete_this_edge( pCEdge1B ) == fmu_kCLEAN_ABORT )
        return fmu_kCLEAN_ABORT;


    //Build 5 new quads
    elmQClass *pCQuadC, *pCQuadD;
    pCQuad0 = new elmQClass( pCCenterNode, pCNodeX, pCNode1, pCNode2,
        pCMaster );
    if( !pCQuad0 ) 
        return fmu_kCLEAN_NO_MEMORY;
    pCQuadA = new elmQClass( pCCenterNode, pCNode0, pCNodeY, pCNodeX,
        pCMaster );
    if( !pCQuadA ) 
        return fmu_kCLEAN_NO_MEMORY;
    pCQuadB = new elmQClass( pCNode0, pCNodeA, pCNodeB, pCNodeY, pCMaster );
    if( !pCQuadB ) 
        return fmu_kCLEAN_NO_MEMORY;
    pCQuadC = new elmQClass( pCNodeY, pCNodeB, pCNodeC, pCNodeD, pCMaster );
    if( !pCQuadC ) 
        return fmu_kCLEAN_NO_MEMORY;
    pCQuadD = new elmQClass( pCNodeX, pCNodeY, pCNodeD, pCNode1, pCMaster );
    if( !pCQuadD ) 
        return fmu_kCLEAN_NO_MEMORY;

    fStatus = fmu_kCLEAN_MESH_CHANGED;

    CSmoothList.append( pCNode0 );
    CSmoothList.append( pCNode1 );
    CSmoothList.append( pCNode2 );
    CSmoothList.append( pCNodeA );
    CSmoothList.append( pCNodeB );
    CSmoothList.append( pCNodeC );
    CSmoothList.append( pCNodeD );
    CSmoothList.append( pCCenterNode );
    smooth_mesh( &CSmoothList );

    if( lDrawCleaning )
    {
        //draw_mesh();
        //pCQuad0->draw_me( dp_kCLR_BLUE );
        //pCQuadA->draw_me( dp_kCLR_GRBLUE );
        //pCQuadB->draw_me( dp_kCLR_LTBLUE );
        //pCQuadC->draw_me( dp_kCLR_GOLD );
        //pCQuadD->draw_me( dp_kCLR_YELLOW );
        // Refresh debug graphics
        //EraseElements();
        pCQuad0->draw_me( dp_kCLR_BLUE );
        pCQuadA->draw_me( dp_kCLR_BLUE );
        pCQuadB->draw_me( dp_kCLR_BLUE );
        pCQuadC->draw_me( dp_kCLR_BLUE );
        pCQuadD->draw_me( dp_kCLR_BLUE );
        //draw_close();
    }

Error:
    reclean( pCNode0 );
    reclean( pCNode1 );
    reclean( pCNodeA );
    reclean( pCNodeX );
    reclean( pCNodeY );
    reclean( pCCenterNode );

    return fStatus;
}



int QuadCleanTool::bndy_valence_clean( nodClass *pCCenterNode,
                                           int iCenterValence,
                                           nodClassDynArray &CNodes,
                                           fmuEdgeDynArray &CEdges,
                                           elmQClassDynArray &CQuads )
{
   
   //Lists aren't rotated in this routine due to the boundary.
   //
   //There are patterns worth checking and others worth ignoring. The
   //ones to ignore may have tria elements that change the valence but
   //not the node count. Here are patterns worth checking:
   //
   //   4   3  2              3   2   1              6 5 4 3 2
   // 5  \    /   1               |                7  \  |  /  1
   //     \  /                    |                    \ | /
   // 6 ___\/___ 0           4 ___|___ 0           8 ___\|/___ 0
   //       :                     :                      :
   //       :                     :                      :
   //  5 valent center node,  4 valent center node  6 valent center node
   //  one virtual edge       one virtual edges     one virtual edge
   //  7 neighbor nodes       5 neighbor nodes      9 neighbor nodes
   //
   //It is also worth checking corresponding patterns that have no virtual
   //edge (boundary is concave) and two more neighbor nodes.
   //There is one exception to this. It is the first case.
   //
   //There are series of "condensed" cases which count up the number of
   //irregular nodes in addition to the center nodes. For a 6 valent
   //center node, one other irregular node is sufficient as this changes
   //a 6 valent to a 5 valent. For 5 valent center nodes, having only on
   //other irregular node still leaves 2 irregular nodes when done.

   int ii;
   int iCallStatus;
   int nIrregular = 0;
   //const int ARRAY_SIZE = 20;
   int node_count = CNodes.length();
   int aiValence[ARRAY_SIZE]={0}, aiHardpt[ARRAY_SIZE]={0};
   double dAngle;

   for( ii = 0; ii < node_count; ++ii )
   {
      if (CNodes[ii])
	  {
		aiValence[ii] = get_valence( CNodes[ii] );
		aiHardpt[ii] = CNodes[ii]->hardpt();
	  }
	  else
		aiValence[ii] = 0;
   }

   if( node_count == 9 &&
       iCenterValence == 5 &&
       aiValence[3] >= 4 &&
       aiValence[4] == 3 &&
       aiValence[5] == 3 &&
       aiValence[6] == 3 &&
       aiValence[7] >= 4 &&
       !aiHardpt[4] &&
       !aiHardpt[5] &&
       !aiHardpt[6] )
   {
      // case 5-0004+3334+ boundary
      CNodes.step(2);
      CEdges.step();
      CQuads.step();
      iCallStatus = bndy_5043334( pCCenterNode, CNodes, CEdges, CQuads );
      if( iCallStatus != fmu_kCLEAN_NO_ACTION ) return iCallStatus;
      //Operation not done, reset the lists
      CNodes.step(-2);
      CEdges.step(-1);
      CQuads.step(-1);
   }

   //Eliminate patterns that are of no interest and that would cause problems
   if( iCenterValence == 3 ) return fmu_kCLEAN_NO_ACTION;
   if( iCenterValence == 4 && node_count < 5 ) return fmu_kCLEAN_NO_ACTION;
   if( iCenterValence == 5 && node_count < 7 ) return fmu_kCLEAN_NO_ACTION;
   if( iCenterValence == 6 && node_count < 9 ) return fmu_kCLEAN_NO_ACTION;

   if( iCenterValence == 4 &&
       aiValence[0] == 4 &&
       aiValence[1] == 3 &&
       aiValence[2] == 4 &&
       aiValence[3] == 4 &&
       aiValence[4] == 4 &&
       !aiHardpt[1] &&
       !aiHardpt[2] &&
       !aiHardpt[3] )
   {
      // case 4-43444 boundary, check for edges to straighten out.
      iCallStatus = bndy_fix_skew( pCCenterNode, CNodes, CQuads, 0 );
      if( iCallStatus != fmu_kCLEAN_NO_ACTION ) return iCallStatus;
   }

   if( iCenterValence == 4 &&
       aiValence[0] == 4 &&
       aiValence[1] == 4 &&
       aiValence[2] == 4 &&
       aiValence[3] == 3 &&
       aiValence[4] == 4 &&
       !aiHardpt[1] &&
       !aiHardpt[2] &&
       !aiHardpt[3] )
   {
      // case 4-44434 boundary, check for edges to straighten out.
      iCallStatus = bndy_fix_skew( pCCenterNode, CNodes, CQuads, 1 );
      if( iCallStatus != fmu_kCLEAN_NO_ACTION ) return iCallStatus;
   }

   // Condensing 7 cases into some simple checks.
   // case 5-0350300 boundary, also 5-0350400, 5-0340300, 5-0450300,
   //                               5-0340400, 5-0450400, 5-0440300
   // All are variations on 5-04-4+04-00 where at least one is not 4
   if( iCenterValence == 5 &&
       aiValence[1] <= 4 &&
       aiValence[2] >= 4 &&
       aiValence[4] <= 4 )
   {
      if( aiValence[1] == 3 ) ++nIrregular;
      if( aiValence[2] == 5 ) ++nIrregular;
      if( aiValence[4] == 3 ) ++nIrregular;
      if( nIrregular > 1 )
      {
         iCallStatus = bndy_rotate_edge( pCCenterNode, CNodes, CEdges,
                                                               CQuads, 1 );
         if( iCallStatus != fmu_kCLEAN_NO_ACTION ) return iCallStatus;
      }
      nIrregular = 0;
   }

   // Condensing another 7 cases into some simple checks.
   // case 5-0030530 boundary, also 5-0040530, 5-0030430, 5-0030540,
   //                               5-0040430, 5-0040540, 5-0030440
   // All are variations on 5-004-04+4-0- where at least one is not 4
   if( iCenterValence == 5 &&
       aiValence[2] <= 4 &&
       aiValence[4] >= 4 &&
       aiValence[5] <= 4 )
   {
      if( aiValence[2] == 3 ) ++nIrregular;
      if( aiValence[4] == 5 ) ++nIrregular;
      if( aiValence[5] == 3 ) ++nIrregular;
      if( nIrregular > 1 )
      {
         CNodes.step(2);
         CEdges.step();
         CQuads.step();
         iCallStatus = bndy_rotate_edge( pCCenterNode, CNodes, CEdges,
                                                               CQuads, 2 );
         if( iCallStatus != fmu_kCLEAN_NO_ACTION ) return iCallStatus;
         //Operation not done, reset the lists
         CNodes.step(-2);
         CEdges.step(-1);
         CQuads.step(-1);
      }
      nIrregular = 0;
   }

   if( iCenterValence == 5 &&
       aiValence[2] <= 4 &&
       aiValence[4] >= 4 &&
       aiValence[5] == 3 &&
       aiValence[6] >= 5 )
   {
      // case 5-004-04+35 boundary
      CNodes.step(2);
      CEdges.step();
      CQuads.step();
      iCallStatus = bndy_rotate_edge( pCCenterNode, CNodes, CEdges,
                                                            CQuads, 2 );
      if( iCallStatus != fmu_kCLEAN_NO_ACTION ) return iCallStatus;
      //Operation not done, reset the lists
      CNodes.step(-2);
      CEdges.step(-1);
      CQuads.step(-1);
   }

   if( iCenterValence == 5 &&
       aiValence[0] >= 5 &&
       aiValence[1] == 3 &&
       aiValence[2] >= 4 &&
       aiValence[4] <= 4 )
   {
      // case 5-534+04-00 boundary
      iCallStatus = bndy_rotate_edge( pCCenterNode, CNodes, CEdges,
                                                            CQuads, 1 );
      if( iCallStatus != fmu_kCLEAN_NO_ACTION ) return iCallStatus;
   }

   if( iCenterValence == 5 &&
       aiValence[0] == 4 &&
       aiValence[1] == 3 &&
       aiValence[2] == 4 &&
       aiValence[3] == 4 &&
       aiValence[4] == 4 &&
       !aiHardpt[1] &&
       !aiHardpt[2] &&
       !aiHardpt[3] )
   {
      // case 5-43444 boundary
      iCallStatus = bndy_543( pCCenterNode, CNodes, CEdges, CQuads );
      if( iCallStatus != fmu_kCLEAN_NO_ACTION ) return iCallStatus;
   }

   if( iCenterValence == 5 &&
       aiValence[2] == 4 &&
       aiValence[3] == 4 &&
       aiValence[4] == 4 &&
       aiValence[5] == 3 &&
       aiValence[6] == 4 &&
       !aiHardpt[3] &&
       !aiHardpt[4] &&
       !aiHardpt[5] )
   {
      // case 5-0044434 boundary
      iCallStatus = bndy_543m( pCCenterNode, CNodes, CEdges, CQuads );
      if( iCallStatus != fmu_kCLEAN_NO_ACTION ) return iCallStatus;
   }

   if( iCenterValence == 5 &&
       aiValence[1] >= 4 &&
       aiValence[2] == 3 &&
       aiValence[3] == 3 &&
       aiValence[4] == 3 &&
       aiValence[5] >= 4 &&
       !aiHardpt[2] &&
       !aiHardpt[3] &&
       !aiHardpt[4] )
   {
      // case 5-04+3334+ boundary
      iCallStatus = bndy_5043334( pCCenterNode, CNodes, CEdges, CQuads );
      if( iCallStatus != fmu_kCLEAN_NO_ACTION ) return iCallStatus;
   }


   if( iCenterValence == 5 &&
       aiValence[0] <= 4 &&
       aiValence[2] >= 4 &&
       aiValence[3] <= 4 )
   {
      // case 5-4-04+4-00 boundary with angle check
      fmuVector CVecC0( pCCenterNode, CNodes.get() );
      fmuVector CVecC2( pCCenterNode, CNodes.next(2) );
      dAngle = CNormal.angle( CVecC0, CVecC2 );
      if( dAngle <= fmu_kBNDY_SKEW_TOL )
      {
         iCallStatus = bndy_rotate_edge( pCCenterNode, CNodes, CEdges,
                                                               CQuads, 0 );
         if( iCallStatus != fmu_kCLEAN_NO_ACTION ) return iCallStatus;
      }
   }

   if( iCenterValence == 5 &&
       aiValence[3] <= 4 &&
       aiValence[4] >= 4 &&
       aiValence[6] <= 4 )
   {
      // case 5-0004-4+4- boundary with angle check
      fmuVector CVecC4( pCCenterNode, CNodes.next(4) );
      fmuVector CVecC6( pCCenterNode, CNodes.next(6) );
      dAngle = CNormal.angle( CVecC4, CVecC6 );
      if( dAngle <= fmu_kBNDY_SKEW_TOL )
      {
         CNodes.step(2);
         CEdges.step();
         CQuads.step();
         iCallStatus = bndy_rotate_edge( pCCenterNode, CNodes, CEdges,
                                                               CQuads, 1 );
         if( iCallStatus != fmu_kCLEAN_NO_ACTION ) return iCallStatus;
         CNodes.step(-2);
         CEdges.step(-1);
         CQuads.step(-1);
      }
   }

   //These cases move a 5-3 transition away from the boundary, but don't
   //appear to improve the mesh. Commented out for now
/**********
   if( iCenterValence == 5 &&
       aiValence[1] == 4 &&
       aiValence[2] == 4 &&
       aiValence[3] == 4 &&
       aiValence[4] == 3 &&
       aiValence[5] == 4 &&
       !aiHardpt[1] &&
       !aiHardpt[2] &&
       !aiHardpt[3] &&
       !aiHardpt[4] &&
       !aiHardpt[5] )
   {
      // case 5-0444340 boundary
      CNodes.step(2);
      CEdges.step();
      CQuads.step();
      iCallStatus = bndy_5043( pCCenterNode, CNodes, CEdges, CQuads, 0 );
      if( iCallStatus != fmu_kCLEAN_NO_ACTION ) return iCallStatus;
      //Operation not done, reset the lists
      CNodes.step(-2);
      CEdges.step(-1);
      CQuads.step(-1);
   }

   if( iCenterValence == 5 &&
       aiValence[1] == 4 &&
       aiValence[2] == 3 &&
       aiValence[3] == 4 &&
       aiValence[4] == 4 &&
       aiValence[5] == 4 &&
       !aiHardpt[1] &&
       !aiHardpt[2] &&
       !aiHardpt[3] &&
       !aiHardpt[4] &&
       !aiHardpt[5] )
   {
      // case 5-0434440 boundary
      iCallStatus = bndy_5043( pCCenterNode, CNodes, CEdges, CQuads, 1 );
      if( iCallStatus != fmu_kCLEAN_NO_ACTION ) return iCallStatus;
   }
**********/

   if( iCenterValence == 6 &&
       aiValence[2] == 3 &&
       aiValence[6] == 3 )
   {
      // case 6-003000300 boundary
      CNodes.step(2);
      CEdges.step();
      CQuads.step();
      iCallStatus = across_open_face( pCCenterNode, CNodes, CEdges, CQuads );
      if( iCallStatus != fmu_kCLEAN_NO_ACTION ) return iCallStatus;
      //Operation not done, reset the lists
      CNodes.step(-2);
      CEdges.step(-1);
      CQuads.step(-1);
   }

   // There are several blocks of condensed checks. They are in order such
   // that central nodes are checked before nodes along the boundary.
   // Condensing 7 cases into some simple checks.
   // case 6-003053000 boundary, also 6-003054000, 6-003043000, 6-004053000,
   //                                 6-003044000, 6-004054000, 6-004043000,
   // All are variations on 6-004-04+4-00 where at least one is not 4
   // and avoiding case 6-000303000
   if( iCenterValence == 6 &&
       aiValence[2] <= 4 &&
       aiValence[4] >= 4 &&
       aiValence[5] <= 4 &&
       aiValence[3] >= 4 )
   {
      if( aiValence[2] == 3 ) ++nIrregular;
      if( aiValence[4] == 5 ) ++nIrregular;
      if( aiValence[5] == 3 ) ++nIrregular;
      if( nIrregular >= 1 )
      {
         CNodes.step(2);
         CEdges.step();
         CQuads.step();
         iCallStatus = bndy_rotate_edge( pCCenterNode, CNodes, CEdges,
                                                               CQuads, 2 );
         if( iCallStatus != fmu_kCLEAN_NO_ACTION ) return iCallStatus;
         CNodes.step(-2);
         CEdges.step(-1);
         CQuads.step(-1);
      }
      nIrregular = 0;
   }

   // Condensing 7 cases into some simple checks.
   // case 6-000350300 boundary, also 6-000350400, 6-000340300, 6-000450300,
   //                                 6-000340400, 6-000450400, 6-000440300
   // All are variations on 6-0004-4+04-00 where at least one is not 4
   // and avoiding case 6-000303000
   if( iCenterValence == 6 &&
       aiValence[3] <= 4 &&
       aiValence[4] >= 4 &&
       aiValence[6] <= 4 &&
       aiValence[5] >= 4 )
   {
      if( aiValence[3] == 3 ) ++nIrregular;
      if( aiValence[4] == 5 ) ++nIrregular;
      if( aiValence[6] == 3 ) ++nIrregular;
      if( nIrregular >= 1 )
      {
         CNodes.step(2);
         CEdges.step();
         CQuads.step();
         iCallStatus = bndy_rotate_edge( pCCenterNode, CNodes, CEdges,
                                                               CQuads, 1 );
         if( iCallStatus != fmu_kCLEAN_NO_ACTION ) return iCallStatus;
         CNodes.step(-2);
         CEdges.step(-1);
         CQuads.step(-1);
      }
      nIrregular = 0;
   }

   // Condensing 7 cases into some simple checks.
   // case 6-035030000 boundary, also 6-035040000, 6-034030000, 6-045030000,
   //                                 6-034040000, 6-045040000, 6-044030000
   // All are variations on 6-04-4+04-0000 where at least one is not 4
   if( iCenterValence == 6 &&
       aiValence[1] <= 4 &&
       aiValence[2] >= 4 &&
       aiValence[4] <= 4 )
   {
      if( aiValence[1] == 3 ) ++nIrregular;
      if( aiValence[2] == 5 ) ++nIrregular;
      if( aiValence[4] == 3 ) ++nIrregular;
      if( nIrregular >= 1 )
      {
         iCallStatus = bndy_rotate_edge( pCCenterNode, CNodes, CEdges,
                                                               CQuads, 1 );
         if( iCallStatus != fmu_kCLEAN_NO_ACTION ) return iCallStatus;
      }
      nIrregular = 0;
   }

   // Condensing 7 cases into some simple checks.
   // case 6-000030530 boundary, also 6-000030540, 6-000030400, 6-000040530,
   //                                 6-000030440, 6-000040540, 6-000040430,
   // All are variations on 6-00004-04+4-0 where at least one is not 4
   if( iCenterValence == 6 &&
       aiValence[4] <= 4 &&
       aiValence[6] >= 4 &&
       aiValence[7] <= 4 )
   {
      if( aiValence[4] == 3 ) ++nIrregular;
      if( aiValence[6] == 5 ) ++nIrregular;
      if( aiValence[7] == 3 ) ++nIrregular;
      if( nIrregular >= 1 )
      {
         CNodes.step(4);
         CEdges.step(2);
         CQuads.step(2);
         iCallStatus = bndy_rotate_edge( pCCenterNode, CNodes, CEdges,
                                                               CQuads, 2 );
         if( iCallStatus != fmu_kCLEAN_NO_ACTION ) return iCallStatus;
         CNodes.step(-4);
         CEdges.step(-2);
         CQuads.step(-2);
      }
      nIrregular = 0;
   }

   // Condensing 7 cases into some simple checks.
   // case 6-305300000 boundary, also 6-305400000, 6-304300000, 6-405300000,
   //                                 6-304400000, 6-405400000, 6-404300000,
   // All are variations on 6-4-04+4-0000 where at least one is not 4
   // and avoiding case 6-000303000
   if( iCenterValence == 6 &&
       aiValence[0] <= 4 &&
       aiValence[2] >= 4 &&
       aiValence[3] <= 4 &&
       aiValence[5] >= 4 )
   {
      if( aiValence[0] == 3 ) ++nIrregular;
      if( aiValence[2] == 5 ) ++nIrregular;
      if( aiValence[3] == 3 ) ++nIrregular;
      if( nIrregular >= 1 )
      {
         iCallStatus = bndy_rotate_edge( pCCenterNode, CNodes, CEdges,
                                                               CQuads, 2 );
         if( iCallStatus != fmu_kCLEAN_NO_ACTION ) return iCallStatus;
      }
      nIrregular = 0;
   }

   // Condensing 7 cases into some simple checks.
   // case 6-000003503 boundary, also 6-000003504, 6-000003403, 6-000004503,
   //                                 6-000003404, 6-000004504, 6-000004403
   // All are variations on 6-000004-4+04- where at least one is not 4
   // and avoiding case 6-000303000
   if( iCenterValence == 6 &&
       aiValence[5] <= 4 &&
       aiValence[6] >= 4 &&
       aiValence[8] <= 4 &&
       aiValence[3] >= 4 )
   {
      if( aiValence[5] == 3 ) ++nIrregular;
      if( aiValence[6] == 5 ) ++nIrregular;
      if( aiValence[8] == 3 ) ++nIrregular;
      if( nIrregular >= 1 )
      {
         CNodes.step(4);
         CEdges.step(2);
         CQuads.step(2);
         iCallStatus = bndy_rotate_edge( pCCenterNode, CNodes, CEdges,
                                                               CQuads, 1 );
         if( iCallStatus != fmu_kCLEAN_NO_ACTION ) return iCallStatus;
         CNodes.step(-4);
         CEdges.step(-2);
         CQuads.step(-2);
      }
      nIrregular = 0;
   }

   if( iCenterValence == 6 &&
       aiValence[1] == 4 &&
       aiValence[2] == 4 &&
       aiValence[3] == 4 &&
       aiValence[4] == 4 &&
       aiValence[5] == 4 &&
       aiValence[6] == 4 &&
       aiValence[7] == 4 &&
       !aiHardpt[2] &&
       !aiHardpt[3] &&
       !aiHardpt[4] &&
       !aiHardpt[5] )
   {
      // case 6-044444440 boundary
      CNodes.step(2);
      CEdges.step();
      CQuads.step();
      iCallStatus = rotate_edge( pCCenterNode, CNodes, CEdges, CQuads, 1 );
      if( iCallStatus != fmu_kCLEAN_NO_ACTION ) return iCallStatus;
      CNodes.step(-2);
      CEdges.step(-1);
      CQuads.step(-1);
   }

   if( iCenterValence == 7 &&
       aiValence[2] == 3 &&
       aiValence[6] == 3 )
   {
      // case 7-00300030000 boundary
      CNodes.step(2);
      CEdges.step();
      CQuads.step();
      iCallStatus = across_open_face( pCCenterNode, CNodes, CEdges, CQuads );
      if( iCallStatus != fmu_kCLEAN_NO_ACTION ) return iCallStatus;
      //Operation not done, reset the lists
      CNodes.step(-2);
      CEdges.step(-1);
      CQuads.step(-1);
   }

   if( iCenterValence == 7 &&
       aiValence[4] == 3 &&
       aiValence[8] == 3 )
   {
      // case 7-00003000300 boundary
      CNodes.step(4);
      CEdges.step(2);
      CQuads.step(2);
      iCallStatus = across_open_face( pCCenterNode, CNodes, CEdges, CQuads );
      if( iCallStatus != fmu_kCLEAN_NO_ACTION ) return iCallStatus;
      //Operation not done, reset the lists
      CNodes.step(-4);
      CEdges.step(-2);
      CQuads.step(-2);
   }

   // Condensing 7 cases into some simple checks.
   // case 7-00305300000 boundary, also 7-00305400000, 7-00304300000,
   //                                   7-00405300000, 7-00304400000,
   //                                   7-00405400000, 7-00404300000
   // All are variations on 7-004-04+4-0000 where at least one is not 4
   // and avoiding case 7-000303000
   if( iCenterValence == 7 &&
       aiValence[2] <= 4 &&
       aiValence[4] >= 4 &&
       aiValence[5] <= 4 &&
       aiValence[3] >= 4 )
   {
      if( aiValence[2] == 3 ) ++nIrregular;
      if( aiValence[4] == 5 ) ++nIrregular;
      if( aiValence[5] == 3 ) ++nIrregular;
      if( nIrregular >= 1 )
      {
         CNodes.step(2);
         CEdges.step();
         CQuads.step();
         iCallStatus = bndy_rotate_edge( pCCenterNode, CNodes, CEdges,
                                                               CQuads, 2 );
         if( iCallStatus != fmu_kCLEAN_NO_ACTION ) return iCallStatus;
         CNodes.step(-2);
         CEdges.step(-1);
         CQuads.step(-1);
      }
      nIrregular = 0;
   }

   // Condensing 7 cases into some simple checks.
   // case 7-00035030000 boundary, also 7-00035040000, 7-00034030000,
   //                                   7-00045030000, 7-00034040000,
   //                                   7-00045040000, 7-00044030000
   // All are variations on 7-0004-4+04-0000 where at least one is not 4
   // and avoiding case 7-00030300000
   if( iCenterValence == 7 &&
       aiValence[3] <= 4 &&
       aiValence[4] >= 4 &&
       aiValence[6] <= 4 &&
       aiValence[5] >= 4 )
   {
      if( aiValence[3] == 3 ) ++nIrregular;
      if( aiValence[4] == 5 ) ++nIrregular;
      if( aiValence[6] == 3 ) ++nIrregular;
      if( nIrregular >= 1 )
      {
         CNodes.step(2);
         CEdges.step();
         CQuads.step();
         iCallStatus = bndy_rotate_edge( pCCenterNode, CNodes, CEdges,
                                                               CQuads, 1 );
         if( iCallStatus != fmu_kCLEAN_NO_ACTION ) return iCallStatus;
         CNodes.step(-2);
         CEdges.step(-1);
         CQuads.step(-1);
      }
      nIrregular = 0;
   }

   // Condensing 7 cases into some simple checks.
   // case 7-00305300000 boundary, also 7-00003054000, 7-00003043000,
   //                                   7-00004053000, 7-00003044000,
   //                                   7-00004054000, 7-00004043000
   // All are variations on 7-00004-04+4-00 where at least one is not 4
   // and avoiding case 7-00000303000
   if( iCenterValence == 7 &&
       aiValence[4] <= 4 &&
       aiValence[6] >= 4 &&
       aiValence[7] <= 4 &&
       aiValence[5] >= 4 )
   {
      if( aiValence[2] == 3 ) ++nIrregular;
      if( aiValence[4] == 5 ) ++nIrregular;
      if( aiValence[5] == 3 ) ++nIrregular;
      if( nIrregular >= 1 )
      {
         CNodes.step(4);
         CEdges.step(2);
         CQuads.step(2);
         iCallStatus = bndy_rotate_edge( pCCenterNode, CNodes, CEdges,
                                                               CQuads, 2 );
         if( iCallStatus != fmu_kCLEAN_NO_ACTION ) return iCallStatus;
         CNodes.step(-4);
         CEdges.step(-2);
         CQuads.step(-2);
      }
      nIrregular = 0;
   }

   // Condensing 7 cases into some simple checks.
   // case 7-00035030000 boundary, also 7-00000350400, 7-00000340300,
   //                                   7-00000450300, 7-00000340400,
   //                                   7-00000450400, 7-00000440300
   // All are variations on 7-000004-4+04-00 where at least one is not 4
   // and avoiding case 7-00000303000
   if( iCenterValence == 7 &&
       aiValence[5] <= 4 &&
       aiValence[6] >= 4 &&
       aiValence[8] <= 4 &&
       aiValence[7] >= 4 )
   {
      if( aiValence[3] == 3 ) ++nIrregular;
      if( aiValence[4] == 5 ) ++nIrregular;
      if( aiValence[6] == 3 ) ++nIrregular;
      if( nIrregular >= 1 )
      {
         CNodes.step(4);
         CEdges.step(2);
         CQuads.step(2);
         iCallStatus = bndy_rotate_edge( pCCenterNode, CNodes, CEdges,
                                                               CQuads, 1 );
         if( iCallStatus != fmu_kCLEAN_NO_ACTION ) return iCallStatus;
         CNodes.step(-4);
         CEdges.step(-2);
         CQuads.step(-2);
      }
      nIrregular = 0;
   }

   if( iCenterValence > 6 &&
       aiValence[2] <= 4 &&
       aiValence[3] <= 4 &&
       aiValence[4] >= 4 &&
       aiValence[6] <= 4 &&
       aiValence[7] <= 4 &&
       !aiHardpt[2] &&
       !aiHardpt[3] &&
       !aiHardpt[4] &&
       !aiHardpt[5] &&
       !aiHardpt[6] &&
       !aiHardpt[7] )
   {
      // case 7-004-4-4+04-4-000 boundary
      //The mirror will be written if encountered
      iCallStatus = bndy_7( pCCenterNode, CNodes, CEdges, CQuads );
      if( iCallStatus != fmu_kCLEAN_NO_ACTION ) return iCallStatus;
   }

   return fmu_kCLEAN_NO_ACTION;
   
}


int QuadCleanTool::bndy_fix_skew( nodClass *pCCenterNode,
                                      nodClassDynArray &CNodes,
                                      elmQClassDynArray &CQuads,
                                      int iSide )
{
   
   //In this case, the change is not made at the boundary (which is
   //all 4 valent, but one row of elements in from the boundary (if any).
   int iCallStatus;
   double dAngle;
   nodClassDynArray CDeleteNodes(10);
   nodClassDynArray CRingNodes(10);
   fmuEdgeDynArray CDeleteEdges(10);
   elmQClassDynArray CDeleteQuads(10);
   nodClass *pCNode0 = CNodes.get();
   nodClass *pCNode1 = CNodes.next();
   nodClass *pCNode2 = CNodes.next(2);
   nodClass *pCNode3 = CNodes.next(3);
   nodClass *pCNode4 = CNodes.next(4);

   if( iSide == 1 )
   {
      fmuVector CVecC2( pCCenterNode, pCNode2 );
      fmuVector CVecC4( pCCenterNode, pCNode4 );
      dAngle = CNormal.angle( CVecC2, CVecC4 );
      if( dAngle > fmu_kBNDY_SKEW_TOL ) return fmu_kCLEAN_NO_ACTION;
   }
   else
   {
      fmuVector CVecC0( pCCenterNode, pCNode0 );
      fmuVector CVecC2( pCCenterNode, pCNode2 );
      dAngle = CNormal.angle( CVecC0, CVecC2 );
      if( dAngle > fmu_kBNDY_SKEW_TOL ) return fmu_kCLEAN_NO_ACTION;
   }

   fmuEdge *pCEdge12 = pCNode1->shared_edge( pCNode2 );
   fmuEdge *pCEdge23 = pCNode2->shared_edge( pCNode3 );
   if( pCEdge12->frozen() ) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge23->frozen() ) return fmu_kCLEAN_NO_ACTION;

   elmQClass *pCQuadA = pCEdge12->other_quad( CQuads.get() );
   nodClass *pCNodeA = pCQuadA->next_node( pCNode1 );
   nodClass *pCNodeB = pCQuadA->opposite_node( pCNode1 );
   if( pCNodeA->hardpt() ) return fmu_kCLEAN_NO_ACTION;
   if( pCNodeB->hardpt() ) return fmu_kCLEAN_NO_ACTION;
   elmQClass *pCQuadB = pCEdge23->other_quad( CQuads.next() );
   nodClass *pCNodeC = pCQuadB->opposite_node( pCNode2 );
   if( pCNodeC->hardpt() ) return fmu_kCLEAN_NO_ACTION;
   fmuEdge *pCEdge2B = pCNode2->shared_edge( pCNodeB );
   if( pCEdge2B->frozen() ) return fmu_kCLEAN_NO_ACTION;
   //Two quads and one edge are to be deleted. No nodes deleted
   CDeleteQuads[0] = pCQuadA;
   CDeleteQuads[1] = pCQuadB;
   CDeleteEdges[0] = pCEdge2B;

   CRingNodes[0] = pCNode1;
   CRingNodes[1] = pCNodeA;
   CRingNodes[2] = pCNodeB;
   CRingNodes[3] = pCNodeC;
   CRingNodes[4] = pCNode3;
   CRingNodes[5] = pCNode2;
   if( iSide == 1 ) CRingNodes.step();

   iCallStatus =  fill_two( CDeleteNodes, CDeleteEdges, CDeleteQuads,
                                                          CRingNodes );
   if( iCallStatus == fmu_kCLEAN_MESH_CHANGED ) return fmu_kCLEAN_BACKSTEP;
   return iCallStatus;
  
}


int QuadCleanTool::bndy_rotate_edge( nodClass *pCCenterNode,
                                         nodClassDynArray &CNodes,
                                         fmuEdgeDynArray &CEdges,
                                         elmQClassDynArray &CQuads,
                                         int iSide )
{
   
   //Check whether permanent nodes will make distorted elements if
   //the change is made
   double dAngle;
   nodClass *pCNode0, *pCNode1, *pCNode2, *pCNode3;
   if( iSide == 1 )
   {
      //Take nodes one step farther around the loop
      pCNode0 = CNodes.next();
      pCNode1 = CNodes.next(2);
      pCNode2 = CNodes.next(3);
      pCNode3 = CNodes.next(4);
   }
   else
   {
      pCNode0 = CNodes.get();
      pCNode1 = CNodes.next(1);
      pCNode2 = CNodes.next(2);
      pCNode3 = CNodes.next(3);
   }

   if( pCNode1->hardpt() )
   {
      if( pCNode0->hardpt() ) return fmu_kCLEAN_NO_ACTION;
      fmuVector CVecC0( pCCenterNode, pCNode0 );
      fmuVector CVecC1( pCCenterNode, pCNode1 );
      if( CVecC1.length() < CVecC0.length() ) return fmu_kCLEAN_NO_ACTION;
      fmuVector CVec01( pCNode0, pCNode1 );
      fmuVector CVec03( pCNode0, pCNode3 );
      dAngle = CNormal.angle( CVec01, CVec03 );
      if( dAngle < fmu_kSMALL_SHAPE_TOL || dAngle > fmu_kSHAPE_ANGLE_TOL )
         return fmu_kCLEAN_NO_ACTION;
   }
   if( pCNode2->hardpt() )
   {
      if( pCNode3->hardpt() ) return fmu_kCLEAN_NO_ACTION;
      fmuVector CVecC3( pCCenterNode, pCNode3 );
      fmuVector CVecC2( pCCenterNode, pCNode2 );
      if( CVecC2.length() < CVecC3.length() ) return fmu_kCLEAN_NO_ACTION;
      fmuVector CVec30( pCNode3, pCNode0 );
      fmuVector CVec32( pCNode3, pCNode2 );
      dAngle = CNormal.angle( CVec30, CVec32 );
      if( dAngle < fmu_kSMALL_SHAPE_TOL || dAngle > fmu_kSHAPE_ANGLE_TOL )
         return fmu_kCLEAN_NO_ACTION;
   }

   //Check if the transition is too great (new elements too skewed) to allow
   //the change to be made
   fmuVector CLeft( pCCenterNode, CNodes.get() );
   fmuVector CCenter( pCCenterNode, CNodes.next(2) );
   fmuVector CRight( pCCenterNode, CNodes.next(4) );
   double dSkewLen = CCenter.length() * 2.0;
   double dLeftLen = CLeft.length();
   double dRightLen = CRight.length();
   if( dSkewLen < dLeftLen ) return fmu_kCLEAN_NO_ACTION;
   if( dSkewLen < dRightLen ) return fmu_kCLEAN_NO_ACTION;

   fmuVector CNew( pCNode0, pCNode3 );
   double dNewLen = CNew.length();
   if( CCenter.length() * 1.5 < dNewLen ) return fmu_kCLEAN_NO_ACTION;

   return rotate_edge( pCCenterNode, CNodes, CEdges, CQuads, iSide );
  
}


int QuadCleanTool::bndy_543( nodClass *pCCenterNode,
                                 nodClassDynArray &CNodes,
                                 fmuEdgeDynArray &CEdges,
                                 elmQClassDynArray &CQuads )
{
   
   double dAngle;
   fmuVector CVecC0( pCCenterNode, CNodes.get() );
   fmuVector CVecC4( pCCenterNode, CNodes.next(4) );
   dAngle = CNormal.angle( CVecC0, CVecC4 );
   if( dAngle < fmu_kRIGHT_LOW_TOL || dAngle > fmu_kRIGHT_HIGH_TOL )
      return fmu_kCLEAN_NO_ACTION;

   return rotate_edge( pCCenterNode, CNodes, CEdges, CQuads, 1 );
   
}


int QuadCleanTool::bndy_543m( nodClass *pCCenterNode,
                                 nodClassDynArray &CNodes,
                                 fmuEdgeDynArray &CEdges,
                                 elmQClassDynArray &CQuads )
{
   
   int iCallStatus;
   double dAngle;
   fmuVector CVecC2( pCCenterNode, CNodes.next(2) );
   fmuVector CVecC6( pCCenterNode, CNodes.next(6) );
   dAngle = CNormal.angle( CVecC2, CVecC6 );
   if( dAngle < fmu_kRIGHT_LOW_TOL || dAngle > fmu_kRIGHT_HIGH_TOL )
      return fmu_kCLEAN_NO_ACTION;

   CNodes.step(2);
   CEdges.step();
   CQuads.step();
   iCallStatus =  rotate_edge( pCCenterNode, CNodes, CEdges, CQuads, 2 );
   if( iCallStatus != fmu_kCLEAN_NO_ACTION ) return iCallStatus;

   //Action not taken, need to reset the list so it is correct for
   //another action.
   CNodes.step(-2);
   CEdges.step(-1);
   CQuads.step(-1);

   return fmu_kCLEAN_NO_ACTION;
  
}


int QuadCleanTool::bndy_5043334( nodClass *pCCenterNode,
                                     nodClassDynArray &CNodes,
                                     fmuEdgeDynArray &CEdges,
                                     elmQClassDynArray &CQuads )
{
   
   int iCallStatus;
   nodClassDynArray CDeleteNodes(10);
   nodClassDynArray CRingNodes(10);
   fmuEdgeDynArray CDeleteEdges(10);
   elmQClassDynArray CDeleteQuads(10);
   nodClass *pCNode0 = CNodes.get();
   nodClass *pCNode1 = CNodes.next();
   nodClass *pCNode2 = CNodes.next(2);
   nodClass *pCNode3 = CNodes.next(3);
   nodClass *pCNode4 = CNodes.next(4);
   nodClass *pCNode5 = CNodes.next(5);
   nodClass *pCNode6 = CNodes.next(6);
   fmuEdge  *pCEdge1 = CEdges.next();
   fmuEdge  *pCEdge2 = CEdges.next(2);
   elmQClass *pCQuad0 = CQuads.get();
   elmQClass *pCQuad1 = CQuads.next();
   elmQClass *pCQuad2 = CQuads.next(2);

   if( pCEdge1->frozen() ) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge2->frozen() ) return fmu_kCLEAN_NO_ACTION;

   fmuEdge *pCEdge12 = pCNode1->shared_edge( pCNode2 );
   if( pCEdge12->frozen() ) return fmu_kCLEAN_NO_ACTION;
   fmuEdge *pCEdge23 = pCNode2->shared_edge( pCNode3 );
   if( pCEdge23->frozen() ) return fmu_kCLEAN_NO_ACTION;
   fmuEdge *pCEdge34 = pCNode3->shared_edge( pCNode4 );
   if( pCEdge34->frozen() ) return fmu_kCLEAN_NO_ACTION;
   fmuEdge *pCEdge45 = pCNode4->shared_edge( pCNode5 );
   if( pCEdge45->frozen() ) return fmu_kCLEAN_NO_ACTION;


   elmQClass *pCQuadA = pCEdge12->other_quad( pCQuad0 );
   nodClass *pCNodeA = pCQuadA->opposite_node( pCNode2 );
   fmuEdge *pCEdge3A = pCNode3->shared_edge( pCNodeA );
   if( pCEdge3A->frozen() ) return fmu_kCLEAN_NO_ACTION;

   elmQClass *pCQuadB = pCEdge34->other_quad( pCQuad1 );

   //Quads to delete
   CDeleteQuads[0] = pCQuad0;
   CDeleteQuads[1] = pCQuad1;
   CDeleteQuads[2] = pCQuad2;
   CDeleteQuads[3] = pCQuadA;
   CDeleteQuads[4] = pCQuadB;
   //Edges to delete
   CDeleteEdges[0] = pCEdge1;
   CDeleteEdges[1] = pCEdge2;
   CDeleteEdges[2] = pCEdge12;
   CDeleteEdges[3] = pCEdge23;
   CDeleteEdges[4] = pCEdge34;
   CDeleteEdges[5] = pCEdge45;
   CDeleteEdges[6] = pCEdge3A;
   //Nodes to delete
   CDeleteNodes[0] = pCNode2;
   CDeleteNodes[1] = pCNode3;
   CDeleteNodes[2] = pCNode4;
   //Ring nodes
   CRingNodes[0] = pCCenterNode;
   CRingNodes[1] = pCNode0;
   CRingNodes[2] = pCNode1;
   CRingNodes[3] = pCNodeA;
   CRingNodes[4] = pCNode5;
   CRingNodes[5] = pCNode6;

   iCallStatus = fill_two( CDeleteNodes, CDeleteEdges, CDeleteQuads,
                                                         CRingNodes );

   return iCallStatus;
   
}


int QuadCleanTool::bndy_5043( nodClass *pCCenterNode,
                                 nodClassDynArray &CNodes,
                                 fmuEdgeDynArray &CEdges,
                                 elmQClassDynArray &CQuads, int iSide )
{
   
   int iCallStatus;
   nodClassDynArray CDeleteNodes(10);
   nodClassDynArray CRingNodes(10);
   fmuEdgeDynArray CDeleteEdges(10);
   elmQClassDynArray CDeleteQuads(10);
   nodClass *pCNode0 = CNodes.get();
   nodClass *pCNode1 = CNodes.next();
   nodClass *pCNode2 = CNodes.next(2);
   nodClass *pCNode3 = CNodes.next(3);
   nodClass *pCNode4 = CNodes.next(4);
   fmuEdge  *pCEdge1 = CEdges.next();
   elmQClass *pCQuad0 = CQuads.get();
   elmQClass *pCQuad1 = CQuads.next();

   if( pCEdge1->frozen() ) return fmu_kCLEAN_NO_ACTION;

   fmuEdge *pCEdge12 = pCNode1->shared_edge( pCNode2 );
   if( pCEdge12->frozen() ) return fmu_kCLEAN_NO_ACTION;
   fmuEdge *pCEdge23 = pCNode2->shared_edge( pCNode3 );
   if( pCEdge23->frozen() ) return fmu_kCLEAN_NO_ACTION;

   elmQClass *pCQuadA = pCEdge12->other_quad( pCQuad0 );
   nodClass *pCNodeA = pCQuadA->opposite_node( pCNode2 );

   //Quads to delete
   CDeleteQuads[0] = pCQuad0;
   CDeleteQuads[1] = pCQuad1;
   CDeleteQuads[2] = pCQuadA;
   //Edges to delete
   CDeleteEdges[0] = pCEdge1;
   CDeleteEdges[1] = pCEdge12;
   CDeleteEdges[2] = pCEdge23;
   //Nodes to delete
   CDeleteNodes[0] = pCNode2;
   //Ring nodes
   CRingNodes[0] = pCNode0;
   CRingNodes[1] = pCNode1;
   CRingNodes[2] = pCNodeA;
   CRingNodes[3] = pCNode3;
   CRingNodes[4] = pCNode4;
   CRingNodes[5] = pCCenterNode;
   if( iSide == 1 ) CRingNodes.step();

   iCallStatus = fill_two( CDeleteNodes, CDeleteEdges, CDeleteQuads,
                                                         CRingNodes );

   return iCallStatus;
}


int QuadCleanTool::bndy_7( nodClass *pCCenterNode,
                               nodClassDynArray &CNodes,
                               fmuEdgeDynArray &CEdges,
                               elmQClassDynArray &CQuads )
{
   
   nodClass *pCNode2 = CNodes.next(2);
   nodClass *pCNode3 = CNodes.next(3);
   nodClass *pCNode4 = CNodes.next(4);
   nodClass *pCNode5 = CNodes.next(5);
   nodClass *pCNode6 = CNodes.next(6);
   nodClass *pCNode7 = CNodes.next(7);
   nodClass *pCNode8 = CNodes.next(8);
   fmuEdge  *pCEdge2 = CEdges.next(2);
   fmuEdge  *pCEdge3 = CEdges.next(3);
   elmQClass *pCQuad1 = CQuads.next();
   elmQClass *pCQuad2 = CQuads.next(2);
   elmQClass *pCQuad3 = CQuads.next(3);

   if( pCEdge2->frozen() ) return fmu_kCLEAN_NO_ACTION;
   if( pCEdge3->frozen() ) return fmu_kCLEAN_NO_ACTION;

   if( lDrawCleaning )
   {
   //draw_open();
   //pCQuad1->draw_me( dp_kCLR_RED );
   //pCQuad2->draw_me( dp_kCLR_GOLD );
   //pCQuad3->draw_me( dp_kCLR_DKGRN );
        // Refresh debug graphics
        //EraseElements();
        pCQuad1->draw_me( dp_kCLR_YELLOW );
        pCQuad2->draw_me( dp_kCLR_YELLOW );
        pCQuad3->draw_me( dp_kCLR_YELLOW );
   //draw_close();
   }

   delete_this_quad( pCQuad1 );
   delete_this_quad( pCQuad2 );
   delete_this_quad( pCQuad3 );

   if( delete_this_edge( pCEdge2 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdge3 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;

   //Three new quads
   pCQuad1 = new elmQClass( pCCenterNode, pCNode2, pCNode7, pCNode8,
                                                                  pCMaster );
   if( !pCQuad1 ) return fmu_kCLEAN_NO_MEMORY;
   pCQuad2 = new elmQClass( pCNode2, pCNode3, pCNode6, pCNode7, pCMaster );
   if( !pCQuad2 ) return fmu_kCLEAN_NO_MEMORY;
   pCQuad3 = new elmQClass( pCNode3, pCNode4, pCNode5, pCNode6, pCMaster );
   if( !pCQuad3 ) return fmu_kCLEAN_NO_MEMORY;

   nodClassDynArray CSmoothList;
   CSmoothList.append( pCNode2 );
   CSmoothList.append( pCNode3 );
   CSmoothList.append( pCNode4 );
   CSmoothList.append( pCNode5 );
   CSmoothList.append( pCNode6 );
   CSmoothList.append( pCNode7 );
   CSmoothList.append( pCNode8 );
   smooth_mesh( &CSmoothList );
   if( lDrawCleaning )
   {
   //draw_mesh();
   //pCQuad1->draw_me( dp_kCLR_RED );
   //pCQuad2->draw_me( dp_kCLR_GOLD );
   //pCQuad3->draw_me( dp_kCLR_DKGRN );
        // Refresh debug graphics
        //EraseElements();
        pCQuad1->draw_me( dp_kCLR_BLUE );
        pCQuad2->draw_me( dp_kCLR_BLUE );
        pCQuad3->draw_me( dp_kCLR_BLUE );
 
   //draw_close();
   }

   return fmu_kCLEAN_MESH_CHANGED;

   
}


int QuadCleanTool::shape_clean()
{
  
   int iDebugFlag = 0;
   int ii, jj;
   int iCallStatus;
   int iStatus = fmu_kCLEAN_NO_ACTION;
   int iInitialLength;
   elmQClass *pCQuad;
   nodClass *pCNode, *pCPrevNode, *pCNextNode, *apCNodes[4];
   double dAngle[4], dAngleSum;

   //Copy the list of quads so that new quads are not checked to avoid
   //infinite loops. When quads are deleted (in delete_this_quad), do not
   //go beyond the current end of the array as intermediate positions
   //will not have initial values.
   CShapeList = (*pCMaster->global_quads());
   iInitialLength = CShapeList.length();

   for( ii = 0; ii < iInitialLength; ++ii )
   {
      pCQuad = CShapeList[ii];
      if( !pCQuad ) continue;
      pCQuad->nodes( apCNodes );

      if( iDebugFlag )
      {
      //draw_mesh();
        // Refresh debug graphics
        //EraseElements();
        pCQuad->draw_me( dp_kCLR_YELLOW );
      //pCQuad->draw_me( dp_kCLR_CYAN );
      //draw_close();
      }

      iCallStatus = fmu_kCLEAN_NO_ACTION;
      //If an angle is bad, use that node for the basis of cleaning the
      //quad. If the angle is very bad, it is likely from a bowtie quad
      //and another angle would be better for the basis of cleaning. There
      //won't be another angle if the quad forms a tight arrowhead.
      dAngleSum = 0.0;
      for( jj = 0; jj < 4; ++jj )
      {
         pCNode = apCNodes[jj];
         pCPrevNode = pCQuad->prev_node( pCNode );
         pCNextNode = pCQuad->next_node( pCNode );
         fmuVector CPrevVec( pCNode, pCPrevNode );
         fmuVector CNextVec( pCNode, pCNextNode );
         dAngle[jj] = CNormal.angle( CNextVec, CPrevVec );
         dAngleSum += dAngle[jj];
      }
      //If the quad is inverted, its angles add up to 4PI, not 2PI
      if( dAngleSum > THREPI )
      {
         iCallStatus = fix_inverted_quad( pCQuad, apCNodes );
         if( iCallStatus == fmu_kCLEAN_NO_MEMORY )
                                              return fmu_kCLEAN_NO_MEMORY;
         if( iCallStatus == fmu_kCLEAN_ABORT )
                                              return fmu_kCLEAN_ABORT;
         if( iCallStatus == fmu_kCLEAN_MESH_CHANGED )
            iStatus = iCallStatus;
      }
      for( jj = 0; jj < 4 && iCallStatus == fmu_kCLEAN_NO_ACTION; ++jj )
      {
         if( dAngle[jj] > fmu_kSHAPE_ANGLE_TOL &&
             dAngle[jj] < fmu_kEXTREME_ANGLE_TOL )
         {
            iCallStatus = fix_bad_angle( pCQuad, apCNodes[jj] );
            if( iCallStatus == fmu_kCLEAN_NO_MEMORY )
                                              return fmu_kCLEAN_NO_MEMORY;
            if( iCallStatus == fmu_kCLEAN_ABORT )
                                              return fmu_kCLEAN_ABORT;
            if( iCallStatus == fmu_kCLEAN_MESH_CHANGED )
               iStatus = iCallStatus;
         }
      }

      if( lDrawCleaning && iCallStatus == fmu_kCLEAN_MESH_CHANGED )
      {
      //draw_mesh();
      //draw_close();
      }
   }
   CShapeList.clear();
   return iStatus;
}


int QuadCleanTool::fix_inverted_quad( elmQClass *pCQuad,
                                           nodClass *apCNodes[4] )
{
   
   int ii, iCallStatus = fmu_kCLEAN_NO_ACTION;
   nodClass *pCNode, *pCPrevNode, *pCNextNode, *pCOppNode;
   nodClass *pCLANode, *pCLBNode, *pCRANode, *pCRBNode;
   elmQClass *pCPrevQuad, *pCNextQuad;
   fmuEdge *pCPrevEdge, *pCNextEdge, *pCPOEdge, *pCNOEdge;
   fmuEdge *pCTLAEdge, *pCTRAEdge, *pCPLBEdge, *pCNRBEdge;
   nodClassDynArray CRingNodes(ARRAY_SIZE);
   nodClassDynArray CNodes(ARRAY_SIZE);
   fmuEdgeDynArray CEdges(ARRAY_SIZE);
   elmQClassDynArray CQuads(ARRAY_SIZE);

   //This routine does not handle bowtie elements. Return these to the
   //caller for a different process.
   pCPrevEdge = apCNodes[0]->shared_edge( apCNodes[1] );
   pCNextEdge = apCNodes[2]->shared_edge( apCNodes[3] );
   if( edges_intersect( pCPrevEdge, pCNextEdge ) ) return fmu_kCLEAN_NO_ACTION;
   pCPrevEdge = apCNodes[0]->shared_edge( apCNodes[3] );
   pCNextEdge = apCNodes[1]->shared_edge( apCNodes[2] );
   if( edges_intersect( pCPrevEdge, pCNextEdge ) ) return fmu_kCLEAN_NO_ACTION;

   //This quad is upside down compared to the surface normal. Squeeze it
   //out by finding a node attached to neighboring elements.
   for( ii = 0; ii < 4; ++ii )
   {
      pCNode = apCNodes[ii];
      pCPrevNode = pCQuad->prev_node( pCNode );
	  if (!pCPrevNode)  return fmu_kCLEAN_NO_ACTION;
      pCNextNode = pCQuad->next_node( pCNode );
	  if (!pCNextNode)  return fmu_kCLEAN_NO_ACTION;
      pCPrevEdge = pCNode->shared_edge( pCPrevNode );
	  if (!pCPrevEdge)  return fmu_kCLEAN_NO_ACTION;
      pCNextEdge = pCNode->shared_edge( pCNextNode );
	  if (!pCNextEdge)  return fmu_kCLEAN_NO_ACTION;
      if( !pCPrevEdge->frozen() &&
          !pCNextEdge->frozen() &&
          ( !pCNode->hardpt() || pCNode->smooth_hardpt() ) )
      {
         //pCLANode may or may not be the same as pCRANode, pCTRAEdge
         //may or may not be the same as pCTLAEdge. This process does
         //not care.
         pCOppNode = pCQuad->opposite_node( pCNode );
		 if (!pCOppNode)  return fmu_kCLEAN_NO_ACTION;
         pCNOEdge = pCNextNode->shared_edge( pCOppNode );
		 if (!pCNOEdge)  return fmu_kCLEAN_NO_ACTION;
         pCPOEdge = pCPrevNode->shared_edge( pCOppNode );
		 if (!pCPOEdge)  return fmu_kCLEAN_NO_ACTION;
         pCPrevQuad = pCPrevEdge->other_quad( pCQuad );
		 if (!pCPrevQuad)  return fmu_kCLEAN_NO_ACTION;
         pCLANode = pCPrevQuad->opposite_node( pCPrevNode );
		 if (!pCLANode)  return fmu_kCLEAN_NO_ACTION;
         pCTLAEdge = pCNode->shared_edge( pCLANode );
		 if (!pCTLAEdge)  return fmu_kCLEAN_NO_ACTION;
         pCLBNode = pCPrevQuad->next_node( pCPrevNode );
		 if (!pCLBNode)  return fmu_kCLEAN_NO_ACTION;
         pCPLBEdge = pCPrevNode->shared_edge( pCLBNode );
		 if (!pCPLBEdge)  return fmu_kCLEAN_NO_ACTION;
         pCNextQuad = pCNextEdge->other_quad( pCQuad );
		 if (!pCNextQuad)  return fmu_kCLEAN_NO_ACTION;
         pCRANode = pCNextQuad->opposite_node( pCNextNode );
		 if (!pCRANode)  return fmu_kCLEAN_NO_ACTION;
         pCTRAEdge = pCNode->shared_edge( pCRANode );
		 if (!pCTRAEdge)  return fmu_kCLEAN_NO_ACTION;
         pCRBNode = pCNextQuad->prev_node( pCNextNode );
		 if (!pCRBNode)  return fmu_kCLEAN_NO_ACTION;
         pCNRBEdge = pCNextNode->shared_edge( pCRBNode );
		 if (!pCNRBEdge)  return fmu_kCLEAN_NO_ACTION;

         if( !pCTLAEdge->frozen() && !pCTRAEdge->frozen() )
         {
            if( edges_intersect( pCTLAEdge, pCNOEdge ) )
            {
               if(edges_intersect( pCPLBEdge, pCNOEdge ) )
               {
                  iCallStatus = split_edge( pCNOEdge );
               }
               else
               {
                  iCallStatus = collapse_inverted( pCQuad, pCNode );
               }
            }
            else if( edges_intersect( pCTRAEdge, pCPOEdge ) )
            {
               if(edges_intersect( pCNRBEdge, pCPOEdge ) )
               {
                  iCallStatus = split_edge( pCPOEdge );
               }
               else
               {
                  iCallStatus = collapse_inverted( pCQuad, pCNode );
               }
            }
            else if( edges_intersect( pCTLAEdge, pCPOEdge ) )
            {
               iCallStatus = collapse_inverted( pCQuad, pCNode );
            }
            else if( edges_intersect( pCTRAEdge, pCNOEdge ) )
            {
               iCallStatus = collapse_inverted( pCQuad, pCNode );
            }
            if( iCallStatus != fmu_kCLEAN_NO_ACTION ) return iCallStatus;
         }
      }
   }
   return fmu_kCLEAN_NO_ACTION;
   
}


int QuadCleanTool::collapse_inverted( elmQClass *pCQuad, nodClass *pCNode )
{
   
   int iDebugFlag = 0;
   int iCallStatus;
   nodClass *pCNextNode;
   fmuEdge *pCNextEdge;
   nodClassDynArray CNodes(ARRAY_SIZE);
   fmuEdgeDynArray CEdges(ARRAY_SIZE);
   elmQClassDynArray CQuads(ARRAY_SIZE);

   if( iDebugFlag )
   {
   //draw_mesh();
        // Refresh debug graphics
        //EraseElements();
        pCQuad->draw_me( dp_kCLR_YELLOW );
      
   //pCQuad->draw_me( dp_kCLR_GRBLUE );
   //draw_close();
   //pCNode->draw_me( dp_kCLR_GOLD );
   }

   pCNextNode = pCQuad->next_node( pCNode );
   pCNextEdge = pCNode->shared_edge( pCNextNode );
   //If smoothing changed the node to hardpoint, then this routine
   //is allowed to change it back.
   pCNode->hardpt( false );
   pCNode->smooth_hardpt( false );
   neighbors( pCNode, CNodes, CEdges, CQuads );
   CNodes.find( pCNextNode );
   CEdges.find( pCNextEdge );
   CQuads.find( pCQuad );
   iCallStatus = across_remove_face( pCNode, CNodes, CEdges, CQuads );
   //Check for and remove 2-valent nodes
   if( iCallStatus == fmu_kCLEAN_MESH_CHANGED )
   {
      iCallStatus = valence_2_check( CNodes );
      return fmu_kCLEAN_MESH_CHANGED;
   }

   return iCallStatus;
  
}


int QuadCleanTool::split_edge( fmuEdge *pCEdge )
{
   
   int iCallStatus;
   elmQClass *pCQuadA, *pCQuadB;
   nodClass *pCNodeA, *pCNodeB;
   nodClassDynArray CRingNodes(ARRAY_SIZE);
   nodClassDynArray CDelNodes(ARRAY_SIZE);
   fmuEdgeDynArray CDelEdges(ARRAY_SIZE);
   elmQClassDynArray CDelQuads(ARRAY_SIZE);

   if( pCEdge->frozen() ) return fmu_kCLEAN_NO_ACTION;

   pCEdge->quads( &pCQuadA, &pCQuadB );
   pCNodeA = pCEdge->start_node();
   pCNodeB = pCEdge->end_node();

   CDelQuads[0] = pCQuadA;
   CDelQuads[1] = pCQuadB;
   CDelEdges[0] = pCEdge;

   if( pCQuadA->prev_node( pCNodeA ) == pCNodeB )
   {
      CRingNodes[0] = pCNodeA;
      CRingNodes[1] = pCQuadA->next_node( pCNodeA );
      CRingNodes[2] = pCQuadA->opposite_node( pCNodeA );
      CRingNodes[3] = pCNodeB;
      CRingNodes[4] = pCQuadB->opposite_node( pCNodeA );
      CRingNodes[5] = pCQuadB->prev_node( pCNodeA );
   }
   else
   {
      CRingNodes[0] = pCNodeA;
      CRingNodes[1] = pCQuadB->next_node( pCNodeA );
      CRingNodes[2] = pCQuadB->opposite_node( pCNodeA );
      CRingNodes[3] = pCNodeB;
      CRingNodes[4] = pCQuadA->opposite_node( pCNodeA );
      CRingNodes[5] = pCQuadA->prev_node( pCNodeA );
   }

   if( get_valence( pCNodeA ) > get_valence( pCNodeB ) ) CRingNodes.step(3);

   iCallStatus =  fill_three( CDelNodes, CDelEdges, CDelQuads, CRingNodes );

   //Check for and remove 2-valent nodes
   if( iCallStatus == fmu_kCLEAN_MESH_CHANGED )
   {
      iCallStatus = valence_2_check( CRingNodes );
      return fmu_kCLEAN_MESH_CHANGED;
   }
   return iCallStatus;
  
}


int QuadCleanTool::fix_bad_angle( elmQClass *pCQuad, nodClass *pCNode )
{
   
   int iDebugFlag = 0;
   nodClass *pCPrevNode, *pCNextNode, *pCOppNode;
   int iPrevValence, iNextValence, iOppValence;
   int iCallStatus;

   pCPrevNode = pCQuad->prev_node( pCNode );
   pCNextNode = pCQuad->next_node( pCNode );
   pCOppNode = pCQuad->opposite_node( pCNode );

   if( iDebugFlag )
   {
   //pCNode->draw_me( dp_kCLR_DKGRN );
   //pCPrevNode->draw_me( dp_kCLR_LTMGNT );
   //pCNextNode->draw_me( dp_kCLR_BLUE );
   //pCOppNode->draw_me( dp_kCLR_YELLOW );
   }

   //Some poorly shaped quads are the best of a bad situation. These
   //all have an opposite node that is on the surface boundary. If so,
   //take a closer look.
   if( pCOppNode->boundary() )
   {
      //If the prev-opp edge is not frozen and the next angle is large,
      //allow the case to proceed. If that edge is not frozen and is
      //considerably longer than the next-opp edge, also allow it to
      //proceed. Also check next-opp edge and prev angle.
      //If both edges frozen, there is no action to take.
      bool lProcess = false;
      fmuEdge *pCEdgePO = pCPrevNode->shared_edge( pCOppNode );
      fmuEdge *pCEdgeNO = pCNextNode->shared_edge( pCOppNode );
      if( pCEdgePO->frozen() && pCEdgeNO->frozen() )
         return fmu_kCLEAN_NO_ACTION;
      fmuVector CVecNO( pCNextNode, pCOppNode );
      fmuVector CVecPO( pCPrevNode, pCOppNode );
      if( !pCEdgePO->frozen() )
      {
         fmuVector CVecNT( pCNextNode, pCNode );
         if( CNormal.angle( CVecNO, CVecNT ) > fmu_kSMALL_SHAPE_TOL )
            lProcess = true;
         if( CVecNO.length() * 3.0 < CVecPO.length() ) 
            lProcess = true;
      }
      if( !pCEdgeNO->frozen() )
      {
         fmuVector CVecPT( pCPrevNode, pCNode );
         if( CNormal.angle( CVecPT, CVecPO ) > fmu_kSMALL_SHAPE_TOL )
            lProcess = true;
         if( CVecPO.length() * 3.0 < CVecNO.length() )
            lProcess = true;
      }
      if( !lProcess ) return fmu_kCLEAN_NO_ACTION;
   }

   //Here is the big list of shape cases.

   // case shape squeeze
   iCallStatus = check_to_squeeze( pCQuad, pCNode, pCNextNode, pCOppNode,
                                                           pCPrevNode );
   if( iCallStatus != fmu_kCLEAN_NO_ACTION ) return iCallStatus;

   //Check for bowtie elements and do a simple fix (regardless of valence)
   fmuEdge *pCEdgeNO = pCNextNode->shared_edge(pCOppNode );
   fmuEdge *pCEdgeTP = pCNode->shared_edge( pCPrevNode );
   if( edges_intersect( pCEdgeNO, pCEdgeTP ) )
   {
      // case shape bowtie next
      iCallStatus = shape_rotate_next( pCQuad, pCNode, pCNextNode,
                                       pCOppNode, pCPrevNode, true );
      if( iCallStatus != fmu_kCLEAN_NO_ACTION ) return iCallStatus;
   }

   fmuEdge *pCEdgePO = pCPrevNode->shared_edge(pCOppNode );
   fmuEdge *pCEdgeTN = pCNode->shared_edge( pCNextNode );
   if( edges_intersect( pCEdgePO, pCEdgeTN ) )
   {
      // case shape bowtie prev
      iCallStatus = shape_rotate_prev( pCQuad, pCNode, pCNextNode,
                                       pCOppNode, pCPrevNode, true );
      if( iCallStatus != fmu_kCLEAN_NO_ACTION ) return iCallStatus;
   }

   iPrevValence = get_valence( pCPrevNode );
   iNextValence = get_valence( pCNextNode );
   iOppValence = get_valence( pCOppNode );

   if( iNextValence >= 5 &&
       iOppValence  >= 5 &&
       iPrevValence == 4 )
   {
      // case shape 554
      iCallStatus = shape_rotate_next( pCQuad, pCNode, pCNextNode,
                                       pCOppNode, pCPrevNode, false );
      if( iCallStatus != fmu_kCLEAN_NO_ACTION ) return iCallStatus;
   }

   if( iNextValence >= 5 &&
       iOppValence  == 4 &&
       iPrevValence == 4 &&
       !pCOppNode->hardpt() )
   {
      // case shape 544
      iCallStatus = shape_544( pCQuad, pCNode, pCNextNode, pCOppNode,
                                                           pCPrevNode );
      if( iCallStatus != fmu_kCLEAN_NO_ACTION ) return iCallStatus;
   }

   if( iNextValence == 4 &&
       iOppValence  >= 5 &&
       iPrevValence >= 5 )
   {
      // case shape 455
      iCallStatus = shape_rotate_prev( pCQuad, pCNode, pCNextNode,
                                       pCOppNode, pCPrevNode, false );
      if( iCallStatus != fmu_kCLEAN_NO_ACTION ) return iCallStatus;
   }

   if( iNextValence == 4 &&
       iOppValence  == 4 &&
       iPrevValence >= 5 &&
       !pCOppNode->hardpt() )
   {
      // case shape 445
      iCallStatus = shape_445( pCQuad, pCNode, pCNextNode, pCOppNode,
                                                           pCPrevNode );
      if( iCallStatus != fmu_kCLEAN_NO_ACTION ) return iCallStatus;
   }

   if( iNextValence == 4 &&
       iOppValence  >= 5 &&
       iPrevValence == 4 &&
       !pCOppNode->hardpt() )
   {
      // case shape 454
      iCallStatus = shape_454( pCQuad, pCNode, pCNextNode );
      if( iCallStatus != fmu_kCLEAN_NO_ACTION ) return iCallStatus;
   }

   if( iNextValence == 4 &&
       iOppValence  == 4 &&
       iPrevValence == 4 && 
       !pCOppNode->hardpt() )
   {
      // case shape 444
      iCallStatus = shape_444( pCQuad, pCNode, pCNextNode, pCOppNode,
                                                           pCPrevNode );
      if( iCallStatus != fmu_kCLEAN_NO_ACTION ) return iCallStatus;
   }

   if( iNextValence >= 4 &&
       iOppValence  == 3 &&
       iPrevValence >= 4 &&
       !pCOppNode->hardpt() )
   {
      // case shape 4+34+
      iCallStatus = shape_434( pCQuad, pCNode, pCNextNode, pCOppNode,
                                                           pCPrevNode );
      if( iCallStatus != fmu_kCLEAN_NO_ACTION ) return iCallStatus;
   }

   //This fix to the mesh assumes that the big angle may result from
   //smoothing. Change the nearby mesh so that smoothing doesn't
   //exert such a pull.
   if( !pCPrevNode->hardpt() )
   {
      iCallStatus = shape_smooth_prev( pCQuad, pCNode, pCPrevNode );
      if( iCallStatus != fmu_kCLEAN_NO_ACTION ) return iCallStatus;
   }
   else if( !pCNextNode->hardpt() )
   {
      iCallStatus = shape_smooth_next( pCQuad, pCNode, pCNextNode );
      if( iCallStatus != fmu_kCLEAN_NO_ACTION ) return iCallStatus;
   }

   // case last chance effort
   iCallStatus = shape_last_chance( pCQuad, pCNode, pCNextNode, pCOppNode,
                                                                pCPrevNode );
   if( iCallStatus != fmu_kCLEAN_NO_ACTION ) return iCallStatus;

   return fmu_kCLEAN_NO_ACTION;
  
}


bool QuadCleanTool::edges_intersect( fmuEdge *pCEdge1, fmuEdge *pCEdge2 )
{
   
   //This routine checks if two edges intersect. It does not compute
   //the actual intersection point. Label the endpoints as A and B from
   //edge 1 and C and D from edge 2. Create vectors AB and CD and unitized
   //versions. Do a dot product of the two vectors to check for parallel.
   //Return false if they are (don't worry about colinear).
   //Cross the AB vector with the normal to get vector N, which is
   //perpendicular to vector AB.
   //Project (dot) the original CD vector onto N and project vector AC onto N.
   //If the ratio of AC projection to CD projection is between zero & 1,
   //then C and D are on opposite sides of AC (since AC and CD are in
   //opposite directions, the negative ratio is used).
   //Project the same two vectors onto the AB vector. Watch out for
   //divide by zero.
   //Use the ratio and the projection of CD to get the distance along AB
   //from the intersection to the end of CD.
   //Add in the distance along AB of the AC vector.
   //Scale this by the length of AB and check for zero to one.
   double dPARALLEL_TOL = 0.99999;
   double dCos;
   fmuVector CvecA( pCEdge1->start_node() );
   fmuVector CvecB( pCEdge1->end_node() );
   fmuVector CvecC( pCEdge2->start_node() );
   fmuVector CvecD( pCEdge2->end_node() );
   fmuVector CvecAB = CvecB - CvecA;
   fmuVector CvecABunit = CvecAB;
   CvecABunit.unitize();
   fmuVector CvecCD = CvecD - CvecC;
   fmuVector CvecCDunit = CvecCD;
   CvecCDunit.unitize();
   dCos = CvecABunit % CvecCDunit;
   if( dCos > dPARALLEL_TOL || dCos < -dPARALLEL_TOL ) return false;

   fmuVector CvecN = CvecABunit * CNormal;
   fmuVector CvecAC = CvecC - CvecA;
   double dACProject = CvecAC % CvecN;
   double dCDProject = CvecCD % CvecN;
   double dCDratio = (-dACProject) / dCDProject;
   if( dCDratio < 0.0 || dCDratio > 1.0 ) return false;
   dACProject = CvecAC % CvecABunit;
   dCDProject = CvecCD % CvecABunit;
   double dABratio = ( dCDratio * dCDProject + dACProject ) /
                                                            CvecAB.length();
   if( dABratio < 0.0 || dABratio > 1.0 ) return false;
   return true;
   
}


int QuadCleanTool::check_to_squeeze( elmQClass *pCQuad, nodClass *pCNode,
                                         nodClass *pCNextNode,
                                         nodClass *pCOppNode,
                                         nodClass *pCPrevNode )
{
   
   int iCallStatus;
   int iPrevValence, iNextValence;
   double dAngle;
   nodClassDynArray CNodes(ARRAY_SIZE);
   fmuEdgeDynArray CEdges(ARRAY_SIZE);
   elmQClassDynArray CQuads(ARRAY_SIZE);

   if( pCNode->hardpt() && pCOppNode->hardpt() ) return fmu_kCLEAN_NO_ACTION;

   iPrevValence = get_valence( pCPrevNode );
   if( pCPrevNode->hardpt() && iPrevValence == 3 )
      return fmu_kCLEAN_NO_ACTION;
   iNextValence = get_valence( pCNextNode );
   if( pCNextNode->hardpt() && iNextValence == 3 )
      return fmu_kCLEAN_NO_ACTION;

   fmuVector CVecNO( pCNextNode, pCOppNode );
   fmuVector CVecNT( pCNextNode, pCNode );
   dAngle = CNormal.angle( CVecNO, CVecNT );
   if( dAngle > fmu_kSMALL_SQUEEZE_ANGLE &&
       dAngle < fmu_kBIG_SQUEEZE_ANGLE ) return fmu_kCLEAN_NO_ACTION;

   fmuVector CVecPO( pCPrevNode, pCOppNode );
   fmuVector CVecPT( pCPrevNode, pCNode );
   dAngle = CNormal.angle( CVecPO, CVecPT );
   if( dAngle > fmu_kSMALL_SQUEEZE_ANGLE  &&
       dAngle < fmu_kBIG_SQUEEZE_ANGLE ) return fmu_kCLEAN_NO_ACTION;

   if( !pCOppNode->hardpt() )
   {
      neighbors( pCOppNode, CNodes, CEdges, CQuads );
      CNodes.find( pCPrevNode );
      CEdges.find( pCPrevNode->shared_edge( pCOppNode ) );
      CQuads.find( pCQuad );
      //Check node after next for bad angles that might result
      nodClass *pCNode3 = CNodes.next(3);

	  // If node not found get out
	  if (!pCNode3)
		  return(fmu_kCLEAN_NO_ACTION);

      if( pCNode3->hardpt() && pCNextNode->hardpt() && pCNode->hardpt() )
      {
         fmuVector CVecN3( pCNextNode, pCNode3 );
         dAngle = CNormal.angle( CVecN3, CVecNT );
         if( dAngle >  fmu_kSHAPE_ANGLE_TOL ) return fmu_kCLEAN_NO_ACTION;
      }
      //Check the node before prev (prev is now at position zero in CNodes)
      nodClass *pCNodeA = CNodes.next(-1);
      if( pCNodeA->hardpt() && pCPrevNode->hardpt() && pCNode->hardpt() )
      {
         fmuVector CVecPA( pCPrevNode, pCNodeA );
         dAngle = CNormal.angle( CVecPT, CVecPA );
         if( dAngle >  fmu_kSHAPE_ANGLE_TOL ) return fmu_kCLEAN_NO_ACTION;
      }
      iCallStatus = across_remove_face( pCOppNode, CNodes, CEdges, CQuads );
      if( iCallStatus == fmu_kCLEAN_MESH_CHANGED )
      {
         //valences at prev and next are now reduced by one but the
         //saved variable was not. The mesh has already been changed so
         //only worry about fatal problems
         if( iPrevValence == 3 )
         {
            iCallStatus = valence_2_clean( pCPrevNode );
            if( iCallStatus == fmu_kCLEAN_NO_MEMORY )
                                         return fmu_kCLEAN_NO_MEMORY;
            if( iCallStatus == fmu_kCLEAN_ABORT )
                                         return fmu_kCLEAN_ABORT;
         }
         if( iNextValence == 3 )
         {
            iCallStatus = valence_2_clean( pCNextNode );
            if( iCallStatus == fmu_kCLEAN_NO_MEMORY )
                                         return fmu_kCLEAN_NO_MEMORY;
            if( iCallStatus == fmu_kCLEAN_ABORT )
                                         return fmu_kCLEAN_ABORT;
         }
         return fmu_kCLEAN_MESH_CHANGED;
      }
   }

   return fmu_kCLEAN_NO_ACTION;

}


int QuadCleanTool::valence_2_check( nodClassDynArray &CRingNodes )
{

   //This routine will have problems if a 2 valent node creates a 2 valent
   //node that is also in the ring.
   int ii;
   int iCallStatus;
   int iStatus = fmu_kCLEAN_NO_ACTION;
   nodClass *pCNode;
   for( ii = 0; ii < 6; ++ii )
   {
      pCNode = CRingNodes[ii];
	  if (!pCNode)
		  return (iStatus);

      if( !pCNode->boundary() )
      {
         if( get_valence( pCNode ) == 2 )
         {
            iCallStatus =  valence_2_clean( pCNode );
            if( iCallStatus == fmu_kCLEAN_ABORT ) return fmu_kCLEAN_ABORT;
            if( iCallStatus == fmu_kCLEAN_NO_MEMORY )
               return fmu_kCLEAN_NO_MEMORY;
            if( iCallStatus == fmu_kCLEAN_MESH_CHANGED )
               iStatus = iCallStatus;
	 }
      }
   }
   return iStatus;

}


int QuadCleanTool::valence_2_clean( nodClass *pCNode )
{

   int iValence, iCallStatus;
   nodClassDynArray CNodes(10);
   fmuEdgeDynArray CEdges(10);
   elmQClassDynArray CQuads(10);
   elmQClass *pCQuad;

   //Can't handle boundary problems here.
   if( pCNode->boundary() ) return fmu_kCLEAN_NO_ACTION;
   //Need a special routine if the problem node is a hardpoint.
   if( pCNode->hardpt() ) return valence_2_hardpt( pCNode );
   iValence = neighbors( pCNode, CNodes, CEdges, CQuads );
   nodClass *pCNode0 = CNodes.get();
   nodClass *pCNode1 = CNodes.next();
   nodClass *pCNode2 = CNodes.next(2);
   nodClass *pCNode3 = CNodes.next(3);
   if(iValence != 2 ) return fmu_kCLEAN_NO_ACTION;
   fmuEdge *edge0 = CEdges.get();
   fmuEdge *edge1 = CEdges.next();
   if( edge0->frozen() ) return fmu_kCLEAN_NO_ACTION;
   if( edge1->frozen() ) return fmu_kCLEAN_NO_ACTION;

   if( lDrawCleaning )
   {
   //draw_open();
   //CQuads.get()->draw_me( dp_kCLR_LTBLUE );
   //CQuads.next()->draw_me( dp_kCLR_GOLD );
        // Refresh debug graphics
        //EraseElements();
        CQuads.get()->draw_me( dp_kCLR_YELLOW );
        CQuads.next()->draw_me( dp_kCLR_YELLOW );
   //draw_close();
   }

   delete_this_quad( CQuads.get() );
   delete_this_quad( CQuads.next() );

   if( delete_this_edge( edge0 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( edge1 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;

   if( delete_this_node( pCNode ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;

   pCQuad = new elmQClass( pCNode0, pCNode1, pCNode2, pCNode3, pCMaster );
   if( !pCQuad ) return fmu_kCLEAN_NO_MEMORY;

   smooth_mesh( &CNodes );
   if( lDrawCleaning )
   {
   //draw_mesh();
   //pCQuad->draw_me( dp_kCLR_LTBLUE );
        // Refresh debug graphics
        //EraseElements();
        pCQuad->draw_me( dp_kCLR_BLUE );
   //draw_close();
   }

   if( !pCNode0->boundary() && get_valence( pCNode0 ) == 2 )
   {
      iCallStatus = valence_2_clean( pCNode0 );
      if( iCallStatus == fmu_kCLEAN_ABORT ) return fmu_kCLEAN_ABORT;
      if( iCallStatus == fmu_kCLEAN_NO_MEMORY ) return fmu_kCLEAN_NO_MEMORY;
   }

   if( !pCNode2->boundary() && get_valence( pCNode2 ) == 2 )
   {
      iCallStatus = valence_2_clean( pCNode2 );
      if( iCallStatus == fmu_kCLEAN_ABORT ) return fmu_kCLEAN_ABORT;
      if( iCallStatus == fmu_kCLEAN_NO_MEMORY ) return fmu_kCLEAN_NO_MEMORY;
   }

   return fmu_kCLEAN_MESH_CHANGED;

}


int QuadCleanTool::valence_2_hardpt( nodClass *pCNode )
{

   int iCallStatus;
   nodClassDynArray CNodes(10);
   nodClassDynArray CRingNodes(10);
   fmuEdgeDynArray CEdges(10);
   elmQClassDynArray CQuads(10);
   elmQClass *pCQuadA;
   nodClass *pCNodeA, *pCNodeB;

   neighbors( pCNode, CNodes, CEdges, CQuads );
   nodClass *pCNode0 = CNodes.get();
   nodClass *pCNode1 = CNodes.next();
   nodClass *pCNode2 = CNodes.next(2);
   nodClass *pCNode3 = CNodes.next(3);
   fmuEdge *pCEdge01 = pCNode0->shared_edge( pCNode1 );
   fmuEdge *pCEdge12 = pCNode1->shared_edge( pCNode2 );
   fmuEdge *pCEdge23 = pCNode2->shared_edge( pCNode3 );
   fmuEdge *pCEdge03 = pCNode0->shared_edge( pCNode3 );

   //If the center node (the hardpoint) is closer to node 3 than node 1 or
   //both outer edges near node 1 are frozen, rotate the orientation.
   fmuVector CVec0( pCNode, pCNode0 );
   fmuVector CVec2( pCNode, pCNode2 );
   if( (CNormal.angle( CVec0, CVec2 ) <= PI &&
       ( !pCEdge23->frozen() || !pCEdge03->frozen() ) ) ||
       ( pCEdge01->frozen() && pCEdge12->frozen() ) )
   {
      CNodes.step(2);
      CEdges.step();
      CQuads.step();
      pCNode0 = CNodes.get();
      pCNode1 = CNodes.next();
      pCNode2 = CNodes.next(2);
      pCEdge01 = pCNode0->shared_edge( pCNode1 );
      pCEdge12 = pCNode1->shared_edge( pCNode2 );
   }

   fmuEdge *pCEdge0 = CEdges.get();
   fmuEdge *pCEdge1 = CEdges.next();
   if( pCEdge01->frozen() &&  pCEdge12->frozen() ) return fmu_kCLEAN_NO_ACTION;
   elmQClass *pCQuad0 = CQuads.get();

   int iSide;
   if( pCEdge01->frozen() ) iSide = 2;
   else if( pCEdge12->frozen() ) iSide = 0;
   else if( edges_intersect( pCEdge1, pCEdge01 ) ) iSide = 0;
   else if( edges_intersect( pCEdge0, pCEdge12 ) ) iSide = 2;
   else
   {
      fmuVector CVec0C( pCNode0, pCNode );
      fmuVector CVec2C( pCNode2, pCNode );
      fmuVector CVec01( pCNode0, pCNode1 );
      fmuVector CVec21( pCNode2, pCNode1 );

      if( CNormal.angle( CVec01, CVec0C ) <
          CNormal.angle( CVec2C, CVec21 ) ) iSide = 0;
      else                                  iSide = 2;
   }

   CQuads.clear();
   CEdges.clear();
   CNodes.clear();
   if( iSide == 0 )
   {
      pCQuadA = pCEdge01->other_quad( pCQuad0 );
      pCNodeA = pCQuadA->next_node( pCNode0 );
      pCNodeB = pCQuadA->opposite_node( pCNode0 );
      //Quads that will be deleted
      CQuads[0] = pCQuad0;
      CQuads[1] = pCQuadA;
      //Edges that will be deleted
      CEdges[0] = pCEdge01;
      //No nodes to delete
      //Nodes of outer ring
      CRingNodes[0] = pCNode;
      CRingNodes[1] = pCNode0;
      CRingNodes[2] = pCNodeA;
      CRingNodes[3] = pCNodeB;
      CRingNodes[4] = pCNode1;
      CRingNodes[5] = pCNode2;

      iCallStatus = fill_two( CNodes, CEdges, CQuads, CRingNodes );
   }
   else
   {
      pCQuadA = pCEdge12->other_quad( pCQuad0 );
      pCNodeA = pCQuadA->prev_node( pCNode2 );
      pCNodeB = pCQuadA->opposite_node( pCNode2 );
      //Quads that will be deleted
      CQuads[0] = pCQuad0;
      CQuads[1] = pCQuadA;
      //Edges that will be deleted
      CEdges[0] = pCEdge12;
      //No nodes to delete
      //Nodes of outer ring
      CRingNodes[0] = pCNode;
      CRingNodes[1] = pCNode0;
      CRingNodes[2] = pCNode1;
      CRingNodes[3] = pCNodeB;
      CRingNodes[4] = pCNodeA;
      CRingNodes[5] = pCNode2;

      iCallStatus = fill_two( CNodes, CEdges, CQuads, CRingNodes );
   }
   return iCallStatus;

}


int QuadCleanTool::shape_rotate_next( elmQClass *pCQuad, nodClass *pCNode,
                                          nodClass *pCNextNode,
                                          nodClass *pCOppNode,
                                          nodClass *pCPrevNode,
                                          bool lForce )
{

   int iCallStatus;
   elmQClassDynArray CQuads;
   nodClassDynArray CNodes, CRingNodes;
   fmuEdgeDynArray CEdges;

   fmuEdge *pCEdgeNO = pCNextNode->shared_edge( pCOppNode );
   if( pCEdgeNO->frozen() ) return fmu_kCLEAN_NO_ACTION;
   elmQClass *pCQuadA = pCEdgeNO->other_quad( pCQuad );
   if (!pCQuadA)return fmu_kCLEAN_NO_ACTION;
   nodClass *pCNodeA = pCQuadA->next_node( pCNextNode );
   nodClass *pCNodeB = pCQuadA->opposite_node( pCNextNode );
   if( !lForce )
   {
      if( get_valence( pCNodeB ) > 4 ) return fmu_kCLEAN_NO_ACTION;
      fmuVector CVecNA( pCNextNode, pCNodeA );
      fmuVector CVecNT( pCNextNode, pCNode );
      if( CNormal.angle( CVecNA, CVecNT ) > fmu_kSHAPE_ANGLE_TOL )
         return fmu_kCLEAN_NO_ACTION;
   }

   //Quads that will be deleted
   CQuads[0] = pCQuad;
   CQuads[1] = pCQuadA;
   //Edges that will be deleted
   CEdges[0] = pCEdgeNO;
   //No nodes to delete
   //Nodes of outer ring
   CRingNodes[0] = pCNode;
   CRingNodes[1] = pCNextNode;
   CRingNodes[2] = pCNodeA;
   CRingNodes[3] = pCNodeB;
   CRingNodes[4] = pCOppNode;
   CRingNodes[5] = pCPrevNode;

   iCallStatus = fill_two( CNodes, CEdges, CQuads, CRingNodes );
   if( iCallStatus != fmu_kCLEAN_MESH_CHANGED ) return iCallStatus;
   if( get_valence( pCOppNode ) == 2 )
   {
      iCallStatus = valence_2_clean( pCOppNode );
      if( iCallStatus == fmu_kCLEAN_ABORT ) return fmu_kCLEAN_ABORT;
      if( iCallStatus == fmu_kCLEAN_NO_MEMORY ) return fmu_kCLEAN_NO_MEMORY;
   }
   if( get_valence( pCNextNode ) == 2 )
   {
      iCallStatus = valence_2_clean( pCNextNode );
      if( iCallStatus == fmu_kCLEAN_ABORT ) return fmu_kCLEAN_ABORT;
      if( iCallStatus == fmu_kCLEAN_NO_MEMORY ) return fmu_kCLEAN_NO_MEMORY;
   }

   return fmu_kCLEAN_MESH_CHANGED;

}


int QuadCleanTool::shape_544( elmQClass *pCQuad, nodClass *pCNode,
                                  nodClass *pCNextNode, nodClass *pCOppNode,
                                  nodClass *pCPrevNode )
{

   int iStatus, iCallStatus;
   elmQClassDynArray CQuads;
   nodClassDynArray CNodes, CRingNodes;
   fmuEdgeDynArray CEdges;

   fmuEdge *pCEdgeNO = pCNextNode->shared_edge( pCOppNode );
   if( pCEdgeNO->frozen() ) return fmu_kCLEAN_NO_ACTION;
   elmQClass *pCQuadA = pCEdgeNO->other_quad( pCQuad );

   // There might be no quad if the other element is a tria.
   if (!pCQuadA) return fmu_kCLEAN_NO_ACTION;
   nodClass *pCNodeA = pCQuadA->next_node( pCNextNode );
   nodClass *pCNodeB = pCQuadA->opposite_node( pCNextNode );

   //An alternate pattern if the fill 3 will only move the bad angle.
   //Do the easy checks first, then the more expensive ones.
   if( pCNode->hardpt() &&
       pCNextNode->hardpt() &&
       !pCPrevNode->hardpt() &&
       !pCNodeA->hardpt() &&
       get_valence( pCPrevNode ) > 3 &&
       get_valence( pCNodeA ) > 3 &&
       get_valence( pCNodeB ) == 3 )
   {
      double dAngle;
      fmuVector CVecNCen( pCNextNode, pCNode );
      fmuVector CVecNA( pCPrevNode, pCNodeA );
      dAngle = CNormal.angle( CVecNA, CVecNCen );
      if( dAngle > fmu_kSHAPE_ANGLE_TOL )
      {
         fmuEdge *pCEdgePO = pCPrevNode->shared_edge( pCOppNode );
         fmuEdge *pCEdgeOB = pCOppNode->shared_edge( pCNodeB );
         fmuEdge *pCEdgeAB = pCNodeA->shared_edge( pCNodeB );
         if( !pCEdgePO->frozen() && !pCEdgeOB->frozen() &&
                                    !pCEdgeAB->frozen() )
         {
            elmQClass *pCQuadB = pCEdgeAB->other_quad( pCQuadA );
            if (!pCQuadB) return fmu_kCLEAN_NO_ACTION;
            elmQClass *pCQuadC = pCEdgeOB->other_quad( pCQuadA );
            if (!pCQuadC) return fmu_kCLEAN_NO_ACTION;
            elmQClass *pCQuadD = pCEdgePO->other_quad( pCQuad );
            if (!pCQuadD) return fmu_kCLEAN_NO_ACTION;
            nodClass *pCNodeC = pCQuadB->opposite_node( pCNodeB );
            nodClass *pCNodeD = pCQuadB->prev_node( pCNodeB );
            nodClass *pCNodeE = pCQuadD->next_node( pCOppNode );
            nodClass *pCNodeF = pCQuadD->opposite_node( pCOppNode );
            fmuEdge *pCEdgeBD = pCNodeB->shared_edge( pCNodeD );
            fmuEdge *pCEdgeOE = pCOppNode->shared_edge( pCNodeE );
            if( !pCEdgeBD->frozen() && !pCEdgeOE->frozen() )
            {
               if( lDrawCleaning )
               {
               //draw_open();
               //pCQuad->draw_me( dp_kCLR_RED );
               //pCQuadA->draw_me( dp_kCLR_LTMGNT );
               //pCQuadB->draw_me( dp_kCLR_BLUE );
               //pCQuadC->draw_me( dp_kCLR_DKGRN );
               //pCQuadD->draw_me( dp_kCLR_PINK );
                    // Refresh debug graphics
                    //EraseElements();
                    pCQuad->draw_me( dp_kCLR_YELLOW );
                    pCQuadA->draw_me( dp_kCLR_YELLOW );
                    pCQuadB->draw_me( dp_kCLR_YELLOW );
                    pCQuadC->draw_me( dp_kCLR_YELLOW );
                    pCQuadD->draw_me( dp_kCLR_YELLOW );
               //draw_close();
               }

               delete_this_quad( pCQuad );
               delete_this_quad( pCQuadA );
               delete_this_quad( pCQuadB );
               delete_this_quad( pCQuadC );
               delete_this_quad( pCQuadD );
               if( delete_this_edge( pCEdgeNO ) == fmu_kCLEAN_ABORT )
                  return fmu_kCLEAN_ABORT;
               if( delete_this_edge( pCEdgePO ) == fmu_kCLEAN_ABORT )
                  return fmu_kCLEAN_ABORT;
               if( delete_this_edge( pCEdgeAB ) == fmu_kCLEAN_ABORT )
                  return fmu_kCLEAN_ABORT;
               if( delete_this_edge( pCEdgeOB ) == fmu_kCLEAN_ABORT )
                  return fmu_kCLEAN_ABORT;
               if( delete_this_edge( pCEdgeBD ) == fmu_kCLEAN_ABORT )
                  return fmu_kCLEAN_ABORT;
               if( delete_this_edge( pCEdgeOE ) == fmu_kCLEAN_ABORT )
                  return fmu_kCLEAN_ABORT;
               if( delete_this_node( pCOppNode ) == fmu_kCLEAN_ABORT )
                  return fmu_kCLEAN_ABORT;
               if( delete_this_node( pCNodeB ) == fmu_kCLEAN_ABORT )
                  return fmu_kCLEAN_ABORT;

               pCQuad =  new elmQClass( pCPrevNode, pCNode, pCNodeE,
                                                pCNodeF, pCMaster ); 
               if( !pCQuad ) return fmu_kCLEAN_NO_MEMORY;
               pCQuadA = new elmQClass( pCNode, pCNextNode, pCNodeD,
                                                pCNodeE, pCMaster ); 
               if( !pCQuadA ) return fmu_kCLEAN_NO_MEMORY;
               pCQuadB = new elmQClass( pCNextNode, pCNodeA, pCNodeC,
                                                pCNodeD, pCMaster ); 
               if( !pCQuadB ) return fmu_kCLEAN_NO_MEMORY;
               nodClassDynArray CSmoothList;
               CSmoothList.append( pCPrevNode );
               CSmoothList.append( pCNodeA );
               CSmoothList.append( pCNodeC );
               CSmoothList.append( pCNodeD );
               CSmoothList.append( pCNodeE );
               CSmoothList.append( pCNodeF );
               smooth_mesh( &CSmoothList );
               if( lDrawCleaning )
               {
               //draw_mesh();
               //pCQuad->draw_me( dp_kCLR_RED );
               //pCQuadA->draw_me( dp_kCLR_DKGRN );
               //pCQuadB->draw_me( dp_kCLR_LTMGNT );
                    // Refresh debug graphics
                    //EraseElements();
                    pCQuad->draw_me( dp_kCLR_BLUE );
                    pCQuadA->draw_me( dp_kCLR_BLUE );
                    pCQuadB->draw_me( dp_kCLR_BLUE );
                  
               //draw_close();
               }

               return fmu_kCLEAN_MESH_CHANGED;

            }
         }
      }
   }

   //The special pattern wasn't done, try the standard fix

   if( pCPrevNode->hardpt() && ( pCOppNode->hardpt() || pCNodeB->hardpt() ) )
   {
      fmuVector CVecOP( pCOppNode, pCPrevNode );
      fmuVector CVecOB( pCOppNode, pCNodeB );
      if( CNormal.angle( CVecOP, CVecOB ) > fmu_kSHAPE_ANGLE_TOL )
         return fmu_kCLEAN_NO_ACTION;
   }

   if( pCNode->hardpt() && pCNextNode->hardpt() && pCNodeA->hardpt() )
   {
      fmuVector CVecNA( pCNextNode, pCNodeA );
      fmuVector CVecNT( pCNextNode, pCNode );
      if( CNormal.angle( CVecNA, CVecNT ) > fmu_kSHAPE_ANGLE_TOL )
         return fmu_kCLEAN_NO_ACTION;
   }

   //Quads that will be deleted
   CQuads[0] = pCQuad;
   CQuads[1] = pCQuadA;
   //Edges that will be deleted
   CEdges[0] = pCEdgeNO;
   //No nodes to delete
   //Nodes of outer ring
   CRingNodes[0] = pCNode;
   CRingNodes[1] = pCNextNode;
   CRingNodes[2] = pCNodeA;
   CRingNodes[3] = pCNodeB;
   CRingNodes[4] = pCOppNode;
   CRingNodes[5] = pCPrevNode;

   if( get_valence( pCNodeA ) > 4 ||
       get_valence( pCNodeB ) == 3 )
      iCallStatus = fill_two( CNodes, CEdges, CQuads, CRingNodes );
   else
      iCallStatus = fill_three( CNodes, CEdges, CQuads, CRingNodes );

   if( iCallStatus == fmu_kCLEAN_MESH_CHANGED )
   {
      iStatus = valence_2_check( CRingNodes );
      if( iStatus == fmu_kCLEAN_ABORT ) return fmu_kCLEAN_ABORT;
      if( iStatus == fmu_kCLEAN_NO_MEMORY )
         return fmu_kCLEAN_NO_MEMORY;
   }

   return iCallStatus;

}


int QuadCleanTool::shape_rotate_prev( elmQClass *pCQuad, nodClass *pCNode,
                                          nodClass *pCNextNode,
                                          nodClass *pCOppNode,
                                          nodClass *pCPrevNode,
                                          bool lForce )
{

   int iCallStatus;
   elmQClassDynArray CQuads;
   nodClassDynArray CNodes, CRingNodes;
   fmuEdgeDynArray CEdges;

   fmuEdge *pCEdgePO = pCPrevNode->shared_edge( pCOppNode );
   if( pCEdgePO->frozen() ) return fmu_kCLEAN_NO_ACTION;
   elmQClass *pCQuadA = pCEdgePO->other_quad( pCQuad );
   if (!pCQuadA) return fmu_kCLEAN_NO_ACTION;
   nodClass *pCNodeA = pCQuadA->prev_node( pCPrevNode );
   nodClass *pCNodeB = pCQuadA->opposite_node( pCPrevNode );
   if( !lForce )
   {
      if( get_valence( pCNodeB ) > 4 ) return fmu_kCLEAN_NO_ACTION;
      fmuVector CVecPT( pCPrevNode, pCNode );
      fmuVector CVecPA( pCPrevNode, pCNodeA );
      if( CNormal.angle( CVecPT, CVecPA ) > fmu_kSHAPE_ANGLE_TOL )
         return fmu_kCLEAN_NO_ACTION;
   }

   //Quads that will be deleted
   CQuads[0] = pCQuad;
   CQuads[1] = pCQuadA;
   //Edges that will be deleted
   CEdges[0] = pCEdgePO;
   //No nodes to delete
   //Nodes of outer ring
   CRingNodes[0] = pCNode;
   CRingNodes[1] = pCNextNode;
   CRingNodes[2] = pCOppNode;
   CRingNodes[3] = pCNodeB;
   CRingNodes[4] = pCNodeA;
   CRingNodes[5] = pCPrevNode;

   iCallStatus = fill_two( CNodes, CEdges, CQuads, CRingNodes );
   if( iCallStatus != fmu_kCLEAN_MESH_CHANGED ) return iCallStatus;
   if( get_valence( pCOppNode ) == 2 )
   {
      iCallStatus = valence_2_clean( pCOppNode );
      if( iCallStatus == fmu_kCLEAN_ABORT ) return fmu_kCLEAN_ABORT;
      if( iCallStatus == fmu_kCLEAN_NO_MEMORY ) return fmu_kCLEAN_NO_MEMORY;
   }
   if( get_valence( pCPrevNode ) == 2 )
   {
      iCallStatus = valence_2_clean( pCPrevNode );
      if( iCallStatus == fmu_kCLEAN_ABORT ) return fmu_kCLEAN_ABORT;
      if( iCallStatus == fmu_kCLEAN_NO_MEMORY ) return fmu_kCLEAN_NO_MEMORY;
   }

   return fmu_kCLEAN_MESH_CHANGED;

}


int QuadCleanTool::shape_445( elmQClass *pCQuad, nodClass *pCNode,
                                  nodClass *pCNextNode, nodClass *pCOppNode,
                                  nodClass *pCPrevNode )
{
    int iStatus, iCallStatus;
    elmQClassDynArray CQuads;
    nodClassDynArray CNodes, CRingNodes;
    fmuEdgeDynArray CEdges;

    fmuEdge *pCEdgePO = pCPrevNode->shared_edge( pCOppNode );
    if( pCEdgePO == NULL )
        return fmu_kCLEAN_NO_ACTION;
    if( pCEdgePO->frozen() )
        return fmu_kCLEAN_NO_ACTION;

    elmQClass *pCQuadA = pCEdgePO->other_quad( pCQuad );
    if( pCQuadA == NULL )
        return fmu_kCLEAN_NO_ACTION;

    nodClass *pCNodeA = pCQuadA->prev_node( pCPrevNode );
    nodClass *pCNodeB = pCQuadA->opposite_node( pCPrevNode );

    //An alternate pattern if the fill 3 will only move the bad angle.
    //Do the easy checks first, then the more expensive ones.
    if( pCNode->hardpt() &&
        pCPrevNode->hardpt() &&
        !pCNextNode->hardpt() &&
        !pCNodeA->hardpt() &&
        get_valence( pCNextNode ) > 3 &&
        get_valence( pCNodeA ) > 3 &&
        get_valence( pCNodeB ) == 3 )
    {
        double dAngle;
        fmuVector CVecPCen( pCPrevNode, pCNode );
        fmuVector CVecPA( pCPrevNode, pCNodeA );
        dAngle = CNormal.angle( CVecPCen, CVecPA );
        if( dAngle > fmu_kSHAPE_ANGLE_TOL )
        {
            fmuEdge *pCEdgeNO = pCNextNode->shared_edge( pCOppNode );
            fmuEdge *pCEdgeOB = pCOppNode->shared_edge( pCNodeB );
            fmuEdge *pCEdgeAB = pCNodeA->shared_edge( pCNodeB );
            if( !pCEdgeNO->frozen() && 
                !pCEdgeOB->frozen() &&
                !pCEdgeAB->frozen() )
            {
                elmQClass *pCQuadB = pCEdgeNO->other_quad( pCQuad );
                elmQClass *pCQuadC = pCEdgeOB->other_quad( pCQuadA );
                elmQClass *pCQuadD = pCEdgeAB->other_quad( pCQuadA );
                nodClass *pCNodeC = pCQuadB->opposite_node( pCOppNode );
                nodClass *pCNodeD = pCQuadB->prev_node( pCOppNode );
                nodClass *pCNodeE = pCQuadC->prev_node( pCNodeB );
                nodClass *pCNodeF = pCQuadD->opposite_node( pCNodeB );
                fmuEdge *pCEdgeOD = pCOppNode->shared_edge( pCNodeD );
                fmuEdge *pCEdgeBE = pCNodeB->shared_edge( pCNodeE );
                if( !pCEdgeOD->frozen() && !pCEdgeBE->frozen() )
                {
                    if( lDrawCleaning )
                    {
                        //draw_open();
                        //pCQuad->draw_me( dp_kCLR_RED );
                        //pCQuadA->draw_me( dp_kCLR_MAGNTA );
                        //pCQuadB->draw_me( dp_kCLR_BLUE );
                        //pCQuadC->draw_me( dp_kCLR_DKGRN );
                        //pCQuadD->draw_me( dp_kCLR_PINK );
                        // Refresh debug graphics
                        //EraseElements();
                        pCQuad->draw_me( dp_kCLR_YELLOW );
                        pCQuadA->draw_me( dp_kCLR_YELLOW );
                        pCQuadB->draw_me( dp_kCLR_YELLOW );
                        pCQuadC->draw_me( dp_kCLR_YELLOW );
                        pCQuadD->draw_me( dp_kCLR_YELLOW );
                        //draw_close();
                    }

                    delete_this_quad( pCQuad );
                    delete_this_quad( pCQuadA );
                    delete_this_quad( pCQuadB );
                    delete_this_quad( pCQuadC );
                    delete_this_quad( pCQuadD );
                    if( delete_this_edge( pCEdgeNO ) == fmu_kCLEAN_ABORT )
                        return fmu_kCLEAN_ABORT;
                    if( delete_this_edge( pCEdgePO ) == fmu_kCLEAN_ABORT )
                        return fmu_kCLEAN_ABORT;
                    if( delete_this_edge( pCEdgeAB ) == fmu_kCLEAN_ABORT )
                        return fmu_kCLEAN_ABORT;
                    if( delete_this_edge( pCEdgeOB ) == fmu_kCLEAN_ABORT )
                        return fmu_kCLEAN_ABORT;
                    if( delete_this_edge( pCEdgeOD ) == fmu_kCLEAN_ABORT )
                        return fmu_kCLEAN_ABORT;
                    if( delete_this_edge( pCEdgeBE ) == fmu_kCLEAN_ABORT )
                        return fmu_kCLEAN_ABORT;
                    if( delete_this_node( pCOppNode ) == fmu_kCLEAN_ABORT )
                        return fmu_kCLEAN_ABORT;
                    if( delete_this_node( pCNodeB ) == fmu_kCLEAN_ABORT )
                        return fmu_kCLEAN_ABORT;

                    pCQuad = new elmQClass( pCNode, pCNextNode, pCNodeC,
                                            pCNodeD, pCMaster ); 
                    if( !pCQuad )
                        return fmu_kCLEAN_NO_MEMORY;

                    pCQuadA = new elmQClass( pCNode, pCNodeD, pCNodeE,
                                            pCPrevNode, pCMaster ); 
                    if( !pCQuadA )
                        return fmu_kCLEAN_NO_MEMORY;

                    pCQuadB = new elmQClass( pCPrevNode, pCNodeE, pCNodeF,
                                            pCNodeA, pCMaster ); 
                    if( !pCQuadB )
                        return fmu_kCLEAN_NO_MEMORY;

                    nodClassDynArray CSmoothList;
                    CSmoothList.append( pCNextNode );
                    CSmoothList.append( pCNodeA );
                    CSmoothList.append( pCNodeC );
                    CSmoothList.append( pCNodeD );
                    CSmoothList.append( pCNodeE );
                    CSmoothList.append( pCNodeF );
                    smooth_mesh( &CSmoothList );
                    if( lDrawCleaning )
                    {
                        //draw_mesh();
                        //pCQuad->draw_me( dp_kCLR_RED );
                        //pCQuadA->draw_me( dp_kCLR_DKGRN );
                        //pCQuadB->draw_me( dp_kCLR_PINK );
                        // Refresh debug graphics
                        //EraseElements();
                        pCQuad->draw_me( dp_kCLR_BLUE );
                        pCQuadA->draw_me( dp_kCLR_BLUE );
                        pCQuadB->draw_me( dp_kCLR_BLUE );
                   
                        //draw_close();
                    }

                    return fmu_kCLEAN_MESH_CHANGED;
                }
            }
        }
    }

    //The special pattern wasn't done, try the standard fix

    if( pCNextNode->hardpt() && ( pCOppNode->hardpt() || pCNodeB->hardpt() ) )
    {
        fmuVector CVecOB( pCOppNode, pCNodeB );
        fmuVector CVecON( pCOppNode, pCNextNode );
        if( CNormal.angle( CVecOB, CVecON ) > fmu_kSHAPE_ANGLE_TOL )
            return fmu_kCLEAN_NO_ACTION;
    }

    if( pCNode->hardpt() && pCPrevNode->hardpt() && pCNodeA->hardpt() )
    {
        fmuVector CVecPT( pCPrevNode, pCNode );
        fmuVector CVecPA( pCPrevNode, pCNodeA );
        if( CNormal.angle( CVecPT, CVecPA ) > fmu_kSHAPE_ANGLE_TOL )
            return fmu_kCLEAN_NO_ACTION;
    }

    //Quads that will be deleted
    CQuads[0] = pCQuad;
    CQuads[1] = pCQuadA;
    //Edges that will be deleted
    CEdges[0] = pCEdgePO;
    //No nodes to delete
    //Nodes of outer ring
    CRingNodes[0] = pCNode;
    CRingNodes[1] = pCNextNode;
    CRingNodes[2] = pCOppNode;
    CRingNodes[3] = pCNodeB;
    CRingNodes[4] = pCNodeA;
    CRingNodes[5] = pCPrevNode;

    if( get_valence( pCNodeA ) > 4 ||
        get_valence( pCNodeB ) == 3 )
        iCallStatus = fill_two( CNodes, CEdges, CQuads, CRingNodes );
    else
        iCallStatus = fill_three( CNodes, CEdges, CQuads, CRingNodes );

    if( iCallStatus == fmu_kCLEAN_MESH_CHANGED )
    {
        iStatus = valence_2_check( CRingNodes );
        if( iStatus == fmu_kCLEAN_ABORT )
            return fmu_kCLEAN_ABORT;
        if( iStatus == fmu_kCLEAN_NO_MEMORY )
            return fmu_kCLEAN_NO_MEMORY;
    }

    return iCallStatus;
}


int QuadCleanTool::shape_454( elmQClass *pCQuad, nodClass *pCNode,
                              nodClass *pCNextNode )
{

   int iIndex, iCallStatus;
   bool lFound;
   nodClass *pCSaveNode;
   nodClassDynArray CNodes(ARRAY_SIZE);
   fmuEdgeDynArray CEdges(ARRAY_SIZE);
   elmQClassDynArray CQuads(ARRAY_SIZE);
   nodClassDynArray *pCBoundary;

   //Need the neighbors. How to get them depends on where this quad is.
   if( pCNode->boundary() )
   {
      lFound = find_on_boundary( pCNode, &pCBoundary, &iIndex, &pCSaveNode );
      if( !lFound ) return fmu_kCLEAN_NO_ACTION;
      boundary_neighbors( pCBoundary, CNodes, CEdges, CQuads );
      //Reset the boundary loop position in case it was important
      pCBoundary->set( iIndex );
   }
   else
   {
      neighbors( pCNode, CNodes, CEdges, CQuads );
   }
   CNodes.find( pCNextNode );
   CEdges.find( pCNode->shared_edge( pCNextNode ) );
   CQuads.find( pCQuad );

   //That's right, use a boundary clean routine - that just happens to
   //do what needs to be done.
   iCallStatus = bndy_flat_5( pCNode, CNodes, CEdges, CQuads );
   return iCallStatus;
   
}


int QuadCleanTool::shape_444( elmQClass *pCQuad, nodClass *pCNode,
                                  nodClass *pCNextNode, nodClass *pCOppNode,
                                  nodClass *pCPrevNode )
{

   //Node designations:
   //          PK
   //          +
   //         / \           PK = Peak
   //        /   \           L = Left
   //     L +     +R         R = Right
   //      / \   / \        LA = LeftA
   //     /   \ /   \       RA = RightA
   //  LA+     +O    +RA     O = Opposite
   //     \   / \   /        P = Previous
   //      \ /   \ /         N = Next
   //       +--+--+          T = This (pCNode)
   //       P  T  N
   //
   //If peak node is 3 or 4 valent, pattern "A" or its mirror is used so that
   //there won't be two 3-valent nodes near each other. This assumes that
   //previous node or next node is not a hardpoint.
   //If peak node is 5 valent, pattern "B" or its mirror is used so that 
   //there won't be a 6 valent node. This requires that previous or next is
   //not a hardpoint and matching LA or RA is not 5 valent.
   int iCallStatus;
   double dAngle, dAngleTPLA, dAngleRANT;
   elmQClassDynArray CQuads;
   nodClassDynArray CNodes, CRingNodes;
   fmuEdgeDynArray CEdges;

   fmuEdge *pCEdgePO = pCPrevNode->shared_edge( pCOppNode );
   if( pCEdgePO->frozen() ) return fmu_kCLEAN_NO_ACTION;
   elmQClass *pCPrevQuad = pCEdgePO->other_quad( pCQuad );
   fmuEdge *pCEdgeNO = pCNextNode->shared_edge( pCOppNode );
   if( pCEdgeNO->frozen() ) return fmu_kCLEAN_NO_ACTION;
   elmQClass *pCNextQuad = pCEdgeNO->other_quad( pCQuad );

   nodClass *pCLeftNode = pCPrevQuad->next_node( pCOppNode );
   nodClass *pCLeftANode = pCPrevQuad->opposite_node( pCOppNode );
   nodClass *pCRightNode = pCNextQuad->prev_node( pCOppNode );
   nodClass *pCRightANode = pCNextQuad->opposite_node( pCOppNode );
   fmuEdge *pCEdgeRO = pCOppNode->shared_edge( pCRightNode );
   if( pCEdgeRO->frozen() ) return fmu_kCLEAN_NO_ACTION;
   fmuEdge *pCEdgeLO = pCOppNode->shared_edge( pCLeftNode );
   if( pCEdgeLO->frozen() ) return fmu_kCLEAN_NO_ACTION;
   elmQClass *pCQuadA = pCEdgeRO->other_quad( pCNextQuad );
   nodClass *pCPeakNode = pCQuadA->opposite_node( pCOppNode );

   int iLeftValence = get_valence( pCLeftNode );
   int iRightValence = get_valence( pCRightNode );
   int iPeakValence = get_valence( pCPeakNode );

   fmuVector CVecPT( pCPrevNode, pCNode );
   fmuVector CVecPLA( pCPrevNode, pCLeftANode );
   dAngleTPLA = CNormal.angle( CVecPT, CVecPLA );

   fmuVector CVecNRA( pCNextNode, pCRightANode );
   fmuVector CVecNT( pCNextNode, pCNode );
   dAngleRANT = CNormal.angle( CVecNRA, CVecNT );

   if( iPeakValence < 5 )
   {
      if( !pCPrevNode->hardpt() && !pCRightANode->hardpt() &&
          dAngleTPLA < fmu_kSHAPE_ANGLE_TOL )
      {
         return shape_444_prev( pCNode, pCNextNode, pCOppNode, pCPrevNode,
                                pCPeakNode, pCRightNode, pCRightANode,
                                pCLeftNode, pCLeftANode,
                                pCEdgePO, pCEdgeNO, pCEdgeLO, pCEdgeRO,
                                pCPrevQuad, pCQuad, pCNextQuad, pCQuadA );
      }

      if( !pCNextNode->hardpt() && !pCLeftANode->hardpt() &&
          dAngleRANT < fmu_kSHAPE_ANGLE_TOL )
      {
         return shape_444_next( pCNode, pCNextNode, pCOppNode, pCPrevNode,
                                pCPeakNode, pCRightNode, pCRightANode,
                                pCLeftNode, pCLeftANode, 
                                pCEdgePO, pCEdgeNO, pCEdgeLO, pCEdgeRO,
                                pCPrevQuad, pCQuad, pCNextQuad, pCQuadA );
      }

      if( iLeftValence != 3 && iRightValence != 3 )
      {
         return shape_444_lowpeak( pCNode, pCNextNode, pCOppNode, pCPrevNode,
                                   pCPeakNode, pCRightNode, pCRightANode,
                                   pCLeftNode, pCLeftANode, 
                                   pCEdgePO, pCEdgeNO, pCEdgeLO, pCEdgeRO,
                                   pCPrevQuad, pCQuad, pCNextQuad, pCQuadA );
      }
   }
   if( iPeakValence > 3 )
   {
      if( !pCPrevNode->hardpt() &&
          get_valence ( pCRightANode ) < 5 &&
          dAngleTPLA < fmu_kSHAPE_ANGLE_TOL )
      {
         //Quads that will be deleted
         CQuads[0] = pCQuad;
         CQuads[1] = pCPrevQuad;
         //Edges that will be deleted
         CEdges[0] = pCEdgePO;
         //No nodes to delete
         //Nodes of outer ring
         CRingNodes[0] = pCNode;
         CRingNodes[1] = pCNextNode;
         CRingNodes[2] = pCOppNode;
         CRingNodes[3] = pCLeftNode;
         CRingNodes[4] = pCLeftANode;
         CRingNodes[5] = pCPrevNode;

         iCallStatus = fill_two( CNodes, CEdges, CQuads, CRingNodes );
         return iCallStatus;
      }
      if( !pCNextNode->hardpt() &&
          get_valence ( pCLeftANode ) < 5 &&
          dAngleRANT < fmu_kSHAPE_ANGLE_TOL )
      {
         //Quads that will be deleted
         CQuads[0] = pCQuad;
         CQuads[1] = pCNextQuad;
         //Edges that will be deleted
         CEdges[0] = pCEdgeNO;
         //No nodes to delete
         //Nodes of outer ring
         CRingNodes[0] = pCNode;
         CRingNodes[1] = pCNextNode;
         CRingNodes[2] = pCRightANode;
         CRingNodes[3] = pCRightNode;
         CRingNodes[4] = pCOppNode;
         CRingNodes[5] = pCPrevNode;

         iCallStatus = fill_two( CNodes, CEdges, CQuads, CRingNodes );
         return iCallStatus;
      }
   }

   fmuVector CVecPO( pCPrevNode, pCOppNode );
   dAngle = CNormal.angle( CVecPT, CVecPO );
   if( iLeftValence == 3 && dAngle < fmu_kRIGHT_LOW_TOL )
   {
      return shape_444_left3( pCNode, pCNextNode, pCOppNode, pCPrevNode,
                                pCLeftNode, pCLeftANode,
                                pCEdgePO, pCPrevQuad, pCQuad );
   }

   fmuVector CVecNO( pCNextNode, pCOppNode );
   dAngle = CNormal.angle( CVecNO, CVecNT );
   if( iRightValence == 3 && dAngle < fmu_kRIGHT_LOW_TOL )
   {
      return shape_444_right3( pCNode, pCNextNode, pCOppNode, pCPrevNode,
                                pCRightNode, pCRightANode, pCEdgeNO,
                                        pCQuad, pCNextQuad );
   }
   return fmu_kCLEAN_NO_ACTION;
   
}


int QuadCleanTool::shape_444_prev( nodClass *pCNode, nodClass *pCNextNode,
                                       nodClass *pCOppNode,
                                       nodClass *pCPrevNode,
                                       nodClass *pCPeakNode,
                                       nodClass *pCRightNode,
                                       nodClass *pCRightANode,
                                       nodClass *pCLeftNode,
                                       nodClass *pCLeftANode,
                                       fmuEdge *pCEdgePO, fmuEdge *pCEdgeNO,
                                       fmuEdge *pCEdgeLO, fmuEdge *pCEdgeRO,
                                       elmQClass *pCPrevQuad,
                                       elmQClass *pCQuad,
                                       elmQClass *pCNextQuad,
                                       elmQClass *pCQuadA )
{
   
   fmuEdge *pCEdgeRRA = 0;
   fmuEdge *pCEdgeRPK = 0;
   elmQClass *pCQuadB = NULL;
   nodClass *pCDelNode = NULL;
   nodClassDynArray CSmoothList;

   int iRightValence = get_valence( pCRightNode );
   if( !pCRightNode->hardpt() && iRightValence == 3 )
   {
      pCEdgeRRA = pCRightNode->shared_edge( pCRightANode );
      pCEdgeRPK = pCRightNode->shared_edge( pCPeakNode );
      if( !pCEdgeRRA->frozen() && !pCEdgeRPK->frozen() )
      {
         pCQuadB = pCEdgeRRA->other_quad( pCNextQuad );
         pCDelNode = pCRightNode;
         pCRightNode = pCQuadB->opposite_node( pCRightNode );
      }
   }
   if( lDrawCleaning )
   {
   //draw_open();
   //pCQuad->draw_me( dp_kCLR_LTBLUE );
   //pCQuadA->draw_me( dp_kCLR_RED );
   //pCPrevQuad->draw_me( dp_kCLR_GOLD );
   //pCNextQuad->draw_me( dp_kCLR_BLUE );
   //if( pCQuadB ) pCQuadB->draw_me( dp_kCLR_CYAN );
                    
       // Refresh debug graphics
       //EraseElements();
       pCQuad->draw_me( dp_kCLR_YELLOW );
       pCQuadA->draw_me( dp_kCLR_YELLOW );
       pCPrevQuad->draw_me( dp_kCLR_YELLOW );
       pCNextQuad->draw_me( dp_kCLR_YELLOW );
                   
   //draw_close();
   }
   delete_this_quad( pCQuad );
   delete_this_quad( pCPrevQuad );
   delete_this_quad( pCNextQuad );
   delete_this_quad( pCQuadA );
   if( pCQuadB ) delete_this_quad( pCQuadB );
   if( delete_this_edge( pCEdgePO ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdgeNO ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdgeLO ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdgeRO ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( pCDelNode )
   {
      if( delete_this_edge( pCEdgeRRA ) == fmu_kCLEAN_ABORT )
         return fmu_kCLEAN_ABORT;
      if( delete_this_edge( pCEdgeRPK ) == fmu_kCLEAN_ABORT )
         return fmu_kCLEAN_ABORT;
      if( delete_this_node( pCDelNode ) == fmu_kCLEAN_ABORT )
         return fmu_kCLEAN_ABORT;
   }
   if( delete_this_node( pCOppNode ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;

   pCQuad     = new elmQClass( pCPrevNode, pCNode, pCLeftNode,
                                           pCLeftANode, pCMaster );
   if( !pCQuad ) return fmu_kCLEAN_NO_MEMORY;
   pCQuadA    = new elmQClass( pCNode, pCNextNode, pCPeakNode,
                                            pCLeftNode, pCMaster );
   if( !pCQuadA ) return fmu_kCLEAN_NO_MEMORY;
   pCNextQuad = new elmQClass( pCNextNode, pCRightANode, pCRightNode,
                                            pCPeakNode, pCMaster );
   if( !pCNextQuad ) return fmu_kCLEAN_NO_MEMORY;

   CSmoothList.append( pCLeftNode );
   CSmoothList.append( pCLeftANode );
   CSmoothList.append( pCRightNode );
   CSmoothList.append( pCRightANode );
   CSmoothList.append( pCPeakNode );
   smooth_mesh( &CSmoothList );
   if( lDrawCleaning )
   {
   //draw_mesh();
   //pCQuad->draw_me( dp_kCLR_LTBLUE );
   //pCQuadA->draw_me( dp_kCLR_RED );
   //pCNextQuad->draw_me( dp_kCLR_BLUE );

        // Refresh debug graphics
        //EraseElements();
        pCQuad->draw_me( dp_kCLR_BLUE );
        pCQuadA->draw_me( dp_kCLR_BLUE );
        pCNextQuad->draw_me( dp_kCLR_BLUE );
                   
   //draw_close();
   }

   return fmu_kCLEAN_MESH_CHANGED;
   
}


int QuadCleanTool::shape_444_next( nodClass *pCNode, nodClass *pCNextNode,
                                       nodClass *pCOppNode,
                                       nodClass *pCPrevNode,
                                       nodClass *pCPeakNode,
                                       nodClass *pCRightNode,
                                       nodClass *pCRightANode,
                                       nodClass *pCLeftNode,
                                       nodClass *pCLeftANode,
                                       fmuEdge *pCEdgePO, fmuEdge *pCEdgeNO,
                                       fmuEdge *pCEdgeLO, fmuEdge *pCEdgeRO,
                                       elmQClass *pCPrevQuad,
                                       elmQClass *pCQuad,
                                       elmQClass *pCNextQuad,
                                       elmQClass *pCQuadA )
{
   
   fmuEdge *pCEdgeLLA = 0;
   fmuEdge *pCEdgeLPK = 0;
   elmQClass *pCQuadB = NULL;
   nodClass *pCDelNode = NULL;
   nodClassDynArray CSmoothList;

   int iLeftValence = get_valence( pCLeftNode );
   if( !pCLeftNode->hardpt() && iLeftValence == 3 )
   {
      pCEdgeLLA = pCLeftNode->shared_edge( pCLeftANode ); 
      pCEdgeLPK = pCLeftNode->shared_edge( pCPeakNode );
      if( !pCEdgeLLA->frozen() && !pCEdgeLPK->frozen() )
      {
         pCQuadB = pCEdgeLLA->other_quad( pCPrevQuad );
         pCDelNode = pCLeftNode;
         pCLeftNode = pCQuadB->opposite_node( pCLeftNode );
      }
   }
   if( lDrawCleaning )
   {
   //draw_open();
   //pCQuad->draw_me( dp_kCLR_LTBLUE );
   //pCQuadA->draw_me( dp_kCLR_RED );
   //pCPrevQuad->draw_me( dp_kCLR_GOLD );
   //pCNextQuad->draw_me( dp_kCLR_BLUE );
   //if( pCQuadB ) pCQuadB->draw_me( dp_kCLR_CYAN );

                    // Refresh debug graphics
                    //EraseElements();
                    pCQuad->draw_me( dp_kCLR_YELLOW );
                    pCQuadA->draw_me( dp_kCLR_YELLOW );
                    pCPrevQuad->draw_me( dp_kCLR_YELLOW );
                    pCNextQuad->draw_me( dp_kCLR_YELLOW );
                    
   //draw_close();
   }
   delete_this_quad( pCQuad );
   delete_this_quad( pCPrevQuad );
   delete_this_quad( pCNextQuad );
   delete_this_quad( pCQuadA );
   if( pCQuadB ) delete_this_quad( pCQuadB );
   if( delete_this_edge( pCEdgePO ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdgeNO ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdgeLO ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdgeRO ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( pCDelNode )
   {
      if( delete_this_edge( pCEdgeLLA ) == fmu_kCLEAN_ABORT )
         return fmu_kCLEAN_ABORT;
      if( delete_this_edge( pCEdgeLPK ) == fmu_kCLEAN_ABORT )
         return fmu_kCLEAN_ABORT;
      if( delete_this_node( pCDelNode ) == fmu_kCLEAN_ABORT )
         return fmu_kCLEAN_ABORT;
   }
   if( delete_this_node( pCOppNode ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;

   pCQuad     = new elmQClass( pCLeftANode, pCPrevNode, pCPeakNode,
                                            pCLeftNode, pCMaster );
   if( !pCQuad ) return fmu_kCLEAN_NO_MEMORY;
   pCQuadA    = new elmQClass( pCPrevNode, pCNode, pCRightNode,
                                            pCPeakNode, pCMaster );
   if( !pCQuadA ) return fmu_kCLEAN_NO_MEMORY;
   pCNextQuad = new elmQClass( pCNode, pCNextNode, pCRightANode,
                                           pCRightNode, pCMaster );
   if( !pCNextQuad ) return fmu_kCLEAN_NO_MEMORY;

   CSmoothList.append( pCLeftNode );
   CSmoothList.append( pCLeftANode );
   CSmoothList.append( pCRightNode );
   CSmoothList.append( pCRightANode );
   CSmoothList.append( pCPeakNode );
   smooth_mesh( &CSmoothList );
   if( lDrawCleaning )
   {
   //draw_mesh();
   //pCQuad->draw_me( dp_kCLR_LTBLUE );
   //pCQuadA->draw_me( dp_kCLR_RED );
   //pCNextQuad->draw_me( dp_kCLR_BLUE );

                    // Refresh debug graphics
                    //EraseElements();
                    pCQuad->draw_me( dp_kCLR_BLUE );
                    pCQuadA->draw_me( dp_kCLR_BLUE );
                    pCNextQuad->draw_me( dp_kCLR_BLUE );
                   
   //draw_close();
   }

   return fmu_kCLEAN_MESH_CHANGED;
   
}


int QuadCleanTool::shape_444_lowpeak( nodClass *pCNode,
                                          nodClass *pCNextNode,
                                          nodClass *pCOppNode,
                                          nodClass *pCPrevNode,
                                          nodClass *pCPeakNode,
                                          nodClass *pCRightNode,
                                          nodClass *pCRightANode,
                                          nodClass *pCLeftNode,
                                          nodClass *pCLeftANode,
                                          fmuEdge *pCEdgePO, fmuEdge *pCEdgeNO,
                                          fmuEdge *pCEdgeLO, fmuEdge *pCEdgeRO,
                                          elmQClass *pCPrevQuad,
                                          elmQClass *pCQuad,
                                          elmQClass *pCNextQuad,
                                          elmQClass *pCQuadA )
{
   
   //A final case for low peak valence
   nodClassDynArray CSmoothList;

   if( lDrawCleaning )
   {
   //draw_mesh();
   //pCQuad->draw_me( dp_kCLR_LTBLUE );
   //pCQuadA->draw_me( dp_kCLR_RED );
   //pCNextQuad->draw_me( dp_kCLR_GOLD );
   //pCPrevQuad->draw_me( dp_kCLR_GRBLUE );

                    // Refresh debug graphics
                    //EraseElements();
                    pCQuad->draw_me( dp_kCLR_YELLOW );
                    pCQuadA->draw_me( dp_kCLR_YELLOW );
                    pCNextQuad->draw_me( dp_kCLR_YELLOW );
                    pCPrevQuad->draw_me( dp_kCLR_YELLOW );
                   
   //draw_close();
   }
   delete_this_quad( pCQuad );
   delete_this_quad( pCPrevQuad );
   delete_this_quad( pCNextQuad );
   delete_this_quad( pCQuadA );
   if( delete_this_edge( pCEdgePO ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdgeNO ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdgeLO ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   if( delete_this_edge( pCEdgeRO ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;

   //Need two nodes on either side of opp node
   nodClass *pCNodeL1 = one_middle_node( pCNode, pCLeftNode );
   nodClass *pCNodeR1 = one_middle_node( pCNode, pCRightNode );

   elmQClass *pCQuadR, *pCQuadL;
   pCQuad     = new elmQClass( pCNode, pCNextNode, pCNodeR1,
                                                    pCOppNode, pCMaster );
   if( !pCQuad ) return fmu_kCLEAN_NO_MEMORY;
   pCNextQuad= new elmQClass( pCNextNode, pCRightANode, pCRightNode,
                                                    pCNodeR1, pCMaster );
   if( !pCNextQuad ) return fmu_kCLEAN_NO_MEMORY;
   pCQuadA    = new elmQClass( pCPrevNode, pCNode, pCOppNode,
                                                    pCNodeL1, pCMaster );
   if( !pCQuadA ) return fmu_kCLEAN_NO_MEMORY;
   pCPrevQuad = new elmQClass( pCPrevNode, pCNodeL1, pCLeftNode,
                                                    pCLeftANode, pCMaster );
   if( !pCPrevQuad ) return fmu_kCLEAN_NO_MEMORY;
   pCQuadR    = new elmQClass( pCOppNode, pCNodeR1, pCRightNode,
                                                    pCPeakNode, pCMaster );
   if( !pCQuadR ) return fmu_kCLEAN_NO_MEMORY;
   pCQuadL    = new elmQClass( pCLeftNode, pCNodeL1, pCOppNode,
                                                    pCPeakNode, pCMaster );
   if( !pCQuadL ) return fmu_kCLEAN_NO_MEMORY;

   CSmoothList.append( pCLeftNode );
   CSmoothList.append( pCLeftANode );
   CSmoothList.append( pCRightNode );
   CSmoothList.append( pCRightANode );
   CSmoothList.append( pCPeakNode );
   CSmoothList.append( pCNodeR1 );
   CSmoothList.append( pCNodeL1 );
   smooth_mesh( &CSmoothList );
   if( lDrawCleaning )
   {
   //draw_mesh();
   //pCQuad->draw_me( dp_kCLR_LTBLUE );
   //pCQuadA->draw_me( dp_kCLR_RED );
   //pCNextQuad->draw_me( dp_kCLR_GOLD );
   //pCPrevQuad->draw_me( dp_kCLR_GRBLUE );
   //pCQuadR->draw_me( dp_kCLR_BLUE );
   //pCQuadL->draw_me( dp_kCLR_CYAN );

                    // Refresh debug graphics
                    //EraseElements();
                    pCQuad->draw_me( dp_kCLR_BLUE );
                    pCQuadA->draw_me( dp_kCLR_BLUE );
                    pCNextQuad->draw_me( dp_kCLR_BLUE );
                    pCPrevQuad->draw_me( dp_kCLR_BLUE );
                    pCQuadR->draw_me( dp_kCLR_BLUE );
                    pCQuadL->draw_me( dp_kCLR_BLUE );

   //draw_close();
   }

   return fmu_kCLEAN_MESH_CHANGED;
   
}



int QuadCleanTool::shape_444_left3( nodClass *pCNode, nodClass *pCNextNode,
                                        nodClass *pCOppNode,
                                        nodClass *pCPrevNode,
                                        nodClass *pCLeftNode,
                                        nodClass *pCLeftANode,
                                        fmuEdge *pCEdgePO,
                                        elmQClass *pCPrevQuad,
                                        elmQClass *pCQuad )
{
   
   nodClassDynArray CSmoothList;

   if( lDrawCleaning )
   {
   //draw_open();
   //pCQuad->draw_me( dp_kCLR_RED );
   //pCPrevQuad->draw_me( dp_kCLR_LTBLUE );

                    // Refresh debug graphics
                    //EraseElements();
                    pCQuad->draw_me( dp_kCLR_YELLOW );
                    pCPrevQuad->draw_me( dp_kCLR_YELLOW );
                   
   //draw_close();
   }

   //Quads to be deleted
   delete_this_quad( pCQuad );
   delete_this_quad( pCPrevQuad );
   if( delete_this_edge( pCEdgePO ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   nodClass *pCNode1, *pCNode2;
   two_middle_nodes(pCLeftANode, pCOppNode, &pCNode1, &pCNode2);
   pCQuad  = new elmQClass( pCNextNode, pCOppNode, pCNode2,
                                        pCNode, pCMaster );
   if( !pCQuad ) return fmu_kCLEAN_NO_MEMORY;
   elmQClass *pCQuad1  = new elmQClass( pCOppNode, pCLeftNode, pCNode1,
                                        pCNode2, pCMaster );
   if( !pCQuad1 ) return fmu_kCLEAN_NO_MEMORY;
   elmQClass *pCQuad2  = new elmQClass( pCNode1, pCLeftNode, pCLeftANode,
                                        pCPrevNode, pCMaster );
   if( !pCQuad2 ) return fmu_kCLEAN_NO_MEMORY;
   elmQClass *pCQuad3  = new elmQClass( pCNode1, pCPrevNode, pCNode,
                                        pCNode2, pCMaster );
   if( !pCQuad3 ) return fmu_kCLEAN_NO_MEMORY;
      
   CSmoothList.clear();
   CSmoothList.append( pCNode );
   CSmoothList.append( pCOppNode );
   CSmoothList.append( pCPrevNode );
   CSmoothList.append( pCNextNode );
   CSmoothList.append( pCNode1 );
   CSmoothList.append( pCNode2 );
   CSmoothList.append( pCLeftNode );
   smooth_mesh( &CSmoothList );
   if( lDrawCleaning )
   {
   //draw_mesh();
   //pCQuad->draw_me( dp_kCLR_LTBLUE );
   //pCQuad1->draw_me( dp_kCLR_RED );
   //pCQuad2->draw_me( dp_kCLR_GOLD );
   //pCQuad3->draw_me( dp_kCLR_BLUE );
                    // Refresh debug graphics
                    //EraseElements();
                    pCQuad->draw_me( dp_kCLR_BLUE );
                    pCQuad1->draw_me( dp_kCLR_BLUE );
                    pCQuad2->draw_me( dp_kCLR_BLUE );
                    pCQuad3->draw_me( dp_kCLR_BLUE );
                  
   //draw_close();
   }

   return fmu_kCLEAN_MESH_CHANGED;

  
}


int QuadCleanTool::shape_444_right3( nodClass *pCNode,
                                         nodClass *pCNextNode,
                                         nodClass *pCOppNode,
                                         nodClass *pCPrevNode,
                                         nodClass *pCRightNode,
                                         nodClass *pCRightANode,
                                         fmuEdge *pCEdgeNO,
                                         elmQClass *pCQuad,
                                         elmQClass *pCNextQuad )
{
   
   nodClassDynArray CSmoothList;

   if( lDrawCleaning )
   {
   //draw_open();
   //pCQuad->draw_me( dp_kCLR_RED );
   //pCNextQuad->draw_me( dp_kCLR_LTBLUE );

                    // Refresh debug graphics
                    //EraseElements();
                    pCQuad->draw_me( dp_kCLR_YELLOW );
                    pCNextQuad->draw_me( dp_kCLR_YELLOW );
                   
   //draw_close();
   }

   //Quads to be deleted
   delete_this_quad( pCQuad );
   delete_this_quad( pCNextQuad );
   if( delete_this_edge( pCEdgeNO ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   nodClass *pCNode1, *pCNode2;
   two_middle_nodes(pCRightANode, pCOppNode, &pCNode1, &pCNode2);
   pCQuad  = new elmQClass( pCPrevNode, pCNode, pCNode2,
                                        pCOppNode, pCMaster );
   if( !pCQuad ) return fmu_kCLEAN_NO_MEMORY;
   elmQClass *pCQuad1  = new elmQClass( pCOppNode, pCNode2, pCNode1,
                                        pCRightNode, pCMaster );
   if( !pCQuad1 ) return fmu_kCLEAN_NO_MEMORY;
   elmQClass *pCQuad2  = new elmQClass( pCNode1, pCNextNode, pCRightANode,
                                        pCRightNode, pCMaster );
   if( !pCQuad2 ) return fmu_kCLEAN_NO_MEMORY;
   elmQClass *pCQuad3  = new elmQClass( pCNode1, pCNode2, pCNode,
                                        pCNextNode, pCMaster );
   if( !pCQuad3 ) return fmu_kCLEAN_NO_MEMORY;

   CSmoothList.clear();
   CSmoothList.append( pCNode );
   CSmoothList.append( pCOppNode );
   CSmoothList.append( pCPrevNode );
   CSmoothList.append( pCNextNode );
   CSmoothList.append( pCNode1 );
   CSmoothList.append( pCNode2 );
   CSmoothList.append( pCRightNode );
   smooth_mesh( &CSmoothList );
   if( lDrawCleaning )
   {
   //draw_mesh();
   //pCQuad->draw_me( dp_kCLR_LTBLUE );
   //pCQuad1->draw_me( dp_kCLR_RED );
   //pCQuad2->draw_me( dp_kCLR_GOLD );
   //pCQuad3->draw_me( dp_kCLR_BLUE );
                    // Refresh debug graphics
                    //EraseElements();
                    pCQuad->draw_me( dp_kCLR_BLUE );
                    pCQuad1->draw_me( dp_kCLR_BLUE );
                    pCQuad2->draw_me( dp_kCLR_BLUE );
                    pCQuad3->draw_me( dp_kCLR_BLUE );
                    
   //draw_close();
   }

   return fmu_kCLEAN_MESH_CHANGED;
  
}


int QuadCleanTool::shape_434( elmQClass *pCQuad, nodClass *pCNode,
                                  nodClass *pCNextNode, nodClass *pCOppNode,
                                  nodClass *pCPrevNode )
{
   
   int iStatus, iCallStatus;
   elmQClassDynArray CQuads;
   nodClassDynArray CNodes, CRingNodes;
   fmuEdgeDynArray CEdges;

   fmuEdge *pCEdgePO = pCPrevNode->shared_edge( pCOppNode );
   if( pCEdgePO->frozen() ) return fmu_kCLEAN_NO_ACTION;
   elmQClass *pCPrevQuad = pCEdgePO->other_quad( pCQuad );
   fmuEdge *pCEdgeNO = pCNextNode->shared_edge( pCOppNode );
   if( pCEdgeNO->frozen() ) return fmu_kCLEAN_NO_ACTION;
   elmQClass *pCNextQuad = pCEdgeNO->other_quad( pCQuad );

   nodClass *pCLeftNode = pCPrevQuad->opposite_node( pCOppNode );
   nodClass *pCRightNode = pCNextQuad->opposite_node( pCOppNode );
   nodClass *pCCommonNode = pCNextQuad->prev_node( pCOppNode );

   fmuVector CVecNR( pCNextNode, pCRightNode );
   fmuVector CVecNT( pCNextNode, pCNode );
   if( CNormal.angle( CVecNR, CVecNT ) > fmu_kSHAPE_ANGLE_TOL )
      return fmu_kCLEAN_NO_ACTION;

   fmuVector CVecPT( pCPrevNode, pCNode );
   fmuVector CVecPL( pCPrevNode, pCLeftNode );
   if( CNormal.angle( CVecPT, CVecPL ) > fmu_kSHAPE_ANGLE_TOL )
      return fmu_kCLEAN_NO_ACTION;

   //Quads that will be deleted
   CQuads[0] = pCQuad;
   CQuads[1] = pCPrevQuad;
   CQuads[2] = pCNextQuad;
   //Edges that will be deleted
   CEdges[0] = pCEdgePO;
   CEdges[1] = pCEdgeNO;
   CEdges[2] = pCOppNode->shared_edge( pCCommonNode );
   //Nodes to delete
   CNodes[0] = pCOppNode;
   //Nodes of outer ring
   CRingNodes[0] = pCNode;
   CRingNodes[1] = pCNextNode;
   CRingNodes[2] = pCRightNode;
   CRingNodes[3] = pCCommonNode;
   CRingNodes[4] = pCLeftNode;
   CRingNodes[5] = pCPrevNode;

   iCallStatus = fill_two( CNodes, CEdges, CQuads, CRingNodes );
   if( iCallStatus == fmu_kCLEAN_MESH_CHANGED )
   {
      iStatus = valence_2_check( CRingNodes );
      if( iStatus == fmu_kCLEAN_ABORT ) return fmu_kCLEAN_ABORT;
      if( iStatus == fmu_kCLEAN_NO_MEMORY )
         return fmu_kCLEAN_NO_MEMORY;
   }
   return iCallStatus;
  
}


int QuadCleanTool::shape_smooth_prev( elmQClass *pCQuad,
                                           nodClass *pCNode,
                                           nodClass *pCPrevNode )
{
   
   int iDebugFlag = 0;
   int iStatus, iCallStatus;
   elmQClassDynArray CQuads;
   nodClassDynArray CNodes, CRingNodes;
   fmuEdgeDynArray CEdges;

   fmuEdge *pCEdgeTP = pCNode->shared_edge( pCPrevNode );
   if( pCEdgeTP->frozen() ) return fmu_kCLEAN_NO_ACTION;
   elmQClass *pCQuadA = pCEdgeTP->other_quad( pCQuad );
   if (pCQuadA == NULL)     return fmu_kCLEAN_NO_ACTION;
   nodClass *pCNodeA = pCQuadA->prev_node( pCNode );
   nodClass *pCNodeB = pCQuadA->opposite_node( pCNode );
   if( pCNodeB->hardpt() ) return fmu_kCLEAN_NO_ACTION;
   if( get_valence( pCNodeB ) < 5 ) return fmu_kCLEAN_NO_ACTION;
   fmuEdge *pCEdgePB = pCPrevNode->shared_edge( pCNodeB );
   if( pCEdgePB->frozen() ) return fmu_kCLEAN_NO_ACTION;
   elmQClass *pCQuadB = pCEdgePB->other_quad( pCQuadA );
   nodClass *pCNodeC = pCQuadB->opposite_node( pCPrevNode );
   nodClass *pCNodeD = pCQuadB->next_node( pCPrevNode );

   if( iDebugFlag )
   {
   //pCNodeA->draw_me( dp_kCLR_GOLD );
   //pCNodeB->draw_me( dp_kCLR_BLUE );
   //pCNodeC->draw_me( dp_kCLR_RED );
   //pCNodeD->draw_me( dp_kCLR_GRBLUE );
   }

   //Quads that will be deleted
   CQuads[0] = pCQuadA;
   CQuads[1] = pCQuadB;
   //Edges that will be deleted
   CEdges[0] = pCEdgePB;
   //No nodes to delete
   //Nodes of outer ring
   CRingNodes[0] = pCPrevNode;
   CRingNodes[1] = pCNodeD;
   CRingNodes[2] = pCNodeC;
   CRingNodes[3] = pCNodeB;
   CRingNodes[4] = pCNodeA;
   CRingNodes[5] = pCNode;

   iCallStatus = fill_three( CNodes, CEdges, CQuads, CRingNodes );
   if( iCallStatus == fmu_kCLEAN_MESH_CHANGED )
   {
      iStatus = valence_2_check( CRingNodes );
      if( iStatus == fmu_kCLEAN_ABORT ) return fmu_kCLEAN_ABORT;
      if( iStatus == fmu_kCLEAN_NO_MEMORY )
         return fmu_kCLEAN_NO_MEMORY;
   }
   return iCallStatus;
  
}


int QuadCleanTool::shape_smooth_next( elmQClass *pCQuad,
                                           nodClass *pCNode,
                                           nodClass *pCNextNode )
{
  
   int iDebugFlag = 0;
   int iStatus, iCallStatus;
   elmQClassDynArray CQuads;
   nodClassDynArray CNodes, CRingNodes;
   fmuEdgeDynArray CEdges;

   fmuEdge *pCEdgeTN = pCNode->shared_edge( pCNextNode );
   if( pCEdgeTN->frozen() ) return fmu_kCLEAN_NO_ACTION;
   elmQClass *pCQuadA = pCEdgeTN->other_quad( pCQuad );
   nodClass *pCNodeA = pCQuadA->next_node( pCNode );
   nodClass *pCNodeB = pCQuadA->opposite_node( pCNode );
   if( pCNodeB->hardpt() ) return fmu_kCLEAN_NO_ACTION;
   if( get_valence( pCNodeB ) < 5 ) return fmu_kCLEAN_NO_ACTION;
   fmuEdge *pCEdgeNB = pCNextNode->shared_edge( pCNodeB );
   if( pCEdgeNB->frozen() ) return fmu_kCLEAN_NO_ACTION;
   elmQClass *pCQuadB = pCEdgeNB->other_quad( pCQuadA );
   nodClass *pCNodeC = pCQuadB->opposite_node( pCNextNode );
   nodClass *pCNodeD = pCQuadB->prev_node( pCNextNode );

   if( iDebugFlag )
   {
   //pCNodeA->draw_me( dp_kCLR_GOLD );
   //pCNodeB->draw_me( dp_kCLR_BLUE );
   //pCNodeC->draw_me( dp_kCLR_RED );
   //pCNodeD->draw_me( dp_kCLR_GRBLUE );
   }

   //Quads that will be deleted
   CQuads[0] = pCQuadA;
   CQuads[1] = pCQuadB;
   //Edges that will be deleted
   CEdges[0] = pCEdgeNB;
   //No nodes to delete
   //Nodes of outer ring
   CRingNodes[0] = pCNextNode;
   CRingNodes[1] = pCNode;
   CRingNodes[2] = pCNodeA;
   CRingNodes[3] = pCNodeB;
   CRingNodes[4] = pCNodeC;
   CRingNodes[5] = pCNodeD;

   iCallStatus = fill_three( CNodes, CEdges, CQuads, CRingNodes );
   if( iCallStatus == fmu_kCLEAN_MESH_CHANGED )
   {
      iStatus = valence_2_check( CRingNodes );
      if( iStatus == fmu_kCLEAN_ABORT ) return fmu_kCLEAN_ABORT;
      if( iStatus == fmu_kCLEAN_NO_MEMORY )
         return fmu_kCLEAN_NO_MEMORY;
   }
   return iCallStatus;
}

int QuadCleanTool::shape_last_chance( elmQClass *pCQuad, nodClass *pCNode,
                                          nodClass *pCNextNode,
                                          nodClass *pCOppNode,
                                          nodClass *pCPrevNode )
{
   
   int iCallStatus, iStatus;
   elmQClass *pCQuadRight = 0, *pCQuadLeft = 0, *pCQuad1 = 0, *pCQuad2 = 0, *pCQuad3 = 0;
   nodClass *pCRANode = 0, *pCRBNode = 0, *pCNodeR = 0;
   nodClass *pCNodeL = 0, *pCLANode = 0, *pCLBNode = 0, *pCNodeC = 0;
   elmQClassDynArray CQuads;
   nodClassDynArray CNodes, CRingNodes;
   nodClassDynArray CSmoothList;
   fmuEdgeDynArray CEdges;

   //There are a few checks for cases not worth handling. Don't bother
   //if both opposite edges are perm (nothing can be done).
   fmuEdge *pCEdgePO = pCPrevNode->shared_edge( pCOppNode );
   fmuEdge *pCEdgeNO = pCNextNode->shared_edge( pCOppNode );
   if( pCEdgePO->frozen() && pCEdgeNO->frozen() ) return fmu_kCLEAN_NO_ACTION;
   //Don't bother if the edges on either side of the node with
   //bad angle are not frozen and if that node is a hardpoint and
   //if any of the other nodes are not hardpoints. The official smoothing
   //will probably resolve it.
   fmuEdge *pCEdgePT = pCPrevNode->shared_edge( pCNode );
   fmuEdge *pCEdgeNT = pCNextNode->shared_edge( pCNode );
   if( !pCEdgePT->frozen() && !pCEdgeNT->frozen() &&
       !pCNode->hardpt() &&
       ( !pCNextNode->hardpt() || !pCOppNode->hardpt() ||
         !pCNextNode->hardpt() ) )            return fmu_kCLEAN_NO_ACTION;

   const double fmu_kANGLE_DIFF = 0.52359878;  //30 degrees
   const double fmu_kLARGE_ANGLE_TOL = 2.3561945; //135 degrees

   fmuVector CVecNO( pCNextNode, pCOppNode );
   fmuVector CVecNT( pCNextNode, pCNode );
   fmuVector CVecPO( pCPrevNode, pCOppNode );
   fmuVector CVecPT( pCPrevNode, pCNode );
   double dAngleN = CNormal.angle( CVecNO, CVecNT );
   double dAngleP = CNormal.angle( CVecPT, CVecPO );
   //If the quad is squeezed at both ends, don't handle it here
   if( dAngleN < fmu_kSMALL_SHAPE_TOL && dAngleP  < fmu_kSMALL_SHAPE_TOL )
      return fmu_kCLEAN_NO_ACTION;
   double dAngleDiff = dAngleN - dAngleP;
   fmuVector CVecOT( pCOppNode, pCNode );
   double dAngleTORB = 0.0, dAngleLBOT = 0.0, dAngleB = 0.0;

   if( !pCEdgeNO->frozen() )
   {
      pCQuadRight = pCEdgeNO->other_quad( pCQuad ); 
	  if( !pCQuadRight ) return fmu_kCLEAN_NO_MEMORY;
	  pCRANode = pCQuadRight->next_node( pCNextNode );
	  //Check whether Two adjacent quads have 3 common nodes
	  if(pCRANode == pCNode) return fmu_kCLEAN_BAD_INPUT;
	  pCRBNode = pCQuadRight->opposite_node( pCNextNode );
	  fmuVector CVecORB( pCOppNode, pCRBNode );
	  dAngleTORB = CNormal.angle( CVecOT, CVecORB );
   }
   if( !pCEdgePO->frozen() )
   {
      pCQuadLeft = pCEdgePO->other_quad( pCQuad );
	  if( !pCQuadLeft ) return fmu_kCLEAN_NO_MEMORY;
	  
	  pCLANode = pCQuadLeft->prev_node( pCPrevNode );
	  //Check whether Two adjacent quads have 3 common nodes
	  if(pCLANode == pCNode) return fmu_kCLEAN_BAD_INPUT;
	  pCLBNode = pCQuadLeft->opposite_node( pCPrevNode );
	  fmuVector CVecOLB( pCOppNode, pCLBNode );
	  dAngleLBOT = CNormal.angle( CVecOLB, CVecOT );
   }

   //Case that uses all 3 quads. Check that:
   //both edges are not frozen, angles are similar,
   //no big angles are made on the far side of the pattern, and
   //the left and right quads do not share an edge.
   if( !pCEdgePO->frozen() && !pCEdgeNO->frozen() &&
       -fmu_kANGLE_DIFF < dAngleDiff && dAngleDiff < fmu_kANGLE_DIFF &&
       dAngleTORB < fmu_kLARGE_ANGLE_TOL &&
       dAngleLBOT < fmu_kLARGE_ANGLE_TOL &&
       pCRBNode != pCLBNode )
   {
      if( lDrawCleaning )
      {
      //draw_open();
      //pCQuadRight->draw_me( dp_kCLR_GRBLUE );
      //pCQuadLeft->draw_me( dp_kCLR_LTBLUE );
      //pCQuad->draw_me( dp_kCLR_DKGRN );

                    // Refresh debug graphics
                    //EraseElements();
                    pCQuadRight->draw_me( dp_kCLR_YELLOW );
                    pCQuadLeft->draw_me( dp_kCLR_YELLOW );
                    pCQuad->draw_me( dp_kCLR_YELLOW );
                   
      //draw_close();
      }

      delete_this_quad( pCQuad );
      delete_this_quad( pCQuadRight );
      delete_this_quad( pCQuadLeft );
      if( delete_this_edge( pCEdgeNO ) == fmu_kCLEAN_ABORT )
         return fmu_kCLEAN_ABORT;
      if( delete_this_edge( pCEdgePO ) == fmu_kCLEAN_ABORT )
         return fmu_kCLEAN_ABORT;

      pCNodeR = one_middle_node( pCNextNode, pCRBNode );
      two_middle_nodes( pCLBNode, pCNodeR, &pCNodeL, &pCNodeC );

      pCQuadLeft  = new elmQClass( pCLANode, pCPrevNode, pCNodeL, pCLBNode,
                                                              pCMaster );
      if( !pCQuadLeft ) return fmu_kCLEAN_NO_MEMORY;
      pCQuad  = new elmQClass( pCPrevNode, pCNode, pCNodeC, pCNodeL,
                                                              pCMaster );
      if( !pCQuad ) return fmu_kCLEAN_NO_MEMORY;
      pCQuad1  = new elmQClass( pCLBNode, pCNodeL, pCNodeC, pCOppNode,
                                                              pCMaster );
      if( !pCQuad1 ) return fmu_kCLEAN_NO_MEMORY;
      pCQuad2  = new elmQClass( pCNode, pCNextNode, pCNodeR, pCNodeC,
                                                              pCMaster );
      if( !pCQuad2 ) return fmu_kCLEAN_NO_MEMORY;
      pCQuad3  = new elmQClass( pCNodeC, pCNodeR, pCRBNode, pCOppNode,
                                                              pCMaster );
      if( !pCQuad3 ) return fmu_kCLEAN_NO_MEMORY;
      pCQuadRight  = new elmQClass( pCNodeR, pCNextNode, pCRANode, pCRBNode,
                                                              pCMaster );
      if( !pCQuadRight ) return fmu_kCLEAN_NO_MEMORY;

      CSmoothList.append( pCNodeR );
      CSmoothList.append( pCNodeL );
      CSmoothList.append( pCNodeC );
      CSmoothList.append( pCNextNode );
      CSmoothList.append( pCOppNode );
      CSmoothList.append( pCPrevNode );
      CSmoothList.append( pCNode );
      CSmoothList.append( pCLANode );
      CSmoothList.append( pCLBNode );
      CSmoothList.append( pCRANode );
      CSmoothList.append( pCRBNode );
      smooth_mesh( &CSmoothList );
      if( lDrawCleaning )
      {
      //draw_mesh();
      //pCQuad->draw_me( dp_kCLR_DKGRN );
      //pCQuadRight->draw_me( dp_kCLR_GRBLUE );
      //pCQuadLeft->draw_me( dp_kCLR_LTBLUE );
      //pCQuad1->draw_me( dp_kCLR_BLUE );
      //pCQuad2->draw_me( dp_kCLR_CYAN );
      //pCQuad3->draw_me( dp_kCLR_BLUE );

            // Refresh debug graphics
            //EraseElements();
            pCQuad->draw_me( dp_kCLR_BLUE );
            pCQuadRight->draw_me( dp_kCLR_BLUE );
            pCQuadLeft->draw_me( dp_kCLR_BLUE );
            pCQuad1->draw_me( dp_kCLR_BLUE );
            pCQuad2->draw_me( dp_kCLR_BLUE );
            pCQuad3->draw_me( dp_kCLR_BLUE );

      //draw_close();
      }

      return fmu_kCLEAN_MESH_CHANGED;
   }

   //Next side edge is not frozen and next angle is smaller (next edge
   //is larger).
   double dLenNO = CVecNO.length();
   double dLenNT = CVecNT.length();
   if( !pCEdgeNO->frozen() &&
       dAngleTORB < fmu_kHIGH_PI_TOL &&
      ( ( ( pCEdgePO->frozen() && dAngleDiff < fmu_kANGLE_DIFF ) ||
        ( !pCEdgePO->frozen() && dAngleDiff <= 0.0 ) ) ||
        dLenNT * 2.0 < dLenNO ) )
   {
      fmuVector CVecNRA( pCNextNode, pCRANode );
      dAngleN = CNormal.angle( CVecNRA, CVecNT );
      if( dAngleN > fmu_kSHAPE_ANGLE_TOL )
      {
         //Check lengths to avoid new elements with bad aspect ratios or
         //bad angles. If bad angles, keep the simpler but bad topology
         //that already exists.
         fmuVector CVecRARB( pCRANode, pCRBNode );
         fmuVector CVecRBO( pCRBNode, pCOppNode );
         if( CVecRBO.length() * 2.5 < CVecRARB.length() )
            return fmu_kCLEAN_NO_ACTION;

         if( lDrawCleaning )
         {
         //draw_open();
         //pCQuad->draw_me( dp_kCLR_DKGRN );
         //pCQuadRight->draw_me( dp_kCLR_GRBLUE );
            // Refresh debug graphics
            //EraseElements();
            pCQuad->draw_me( dp_kCLR_YELLOW );
            pCQuadRight->draw_me( dp_kCLR_YELLOW );
            
         //draw_close();
         }

         delete_this_quad( pCQuad );
         delete_this_quad( pCQuadRight );
         if( delete_this_edge( pCEdgeNO ) == fmu_kCLEAN_ABORT )
            return fmu_kCLEAN_ABORT;

         pCNodeL = one_middle_node( pCOppNode, pCNode );
         pCNodeR = one_middle_node( pCRBNode, pCNextNode );

         pCQuad  = new elmQClass( pCPrevNode, pCNode, pCNodeL,
                                           pCOppNode, pCMaster );
         if( !pCQuad ) return fmu_kCLEAN_NO_MEMORY;
         pCQuad1  = new elmQClass( pCOppNode, pCNodeL, pCNodeR,
                                           pCRBNode, pCMaster );
         if( !pCQuad1 ) return fmu_kCLEAN_NO_MEMORY;
         pCQuad2  = new elmQClass( pCNodeR, pCNextNode, pCRANode,
                                           pCRBNode, pCMaster );
         if( !pCQuad2 ) return fmu_kCLEAN_NO_MEMORY;
         pCQuad3  = new elmQClass( pCNode, pCNextNode, pCNodeR,
                                           pCNodeL, pCMaster );
         if( !pCQuad3 ) return fmu_kCLEAN_NO_MEMORY;

         CSmoothList.append( pCNodeR );
         CSmoothList.append( pCNodeL );
         CSmoothList.append( pCNode );
         CSmoothList.append( pCNextNode );
         CSmoothList.append( pCPrevNode );
         CSmoothList.append( pCRANode );
         CSmoothList.append( pCRBNode );
         CSmoothList.append( pCOppNode );
         smooth_mesh( &CSmoothList );
         if( lDrawCleaning )
         {
         //draw_mesh();
         //pCQuad->draw_me( dp_kCLR_GRBLUE );
         //pCQuad1->draw_me( dp_kCLR_DKGRN );
         //pCQuad2->draw_me( dp_kCLR_LTBLUE );
         //pCQuad3->draw_me( dp_kCLR_BLUE );
            // Refresh debug graphics
            //EraseElements();
            pCQuad->draw_me( dp_kCLR_BLUE );
            pCQuad1->draw_me( dp_kCLR_BLUE );
            pCQuad2->draw_me( dp_kCLR_BLUE );
            pCQuad3->draw_me( dp_kCLR_BLUE );

         //draw_close();
         }

         return fmu_kCLEAN_MESH_CHANGED;

      }
      else
      {

         //A tiny angle that any action will make worse.
         if( pCNextNode->hardpt() && pCRANode->hardpt() && pCRBNode->hardpt() )
         {
            fmuVector CVecRBRA( pCRBNode, pCRANode );
            dAngleB = CNormal.angle( CVecRBRA, CVecNRA );
            if( dAngleB < fmu_kSMALL_SHAPE_TOL ) return fmu_kCLEAN_NO_ACTION;
         }

         //Quads that will be deleted
         CQuads[0] = pCQuad;
         CQuads[1] = pCQuadRight;
         //Edges that will be deleted
         CEdges[0] = pCEdgeNO;
         //No nodes to delete
         //Nodes of outer ring
         CRingNodes[0] = pCNode;
         CRingNodes[1] = pCNextNode;
         CRingNodes[2] = pCRANode;
         CRingNodes[3] = pCRBNode;
         CRingNodes[4] = pCOppNode;
         CRingNodes[5] = pCPrevNode;

         if( !pCPrevNode->hardpt() && !pCRANode->hardpt() &&
             get_valence( pCPrevNode ) == 3 &&
             get_valence( pCRANode ) <= 4 )
         {
            //From this pattern of 3 valent nodes, the bad angle is not
            //split by an edge. Instead the "sides" are pulled in.
            CRingNodes.step(-1);
            iCallStatus = fill_two( CNodes, CEdges, CQuads, CRingNodes );
            if( iCallStatus == fmu_kCLEAN_MESH_CHANGED )
            {
               iStatus = valence_2_check( CRingNodes );
               if( iStatus == fmu_kCLEAN_ABORT ) return fmu_kCLEAN_ABORT;
               if( iStatus == fmu_kCLEAN_NO_MEMORY )
                  return fmu_kCLEAN_NO_MEMORY;
            }
         }
         else
         {
            iCallStatus = fill_three( CNodes, CEdges, CQuads, CRingNodes );
            if( iCallStatus == fmu_kCLEAN_MESH_CHANGED )
            {
               iStatus = valence_2_check( CRingNodes );
               if( iStatus == fmu_kCLEAN_ABORT ) return fmu_kCLEAN_ABORT;
               if( iStatus == fmu_kCLEAN_NO_MEMORY )
                  return fmu_kCLEAN_NO_MEMORY;
            }
         }
         return iCallStatus;
      }
   }

   //Prev side edge is not frozen and prev angle is smaller (prev edge
   //is larger).
   double dLenPO = CVecPO.length();
   double dLenPT = CVecPT.length();
   if( !pCEdgePO->frozen() &&
       dAngleLBOT < fmu_kHIGH_PI_TOL &&
      ( ( ( pCEdgeNO->frozen() && dAngleDiff > -fmu_kANGLE_DIFF ) ||
        ( !pCEdgeNO->frozen() && dAngleDiff > 0.0 ) ) ||
        dLenPT * 2.0 < dLenPO ) )
   {
      fmuVector CVecPLA( pCPrevNode, pCLANode );
      dAngleP = CNormal.angle( CVecPT, CVecPLA );
      if( dAngleP > fmu_kSHAPE_ANGLE_TOL )
      {
         //Check lengths to avoid new elements with bad aspect ratios or
         //bad angles. If bad angles, keep the simpler but bad topology
         //that already exists.
         fmuVector CVecLALB( pCLANode, pCLBNode );
         fmuVector CVecLBO( pCLBNode, pCOppNode );
         if( CVecLBO.length() * 2.5 < CVecLALB.length() )
            return fmu_kCLEAN_NO_ACTION;

         if( lDrawCleaning )
         {
         //draw_open();
         //pCQuad->draw_me( dp_kCLR_DKGRN );
         //pCQuadLeft->draw_me( dp_kCLR_GRBLUE );
            // Refresh debug graphics
            //EraseElements();
            pCQuad->draw_me( dp_kCLR_YELLOW );
            pCQuadLeft->draw_me( dp_kCLR_YELLOW );
            
         //draw_close();
         }

         delete_this_quad( pCQuad );
         delete_this_quad( pCQuadLeft );
         if( delete_this_edge( pCEdgePO ) == fmu_kCLEAN_ABORT )
            return fmu_kCLEAN_ABORT;

         pCNodeL = one_middle_node( pCLBNode, pCPrevNode );
         pCNodeR = one_middle_node( pCOppNode, pCNode );

         pCQuad  = new elmQClass( pCLANode, pCPrevNode, pCNodeL,
                                           pCLBNode, pCMaster );
         if( !pCQuad ) return fmu_kCLEAN_NO_MEMORY;
         pCQuad1  = new elmQClass( pCLBNode, pCNodeL, pCNodeR,
                                           pCOppNode, pCMaster );
         if( !pCQuad1 ) return fmu_kCLEAN_NO_MEMORY;
         pCQuad2  = new elmQClass( pCNodeR, pCNode, pCNextNode,
                                           pCOppNode, pCMaster );
         if( !pCQuad2 ) return fmu_kCLEAN_NO_MEMORY;
         pCQuad3  = new elmQClass( pCPrevNode, pCNode, pCNodeR,
                                           pCNodeL, pCMaster );
         if( !pCQuad3 ) return fmu_kCLEAN_NO_MEMORY;

         CSmoothList.append( pCNodeR );
         CSmoothList.append( pCNodeL );
         CSmoothList.append( pCNode );
         CSmoothList.append( pCNextNode );
         CSmoothList.append( pCPrevNode );
         CSmoothList.append( pCLANode );
         CSmoothList.append( pCLBNode );
         CSmoothList.append( pCOppNode );
         smooth_mesh( &CSmoothList );
         if( lDrawCleaning )
         {
         //draw_mesh();
         //pCQuad->draw_me( dp_kCLR_GRBLUE );
         //pCQuad1->draw_me( dp_kCLR_DKGRN );
         //pCQuad2->draw_me( dp_kCLR_LTBLUE );
         //pCQuad3->draw_me( dp_kCLR_BLUE );
            // Refresh debug graphics
            //EraseElements();
            pCQuad->draw_me( dp_kCLR_BLUE );
            pCQuad1->draw_me( dp_kCLR_BLUE );
            pCQuad2->draw_me( dp_kCLR_BLUE );
            pCQuad3->draw_me( dp_kCLR_BLUE );
         //draw_close();
         }

         return fmu_kCLEAN_MESH_CHANGED;

      }
      else
      {
         //A tiny angle that any action will make worse.
         if( pCPrevNode->hardpt() && pCLANode->hardpt() && pCLBNode->hardpt() )
         {
            fmuVector CVecLBLA( pCLBNode, pCLANode );
            dAngleB = CNormal.angle( CVecPLA, CVecLBLA );
            if( dAngleB < fmu_kSMALL_SHAPE_TOL ) return fmu_kCLEAN_NO_ACTION;
         }

         //Quads that will be deleted
         CQuads[0] = pCQuad;
         CQuads[1] = pCQuadLeft;
         //Edges that will be deleted
         CEdges[0] = pCEdgePO;
         //No nodes to delete
         //Nodes of outer ring
         CRingNodes[0] = pCNode;
         CRingNodes[1] = pCNextNode;
         CRingNodes[2] = pCOppNode;
         CRingNodes[3] = pCLBNode;
         CRingNodes[4] = pCLANode;
         CRingNodes[5] = pCPrevNode;

         if( !pCNextNode->hardpt() && !pCLANode->hardpt() &&
             get_valence( pCNextNode ) == 3 &&
             get_valence( pCLANode ) <= 4 )
         {
            //From this pattern of 3 valent nodes, the bad angle is not
            //split by an edge. Instead the "sides" are pulled in.
            CRingNodes.step();
            iCallStatus = fill_two( CNodes, CEdges, CQuads, CRingNodes );
            if( iCallStatus == fmu_kCLEAN_MESH_CHANGED )
            {
               iStatus = valence_2_check( CRingNodes );
               if( iStatus == fmu_kCLEAN_ABORT ) return fmu_kCLEAN_ABORT;
               if( iStatus == fmu_kCLEAN_NO_MEMORY )
                  return fmu_kCLEAN_NO_MEMORY;
            }
         }
         else
         {
            iCallStatus = fill_three( CNodes, CEdges, CQuads, CRingNodes );
            if( iCallStatus == fmu_kCLEAN_MESH_CHANGED )
            {
               iStatus = valence_2_check( CRingNodes );
               if( iStatus == fmu_kCLEAN_ABORT ) return fmu_kCLEAN_ABORT;
               if( iStatus == fmu_kCLEAN_NO_MEMORY )
                  return fmu_kCLEAN_NO_MEMORY;
            }
         }
         return iCallStatus;
      }
   }

   return fmu_kCLEAN_NO_ACTION;
   }


int QuadCleanTool::fill_choice( nodClassDynArray &CDelNodes,
                                    fmuEdgeDynArray &CDelEdges,
                                    elmQClassDynArray &CDelQuads,
                                    nodClassDynArray &CRingNodes )
{
  
   int iDebugFlag = 0;
   const int kMUST_HAVE_EDGE = 1;
   const int kNEUTRAL_EDGE = 2;
   const int kCANT_HAVE_EDGE = 3;
   const double fmu_kRIGHT_ANGLE = 1.5707963;
   int ii, iCallStatus;
   int iValence, nDelEdges;
   int aiClassify[6];
   int aiHardpt[6];
   double dAngle;

   nodClass *pCNode = 0, *pCPrev = 0, *pCNext = 0, *pCStartNode  = 0;

   if( CRingNodes.length() != 6 ) return fmu_kCLEAN_NO_ACTION;

   //Set the ring nodes internal iterator to a known position.
   CRingNodes.set( 0 );

   //No need to put in exactly what was taken out. Remember some details
   nDelEdges = CDelEdges.length();
   if( nDelEdges > 0 ) pCStartNode = CDelEdges[0]->start_node();

   //From the ring nodes, analyze which nodes get edges attached to them.
   //This analysis takes into account that the old stuff isn't deleted yet.
   for( ii = 0; ii < 6; ++ii )
   {
      pCNode = CRingNodes.get();
      aiHardpt[ii] = pCNode->hardpt();
      iValence = get_valence( pCNode );
      if( iValence > 5 ) aiClassify[ii] = kCANT_HAVE_EDGE;
      else if( iValence == 3) aiClassify[ii] = kMUST_HAVE_EDGE;
      else
      {
         pCPrev = CRingNodes.prev();
         pCNext = CRingNodes.next();
         fmuVector CToPrev( pCNode, pCPrev );
         fmuVector CToNext( pCNode, pCNext );
         dAngle = CNormal.angle( CToNext, CToPrev );
         if( pCNode->hardpt() )
         {
            //Angle of 145 degrees
            if( dAngle > 2.53 ) aiClassify[ii] = kMUST_HAVE_EDGE;
            //Angle of 110 degrees
            else if( dAngle < 1.9198622 ) aiClassify[ii] = kCANT_HAVE_EDGE;
            else aiClassify[ii] = kNEUTRAL_EDGE;
         }
         else
         {
            //Angle of 190 degrees - different as smoothing can adjust it.
            if( dAngle > 3.31612556 ) aiClassify[ii] = kMUST_HAVE_EDGE;
            //Angle of 110 degrees
            else if( dAngle < 1.9198622 ) aiClassify[ii] = kCANT_HAVE_EDGE;
            else aiClassify[ii] = kNEUTRAL_EDGE;
         }
      }
      if( iDebugFlag )
      {
      //if( aiClassify[ii] == kMUST_HAVE_EDGE )
      //pCNode->draw_me( dp_kCLR_BLUE );
      //else if( aiClassify[ii] == kCANT_HAVE_EDGE )
      //pCNode->draw_me( dp_kCLR_RED );
      //else
      //pCNode->draw_me( dp_kCLR_YELLOW );
      }
      CRingNodes.step();
   }

   CRingNodes.set( 0 );

   //Determine if two quads can fill the ring
   for( ii = 0; ii < 6; ++ii )
   {
      if( aiClassify[0] == kMUST_HAVE_EDGE &&
          aiClassify[1] != kMUST_HAVE_EDGE &&
          aiClassify[2] != kMUST_HAVE_EDGE &&
          aiClassify[3] != kCANT_HAVE_EDGE &&
          aiClassify[4] != kMUST_HAVE_EDGE &&
          aiClassify[5] != kMUST_HAVE_EDGE &&
          ( nDelEdges != 1 ||
            ( pCStartNode != CRingNodes.get() &&
              pCStartNode != CRingNodes.next(3) ) ) )
      {
         iCallStatus = fill_two( CDelNodes, CDelEdges, CDelQuads,
                                                       CRingNodes );
         return iCallStatus;
      }
      if( aiClassify[0] != kCANT_HAVE_EDGE &&
          aiClassify[1] == kCANT_HAVE_EDGE &&
          aiClassify[2] == kCANT_HAVE_EDGE &&
          aiClassify[3] != kCANT_HAVE_EDGE &&
          aiClassify[4] != kMUST_HAVE_EDGE &&
          aiClassify[5] != kMUST_HAVE_EDGE &&
          ( nDelEdges != 1 ||
            ( pCStartNode != CRingNodes.get() &&
              pCStartNode != CRingNodes.next(3) ) ) )
      {
         iCallStatus = fill_two( CDelNodes, CDelEdges, CDelQuads,
                                                       CRingNodes );
         return iCallStatus;
      }
      if( aiClassify[0] != kCANT_HAVE_EDGE &&
          aiClassify[1] != kMUST_HAVE_EDGE &&
          aiClassify[2] != kMUST_HAVE_EDGE &&
          aiClassify[3] != kCANT_HAVE_EDGE &&
          aiClassify[4] != kMUST_HAVE_EDGE &&
          aiClassify[5] != kMUST_HAVE_EDGE &&
          aiHardpt[1] && aiHardpt[2] &&
          ( nDelEdges != 1 ||
            ( pCStartNode != CRingNodes.get() &&
              pCStartNode != CRingNodes.next(3) ) ) )
      {
         iCallStatus = fill_two( CDelNodes, CDelEdges, CDelQuads,
                                                       CRingNodes );
         return iCallStatus;
      }
      CRingNodes.step();
      rotate_list( aiClassify, 6, 1 );
      rotate_list( aiHardpt, 6, 1 );
   }

   double adAngles[6], dAngleDiff, dAngleDiffa, dAngleDiffb;
   for( ii = 0; ii < 6; ++ii )
   {
      fmuVector CVecTN( CRingNodes.get(), CRingNodes.next() );
      fmuVector CVecTP( CRingNodes.get(), CRingNodes.prev() );
      adAngles[ii] = CNormal.angle( CVecTN, CVecTP );
      adAngles[ii] = fabs( fmu_kRIGHT_ANGLE - adAngles[ii] );
      CRingNodes.step();
   }
   for( ii = 0; ii < 3; ++ii )
   {
      dAngleDiff  = adAngles[1] + adAngles[2] + adAngles[4] + adAngles[5];
      dAngleDiffa = adAngles[0] + adAngles[1] + adAngles[3] + adAngles[4];
      dAngleDiffb = adAngles[0] + adAngles[2] + adAngles[3] + adAngles[5];
      if( aiClassify[0] <= kNEUTRAL_EDGE &&   //Same as neutral or must
          aiClassify[1] >= kNEUTRAL_EDGE &&   //Same as neutral or cant
          aiClassify[2] >= kNEUTRAL_EDGE &&
          aiClassify[3] <= kNEUTRAL_EDGE &&
          aiClassify[4] >= kNEUTRAL_EDGE &&
          aiClassify[5] >= kNEUTRAL_EDGE &&
          dAngleDiff < dAngleDiffa && dAngleDiff < dAngleDiffb &&
          ( nDelEdges != 1 ||
            ( pCStartNode != CRingNodes.get() &&
              pCStartNode != CRingNodes.next(3) ) ) )
      {
         iCallStatus = fill_two( CDelNodes, CDelEdges, CDelQuads,
                                                       CRingNodes );
         return iCallStatus;
      }
      CRingNodes.step();
      rotate_list( aiClassify, 6, 1 );
      rotate_list( aiHardpt, 6, 1 );
      rotate_doubles( adAngles, 6, 1 );
   }

   //Determine if three quads can fill the ring, two choices available
   for( ii = 0; ii < 2; ++ii )
   {
      if( aiClassify[0] <= kNEUTRAL_EDGE &&
          aiClassify[1] >= kNEUTRAL_EDGE &&
          aiClassify[2] <= kNEUTRAL_EDGE &&
          aiClassify[3] >= kNEUTRAL_EDGE &&
          aiClassify[4] <= kNEUTRAL_EDGE &&
          aiClassify[5] >= kNEUTRAL_EDGE   )
      {
         iCallStatus = fill_three( CDelNodes, CDelEdges, CDelQuads,
                                                         CRingNodes );
         return iCallStatus;
      }
      CRingNodes.step();
      rotate_list( aiClassify, 6, 1 );
   }

   //Since filling the space with 4 will probably only return us to where
   //we started, try an alternate maneuver of including neighboring quads
   //in the process.
   for( ii = 0; ii < 6; ++ii )
   {
      if( aiHardpt[0] && aiClassify[0] == kCANT_HAVE_EDGE &&
          aiHardpt[1] && aiClassify[1] <= kNEUTRAL_EDGE   &&
          aiHardpt[2] && aiClassify[2] <= kNEUTRAL_EDGE   &&
          aiHardpt[3] && aiClassify[3] == kCANT_HAVE_EDGE &&
          !aiHardpt[4] &&
          !aiHardpt[5]  )
      {
         if( aiClassify[4] <= kNEUTRAL_EDGE &&
             aiClassify[5] == kCANT_HAVE_EDGE )
         {
            iCallStatus = fill_choice_side( CDelNodes, CDelEdges, CDelQuads,
                                                                CRingNodes );
            return iCallStatus;
         }

         if( aiClassify[4] == kCANT_HAVE_EDGE &&
             aiClassify[5] <= kNEUTRAL_EDGE )
         {
            iCallStatus = fill_choice_sidem( CDelNodes, CDelEdges, CDelQuads,
                                                                CRingNodes );
            return iCallStatus;
         }

         if( aiClassify[4] <= kNEUTRAL_EDGE &&
             aiClassify[5] <= kNEUTRAL_EDGE )
         {
            //Base the decision on which angle is smaller
            fmuVector CVec43( CRingNodes.next(4), CRingNodes.next(3) );
            CVec43.unitize();
            fmuVector CVec45( CRingNodes.next(4), CRingNodes.next(5) );
            CVec45.unitize();
            fmuVector CVec50( CRingNodes.next(5), CRingNodes.get() );
            CVec50.unitize();
            double dCos345 = CVec43 % CVec45;
            double dCos450 = ( CVec45 * -1.0 ) % CVec50;
            if( dCos345 < dCos450 )
            {
               iCallStatus = fill_choice_side( CDelNodes, CDelEdges,
                                               CDelQuads, CRingNodes );
               return iCallStatus;
            }
            else
            {
               iCallStatus = fill_choice_sidem( CDelNodes, CDelEdges,
                                                CDelQuads, CRingNodes );
               return iCallStatus;
            }
         }


      }

      //If necessary, other patterns can be added here

      CRingNodes.step();
      rotate_list( aiClassify, 6, 1 );
      rotate_list( aiHardpt, 6, 1 );
   }

   for( ii = 0; ii < 6; ++ii )
   {
      if( aiClassify[0] == kCANT_HAVE_EDGE &&
          aiClassify[1] <= kNEUTRAL_EDGE   &&
          aiClassify[2] <= kNEUTRAL_EDGE   &&
          aiClassify[3] == kCANT_HAVE_EDGE &&
          aiClassify[4] == kMUST_HAVE_EDGE &&
          aiClassify[5] == kMUST_HAVE_EDGE )
      {
         //Squeeze out the mesh (when a case is encountered)
         //might also encounter a case where both nodes are "cant"
         iCallStatus = fill_choice_squeeze( CDelNodes, CDelEdges,
                                            CDelQuads, CRingNodes );
         return iCallStatus;
      }
      CRingNodes.step();
      rotate_list( aiClassify, 6, 1 );
      rotate_list( aiHardpt, 6, 1 );
   }

   return fmu_kCLEAN_NO_ACTION;
   
}


int QuadCleanTool::fill_two( nodClassDynArray &CDelNodes,
                                 fmuEdgeDynArray &CDelEdges,
                                 elmQClassDynArray &CDelQuads,
                                 nodClassDynArray &CRingNodes )
{
   
   int ii, nCount;
   int iDebugFlag = 0;

   if( CRingNodes.length() != 6 ) return fmu_kCLEAN_NO_ACTION;

   if( iDebugFlag )
   {
   //draw_mesh();
   //draw_close();
   //CRingNodes[0]->draw_me( dp_kCLR_YELLOW );
   //CRingNodes[1]->draw_me( dp_kCLR_CYAN );
   //CRingNodes[2]->draw_me( dp_kCLR_RED );
   //CRingNodes[3]->draw_me( dp_kCLR_DKGRN );
   //CRingNodes[4]->draw_me( dp_kCLR_LTBLUE );
   //CRingNodes[5]->draw_me( dp_kCLR_MAGNTA );
   }

   nCount = CDelQuads.length();
   if( lDrawCleaning )
   {
   //draw_open();
    // Refresh debug graphics
    //EraseElements();
          
   for( ii = 0; ii < nCount; ++ii )
   {
        CDelQuads[ii]->draw_me( dp_kCLR_YELLOW );
   }
   //draw_close();
   }
   //Delete everything in the deletion lists
   for( ii = 0; ii < nCount; ++ii )
   {
      delete_this_quad( CDelQuads[ii] );
   }
   nCount = CDelEdges.length();
   for( ii = 0; ii < nCount; ++ii )
   {
      if( delete_this_edge( CDelEdges[ii] ) == fmu_kCLEAN_ABORT )
         return fmu_kCLEAN_ABORT;
   }
   nCount = CDelNodes.length();
   for( ii = 0; ii < nCount; ++ii )
   {
      if( delete_this_node( CDelNodes[ii] ) == fmu_kCLEAN_ABORT )
         return fmu_kCLEAN_ABORT;
   }

   //Make two new CQuads out of the ring of nodes
   elmQClass *pCQuad1 = new elmQClass( CRingNodes.get(), CRingNodes.next(),
                                       CRingNodes.next(2), CRingNodes.next(3),
                                       pCMaster );
   if( !pCQuad1 ) return fmu_kCLEAN_NO_MEMORY;
   elmQClass *pCQuad2 = new elmQClass( CRingNodes.next(3), CRingNodes.next(4),
                                       CRingNodes.next(5), CRingNodes.get(),
                                       pCMaster );
   if( !pCQuad2 ) return fmu_kCLEAN_NO_MEMORY;

   smooth_mesh( &CRingNodes );
   if( lDrawCleaning )
   {
   //draw_mesh();
   //pCQuad1->draw_me( dp_kCLR_YELLOW );
   //pCQuad2->draw_me( dp_kCLR_GOLD );
         // Refresh debug graphics
        //EraseElements();
        pCQuad1->draw_me( dp_kCLR_BLUE );
        pCQuad2->draw_me( dp_kCLR_BLUE );
   //draw_close();
   }


   for( ii = 0; ii < CRingNodes.length(); ++ii ) reclean( CRingNodes[ii] );

   return fmu_kCLEAN_MESH_CHANGED;
  
}


int QuadCleanTool::fill_three( nodClassDynArray &CDelNodes,
                               fmuEdgeDynArray  &CDelEdges,
                               elmQClassDynArray &CDelQuads,
                               nodClassDynArray &CRingNodes )
{
    int nCount;
    int iDebugFlag = 0;
    double dQ = 0.0; 
    bool bLockNode = false;

    if( CRingNodes.length() != 6 )
        return fmu_kCLEAN_NO_ACTION;

    if( iDebugFlag )
    {
        //draw_mesh();
        //draw_close();
        //CRingNodes[0]->draw_me( dp_kCLR_YELLOW );
        //CRingNodes[1]->draw_me( dp_kCLR_CYAN );
        //CRingNodes[2]->draw_me( dp_kCLR_RED );
        //CRingNodes[3]->draw_me( dp_kCLR_DKGRN );
        //CRingNodes[4]->draw_me( dp_kCLR_LTBLUE );
        //CRingNodes[5]->draw_me( dp_kCLR_MAGNTA );
    }

    nCount = CDelQuads.length();
    if( lDrawCleaning )
    {
        //draw_open();
        // Refresh
        //EraseElements();
        for( int ii = 0; ii < nCount; ++ii )
        {
        //   CDelQuads[ii]->draw_me( dp_kCLR_RED + ii );
            CDelQuads[ii]->draw_me( dp_kCLR_YELLOW );
        }
        //draw_close();
    }

    //Delete everything in the deletion lists
    for( int ii = 0; ii < nCount; ++ii )
    {
        delete_this_quad( CDelQuads[ii] );
        CDelQuads[ii] = NULL;
    }

    nCount = CDelEdges.length();
    for( int ii = 0; ii < nCount; ++ii )
    {
        int ret = delete_this_edge( CDelEdges[ii] );
        CDelEdges[ii] = NULL;
        if( ret == fmu_kCLEAN_ABORT )
            return fmu_kCLEAN_ABORT;
    }

    nCount = CDelNodes.length();
    for( int ii = 0; ii < nCount; ++ii )
    {
        int ret = delete_this_node( CDelNodes[ii] );
        CDelNodes[ii] = NULL;
        if( ret == fmu_kCLEAN_ABORT )
            return fmu_kCLEAN_ABORT;
    }

    //Make a node in the middle of the ring
    nodClass *pCNewNode = one_middle_node( CRingNodes.get(), CRingNodes.next(3) );

    // Get new node coordinates
    double dNewX=0.0, dNewY=0.0, dNewX1=0.0, dNewY1=0.0;

    dNewX = pCNewNode->node_x();
    dNewY = pCNewNode->node_y();

    //Make three new CQuads out of the ring of nodes
    elmQClass *pCQuad1 = new elmQClass( CRingNodes.get(), CRingNodes.next(),
                                        CRingNodes.next(2), pCNewNode, pCMaster );
    if( !pCQuad1 )
        return fmu_kCLEAN_NO_MEMORY;

#if NM_TEST_QUAD_QUALITY

    dQ = test_quad_quality(pCQuad1);
#endif

    elmQClass *pCQuad2 = new elmQClass( CRingNodes.next(2), CRingNodes.next(3),
                                        CRingNodes.next(4), pCNewNode, pCMaster );
    if( !pCQuad2 )
        return fmu_kCLEAN_NO_MEMORY;

#if NM_TEST_QUAD_QUALITY

    dQ = test_quad_quality( pCQuad2);

#endif

    elmQClass *pCQuad3 = new elmQClass( CRingNodes.next(4), CRingNodes.next(5),
                                        CRingNodes.next(6), pCNewNode, pCMaster );
    if( !pCQuad3 )
        return fmu_kCLEAN_NO_MEMORY;

#if NM_TEST_QUAD_QUALITY

    dQ = test_quad_quality( pCQuad3);

#endif

    smooth_mesh( &CRingNodes );

    dNewX1 = pCNewNode->node_x();
    dNewY1 = pCNewNode->node_y();

    dQ = test_quad_quality( pCQuad1 );
    if (dQ < 0.1)
        bLockNode = true;

    dQ = test_quad_quality( pCQuad2 );
    if (dQ < 0.1)
        bLockNode = true;

     dQ = test_quad_quality( pCQuad3);

     if (dQ < 0.1)
        bLockNode = true;

    // Undo node smoothing lock node from moving...
    if (bLockNode)
    {
        pCNewNode->set(dNewX, dNewY, 0.0);
        pCNewNode->lock_node();
    }

    if( lDrawCleaning )
    {
        //draw_mesh();
        //pCQuad1->draw_me( dp_kCLR_RED );
        //pCQuad2->draw_me( dp_kCLR_MAGNTA );
        //pCQuad3->draw_me( dp_kCLR_PINK );
         // Refresh debug graphics
        //EraseElements();
        pCQuad1->draw_me( dp_kCLR_BLUE );
        pCQuad2->draw_me( dp_kCLR_BLUE );
        pCQuad3->draw_me( dp_kCLR_BLUE );
        //draw_close();
    }

    for( int ii = 0; ii < CRingNodes.length(); ++ii )
        reclean( CRingNodes[ii] );

    return fmu_kCLEAN_MESH_CHANGED;
}


int QuadCleanTool::fill_choice_side( nodClassDynArray &CDelNodes,
                                         fmuEdgeDynArray &CDelEdges,
                                         elmQClassDynArray &CDelQuads,
                                         nodClassDynArray &CRingNodes )
{
  
   int ii, nCount;

   if( CRingNodes.length() != 6 ) return fmu_kCLEAN_NO_ACTION;

   //The ringnodes list has an assumed orientation
   nodClass *pCNode0 = CRingNodes.get();
   nodClass *pCNode1 = CRingNodes.next();
   nodClass *pCNode2 = CRingNodes.next(2);
   nodClass *pCNode3 = CRingNodes.next(3);
   nodClass *pCNode4 = CRingNodes.next(4);
   nodClass *pCNode5 = CRingNodes.next(5);
   if( get_valence( pCNode3 ) <=4 && pCNode3->hardpt() )
      return fmu_kCLEAN_NO_ACTION;
   fmuEdge *pCEdge34 = pCNode3->shared_edge( pCNode4 );
   if( pCEdge34->frozen() ) return fmu_kCLEAN_NO_ACTION;
   fmuEdge *pCEdge45 = pCNode4->shared_edge( pCNode5 );
   if( pCEdge45->frozen() ) return fmu_kCLEAN_NO_ACTION;

   //Orientation of CDelQuads does not matche orientation of CRingNodes
   elmQClass *pCQuadA, *pCQuadX, *pCQuadY;
   pCEdge34->quads( &pCQuadX, &pCQuadY );
   if( CDelQuads.find( pCQuadX ) > -1 ) pCQuadA = pCQuadY;
   else                                 pCQuadA = pCQuadX;

   nodClass *pCNodeA = pCQuadA->next_node( pCNode3 );
   nodClass *pCNodeB = pCQuadA->opposite_node( pCNode3 );

   if( pCNode2->hardpt() && pCNode3->hardpt() && pCNodeA->hardpt() )
   {
      fmuVector CVec3A( pCNode3, pCNodeA );
      fmuVector CVec32( pCNode3, pCNode2 );
      if( CNormal.angle( CVec3A, CVec32 ) > fmu_kSHAPE_ANGLE_TOL )
         return fmu_kCLEAN_NO_ACTION;
   }

   nCount = CDelQuads.length();
   if( lDrawCleaning )
   {
   //draw_open();
   // Refresh debug graphics
   //EraseElements();
   for( ii = 0; ii < nCount; ++ii )
   {
   //CDelQuads[ii]->draw_me( dp_kCLR_YELLOW + ii );
       CDelQuads[ii]->draw_me( dp_kCLR_YELLOW );
   }
   //pCQuadA->draw_me( dp_kCLR_BLUE );
      pCQuadA->draw_me( dp_kCLR_YELLOW );
       
   //draw_close();
   }

   //Delete everything in the deletion lists, plus an extra edge & quad
   for( ii = 0; ii < nCount; ++ii )
   {
      delete_this_quad( CDelQuads[ii] );
   }
   delete_this_quad( pCQuadA );
   nCount = CDelEdges.length();
   for( ii = 0; ii < nCount; ++ii )
   {
      if( delete_this_edge( CDelEdges[ii] ) == fmu_kCLEAN_ABORT )
         return fmu_kCLEAN_ABORT;
   }
   if( delete_this_edge( pCEdge34 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   nCount = CDelNodes.length();
   for( ii = 0; ii < nCount; ++ii )
   {
      if( delete_this_node( CDelNodes[ii] ) == fmu_kCLEAN_ABORT )
         return fmu_kCLEAN_ABORT;
   }

   //Make three new CQuads
   elmQClass *pCQuad1 = new elmQClass( pCNode0, pCNode1, pCNode4, pCNode5,
                                       pCMaster );
   if( !pCQuad1 ) return fmu_kCLEAN_NO_MEMORY;
   elmQClass *pCQuad2 = new elmQClass( pCNode1, pCNode2, pCNodeB, pCNode4,
                                       pCMaster );
   if( !pCQuad2 ) return fmu_kCLEAN_NO_MEMORY;
   elmQClass *pCQuad3 = new elmQClass( pCNode2, pCNode3, pCNodeA, pCNodeB,
                                       pCMaster );
   if( !pCQuad3 ) return fmu_kCLEAN_NO_MEMORY;

   if( lDrawCleaning )
   {
   //nodClassDynArray CSmoothList;
   //CSmoothList.append( pCNodeA );
   //CSmoothList.append( pCNodeB );
   //for( ii = 0; ii < CRingNodes.length(); ++ii )
   //CSmoothList.append( CRingNodes[ii] );
   //smooth_mesh( &CSmoothList );
        nodClassDynArray CSmoothList;
        CSmoothList.append( pCNodeA );
        CSmoothList.append( pCNodeB );
        for( ii = 0; ii < CRingNodes.length(); ++ii )
            CSmoothList.append( CRingNodes[ii] );
        smooth_mesh( &CSmoothList );
   //draw_mesh();
   //pCQuad1->draw_me( dp_kCLR_YELLOW );
   //pCQuad2->draw_me( dp_kCLR_GOLD );
   //pCQuad3->draw_me( dp_kCLR_BLUE );
        // Refresh debug graphics
        //EraseElements();
        pCQuad1->draw_me( dp_kCLR_BLUE );
        pCQuad2->draw_me( dp_kCLR_BLUE );
        pCQuad3->draw_me( dp_kCLR_BLUE );
   //draw_close();
   }


   for( ii = 0; ii < CRingNodes.length(); ++ii ) reclean( CRingNodes[ii] );
   reclean( pCNodeB );

   return fmu_kCLEAN_MESH_CHANGED;
  
}


int QuadCleanTool::fill_choice_sidem( nodClassDynArray &CDelNodes,
                                          fmuEdgeDynArray &CDelEdges,
                                          elmQClassDynArray &CDelQuads,
                                          nodClassDynArray &CRingNodes )
{
   
   int ii, nCount;

   if( CRingNodes.length() != 6 ) return fmu_kCLEAN_NO_ACTION;

   //The ringnodes list has an assumed orientation
   nodClass *pCNode0 = CRingNodes.get();
   nodClass *pCNode1 = CRingNodes.next();
   nodClass *pCNode2 = CRingNodes.next(2);
   nodClass *pCNode3 = CRingNodes.next(3);
   nodClass *pCNode4 = CRingNodes.next(4);
   nodClass *pCNode5 = CRingNodes.next(5);
   fmuEdge *pCEdge45 = pCNode4->shared_edge( pCNode5 );
   if( pCEdge45->frozen() ) return fmu_kCLEAN_NO_ACTION;
   fmuEdge *pCEdge50 = pCNode5->shared_edge( pCNode0 );
   if( pCEdge50->frozen() ) return fmu_kCLEAN_NO_ACTION;

   //Orientation of CDelQuads does not matche orientation of CRingNodes
   elmQClass *pCQuadA, *pCQuadX, *pCQuadY;
   pCEdge50->quads( &pCQuadX, &pCQuadY );
   if( CDelQuads.find( pCQuadX ) > -1 ) pCQuadA = pCQuadY;
   else                                 pCQuadA = pCQuadX;

   nodClass *pCNodeA = pCQuadA->next_node( pCNode5 );
   nodClass *pCNodeB = pCQuadA->opposite_node( pCNode5 );

   if( pCNode0->hardpt() && pCNode1->hardpt() && pCNodeB->hardpt() )
   {
      fmuVector CVec01( pCNode0, pCNode1 );
      fmuVector CVec0B( pCNode0, pCNodeB );
      if( CNormal.angle( CVec01, CVec0B ) > fmu_kSHAPE_ANGLE_TOL )
         return fmu_kCLEAN_NO_ACTION;
   }

   nCount = CDelQuads.length();
   if( lDrawCleaning )
   {
   //draw_open();
        // Refresh debug graphics
        //EraseElements();
        for( ii = 0; ii < nCount; ++ii )
        {
   //CDelQuads[ii]->draw_me( dp_kCLR_YELLOW + ii );
            CDelQuads[ii]->draw_me(  dp_kCLR_YELLOW );
        }
   //pCQuadA->draw_me( dp_kCLR_BLUE );
        pCQuadA->draw_me( dp_kCLR_YELLOW );
   //draw_close();
   }

   //Delete everything in the deletion lists, plus an extra edge & quad
   for( ii = 0; ii < nCount; ++ii )
   {
      delete_this_quad( CDelQuads[ii] );
   }
   delete_this_quad( pCQuadA );
   nCount = CDelEdges.length();
   for( ii = 0; ii < nCount; ++ii )
   {
      if( delete_this_edge( CDelEdges[ii] ) == fmu_kCLEAN_ABORT )
         return fmu_kCLEAN_ABORT;
   }
   if( delete_this_edge( pCEdge50 ) == fmu_kCLEAN_ABORT )
      return fmu_kCLEAN_ABORT;
   nCount = CDelNodes.length();
   for( ii = 0; ii < nCount; ++ii )
   {
      if( delete_this_node( CDelNodes[ii] ) == fmu_kCLEAN_ABORT )
         return fmu_kCLEAN_ABORT;
   }

   //Make three new CQuads
   elmQClass *pCQuad1 = new elmQClass( pCNode0, pCNode1, pCNodeA, pCNodeB,
                                       pCMaster );
   if( !pCQuad1 ) return fmu_kCLEAN_NO_MEMORY;
   elmQClass *pCQuad2 = new elmQClass( pCNode1, pCNode2, pCNode5, pCNodeA,
                                       pCMaster );
   if( !pCQuad2 ) return fmu_kCLEAN_NO_MEMORY;
   elmQClass *pCQuad3 = new elmQClass( pCNode2, pCNode3, pCNode4, pCNode5,
                                       pCMaster );
   if( !pCQuad3 ) return fmu_kCLEAN_NO_MEMORY;

   if( lDrawCleaning )
   {
   //nodClassDynArray CSmoothList;
   //CSmoothList.append( pCNodeA );
   //CSmoothList.append( pCNodeB );
   //for( ii = 0; ii < CRingNodes.length(); ++ii )
   //CSmoothList.append( CRingNodes[ii] );
   //smooth_mesh( &CSmoothList );
   //draw_mesh();
    // Refresh debug graphics
    //EraseElements();
   //pCQuad1->draw_me( dp_kCLR_YELLOW );
   //pCQuad2->draw_me( dp_kCLR_GOLD );
   //pCQuad3->draw_me( dp_kCLR_BLUE );
        nodClassDynArray CSmoothList;
        CSmoothList.append( pCNodeA );
        CSmoothList.append( pCNodeB );
        for( ii = 0; ii < CRingNodes.length(); ++ii )
            CSmoothList.append( CRingNodes[ii] );
        smooth_mesh( &CSmoothList );
   //draw_mesh();
        pCQuad1->draw_me( dp_kCLR_BLUE );
        pCQuad2->draw_me( dp_kCLR_BLUE );
        pCQuad3->draw_me( dp_kCLR_BLUE );
   //draw_close();
   }


   for( ii = 0; ii < CRingNodes.length(); ++ii ) reclean( CRingNodes[ii] );
   reclean( pCNodeB );

   return fmu_kCLEAN_MESH_CHANGED;
   
}


int QuadCleanTool::fill_choice_squeeze( nodClassDynArray &CDelNodes,
                                            fmuEdgeDynArray &CDelEdges,
                                            elmQClassDynArray &CDelQuads,
                                            nodClassDynArray &CRingNodes )
{
   
   int ii, jj, kk;
   int iValence;
   elmQClassDynArray COuterQuads(10);
   elmQClassDynArray CQuads(10);
   elmQClassDynArray CNewQuads(10);
   nodClassDynArray CNodes(10);
   fmuEdgeDynArray CEdges(10);
   //This routine squeezes out a double three pattern by merging nodes
   //1 with 5 and 2 with 4. Which node is kept depends on what side has
   //frozen edges. If both do, don't merge. The merge is accomplished by
   //deleting all quads in the interior of the ring, then deleting and
   //rebuilding the quads along one side.
   nodClass *pCNode0 = CRingNodes.get();
   nodClass *pCNode1 = CRingNodes.next();
   nodClass *pCNode2 = CRingNodes.next(2);
   nodClass *pCNode3 = CRingNodes.next(3);
   fmuEdge *pCEdge01 = pCNode0->shared_edge( pCNode1 );
   fmuEdge *pCEdge12 = pCNode1->shared_edge( pCNode2 );
   fmuEdge *pCEdge23 = pCNode2->shared_edge( pCNode3 );
   if( pCEdge01->frozen() || pCEdge12->frozen() || pCEdge23->frozen() ||
       pCNode1->hardpt() || pCNode2->hardpt() )
   {
      //Rotate the list halfway and check the other side for frozen edges
      CRingNodes.step(3);
      pCNode0 = CRingNodes.get();
      pCNode1 = CRingNodes.next();
      pCNode2 = CRingNodes.next(2);
      pCNode3 = CRingNodes.next(3);
      pCEdge01 = pCNode0->shared_edge( pCNode1 );
      pCEdge12 = pCNode1->shared_edge( pCNode2 );
      pCEdge23 = pCNode2->shared_edge( pCNode3 );
      if( pCEdge01->frozen() || pCEdge12->frozen() || pCEdge23->frozen() ||
          pCNode1->hardpt() || pCNode2->hardpt() )
         return fmu_kCLEAN_NO_ACTION;
   }
   nodClass *pCNode4 = CRingNodes.next(4);
   nodClass *pCNode5 = CRingNodes.next(5);

   nodClassDynArray CSmoothList;
   CSmoothList.append( pCNode0 );
   CSmoothList.append( pCNode3 );
   CSmoothList.append( pCNode4 );
   CSmoothList.append( pCNode5 );

   elmQClass *pCQuadA, *pCQuadB;
   pCEdge01->quads( &pCQuadA, &pCQuadB );
   if( CDelQuads.find( pCQuadA ) > -1 ) pCQuadA = pCQuadB;
   //Check whether a bad angle will be created
   nodClass *pCNodeX = pCQuadA->next_node( pCNode0 );
   if( pCNode0->hardpt() && pCNode5->hardpt() && pCNodeX->hardpt() )
   {
      fmuVector CVec0X( pCNode0, pCNodeX );
      fmuVector CVec05( pCNode0, pCNode5 );
      if( CNormal.angle( CVec0X, CVec05 ) > fmu_kSHAPE_ANGLE_TOL )
         return fmu_kCLEAN_NO_ACTION;
   }

   //Don't do the squeeze if the valence is too low or too high. Do a
   //valence pattern instead. The rest of the processing allows for 4
   //and 5 valent nodes.
   iValence = neighbors( pCNode1, CNodes, CEdges, CQuads );
   if( iValence < 4 || iValence > 5 ) return fmu_kCLEAN_NO_ACTION;
   CNodes.find( pCNode0 );
   CEdges.find( pCEdge01 );
   CQuads.find( pCQuadA );
   fmuEdge *pCNextEdge = NULL;
   fmuEdge *pCEdge = 0;
   elmQClass *pCNextQuad = 0;
   for( ii = 0; pCNextEdge != pCEdge12; ++ii )
   {
      pCNextEdge = CEdges.next(ii+1);
      if( pCNextEdge->frozen() )return fmu_kCLEAN_NO_ACTION;
      if( pCNextEdge != pCEdge12 ) CDelEdges.append( pCNextEdge );
      pCNextQuad = CQuads.next(ii);
      COuterQuads.append( pCNextQuad );
      CSmoothList.append( CNodes.next(ii*2+1));
   }

   CNodes.clear();
   CEdges.clear();
   CQuads.clear();
   iValence = neighbors( pCNode2, CNodes, CEdges, CQuads );
   if( iValence < 4 || iValence > 5 ) return fmu_kCLEAN_NO_ACTION;
   CNodes.find( pCNode1 );
   CNodes.step(2);
   CEdges.find( pCEdge12 );
   CEdges.step();
   CQuads.find( pCNextQuad );
   CQuads.step();
   pCNextEdge = NULL;
   for( ii = 0; pCNextEdge != pCEdge23; ++ii )
   {
      pCEdge = CEdges.next(ii);
      if( pCEdge->frozen() )return fmu_kCLEAN_NO_ACTION;
      CDelEdges.append( pCEdge );
      pCNextQuad = CQuads.next(ii);
      COuterQuads.append( pCNextQuad );
      pCNextEdge = CEdges.next(ii+1);
      CSmoothList.append( CNodes.next(ii*2+1));
   }
   //Check whether a bad angle will be created
   pCNodeX = pCNextQuad->prev_node( pCNode3 );
   if( pCNode4->hardpt() && pCNode3->hardpt() && pCNodeX->hardpt() )
   {
      fmuVector CVec34( pCNode3, pCNode4 );
      fmuVector CVec3X( pCNode3, pCNodeX );
      if( CNormal.angle( CVec34, CVec3X ) > fmu_kSHAPE_ANGLE_TOL )
         return fmu_kCLEAN_NO_ACTION;
   }
   

   CDelEdges.append( pCEdge01 );
   CDelEdges.append( pCEdge12 );
   CDelEdges.append( pCEdge23 );
   CDelNodes.append( pCNode1 );
   CDelNodes.append( pCNode2 );

   if( lDrawCleaning )
   {
   //draw_open();
   //for( ii = 0; ii < CDelQuads.length(); ++ii )
   //{
   //CDelQuads[ii]->draw_me( dp_kCLR_BLUE + ii );
   //}
   //for( ii = 0; ii < COuterQuads.length(); ++ii )
   //{
   //COuterQuads[ii]->draw_me( dp_kCLR_YELLOW + ii );
   //}
    // Refresh debug graphics
    //EraseElements();

    for( ii = 0; ii < CDelQuads.length(); ++ii )
    {
        CDelQuads[ii]->draw_me( dp_kCLR_YELLOW );
    }
    for( ii = 0; ii < COuterQuads.length(); ++ii )
    {
        COuterQuads[ii]->draw_me( dp_kCLR_YELLOW );
    }
   //draw_close();
   }

   //The outer quads are used as a guide for creating the new quads so
   //make a table of the nodes in the new quads before deleting the old ones.
   int nNewNodes = 0;
   nodClass *apNodeList[28];
   kk = 0;
   for( ii = 0; ii < COuterQuads.length(); ++ii )
   {
      pCNextQuad = COuterQuads[ii];
      pCNextQuad->nodes( &(apNodeList[kk]) );
      for( jj = 0; jj < 4; ++jj)
      {
         if( apNodeList[kk+jj] == pCNode1 )      apNodeList[kk+jj] = pCNode5;
         else if( apNodeList[kk+jj] == pCNode2 ) apNodeList[kk+jj] = pCNode4;
      }
      kk += 4;
      ++nNewNodes;
   }
   //Delete everything in the deletion lists
   int nCount;
   nCount = CDelQuads.length();
   for( ii = 0; ii < nCount; ++ii )
   {
      delete_this_quad( CDelQuads[ii] );
   }
   nCount = COuterQuads.length();
   for( ii = 0; ii < nCount; ++ii )
   {
      delete_this_quad( COuterQuads[ii] );
   }
   nCount = CDelEdges.length();
   for( ii = 0; ii < nCount; ++ii )
   {
      if( delete_this_edge( CDelEdges[ii] ) == fmu_kCLEAN_ABORT )
         return fmu_kCLEAN_ABORT;
   }
   nCount = CDelNodes.length();
   for( ii = 0; ii < nCount; ++ii )
   {
      if( delete_this_node( CDelNodes[ii] ) == fmu_kCLEAN_ABORT )
         return fmu_kCLEAN_ABORT;
   }

   kk =0;
   for( ii = 0; ii < nNewNodes; ++ii)
   {
      CNewQuads[ii] = new elmQClass( apNodeList[kk], apNodeList[kk+1],
                            apNodeList[kk+2], apNodeList[kk+3], pCMaster );
      if( !CNewQuads[ii] ) return fmu_kCLEAN_NO_MEMORY;
      kk+= 4;
   }

   for( ii = 0; ii < CSmoothList.length(); ++ii ) reclean( CSmoothList[ii] );

   if( lDrawCleaning )
   {
        smooth_mesh( &CSmoothList );
   
        for( ii = 0; ii < CNewQuads.length(); ++ii )
        {
            CNewQuads[ii]->draw_me(  dp_kCLR_YELLOW );
        }
   //draw_close();
   }

   return fmu_kCLEAN_MESH_CHANGED;

  
}


void QuadCleanTool::two_middle_nodes( nodClass* pCNode1, nodClass* pCNode2,
                                      nodClass** ppCNodeA,
                                      nodClass** ppCNodeB )
{
  
   fmuVector CDirectionVec( pCNode1, pCNode2 );
   CDirectionVec /= 3.0;
   fmuVector CVec1( pCNode1 );
   fmuVector CNodeCoords = CVec1 + CDirectionVec;
   *ppCNodeA = new nodClass( CNodeCoords, pCMaster );
   CNodeCoords = CVec1 + ( CDirectionVec * 2.0 );
   *ppCNodeB = new nodClass( CNodeCoords, pCMaster );
  
}

void QuadCleanTool::get_two_middle_node_location( nodClass* pCNode1, nodClass* pCNode2,
                                                  double adNodeX[], double adNodeY[] )
{
  
   fmuVector CDirectionVec( pCNode1, pCNode2 );
   CDirectionVec /= 3.0;
   fmuVector CVec1( pCNode1 );
   fmuVector CNodeCoords = CVec1 + CDirectionVec;
   adNodeX[0]=CNodeCoords.x();
   adNodeX[1]=CNodeCoords.y();
   CNodeCoords = CVec1 + ( CDirectionVec * 2.0 );
   adNodeY[0]=CNodeCoords.x();
   adNodeY[1]=CNodeCoords.y();
}

void QuadCleanTool::delete_this_quad( elmQClass *pCQuad )
{
    if(!pCQuad)
        return;

    int iPosition = pCQuad->id();
    (*pCMaster->global_quads())[iPosition] = NULL;

    //Don't extend the shape list.
    if(iPosition < CShapeList.length())
        CShapeList[iPosition] = NULL;

    delete pCQuad;
}


int QuadCleanTool::delete_this_edge( fmuEdge *pCEdge )
{
    if(!pCEdge)
        return fmu_kCLEAN_NO_ACTION;

    if( pCEdge->number_elems() != 0 )
        return fmu_kCLEAN_ABORT;

    int iPosition = pCEdge->id();
    (*pCMaster->global_edges())[iPosition] = NULL;
    delete pCEdge;

    return fmu_kCLEAN_MESH_CHANGED;
}


int QuadCleanTool::delete_this_node( nodClass *pCNode )
{
    if( !pCNode )
        return fmu_kCLEAN_NO_ACTION;
    if( pCNode->number_edges() != 0 )
        return fmu_kCLEAN_ABORT;

    int iPosition = pCMaster->global_interior()->find( pCNode );
    if( iPosition > -1 )
        (*pCMaster->global_interior())[iPosition] = NULL;

    iPosition = CCleaningList.find( pCNode );
    if( iPosition > -1 )
        CCleaningList[iPosition] = NULL;

    delete pCNode;

    return fmu_kCLEAN_MESH_CHANGED;
}


void QuadCleanTool::reclean( nodClass *pCNode )
{
   
   int ii, iValence;
   fmuEdgeDynArray CEdges(ARRAY_SIZE);
   fmuEdge *pCEdge;
   nodClass *pCOtherNode;

   if (!pCNode)
	   return;

   //Assume that nodes currently in the cleaning list have the flag set. As
   //nodes are put back in the list, set the flag again. Asking the node if
   //it is in the list is not ideal, but is faster than the alternative.
   if( !pCNode->boundary() )
   {
      if( !(pCNode->in_clean_list()) )
      {
         CCleaningList.append( pCNode );
         pCNode->in_clean_list(true);
      }
   }

   iValence = pCNode->number_edges();
   pCNode->edge_list( CEdges );
   for( ii = 0; ii < iValence; ++ii )
   {
      pCEdge = CEdges[ii];
      pCOtherNode = pCEdge->other_node( pCNode );
      if( pCOtherNode->boundary() ) continue;
      if( !(pCOtherNode->in_clean_list()) )
      {
         CCleaningList.append( pCOtherNode );
         pCOtherNode->in_clean_list(true);
      }
   }
 
}


//  quad mesh cleaning Member Function Descriptor Block -----
//
// QuadCleanTool::AnalyzeMesh
// class - QuadCleanTool
// 
//
// Description:
//
// Perform a cleaning analysis on the mesh.
//
// Access:
//
// AnalyzeMesh( pCMeshMaster )
//
// Input Parameters:
//
// pCMeshMaster         fmuMaster*          Lists of nodes, edges, elements
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
// The mesh analysis is appended to the end of file "mesh.clean". The
// analysis is given a header to show differentiate the report for each
// surface in case more than one surface is meshed at a time.
//
// The report has 3 sections. The first lists valence patterns that probably
// can be cleaned up. These patterns have 3 or more irregular (non-4 valent)
// nodes. All nodes with a valence of 2 or less or 6 or more are reported.
//
// The second part lists bad boundary conditions: (1) too few or too namy
// edges into a boundary node based on the angle at the node, and (2)
// boundary valence patterns that probably be cleaned.
//
// The third part lists shape problems. It includes nodes that are too big
// and are too small.
//
// ------------------------------------------
void QuadCleanTool::AnalyzeMesh( fmuMaster *pCMeshMaster )
{
  
   int iDebugFlag = 0;
   int ii, jj, kk;
   int iValence, iCheckValence, iSurfaceValence, iBoundaryValence;
   int iNumberNodes;
   int iIrregCount;
   int aiValenceList[10];
   char cMessage[80];
   double dAngle, dDegrees;
   nodClass *pCNode, *pCPrevNode, *pCNextNode, *apCNodes[4];
   elmQClass *pCQuad;
   nodClassDynArray CNodes(ARRAY_SIZE);
   fmuEdgeDynArray CEdges(ARRAY_SIZE);
   elmQClassDynArray CQuads(ARRAY_SIZE);

   pCMaster = pCMeshMaster;

   FILE *pzReport;

   pzReport = fopen( "d:\\workdir\\sam\\mesh.clean", "a" );

   fprintf( pzReport, "--------------------------------------------------\n");
   fprintf( pzReport, "                 Mesh analysis\n");
   fprintf( pzReport, "--------------------------------------------------\n");

   //Valence analysis
   for( ii = 0; ii < pCMaster->global_interior()->length(); ++ii )
   {
      CNodes.clear();
      CEdges.clear();
      CQuads.clear();
      pCNode = (*pCMaster->global_interior())[ii];
      if( !pCNode ) continue;
      //Nodes around trias are still in the interior list but cannot be
      //analyzed that way.
      if( pCNode->boundary() ) continue;

      iValence = get_valence( pCNode );

      if( iValence < 3 || iValence > 5 )
      {
         sprintf( cMessage, "%1d valent node", iValence );
         print_node( pzReport, pCNode, cMessage );
      }
      else
      {
         iCheckValence = neighbors( pCNode, CNodes, CEdges, CQuads );
         if( iCheckValence == 0 )
            print_node( pzReport, pCNode, "Disassociated from mesh");
         else
         {
            if( iDebugFlag )
            {
            //draw_mesh();
            //draw_close();
            //for( jj = 0; jj < iValence * 2; ++jj )
            //{
            //nodClass *pCNode = CNodes[jj];
            //if( pCNode->boundary() )
            //pCNode->draw_me( dp_kCLR_BLUE );
            //else if( pCNode->hardpt() )
            //pCNode->draw_me( dp_kCLR_RED );
            //else
            //pCNode->draw_me( dp_kCLR_GOLD );
            //}
            }
            if( iValence == 4 ) iIrregCount = 0;
            else                iIrregCount = 1;
            for( jj = 0; jj < iValence * 2; ++jj )
            {
               aiValenceList[jj] = get_valence( CNodes[jj] );
               if( aiValenceList[jj]  != 4 ) ++iIrregCount;
            }
            if( iIrregCount > 2 )
            {
               iNumberNodes = iValence * 2;
               for( jj = 0;
                    jj < iNumberNodes && aiValenceList[jj] == 4 &&
                                         aiValenceList[jj+1] == 4; jj+=2);
               if(jj > 0 ) rotate_list( aiValenceList, iNumberNodes, jj );
               sprintf( cMessage, "Valence score  %1d, Pattern %1d-",
                        iIrregCount, iValence );
               for( jj = 0; jj < iNumberNodes; ++jj )
               {
                  sprintf( &(cMessage[ 28+jj ]), "%1d", aiValenceList[jj]);
               }
               print_node( pzReport, pCNode, cMessage );
            }
         }
      }
   }

   //Boundary analysis
   nodClassDynArrayDynArray *boundaries = pCMaster->global_boundary();
   for( ii = 0; ii < boundaries->length(); ++ii )
   {
      nodClassDynArray *boundary = (*boundaries)[ii];
      for( kk = 0; kk < boundary->length(); ++kk )
      {
         CNodes.clear();
         CEdges.clear();
         CQuads.clear();
         boundary->set( kk );
         pCNode = boundary->get();
         pCPrevNode = boundary->prev();
         pCNextNode = boundary->next();
         iSurfaceValence = pCNode->number_edges() - 2;
         fmuVector CPrevVec( pCNode, pCPrevNode );
         fmuVector CNextVec( pCNode, pCNextNode );
         dAngle = CNormal.angle( CNextVec, CPrevVec );
         if( dAngle < fmu_kTWO_VALENT_ANGLE )
         {
            if( iSurfaceValence > 1 )
            {
               sprintf( cMessage, "valence of %d in small angle",
                                                      iSurfaceValence );
               print_node( pzReport, pCNode, cMessage );
            }
         }
         else if( dAngle < fmu_kONE_VALENT_ANGLE )
         {
            if ( iSurfaceValence < 1 || iSurfaceValence > 2  )
            {
               sprintf( cMessage, "valence of %d in flat angle",
                                                      iSurfaceValence );
               print_node( pzReport, pCNode, cMessage );
            }
         }
         else if( iSurfaceValence < 2 || iSurfaceValence > 4  )
         {
            sprintf( cMessage, "valence of %d in big angle",
                                                      iSurfaceValence );
            print_node( pzReport, pCNode, cMessage );
         }
   
         iBoundaryValence = boundary_neighbors( boundary,
                                                CNodes, CEdges, CQuads );
 
         if( iBoundaryValence == 0 )
            print_node( pzReport, pCNode, "Disassociated from mesh");
         else
         {
            iValence = get_valence( pCNode );
            if( iValence == 4 ) iIrregCount = 0;
            else                iIrregCount = 1;
            iNumberNodes = CNodes.length();
            for( jj = 0; jj < iNumberNodes; ++jj )
            {
               aiValenceList[jj] = get_valence( CNodes[jj] );
               if( aiValenceList[jj]  != 4 ) ++iIrregCount;
            }
            if( iIrregCount > 2 )
            {
               sprintf( cMessage, "Boundary score %1d, Pattern %1d-",
                              iIrregCount, iValence );
               for( jj = 0; jj < iNumberNodes; ++jj )
               {
                  sprintf( &(cMessage[ 28+jj ]), "%1d", aiValenceList[jj]);
               }
               print_node( pzReport, pCNode, cMessage );
            }
         }
         boundary->step();
      }
   }

   //Shape analysis
   for( ii = 0; ii < pCMaster->global_quads()->length(); ++ii )
   {
      pCQuad = (*pCMaster->global_quads())[ii];
      if( !pCQuad ) continue;
      pCQuad->nodes( apCNodes );
      for( jj = 0; jj < 4; ++jj )
      {
         pCNode = apCNodes[jj];
         pCPrevNode = pCQuad->prev_node( pCNode );
         pCNextNode = pCQuad->next_node( pCNode );
         fmuVector CPrevVec( pCNode, pCPrevNode );
         fmuVector CNextVec( pCNode, pCNextNode );
         dAngle = CNormal.angle( CNextVec, CPrevVec );
         dDegrees = dAngle / PI * 180.0;
         if( dDegrees < 30.0 || dDegrees > 160.0 )
         {
            sprintf( cMessage, " bad angle of %6.1f", dDegrees );
            print_node( pzReport, pCNode, cMessage );
         }
      }
   }

   fclose( pzReport );
   
}


int QuadCleanTool::neighbors( nodClass *pCNode, nodClassDynArray &CNodes,
                                  fmuEdgeDynArray &CEdges,
                                  elmQClassDynArray &CQuads )
{
   
   int ii=0, jj=0;
   int iValence=0;
   fmuEdge *pCNextEdge=0;
   nodClass *pCNextNode=0;
   elmQClass *pCNextQuad=0, *pCQuad1=0, *pCQuad2=0;

   iValence = pCNode->number_edges();
   // This should "never" happen in a good mesh. It is here in case the
   // mesh is screwed up. It means "don't clean this node". Other anomalies
   // are left to cause a segfault as it is likely a cleaning programmer
   // error (not checking something).
   if( iValence == 0 ) return 0;

   pCNextEdge = pCNode->first_edge();
   pCNextNode = pCNextEdge->other_node( pCNode );
   pCNextEdge->quads( &pCQuad1, &pCQuad2 );

   if (!pCQuad1 || !pCQuad2)
       return 0;

   if( pCQuad1->next_node( pCNode ) == pCNextNode )
      pCNextQuad = pCQuad1;
   else
   {
      if( pCQuad2->next_node( pCNode ) == pCNextNode )
         pCNextQuad = pCQuad2;
      else
         return 0;
   }

   jj = 0;
   for( ii = 0; ii < iValence; ++ii )
   {
      CNodes[jj] = pCNextNode;
      CNodes[jj+1] = pCNextQuad->opposite_node( pCNode );
      CEdges[ii] = pCNextEdge;
      CQuads[ii] = pCNextQuad;
      pCNextNode = pCNextQuad->prev_node( pCNode );
      if( !pCNextNode )
         return 0;
      pCNextEdge = pCNode->shared_edge( pCNextNode );
      if( !pCNextEdge )
         return 0;
      pCNextQuad = pCNextEdge->other_quad( pCNextQuad );
      if( !pCNextQuad )
         return 0;
      jj += 2;
   }

   return iValence;
  
}


int QuadCleanTool::boundary_neighbors( nodClassDynArray *pCBoundary,
                                       nodClassDynArray &CNodes,
                                       fmuEdgeDynArray &CEdges,
                                       elmQClassDynArray &CQuads )
{
    if( pCBoundary == NULL )
        return 0;

    if( pCBoundary->empty() )
        return 0;

    int iDebugFlag = 0;
    int ii, jj;

    elmQClass *pCNextQuad, *pCQuad1, *pCQuad2;
    fmuEdge *pCLastEdge;

    nodClass *pCNode = pCBoundary->get();
    nodClass *pCPrevNode = pCBoundary->prev();
    nodClass *pCNextNode = pCBoundary->next();
    fmuEdge *pCPrevEdge = pCNode->shared_edge( pCPrevNode );
    fmuEdge *pCNextEdge = pCNode->shared_edge( pCNextNode );

    if( iDebugFlag )
    {
        //draw_mesh();
        //pCPrevEdge->draw_me( dp_kCLR_RED );
        //pCNextEdge->draw_me( dp_kCLR_CYAN );
        //draw_close();
    }

    //These lists are sometimes used consecutively for different nodes.
    //Force them to be empty.
    CNodes.clear();
    CEdges.clear();
    CQuads.clear();

    //This edge is along the boundary so there is no second quad. There
    //might be a tria which this boundary is isolating from the quad mesh.
    pCNextEdge->quads( &pCQuad1, &pCQuad2 );
    if(!pCQuad1)
        return -1;
    if( pCQuad1->type() == SAM_TRIA )
    {
        if( !pCQuad2 )
            return 0;
        pCQuad1 = pCQuad2;
    }

    if( pCQuad1->next_node( pCNode ) == pCNextNode )
        pCNextQuad = pCQuad1;
    else
        return 0;

    if( iDebugFlag )
    {
        //draw_open();
        // Refresh debug graphics
        //EraseElements();
        //pCNextQuad->draw_me( dp_kCLR_BLUE );
        pCNextQuad->draw_me( dp_kCLR_YELLOW );
        //draw_close();
    }
    jj = 0;
    int iValence = 0;
    for( ii = 0; pCNextEdge != pCPrevEdge; ++ii )
    {
        CNodes[jj] = pCNextNode;
        CNodes[jj+1] = pCNextQuad->opposite_node( pCNode );
        CEdges[ii] = pCNextEdge;
        CQuads[ii] = pCNextQuad;
        pCNextNode = pCNextQuad->prev_node( pCNode );
        if( !pCNextNode )
            return 0;
        pCLastEdge = pCNextEdge;
        pCNextEdge = pCNode->shared_edge( pCNextNode );
        if( !pCNextEdge || pCNextEdge == pCLastEdge )
            return 0;
        pCNextQuad = pCNextEdge->other_quad( pCNextQuad );
        //Next quad doesn't matter if the next edge is the last one.
        if( pCNextEdge != pCPrevEdge && !pCNextQuad )
            return 0;
        jj += 2;
        ++iValence;
    }

    CNodes[jj] = pCPrevNode;
    CEdges[iValence] = pCPrevEdge;
    return iValence;
}


double QuadCleanTool::angle_on_boundary ( nodClassDynArray *boundary )
{
   

   int iDebugFlag = 0;
   nodClass *pCNode, *pCPrevNode, *pCNextNode;
   double dAngle;
   int iBndyIndex, iHalfLoop;

   pCNode = boundary->get();
   pCPrevNode = boundary->prev();
   pCNextNode = boundary->next();
   if( iDebugFlag )
   {
   //draw_mesh();
   //draw_close();
   //pCNode->draw_me( dp_kCLR_MAGNTA );
   //pCPrevNode->draw_me( dp_kCLR_YELLOW );
   //pCNextNode->draw_me( dp_kCLR_LTBLUE );
   }

   fmuVector CPrevVec( pCNode, pCPrevNode );
   fmuVector CNextVec( pCNode, pCNextNode );
   dAngle = CNormal.angle( CNextVec, CPrevVec );

   if( dAngle < fmu_kSMALL_BNDY_ANGLE_TOL || 
       dAngle > fmu_kBIG_BNDY_ANGLE_TOL )
   {
      iBndyIndex=1;
      iHalfLoop = boundary->length() / 2;
      while( dAngle < fmu_kSMALL_BNDY_ANGLE_TOL || 
             dAngle > fmu_kBIG_BNDY_ANGLE_TOL )
      {
         pCPrevNode = boundary->prev(iBndyIndex);
         pCNextNode = boundary->next(iBndyIndex);

         //Give up if wrapped around the loop
         if( !pCPrevNode || !pCNextNode ) break;
         if( pCPrevNode == pCNextNode ) break;
         if( iBndyIndex > iHalfLoop ) break;

         if( iDebugFlag )
         {
         //draw_mesh();
         //draw_close();
         //pCNode->draw_me( dp_kCLR_MAGNTA );
         //pCPrevNode->draw_me( dp_kCLR_YELLOW );
         //pCNextNode->draw_me( dp_kCLR_LTBLUE );
         }

         fmuVector CPrevVec( pCNode, pCPrevNode );
         fmuVector CNextVec( pCNode, pCNextNode );
         dAngle = CNormal.angle( CNextVec, CPrevVec );
         iBndyIndex++;
      }
      if(dAngle > PI) 
         dAngle = TWOPI;
      else
         dAngle = 0;
   }

   return dAngle;
  
}


int QuadCleanTool::get_valence( nodClass *pCNode )
{
   
   if (!pCNode)
	   return(0);

   if( pCNode->boundary() )
   {
      int iValence, iIndex;
      double dAngle;
      nodClass *pCSaveNode;
      nodClassDynArray *pCLoop;
      iValence = pCNode->number_edges();
      if( find_on_boundary( pCNode, &pCLoop, &iIndex, &pCSaveNode ) )
      {
         fmuVector CToPrev( pCNode, pCLoop->prev() );
         fmuVector CToNext( pCNode, pCLoop->next() );
         dAngle = CNormal.angle( CToNext, CToPrev );
         //If angle is close to 2PI, treat the two boundary edges as one edge
         if( dAngle > fmu_kZERO_VALENT_ANGLE ) iValence-=1;
         //If angle is small, assume two outer edges
         else if( dAngle < fmu_kTWO_VALENT_ANGLE ) iValence += 2;
         //If angle is moderate, assume one outer edge
         else if( dAngle < fmu_kONE_VALENT_ANGLE ) ++iValence;
      }
      pCLoop->set( iIndex );
      return iValence;
   }
   else
      return pCNode->number_edges();
  
}


bool QuadCleanTool::find_on_boundary( nodClass *pCSearchNode,
                                         nodClassDynArray **ppCLoop,
                                         int *piIndex,
                                         nodClass **ppCSaveNode )
{
 
   //The index returned is the loop position before the search was done.
   int ii, iPosition;
   nodClassDynArrayDynArray *pCboundaries = pCMaster->global_boundary();
   for( ii = 0; ii < pCboundaries->length(); ++ii )
   {
      *ppCLoop = (*pCboundaries)[ii];
      if((*ppCLoop)->length() == 0) continue;

      *ppCSaveNode = (*ppCLoop)->get();
      //Just a way to get the index
      *piIndex = (*ppCLoop)->step(0);
      //If at the current position, don't search as there might be
      //duplicates due to loops around triangles
      if( *ppCSaveNode == pCSearchNode )
      {
         return true;
      }
      iPosition =(*ppCLoop)->find_near( pCSearchNode );
      if( iPosition >= 0 )
      {
         return true;
      }
   }
   return false;
   
}


void QuadCleanTool::rotate_list( int *aiList, int nCount,
                                                  int iRotate )
{

   int aiHolder[ARRAY_SIZE];
   int ii, jj;
   jj = iRotate;
   for( ii = 0; ii < nCount; ++ii)
   {
      aiHolder[ii] = aiList[jj];
      ++jj;
      if( jj == nCount) jj = 0;
   }
   for( ii = 0; ii < nCount; ++ii)
      aiList[ii] = aiHolder[ii];

}


void QuadCleanTool::rotate_doubles( double *adList,
                                    int nCount,
                                    int iRotate )
{

   double adHolder[ARRAY_SIZE];
   int ii, jj;
   jj = iRotate;
   for( ii = 0; ii < nCount; ++ii)
   {
      adHolder[ii] = adList[jj];
      ++jj;
      if( jj == nCount) jj = 0;
   }
   for( ii = 0; ii < nCount; ++ii)
      adList[ii] = adHolder[ii];

}


void QuadCleanTool::print_node( FILE *pzReport, nodClass *pCNode, const char *pcMessage )
{

      //int iUnits = ULEN;
      //int iType = 4;
      //int nVal = 3;
   //double adCoords[3], adUserCoords[3];
   //adCoords[0] = pCNode->node_x();
   //adCoords[1] = pCNode->node_y();
   //adCoords[2] = pCNode->node_z();

   //CALL_FORTRAN(untstd) ( &iUnits, &iType, &nVal, (float*) adCoords, 
   //                                               (float*) adUserCoords );

   //fprintf( pzReport, "Node %4d at %8.3f, %8.3f, %8.3f: %s \n", pCNode->id(),
   //         adUserCoords[0], adUserCoords[1], adUserCoords[2], pcMessage );
   
}


double QuadCleanTool::edge_ratio( fmuEdge *pCEdge )
{
   
   double dSX, dSY, dSZ, dEX, dEY, dEZ, dDiffX, dDiffY, dDiffZ;
   double dX, dY, dZ, dLen, dSize, dRatio;
   nodClass *pCNode = pCEdge->start_node();
   dSX = pCNode->node_x();
   dSY = pCNode->node_y();
   dSZ = pCNode->node_z();
   pCNode = pCEdge->end_node();
   dEX = pCNode->node_x();
   dEY = pCNode->node_y();
   dEZ = pCNode->node_z();

   dX = ( dSX + dEX ) / 2.0;
   dY = ( dSY + dEY ) / 2.0;
   dZ = ( dSZ + dEZ ) / 2.0;
   dDiffX = dSX - dEX;
   dDiffY = dSY - dEY;
   dDiffZ = dSZ - dEZ;
   dLen = sqrt( dDiffX * dDiffX + dDiffY * dDiffY + dDiffZ * dDiffZ );

   dSize = get_size( dX, dY, dZ );

   dRatio = dLen / dSize;

   return dRatio;
   
}

// Erase the temp graphics
void QuadCleanTool::EraseElements()
{
     if ( !CMeshDebugOptions::debugCleaner()) 
        return;
     
MeshingDebugger::remove2dDisplay();
}


void QuadCleanTool::draw_mesh()
{
   
      //if( !lDrawSize )
      //{
      //CleanTool::draw_mesh();
      //return;
      //}

      //int  ii;
      //fmuEdge *pCEdge;
      //double dRatio;

   //Erase what is there
   //gg3pger();

      //draw_open();

      //for( ii = 0; ii < pCMaster->global_edges()->length(); ++ii )
      //{
      //pCEdge = (*pCMaster->global_edges())[ii];
      //if( pCEdge)
      //{
      //dRatio = edge_ratio( pCEdge );
      //if( dRatio < 0.1666666 )
      //pCEdge->draw_me( dp_kCLR_RED );
      //else if( dRatio < 0.20 )
      //pCEdge->draw_me( dp_kCLR_MAGNTA );
      //else if( dRatio < 0.25 )
      //pCEdge->draw_me( dp_kCLR_LTMGNT );
      //else if( dRatio < 0.33333333 )
      //pCEdge->draw_me( dp_kCLR_PINK );
      //else if( dRatio < 0.50 )
      //pCEdge->draw_me( dp_kCLR_CYAN );
      //else if( dRatio < 0.75 )
      //pCEdge->draw_me( dp_kCLR_ORANGE );
      //else if( dRatio < 0.999 )
      //pCEdge->draw_me( dp_kCLR_GOLD );
      //else if( dRatio < 1.001 )
      //pCEdge->draw_me( dp_kCLR_WHITE );
      //else if( dRatio < 1.33333 )
      //pCEdge->draw_me( dp_kCLR_BLUE );
      //else if( dRatio < 1.66666 )
      //pCEdge->draw_me( dp_kCLR_DKGRN );
      //else if( dRatio < 2.0 )
      //pCEdge->draw_me( dp_kCLR_DKOLIV );
      //else if( dRatio < 2.33333 )
      //pCEdge->draw_me( dp_kCLR_LTBLUE );
      //else if( dRatio < 2.66666 )
      //pCEdge->draw_me( dp_kCLR_GRBLUE );
      //else
      //pCEdge->draw_me( dp_kCLR_BLUE );
      //}
      //}
  
}

void QuadCleanTool::load_size()
{
   
   //This routine uses structures instead of classes for speed
   int ii, jj;
   fmuEdgeDynArray *CEdges;
   fmuEdge *pCEdge;
   CEdges = pCMaster->global_edges();
   nSize = 0;
   for( ii = 0; ii < CEdges->length(); ++ii )
   {
      pCEdge = (*CEdges)[ii];
      if( pCEdge )
         if( pCEdge->frozen() ) ++nSize;
   }

   zSize.adX = (double*) malloc( sizeof(double) * nSize );
   zSize.adY = (double*) malloc( sizeof(double) * nSize );
   zSize.adZ = (double*) malloc( sizeof(double) * nSize );
   zSize.adLen = (double*) malloc( sizeof(double) * nSize );

   jj = 0;
   for( ii = 0; ii < CEdges->length(); ++ii )
   {
      pCEdge = (*CEdges)[ii];
      if( !pCEdge ) continue;
      if( !pCEdge->frozen() ) continue;

      double dSX, dSY, dSZ, dEX, dEY, dEZ, dDiffX, dDiffY, dDiffZ;
      nodClass *pCNode = pCEdge->start_node();
      dSX = pCNode->node_x();
      dSY = pCNode->node_y();
      dSZ = pCNode->node_z();
      pCNode = pCEdge->end_node();
      dEX = pCNode->node_x();
      dEY = pCNode->node_y();
      dEZ = pCNode->node_z();

      zSize.adX[jj] = ( dSX + dEX ) / 2.0;
      zSize.adY[jj] = ( dSY + dEY ) / 2.0;
      zSize.adZ[jj] = ( dSZ + dEZ ) / 2.0;
      dDiffX = dSX - dEX;
      dDiffY = dSY - dEY;
      dDiffZ = dSZ - dEZ;
      zSize.adLen[jj] =
                  sqrt( dDiffX * dDiffX + dDiffY * dDiffY + dDiffZ * dDiffZ );
      ++jj;
   }
  
}


double QuadCleanTool::get_size( double dX,
                                          double dY,
                                          double dZ )
{
  
   int iRc;
   double dGrading,  dDefGrading = 1.0;

   dGrading = sam_ComputeGradingValue( this->pSMGradeFunc,dX, dY, dDefGrading, &iRc );

   if( iRc == 0 && dGrading < 1.0 ) return dGrading;

   int ii;
   double dDiffX, dDiffY, dDiffZ, dLen;
   double dSumNumerator = 0.0;
   double dSumDenominator = 0.0;
   for( ii = 0; ii < nSize; ++ii )
   {
      dDiffX = zSize.adX[ii] - dX;
      dDiffY = zSize.adY[ii] - dY;
      dDiffZ = zSize.adZ[ii] - dZ;
      dLen = sqrt( dDiffX * dDiffX + dDiffY * dDiffY + dDiffZ * dDiffZ );
      if( dLen < 0.0000001 ) return zSize.adLen[ii];
      dSumNumerator += zSize.adLen[ii] / dLen;
      dSumDenominator += 1.0 / dLen;
   }

   return dSumNumerator / dSumDenominator;
  
}


void QuadCleanTool::unload_size()
{
   
   free (zSize.adX);
   free (zSize.adY);
   free (zSize.adZ);
   free (zSize.adLen);
  
}


void QuadCleanTool::draw_hunting( nodClass* pCCurrentNode, int index )
{
   
   //Node passed in as the current position of CleaningList has just
   //been set to NULL.
      //int ii;
      //nodClass *pCNode;

      //draw_mesh();
      //draw_close();
      //pCCurrentNode->draw_me( dp_kCLR_YELLOW );

      //for( ii = index + 1; ii < CCleaningList.length(); ++ii )
      //{
      //pCNode = CCleaningList[ii];
      //if( !pCNode ) continue;
      //pCNode->draw_me( dp_kCLR_RED );
      //}
  
}


/*============================================================================================
test_quad_quality - Tests the quality of a quad element.
                    Given the 2D coordinates of the 4 nodes, this function
					computes the Jacobian Ratio = Min. Jacobian Det/Max Jac. Det.
					The min/max Jacobian determinants are from the ones
					computed at the four Gauss points. 

					A Ratio os 1.0 measures a perfect quad, a ratio of zero
					measures a quad that is squished to a straight line.
					Negative values would indicate degenerate quads like
					arrow-head quads etc.

                    Input - 2d coordinates of the 4 nodes in CCW order
					        dx1,dy1,dx2,dy2,dx3,dy3,dx4,dy4

				    Returns Jacobian Ratio

=============================================================================================*/
double QuadCleanTool::test_quad_quality (	double dx1, double dy1,
											double dx2, double dy2,
											double dx3, double dy3,
											double dx4, double dy4)
{
   
	double dQ=-1.0;
	double adQuad[8]={0.0};
	double adDetJac[4]  = {0.0, 0.0, 0.0, 0.0}, dMinJac =0.0, dMaxJac=0.0, dArea=0.0;

 
	// Get the coordinates of the would-be-quad
    adQuad[0] = dx1;
    adQuad[2] = dx2;
    adQuad[4] = dx3;
    adQuad[6] = dx4;

    adQuad[1] = dy1;
    adQuad[3] = dy2;
    adQuad[5] = dy3;
    adQuad[7] = dy4;

	// Compute the Jacobian Ratio
	sam_CalculateQuad4Jacobian (adQuad, 0, adDetJac, &dMinJac, &dMaxJac, &dArea, &dQ);
    return(dQ);
}

double QuadCleanTool::test_quad_quality (elmQClass *pCQuad)
{
    nodClass *apCNodes[4]= {0};
    double dQ=0.0, dx1=0.0, dy1=0.0, dx2=0.0, dy2=0.0, dx3=0.0, dy3=0.0, dx4=0.0, dy4=0.0;

    pCQuad->nodes(apCNodes);
    dx1 = apCNodes[0]->node_x();
    dy1 = apCNodes[0]->node_y();
    dx2 = apCNodes[1]->node_x();
    dy2 = apCNodes[1]->node_y();
    dx3 = apCNodes[2]->node_x();
    dy3 = apCNodes[2]->node_y();
    dx4 = apCNodes[3]->node_x();
    dy4 = apCNodes[3]->node_y();

    dQ = test_quad_quality (dx1, dy1, dx2, dy2, dx3, dy3, dx4, dy4);
    return(dQ);
}

