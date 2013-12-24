/*

*/

#include "cleaner/fmutclean.hxx"
#include "cleaner/fmuedgar.hxx"
#include "cleaner/elmtclass.hxx"
#include "cleaner/elmqclass.hxx"
#include "cleaner/nodclass.hxx"
#include "cleaner/fmuedge.hxx"
#include "cleaner/fmumaster.hxx"
#include "cleaner/fmu.x"
#include "quality_check/sam_elemquality.hxx"
#include <cmath>

//============================================================================================================================
//File name:			fmuqclean.cxx
//
//  mixed element mesh cleaning Member Function Descriptor Block -----
//
// TriaCleanTool::constructor
// class - QuadCleanTool
// 
//
// Description:
//
// constructor for the class.
//
// Access:
//
// TriaMeshClean()
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
// The constructor only initializes some basic values. While the
// QuadCleanTool cleans a quad only mesh, the TriaCleanTool does
// minimal work on a mixed element mesh and then isolates the
// trias behind frozen edges for input into QuadCleanTool.
//
///


const double fmu_kBIG_ANGLE_TOL = 2.7925268;   // 160 degrees

TriaCleanTool::~TriaCleanTool()
{
}

TriaCleanTool::TriaCleanTool()
{
}

int TriaCleanTool::CleanMesh( fmuMaster *pCMeshMaster,
                              fmuVector &CPlaneNormal )
{
    int iDebugFlag = 0;
    int ii, kk, nQuads;
    int iStatus = fmu_kCLEAN_NO_ACTION;
    int iCallStatus = 0;
    int iCaseCount;
    bool lMeshChanged;
    nodClass *apCNodes[3], *apCQNodes[3];
    fmuEdge *apCEdges[3];
    elmQClass *apCQuads[3];
    elmTClass *pCTria = 0;
    int iInitialLength;
    bool fAllowTrias = false;

    pCMaster = pCMeshMaster;

    CNormal = CPlaneNormal;

    smooth_mesh( pCMaster->global_interior() );

    if( lDrawCleaning )
    {
        //draw_mesh();
        //draw_close();
    }

    fAllowTrias = pCMaster->triqua_mesh_flag();
    CTriaList.clear();
    CTriaList = (*pCMaster->global_trias() );

    //First pass loop, check for best shape between a tria and a quad neighbor
    for( kk = 0; kk < CTriaList.length(); ++kk )
    {
        pCTria = CTriaList[kk];
        if( !pCTria )
            continue;
        if( pCTria->frozen() )
            continue;

        if( iDebugFlag )
        {
            //draw_mesh();
            //pCTria->draw_me( dp_kCLR_GREEN );
            //draw_close();
        }

        pCTria->nodes( apCNodes );

        nQuads = neighbor_quads( pCTria, apCNodes, apCQuads, apCEdges, apCQNodes );
        if(nQuads == 3 && fAllowTrias)
            continue;

        for( ii = 0; ii < nQuads; ++ii )
        {
            iCallStatus = reshape_quad( pCTria, apCQuads[ii], apCEdges[ii],
                apCQNodes[ii] );
            if( iCallStatus != fmu_kCLEAN_NO_ACTION )
                break;
        }
        if( iCallStatus == fmu_kCLEAN_MESH_CHANGED )
        {
            iStatus = true;
            if( lDrawCleaning )
            {
                //draw_mesh();
                //draw_close();
            }
        }
    }

    smooth_mesh( pCMaster->global_interior() );

    //Master loop. Execute this no more than 3 times.
    lMeshChanged = true;
    for( iCaseCount = 0; iCaseCount < 3 && lMeshChanged; ++iCaseCount )
    {
        lMeshChanged = false;
        CTriaList.clear();
        CTriaList = (*pCMaster->global_trias() );

        //Check every tria for possible combinations with neighbors
        iInitialLength = CTriaList.length();
        for( ii = 0; ii < iInitialLength; ++ii )
        {
            pCTria = CTriaList[ii];
            if( !pCTria )
                continue;
            if( pCTria->frozen() )
                continue;

            if( iDebugFlag )
            {
                //draw_mesh();
                //pCTria->draw_me( dp_kCLR_GREEN );
                //draw_close();
            }

            pCTria->nodes( apCNodes );
            nQuads = neighbor_quads( pCTria, apCNodes, apCQuads, apCEdges, apCQNodes );
            if(nQuads == 3 && fAllowTrias)
                continue;

            iCallStatus = clean_tria( pCTria );
            if( iCallStatus == fmu_kCLEAN_MESH_CHANGED )
            {
                lMeshChanged = true;
                if( lDrawCleaning )
                {
                    //draw_mesh();
                    //draw_close();
                }
            }
        }

        smooth_mesh( pCMaster->global_interior() );

        if( lMeshChanged )
            iStatus = fmu_kCLEAN_MESH_CHANGED;
    }

    isolate_trias();

    return iStatus;
}


int TriaCleanTool::clean_tria( elmTClass *pCTria )
{
    // Idiot-proofing - NM
    if (!pCTria)
        return fmu_kCLEAN_NO_ACTION;

    int ii, jj;
    int iCallStatus;
    int nQuads;
    nodClass *apCNodes[3], *apCQNodes[3];
    fmuEdge *pCTestEdge = 0, *apCEdges[3];
    elmQClass *apCQuads[3];

    pCTria->nodes( apCNodes );

    //Check for two trias before checking for tria along a quad
    for( ii = 0; ii < 3; ++ii )
    {
        jj = ii + 1;
        if( jj == 3 ) jj = 0;

        pCTestEdge = apCNodes[ii]->shared_edge( apCNodes[jj] );

        // If not found continue - NM
        if (!pCTestEdge)
            continue;

        if( pCTestEdge->frozen() ) continue;

        elmBase *pCOtherElem = 0;
        pCOtherElem = pCTestEdge->other_elem( pCTria );

        if (pCOtherElem)
        {
            if( pCOtherElem->type() == SAM_TRIA )
            {
                elmTClass *pCOtherTria = (elmTClass*) pCOtherElem;
                iCallStatus = double_trias( pCTria, pCOtherTria, pCTestEdge,
                    apCNodes[ii] );
                if( iCallStatus == fmu_kCLEAN_MESH_CHANGED )
                    return fmu_kCLEAN_MESH_CHANGED;
            }
        }
    }

    nQuads = neighbor_quads( pCTria, apCNodes, apCQuads, apCEdges, apCQNodes );

    //Check for all combinations of more than one tria & one quad after
    //checking for proper shape in a single tria-quad pair.
    for( ii = 0; ii < nQuads; ++ii )
    {
        iCallStatus = reshape_quad( pCTria, apCQuads[ii], apCEdges[ii],
            apCQNodes[ii] );
        if( iCallStatus == fmu_kCLEAN_MESH_CHANGED )
            return fmu_kCLEAN_MESH_CHANGED;
    }
    for( ii = 0; ii < nQuads; ++ii )
    {
        iCallStatus = tria_and_quad( pCTria, apCQuads[ii], apCEdges[ii],
            apCQNodes[ii] );
        if( iCallStatus == fmu_kCLEAN_MESH_CHANGED )
            return fmu_kCLEAN_MESH_CHANGED;
    }

    return fmu_kCLEAN_NO_ACTION;
}


int TriaCleanTool::neighbor_quads( elmTClass *pCTria, 
                                  nodClass *apCNodes[3],
                                  elmQClass *apCQuads[3],
                                  fmuEdge *apCEdges[3],
                                  nodClass *apCQNodes[3] )
{
    // Idiot-proofing - NM
    if (!pCTria)
        return fmu_kCLEAN_NO_ACTION;

    fmuEdge *pCTestEdge = 0;
    int ii, jj, nQuads = 0;

    for( ii = 0; ii < 3; ++ii )
    {
        jj = ii + 1;
        if( jj == 3 ) jj = 0;

        pCTestEdge = apCNodes[ii]->shared_edge( apCNodes[jj] );

        // If not found continue - NM
        if (!pCTestEdge)
            continue;

        if( pCTestEdge->frozen() )
            continue;

        elmBase *pCOtherElem = pCTestEdge->other_elem( pCTria );

        if (pCOtherElem)
        {
            if( pCOtherElem->type() == SAM_QUAD )
            {
                elmQClass *pCQuad = (elmQClass*) pCOtherElem;
                apCQuads[nQuads] = pCQuad;
                apCEdges[nQuads] = pCTestEdge;
                apCQNodes[nQuads] = apCNodes[ii];
                ++nQuads;
            }
        }
    }
    return nQuads;
}


int TriaCleanTool::double_trias( elmTClass *pCTria1,
                                elmTClass *pCTria2,
                                fmuEdge *pCEdge,
                                nodClass *pCNode  )
{
    double dQ = -1.0;
	
    // Idiot-proofing - NM
    if (!pCTria1 ||
        !pCTria2 ||
        !pCNode	||
        !pCEdge)
        return fmu_kCLEAN_NO_ACTION;

    if( pCTria2->frozen() )
        return fmu_kCLEAN_NO_ACTION;

    nodClassDynArray CSmoothList;
    nodClass *pCNextNode = 0, *pCOppNode = 0, *pCPrevNode = 0;
    pCNextNode = pCTria2->next_node( pCNode );
    if (pCNextNode == NULL)
        return fmu_kCLEAN_NO_ACTION;
    pCOppNode  = pCTria1->next_node( pCNode );
    if (pCOppNode == NULL)
        return fmu_kCLEAN_NO_ACTION;
    pCPrevNode = pCTria1->prev_node( pCNode );
    if (pCPrevNode == NULL)
        return fmu_kCLEAN_NO_ACTION;

    //If opposite nodes of a future quad angle are hardpoints, any flat
    //angles will not be improved with smoothing. It is best to attemp
    // two good triagles instead of making a bad quad.
    fmuVector CVecTP( pCNode, pCPrevNode );
    fmuVector CVecTN( pCNode, pCNextNode );
    fmuVector CVecNO( pCNextNode, pCOppNode );
    fmuVector CVecOP( pCOppNode, pCPrevNode );
    if( ( pCNode->hardpt() && pCOppNode->hardpt() ) ||
        ( pCPrevNode->hardpt() && pCNextNode->hardpt() ) )
    {
        if( CNormal.angle( CVecTN, CVecTP ) > fmu_kBIG_ANGLE_TOL )
            return two_better_trias( pCTria1, pCTria2, pCEdge, pCNode,
                                     pCNextNode, pCOppNode, pCPrevNode );

        if( CNormal.angle( CVecNO, -CVecTN ) > fmu_kBIG_ANGLE_TOL )
            return two_better_trias( pCTria1, pCTria2, pCEdge, pCNextNode,
                                     pCOppNode, pCPrevNode, pCNode );

        if( CNormal.angle( CVecOP, -CVecNO ) > fmu_kBIG_ANGLE_TOL )
            return two_better_trias( pCTria1, pCTria2, pCEdge, pCOppNode,
                                     pCPrevNode, pCNode, pCNextNode );

        if( CNormal.angle( CVecTP, CVecOP ) > fmu_kBIG_ANGLE_TOL )
            return two_better_trias( pCTria1, pCTria2, pCEdge, pCPrevNode,
                                     pCNode, pCNextNode, pCOppNode );
    }

    if( lDrawCleaning )
    {
        //draw_open();
        //pCTria1->draw_me( dp_kCLR_RED );
        //pCTria2->draw_me( dp_kCLR_YELLOW );
        //draw_close();
    }

	// Is the would-be-quad good ?
	// Compute Jac ratio of the would-be elements
	dQ = test_quad_quality( pCNode->node_x(),  pCNode->node_y(),
                            pCNextNode->node_x(), pCNextNode->node_y(),
                            pCOppNode->node_x(),  pCOppNode->node_y(),
                            pCPrevNode->node_x(), pCPrevNode->node_y());
    if (dQ < 0.1)
       return fmu_kCLEAN_NO_ACTION;
		
		

    delete_this_tria( pCTria1 );
    delete_this_tria( pCTria2 );

    delete_this_edge( pCEdge );

    elmQClass *pCQuad = 0;
    pCQuad = new elmQClass( pCNode, pCNextNode, pCOppNode, pCPrevNode, pCMaster );

    // Added abort condition - NM
    if (!pCQuad)
        return fmu_kCLEAN_ABORT;

    CSmoothList.append( pCNode );
    CSmoothList.append( pCNextNode );
    CSmoothList.append( pCOppNode );
    CSmoothList.append( pCPrevNode );
    smooth_mesh( &CSmoothList );
    if( lDrawCleaning )
    {
        //draw_mesh();
        //pCQuad->draw_me( dp_kCLR_RED );
        //draw_close();
    }

    return fmu_kCLEAN_MESH_CHANGED;
}


int TriaCleanTool::two_better_trias( elmTClass *pCTria1,
                                    elmTClass *pCTria2,
                                    fmuEdge *pCEdge, nodClass *pCNode,
                                    nodClass *pCNextNode,
                                    nodClass *pCOppNode,
                                    nodClass *pCPrevNode )
{
    // Idiot-proofing - NM
    if (!pCTria1 ||
        !pCTria2)
        return fmu_kCLEAN_NO_ACTION;

    if (pCNextNode == NULL ||
        pCPrevNode == NULL)
        return fmu_kCLEAN_NO_ACTION;

	nodClass *apCTriaNodes[3] = {0};
	double dTriaAreaOrig = 0.0, dTriaArea = 0.0;
	double dAreaTol=kdfToler;

    //If there is an edge between this and the opposite node, do nothing.
    //Otherwise delete the two trias and make two more trias.

    if(!pCNode || !pCOppNode )
        return fmu_kCLEAN_NO_ACTION;

    fmuEdge *pCCrossEdge = pCNode->shared_edge( pCOppNode );
    if( pCCrossEdge )
        return fmu_kCLEAN_NO_ACTION;

	
	// Check quality of would-be trias 
	// Compute area of input Trias
	pCTria1->nodes( apCTriaNodes );
	dTriaAreaOrig = get_tria_area (apCTriaNodes[0]->node_x(),  apCTriaNodes[0]->node_y(),
                                   apCTriaNodes[1]->node_x(),  apCTriaNodes[1]->node_y(),
                                   apCTriaNodes[2]->node_x(),  apCTriaNodes[2]->node_y());
	dTriaArea     = get_tria_area (pCNode->node_x(), pCNode->node_y(), 
		                           pCNextNode->node_x(), pCNextNode->node_y(),
								   pCOppNode->node_x(), pCOppNode->node_y());

   
	// The original tria and the would-be-tria must have the same sense. Their areas must not change signs.
	// Furthermore the tria must not be collapsed.
	if (fabs(dTriaArea) < dAreaTol)
		return fmu_kCLEAN_NO_ACTION;
	else if ((dTriaAreaOrig/dTriaArea) < dAreaTol)
		return fmu_kCLEAN_NO_ACTION;

	pCTria2->nodes( apCTriaNodes );
	dTriaAreaOrig = get_tria_area (apCTriaNodes[0]->node_x(),  apCTriaNodes[0]->node_y(),
                                   apCTriaNodes[1]->node_x(),  apCTriaNodes[1]->node_y(),
                                   apCTriaNodes[2]->node_x(),  apCTriaNodes[2]->node_y());

	dTriaArea     = get_tria_area (pCNode->node_x(), pCNode->node_y(), 
		                           pCOppNode->node_x(), pCOppNode->node_y(),
								   pCPrevNode->node_x(), pCPrevNode->node_y());

   
	// The original tria and the would-be-tria must have the same sense. Their areas must not change signs.
	// Furthermore the tria must not be collapsed.
	if (fabs(dTriaArea) < dAreaTol)
		return fmu_kCLEAN_NO_ACTION;
	else if ((dTriaAreaOrig/dTriaArea) < dAreaTol)
		return fmu_kCLEAN_NO_ACTION;

    if( lDrawCleaning )
    {
        //draw_open();
        //pCTria1->draw_me( dp_kCLR_RED );
        //pCTria2->draw_me( dp_kCLR_ORANGE );
        //draw_close();
    }

    delete_this_tria( pCTria1 );
    delete_this_tria( pCTria2 );

    delete_this_edge( pCEdge );

    pCTria1 = new elmTClass( pCNode, pCNextNode, pCOppNode, pCMaster );
    pCTria2 = new elmTClass( pCNode, pCOppNode, pCPrevNode, pCMaster );

    if (!pCTria1 ||
        !pCTria2)
        return fmu_kCLEAN_ABORT;

    nodClassDynArray CSmoothList;
    CSmoothList.append( pCNode );
    CSmoothList.append( pCNextNode );
    CSmoothList.append( pCOppNode );
    CSmoothList.append( pCPrevNode );
    smooth_mesh( &CSmoothList );

    if( lDrawCleaning )
    {
        //draw_mesh();
        //pCTria1->draw_me( dp_kCLR_RED );
        //pCTria2->draw_me( dp_kCLR_ORANGE );
        //draw_close();
    }

    return fmu_kCLEAN_MESH_CHANGED;
}


int TriaCleanTool::tria_and_quad( elmTClass *pCTria,
                                 elmQClass *pCQuad,
                                 fmuEdge *pCEdge,
                                 nodClass *pCNode  )
{
    // Check for NULL objects - NM
    if (pCTria == NULL ||
        pCQuad == NULL ||
        pCEdge == NULL ||
        pCNode == NULL)
        return fmu_kCLEAN_NO_ACTION;

    int iCallStatus;
    //See if this tria with this quad fits into one of the standard patterns.
    //Start by identifying nodes around the tria and the quad. Work outward
    //from there.

    nodClass *pCEdgeNode = 0;
    pCEdgeNode = pCEdge->other_node( pCNode );
    if (!pCEdgeNode)
return fmu_kCLEAN_NO_ACTION;

    nodClass *pCTriaPrevNode = 0;
    pCTriaPrevNode = pCTria->prev_node( pCNode );
    if (!pCTriaPrevNode)
return fmu_kCLEAN_NO_ACTION;

    nodClass *pCQuadNextNode = 0;
    pCQuadNextNode = pCQuad->next_node( pCNode );
    if (!pCQuadNextNode)
return fmu_kCLEAN_NO_ACTION;

    nodClass *pCQuadOppNode = 0;
    pCQuadOppNode = pCQuad->opposite_node( pCNode );
    if (!pCQuadOppNode)
return fmu_kCLEAN_NO_ACTION;

    //Check if there is another tria that is vertex-to-vertex with the
    //original.
    elmBase *pCPrevElem = NULL;
    fmuEdge *pCPrevEdge = pCEdgeNode->shared_edge( pCQuadOppNode );
    if( !pCPrevEdge->frozen() )
    {
        pCPrevElem = pCPrevEdge->other_elem( pCQuad );
        if( pCPrevElem->type() == SAM_TRIA )
        {
            elmTClass *pCPrevTria = (elmTClass*) pCPrevElem;
            iCallStatus = two_of_each( pCTria, pCQuad, pCEdge, pCPrevTria, 
                                        pCPrevEdge, pCNode, pCEdgeNode,
                                        pCTriaPrevNode, pCQuadNextNode,
                                        pCQuadOppNode );
            if( iCallStatus == fmu_kCLEAN_MESH_CHANGED )
                return fmu_kCLEAN_MESH_CHANGED;
        }
    }

    //Check if there is another tria on the opposite side of the quad
    //from the first tria.
    elmBase *pCOppElem = NULL;
    fmuEdge *pCOppEdge = pCQuadNextNode->shared_edge( pCQuadOppNode );
    if( !pCOppEdge->frozen() )
    {
        pCOppElem = pCOppEdge->other_elem( pCQuad );
        if( pCOppElem->type() == SAM_TRIA )
        {
            elmTClass *pCOppTria = (elmTClass*) pCOppElem;
            iCallStatus = tria_quad_tria( pCTria, pCQuad, pCEdge, pCOppTria,
                                            pCOppEdge, pCNode, pCEdgeNode,
                                            pCTriaPrevNode, pCQuadNextNode,
                                            pCQuadOppNode );
            if( iCallStatus == fmu_kCLEAN_MESH_CHANGED )
                return fmu_kCLEAN_MESH_CHANGED;
        }
    }

    //One of the secondary elements could be a quad.
    if( pCPrevElem && pCPrevElem->type() == SAM_QUAD )
    {
        elmQClass *pCPrevQuad = (elmQClass*) pCPrevElem;
        iCallStatus = tria_quad_quad( pCTria, pCQuad, pCEdge, pCPrevQuad,
                                        pCPrevEdge, pCNode, pCEdgeNode,
                                        pCTriaPrevNode, pCQuadNextNode,
                                        pCQuadOppNode );
        if( iCallStatus == fmu_kCLEAN_MESH_CHANGED )
            return fmu_kCLEAN_MESH_CHANGED;
    }
    return fmu_kCLEAN_NO_ACTION;
}


int TriaCleanTool::two_of_each( elmTClass *pCTria, elmQClass *pCQuad,
                               fmuEdge *pCEdge, elmTClass *pCAltTria,
                               fmuEdge *pCPrevEdge, nodClass *pCNode,
                               nodClass *pCEdgeNode,
                               nodClass *pCTriaPrevNode,
                               nodClass *pCQuadNextNode,
                               nodClass *pCQuadOppNode )
{
    // Check for NULL objects - NM
    if (!pCAltTria)
        return fmu_kCLEAN_NO_ACTION;
    if( pCAltTria->frozen() )
        return fmu_kCLEAN_NO_ACTION;

    if (!pCTriaPrevNode)
        return fmu_kCLEAN_NO_ACTION;

    //There must be a second quad that shares an edge with the two
    //triangles. If not, this situation is appropriate for a different
    //cleaning routine.
    fmuEdge *pCEdge1 = pCEdgeNode->shared_edge( pCTriaPrevNode );
    if (!pCEdge1)
        return fmu_kCLEAN_NO_ACTION; // Check for valid edge - NM

    if( pCEdge1->frozen() )
    {
        return nose_to_nose( pCTria, pCQuad, pCEdge, pCAltTria, 
                                pCPrevEdge, pCNode, pCEdgeNode,
                                pCTriaPrevNode, pCQuadNextNode, pCQuadOppNode );
    }

    elmBase *pCOuterElem = pCEdge1->other_elem( pCTria );
    if (!pCOuterElem)
        return fmu_kCLEAN_NO_ACTION;

    if( pCOuterElem->type() != SAM_QUAD )
    {
        return nose_to_nose( pCTria, pCQuad, pCEdge, pCAltTria, 
            pCPrevEdge, pCNode, pCEdgeNode,
            pCTriaPrevNode, pCQuadNextNode, pCQuadOppNode );
    }
    elmQClass *pCOuterQuad = dynamic_cast< elmQClass* >( pCOuterElem );

    if (!pCEdgeNode)
        return fmu_kCLEAN_NO_ACTION;

    nodClass *pCSharedNode = pCOuterQuad->next_node( pCEdgeNode );
    if( pCSharedNode != pCAltTria->prev_node( pCEdgeNode ) )
    {
        return nose_to_nose( pCTria, pCQuad, pCEdge, pCAltTria, 
                            pCPrevEdge, pCNode, pCEdgeNode,
                            pCTriaPrevNode, pCQuadNextNode, pCQuadOppNode );
    }

    if (!pCSharedNode)
        return fmu_kCLEAN_NO_ACTION;
    fmuEdge *pCEdge2 = pCEdgeNode->shared_edge( pCSharedNode );
    if (!pCEdge2)
        return fmu_kCLEAN_NO_ACTION;

    if( pCEdge2->frozen() )
    {
        return nose_to_nose( pCTria, pCQuad, pCEdge, pCAltTria, 
            pCPrevEdge, pCNode, pCEdgeNode,
            pCTriaPrevNode, pCQuadNextNode, pCQuadOppNode );
    }

    if (!pCEdgeNode)
        return fmu_kCLEAN_NO_ACTION;
    if( pCEdgeNode->hardpt() )
        return fmu_kCLEAN_NO_ACTION;

    nodClass *pCOuterNode = pCOuterQuad->opposite_node( pCEdgeNode );
    if (!pCOuterNode)
        return fmu_kCLEAN_NO_ACTION;

    if( lDrawCleaning )
    {
        //draw_open();
        //pCQuad->draw_me( dp_kCLR_GREEN );
        //pCOuterQuad->draw_me( dp_kCLR_BLUE );
        //pCTria->draw_me( dp_kCLR_YELLOW );
        //pCAltTria->draw_me( dp_kCLR_RED );
        //draw_close();
    }

    delete_this_tria( pCTria );
    delete_this_tria( pCAltTria );
    delete_this_quad( pCQuad );
    delete_this_quad( pCOuterQuad );
    delete_this_edge( pCEdge );
    delete_this_edge( pCPrevEdge );
    delete_this_edge( pCEdge1 );
    delete_this_edge( pCEdge2 );
    delete_this_node( pCEdgeNode );

    nodClassDynArray CRingNodes;
    CRingNodes[0] = pCNode;
    CRingNodes[1] = pCQuadNextNode;
    CRingNodes[2] = pCQuadOppNode;
    CRingNodes[3] = pCSharedNode;
    CRingNodes[4] = pCOuterNode;
    CRingNodes[5] = pCTriaPrevNode;

    fill_6( CRingNodes );

    return fmu_kCLEAN_MESH_CHANGED;
}


int TriaCleanTool::nose_to_nose( elmTClass *pCTria, elmQClass *pCQuad,
                                fmuEdge *pCEdge, elmTClass *pCAltTria,
                                fmuEdge *pCPrevEdge, nodClass *pCNode,
                                nodClass *pCEdgeNode,
                                nodClass *pCTriaPrevNode,
                                nodClass *pCQuadNextNode,
                                nodClass *pCQuadOppNode )
{
    double dQ = 1.0;
	nodClassDynArray CSmoothList;
	elmQClass *pCQuadA=0, *pCQuadB=0;
	
    if (!pCTria)
        return fmu_kCLEAN_NO_ACTION;

    if (!pCEdgeNode)
        return fmu_kCLEAN_NO_ACTION;

    nodClass *pCOuterNode = pCAltTria->prev_node( pCEdgeNode );
    if (!pCOuterNode)
        return fmu_kCLEAN_NO_ACTION;

    if (!pCNode)
        return fmu_kCLEAN_NO_ACTION;

    //Don't proceed with the change if the new quads will be poorly shaped.
    if( pCEdgeNode->hardpt() && pCNode->hardpt() && pCOuterNode->hardpt() )
    {
        fmuVector CVecENN( pCEdgeNode, pCNode );
        fmuVector CVecENON( pCEdgeNode, pCOuterNode );
        if( CNormal.angle( CVecENN, CVecENON ) > PI )
            return fmu_kCLEAN_NO_ACTION;
    }

    if (!pCTriaPrevNode)
        return fmu_kCLEAN_NO_ACTION;

    if (!pCQuadOppNode)
        return fmu_kCLEAN_NO_ACTION;

    if( pCEdgeNode->hardpt() && 
        pCTriaPrevNode->hardpt() &&
        pCQuadOppNode->hardpt() )
    {
        fmuVector CVecENTP( pCEdgeNode, pCTriaPrevNode );
        fmuVector CVecENQO( pCEdgeNode, pCQuadOppNode );
        if( CNormal.angle( CVecENTP, CVecENQO ) > PI )
            return fmu_kCLEAN_NO_ACTION;
    }

    if( lDrawCleaning )
    {
        //draw_open();
        //pCTria->draw_me( dp_kCLR_GOLD );
        //pCAltTria->draw_me( dp_kCLR_YELLOW );
        //pCQuad->draw_me( dp_kCLR_BLUE );
        //draw_close();
    }

    nodClass *pCNewNode = one_middle_node( pCEdgeNode, pCQuadNextNode );

	// Are these 3 would-be-quads good ? If not skip it.
	// Compute Jac ratio of the would-be elements
	dQ = test_quad_quality( pCEdgeNode->node_x(),		pCEdgeNode->node_y(),
                            pCTriaPrevNode->node_x(),	pCTriaPrevNode->node_y(),
                            pCNode->node_x(),			pCNode->node_y(),
                            pCNewNode->node_x(),		pCNewNode->node_y());
    if (dQ < 0.1)
		goto Error;

	dQ = test_quad_quality( pCNewNode->node_x(),		pCNewNode->node_y(),
                            pCNode->node_x(),			pCNode->node_y(),
                            pCQuadNextNode->node_x(),	pCQuadNextNode->node_y(),
                            pCQuadOppNode->node_x(),	pCQuadOppNode->node_y());
    if (dQ < 0.1)
       goto Error;

	dQ = test_quad_quality( pCNewNode->node_x(),		pCNewNode->node_y(),
                            pCQuadOppNode->node_x(),	pCQuadOppNode->node_y(),
                            pCOuterNode->node_x(),		pCOuterNode->node_y(),
                            pCEdgeNode->node_x(),		pCEdgeNode->node_y());
    if (dQ < 0.1)
       goto Error;


    delete_this_tria( pCTria );
    delete_this_tria( pCAltTria );
    delete_this_quad( pCQuad );
    delete_this_edge( pCEdge );
    delete_this_edge( pCPrevEdge );

    

    pCQuad = new elmQClass( pCEdgeNode, pCTriaPrevNode, pCNode, pCNewNode, pCMaster );
    pCQuadA = new elmQClass( pCNewNode, pCNode, pCQuadNextNode, pCQuadOppNode, pCMaster );
    pCQuadB = new elmQClass( pCNewNode, pCQuadOppNode, pCOuterNode, pCEdgeNode, pCMaster );

    
    CSmoothList.append( pCNode );
    CSmoothList.append( pCEdgeNode );
    CSmoothList.append( pCTriaPrevNode );
    CSmoothList.append( pCQuadNextNode );
    CSmoothList.append( pCQuadOppNode );
    CSmoothList.append( pCOuterNode );
    smooth_mesh( &CSmoothList );

    if( lDrawCleaning )
    {
        //draw_mesh();
        //pCQuad->draw_me( dp_kCLR_GOLD );
        //pCQuadA->draw_me( dp_kCLR_BLUE );
        //pCQuadB->draw_me( dp_kCLR_YELLOW );
        //draw_close();
    }

    return fmu_kCLEAN_MESH_CHANGED;

Error:
	if (pCNewNode)
		delete_this_node(pCNewNode); 
	return fmu_kCLEAN_NO_ACTION;

}


int TriaCleanTool::tria_quad_tria( elmTClass *pCTria, elmQClass *pCQuad,
                                  fmuEdge *pCEdge, elmTClass *pCAltTria,
                                  fmuEdge *pCOppEdge, nodClass *pCNode,
                                  nodClass *pCEdgeNode,
                                  nodClass *pCTriaPrevNode,
                                  nodClass *pCQuadNextNode,
                                  nodClass *pCQuadOppNode )
{
    if( pCTria == NULL )
        return fmu_kCLEAN_NO_ACTION;

    if( pCAltTria == NULL )
        return fmu_kCLEAN_NO_ACTION;

    if( pCAltTria->frozen() )
        return fmu_kCLEAN_NO_ACTION;

    if( pCQuadNextNode == NULL )
        return fmu_kCLEAN_NO_ACTION;

    nodClass *pCOuterNode = pCAltTria->next_node( pCQuadNextNode );
    if( pCOuterNode == NULL )
        return fmu_kCLEAN_NO_ACTION;

    nodClassDynArray CRingNodes;
    CRingNodes[0] = pCNode;
    CRingNodes[1] = pCQuadNextNode;
    CRingNodes[2] = pCOuterNode;
    CRingNodes[3] = pCQuadOppNode;
    CRingNodes[4] = pCEdgeNode;
    CRingNodes[5] = pCTriaPrevNode;

    int ii, iRet = 0;
    int nHardpt = 0;

    //If there are lots of hardpoints or boundary nodes here, the two
    //trias and the single quad are probabaly the best shape possible.
    for( ii = 0; ii < 6; ++ii)
    {
        if( CRingNodes[ii] == NULL )
            return fmu_kCLEAN_NO_ACTION;
        if( CRingNodes[ii]->hardpt() )
            ++nHardpt;
    }

    if( nHardpt > 3 )
        return fmu_kCLEAN_NO_ACTION;

    if( lDrawCleaning )
    {
        //draw_open();
        //pCTria->draw_me( dp_kCLR_MAGNTA );
        //pCAltTria->draw_me( dp_kCLR_PINK );
        //pCQuad->draw_me( dp_kCLR_BLUE );
        //draw_close();
    }

    iRet = check_fill_6( CRingNodes );
    if (iRet != 0)
        return fmu_kCLEAN_NO_ACTION;

    delete_this_tria( pCTria );
    delete_this_tria( pCAltTria );

    if( pCQuad == NULL )
        return fmu_kCLEAN_NO_ACTION;

    delete_this_quad( pCQuad );

    if( pCEdge == NULL )
        return fmu_kCLEAN_NO_ACTION;

    delete_this_edge( pCEdge );

    if( pCOppEdge == NULL )
        return fmu_kCLEAN_NO_ACTION;

    delete_this_edge( pCOppEdge );

    fill_6( CRingNodes );

    return fmu_kCLEAN_MESH_CHANGED;
}


int TriaCleanTool::tria_quad_quad( elmTClass *pCTria, elmQClass *pCQuad,
                                  fmuEdge *pCEdge, elmQClass *pCAltQuad,
                                  fmuEdge *pCPrevEdge, nodClass *pCNode,
                                  nodClass *pCEdgeNode,
                                  nodClass *pCTriaPrevNode,
                                  nodClass *pCQuadNextNode,
                                  nodClass *pCQuadOppNode )
{

    int iCallStatus;
    nodClassDynArray CRingNodes;

    nodClass *pCOuterNode = pCAltQuad->opposite_node( pCEdgeNode );
    nodClass *pCAltNode = pCAltQuad->prev_node( pCEdgeNode );
    fmuEdge *pCTriaEdge = pCTriaPrevNode->shared_edge( pCEdgeNode );
    //If things are not quite to your liking here, try an alternative.
    if( pCAltNode != pCTriaPrevNode )
    {
        return tria_hourglass( pCTria, pCQuad, pCEdge, pCAltQuad,
            pCPrevEdge, pCNode, pCEdgeNode,
            pCTriaPrevNode, pCQuadNextNode,
            pCQuadOppNode, pCOuterNode, pCAltNode );
    }

    // Check for NULL pointer - NM
    if (!pCEdgeNode) return fmu_kCLEAN_NO_ACTION;
    if( pCEdgeNode->hardpt() ) return fmu_kCLEAN_NO_ACTION;

    if (!pCTriaEdge) return fmu_kCLEAN_NO_ACTION;
    // Check for NULL pointer - NM
    if( pCTriaEdge->frozen() ) return fmu_kCLEAN_NO_ACTION;

    if( lDrawCleaning )
    {
        //draw_open();
        //pCTria->draw_me( dp_kCLR_DKGRN );
        //pCAltQuad->draw_me( dp_kCLR_PINK );
        //pCQuad->draw_me( dp_kCLR_GRBLUE );
        //draw_close();
    }

    delete_this_tria( pCTria );
    delete_this_quad( pCQuad );
    delete_this_quad( pCAltQuad );
    delete_this_edge( pCEdge );
    delete_this_edge( pCPrevEdge );
    delete_this_edge( pCTriaEdge );
    delete_this_node( pCEdgeNode );

    CRingNodes[0] = pCNode;
    CRingNodes[1] = pCQuadNextNode;
    CRingNodes[2] = pCQuadOppNode;
    CRingNodes[3] = pCOuterNode;
    CRingNodes[4] = pCTriaPrevNode;
    iCallStatus = fill_5( CRingNodes, NULL, NULL, NULL );

    return iCallStatus;

}


int TriaCleanTool::tria_hourglass( elmTClass *pCTria, elmQClass *pCQuad,
                                  fmuEdge *pCEdge, elmQClass *pCAltQuad,
                                  fmuEdge *pCPrevEdge, nodClass *pCNode,
                                  nodClass *pCEdgeNode,
                                  nodClass *pCTriaPrevNode,
                                  nodClass *pCQuadNextNode,
                                  nodClass *pCQuadOppNode,
                                  nodClass *pCOuterNode,
                                  nodClass *pCAltNode )
{
    nodClassDynArray CSmoothList;
    elmQClass *pCQuadA = 0, *pCQuadB = 0;
    fmuEdge *pCAltEdge = 0;

    pCAltEdge = pCEdgeNode->shared_edge( pCAltNode );
    // Check for NULL pointer - NM
    if (!pCAltEdge) return fmu_kCLEAN_NO_ACTION;
    if( pCAltEdge->frozen() ) return fmu_kCLEAN_NO_ACTION;

    elmBase *pCLastElem = 0;
    pCLastElem = pCAltEdge->other_elem( pCAltQuad );
    // Check for NULL pointer - NM
    if (!pCLastElem) return fmu_kCLEAN_NO_ACTION;
    if( pCLastElem->type() != SAM_TRIA ) return fmu_kCLEAN_NO_ACTION;

    elmTClass *pCAltTria = 0;
    pCAltTria = (elmTClass*) pCLastElem;
    // Check for NULL pointer - NM
    if (!pCAltTria) return fmu_kCLEAN_NO_ACTION;
    if( pCAltTria->frozen() ) return fmu_kCLEAN_NO_ACTION;

    nodClass *pCAltTriaNode = 0;
    pCAltTriaNode = pCAltTria->next_node( pCAltNode );
    if( pCTriaPrevNode == pCAltTriaNode ) return fmu_kCLEAN_NO_ACTION;

    if( lDrawCleaning )
    {
        //draw_open();
        //pCTria->draw_me( dp_kCLR_RED );
        //pCQuad->draw_me( dp_kCLR_LTBLUE );
        //pCAltQuad->draw_me( dp_kCLR_PINK );
        //pCAltTria->draw_me( dp_kCLR_CYAN );
        //draw_close();
    }

	double dQ = -1.0;
	nodClass *pCNewNode = 0;
    pCNewNode = one_middle_node( pCEdgeNode, pCQuadOppNode );

	// Check the quality of all would-be-quads. If any one has a poor Jac Ratio
	// we skip this template.
	dQ = test_quad_quality( pCEdgeNode->node_x(),		pCEdgeNode->node_y(),
                            pCTriaPrevNode->node_x(),	pCTriaPrevNode->node_y(),
                            pCNode->node_x(),			pCNode->node_y(),
                            pCNewNode->node_x(),		pCNewNode->node_y());
    if (dQ < 0.1)
       goto Error;

	dQ = test_quad_quality( pCNewNode->node_x(),		pCNewNode->node_y(),
                            pCNode->node_x(),			pCNode->node_y(),
							pCQuadNextNode->node_x(),	pCQuadNextNode->node_y(),
                            pCQuadOppNode->node_x(),	pCQuadOppNode->node_y());
    if (dQ < 0.1)
       goto Error;

	dQ = test_quad_quality( pCAltTriaNode->node_x(),	pCAltTriaNode->node_y(),
                            pCEdgeNode->node_x(),		pCEdgeNode->node_y(),
                            pCNewNode->node_x(),		pCNewNode->node_y(),
                            pCAltNode->node_x(),		pCAltNode->node_y());
    if (dQ < 0.1)
       goto Error;

	dQ = test_quad_quality( pCAltNode->node_x(),		pCAltNode->node_y(),
							pCNewNode->node_x(),		pCNewNode->node_y(),
							pCQuadOppNode->node_x(),	pCQuadOppNode->node_y(),
                            pCOuterNode->node_x(),		pCOuterNode->node_y()
                          );
    if (dQ < 0.1)
       goto Error;

    delete_this_tria( pCTria );
    delete_this_tria( pCAltTria );
    delete_this_quad( pCQuad );
    delete_this_quad( pCAltQuad );
    delete_this_edge( pCEdge );
    delete_this_edge( pCPrevEdge );
    delete_this_edge( pCAltEdge );


    pCQuadA = new elmQClass( pCEdgeNode, pCTriaPrevNode, pCNode, pCNewNode,
        pCMaster );
    pCQuad = new elmQClass( pCNewNode, pCNode, pCQuadNextNode, pCQuadOppNode,
        pCMaster );
    pCAltQuad = new elmQClass( pCAltTriaNode, pCEdgeNode, pCNewNode,
        pCAltNode, pCMaster );
    pCQuadB = new elmQClass( pCAltNode, pCNewNode, pCQuadOppNode,
        pCOuterNode, pCMaster );

    CSmoothList.append( pCNode );
    CSmoothList.append( pCEdgeNode );
    CSmoothList.append( pCTriaPrevNode );
    CSmoothList.append( pCQuadNextNode );
    CSmoothList.append( pCQuadOppNode );
    CSmoothList.append( pCOuterNode );
    CSmoothList.append( pCAltNode );
    CSmoothList.append( pCAltTriaNode );
    smooth_mesh( &CSmoothList );
    if( lDrawCleaning )
    {
        //draw_mesh();
        //pCQuadA->draw_me( dp_kCLR_RED );
        //pCQuad->draw_me( dp_kCLR_LTBLUE );
        //pCAltQuad->draw_me( dp_kCLR_PINK );
        //pCQuadB->draw_me( dp_kCLR_CYAN );
        //draw_close();
    }

    return fmu_kCLEAN_MESH_CHANGED;

Error:
	if (pCNewNode)
		delete_this_node(pCNewNode); 
	return fmu_kCLEAN_NO_ACTION;
}


int TriaCleanTool::reshape_quad( elmTClass *pCTria, elmQClass *pCQuad,
                                fmuEdge *pCEdge, nodClass *pCNode )
{
    int iCallStatus;
    nodClassDynArray CRingNodes;

    if (!pCNode)
        return fmu_kCLEAN_NO_ACTION; //sanity check - NM
    CRingNodes[0] = pCNode;
    if (!pCQuad)
        return fmu_kCLEAN_NO_ACTION; //sanity check - NM
    CRingNodes[1] = pCQuad->next_node( pCNode );
    CRingNodes[2] = pCQuad->opposite_node( pCNode );
    if (!pCEdge)
        return fmu_kCLEAN_NO_ACTION; //sanity check - NM
    CRingNodes[3] = pCEdge->other_node( pCNode );
    if (!pCTria)
        return fmu_kCLEAN_NO_ACTION; //sanity check - NM
    CRingNodes[4] = pCTria->prev_node( pCNode );
    iCallStatus = fill_5( CRingNodes, pCTria, pCQuad, pCEdge );

    return iCallStatus;
}


void TriaCleanTool::fill_6( nodClassDynArray &CRingNodes )
{
    //Cut across the loop from the node with the biggest angle making 2 quads.
    //This may not be the best pattern, but quad cleaning will improve it.
    int ii;
    double dAngle, dBigAngle = 0.0, dQ=1.0;
    nodClass *pCFlatNode = 0, *pCNode = 0;
    elmQClass *pCQuadA = 0, *pCQuadB = 0;
    for( ii = 0; ii < 6; ++ii )
    {
        pCNode = CRingNodes.get();
        if (!pCNode) return; // sanity check - NM

        fmuVector CToPrev( pCNode, CRingNodes.prev() );
        fmuVector CToNext( pCNode, CRingNodes.next() );
        dAngle = CNormal.angle( CToNext, CToPrev );
        if( dAngle > dBigAngle )
        {
            dBigAngle = dAngle;
            pCFlatNode = pCNode;
        }
    }

    CRingNodes.find( pCFlatNode );

    pCQuadA = new elmQClass( pCFlatNode, CRingNodes.next(), CRingNodes.next(2),
        CRingNodes.next(3), pCMaster );

    pCQuadB = new elmQClass( pCFlatNode, CRingNodes.next(3), CRingNodes.next(4),
        CRingNodes.next(5), pCMaster );

    smooth_mesh( &CRingNodes );
    if( lDrawCleaning )
    {
        //draw_mesh();
        //pCQuadA->draw_me( dp_kCLR_RED );
        //pCQuadB->draw_me( dp_kCLR_MAGNTA );
        //draw_close();
    }
}

// Check if the fill_6 template will work. 
// Mostly perform quality check on the would-be-quads
int TriaCleanTool::check_fill_6( nodClassDynArray &CRingNodes )
{
    //Cut across the loop from the node with the biggest angle making 2 quads.
    //This may not be the best pattern, but quad cleaning will improve it.
    int ii;
    double dAngle, dBigAngle = 0.0, dQ=1.0;
    nodClass *pCFlatNode = 0, *pCNode = 0;
    for( ii = 0; ii < 6; ++ii )
    {
        pCNode = CRingNodes.get();
        if (!pCNode) return(1); // sanity check - NM

        fmuVector CToPrev( pCNode, CRingNodes.prev() );
        fmuVector CToNext( pCNode, CRingNodes.next() );
        dAngle = CNormal.angle( CToNext, CToPrev );
        if( dAngle > dBigAngle )
        {
            dBigAngle = dAngle;
            pCFlatNode = pCNode;
        }
    }

    CRingNodes.find( pCFlatNode );

    // Compute Jac ratio of the would-be elements
    dQ = test_quad_quality(pCFlatNode->node_x(), pCFlatNode->node_y(),
        CRingNodes.next()->node_x(), CRingNodes.next()->node_y(),
        CRingNodes.next(2)->node_x(), CRingNodes.next(2)->node_y(),
        CRingNodes.next(3)->node_x(), CRingNodes.next(3)->node_y());

    if (dQ < 0.1)
        return(1);

    // Compute Jac ratio of the would-be elements
    dQ = test_quad_quality(pCFlatNode->node_x(), pCFlatNode->node_y(),
        CRingNodes.next(3)->node_x(), CRingNodes.next(3)->node_y(),
        CRingNodes.next(4)->node_x(), CRingNodes.next(4)->node_y(),
        CRingNodes.next(5)->node_x(), CRingNodes.next(5)->node_y());
    if (dQ < 0.1)
        return(1);

    return(0);
}


int TriaCleanTool::fill_5( nodClassDynArray &CRingNodes,
                          elmTClass *pCTria, elmQClass *pCQuad,
                          fmuEdge *pCEdge )
{
    int ii, jj, kk, mm, nn, ita, itb, itc;
    int iBestIndex;
    double adAngleDiffs[30];
    double dBestSum, dSumDiffs;
    double dFullAngle, dLeftAngle, dRightAngle;
    double kTRIA_ANGLE = PI / 3.0;
    double kRIGHT_ANGLE = PI / 2.0;
    double kQUAD_FACTOR = 0.52359878;
    double kTRIA_FACTOR = 0.34906585;
	double dQ = 0.0, dTriaAreaOrig = 0.0, dTriaArea=0.0, dAreaTol=kdfToler;

	nodClass *pCNode = 0, *apCTriaNodes[3]={0};
    //Align internal and external iterators
    CRingNodes.find( CRingNodes[0] );
    //adAngleDiffs has difference between angle at 90 degrees. 
    //[0-4] is for at node, [5-9]is for from prev to after next
    //[10-14] is for from before prev to next.
    //[15-19] is angle at node diff from 60 degrees,
    //[20-24] is for prev to before prev, [25-29] is after next to next.
    for( ii = 0; ii < 5; ++ii )
    {
        pCNode = CRingNodes.get();
        if (!pCNode) return fmu_kCLEAN_NO_ACTION;

        fmuVector CToPrev( pCNode, CRingNodes.prev() );
        fmuVector CToNext( pCNode, CRingNodes.next() );
        fmuVector CToPrev2( pCNode, CRingNodes.prev(2) );
        fmuVector CToNext2( pCNode, CRingNodes.next(2) );
        dFullAngle = CNormal.angle( CToNext, CToPrev );
        dLeftAngle = CNormal.angle( CToNext2, CToPrev );
        dRightAngle = CNormal.angle( CToNext, CToPrev2 );
        adAngleDiffs[ii] = fabs( dFullAngle - kRIGHT_ANGLE );
        adAngleDiffs[ii+5] = fabs( dLeftAngle - kRIGHT_ANGLE );
        adAngleDiffs[ii+10] = fabs( dRightAngle - kRIGHT_ANGLE );
        adAngleDiffs[ii+15] = fabs( dFullAngle - kTRIA_ANGLE );
        //Check for concave angles and mark as "impossible"
        if( dRightAngle < dFullAngle )
            adAngleDiffs[ii+20] = fabs( (dFullAngle-dRightAngle) - kTRIA_ANGLE );
        else
            adAngleDiffs[ii+20] = TWOPI;
        if( dLeftAngle < dFullAngle )
            adAngleDiffs[ii+25] = fabs( (dFullAngle-dLeftAngle) - kTRIA_ANGLE );
        else
            adAngleDiffs[ii+25] = TWOPI;
        CRingNodes.step();
    }
    //Penalize large angle errors. Errors less than 30 degrees for quads (20
    //for trias have no penalty. Errors between 60 and 90 degrees for quads
    //(40 to 60 for trias) have twice the penalty as errors between 30 and 60.
    for( ii = 0; ii < 15; ++ii )
    {
        if( adAngleDiffs[ii] > kQUAD_FACTOR )
            adAngleDiffs[ii] *= pow( 2.0, ( adAngleDiffs[ii] - kQUAD_FACTOR ) );
    }
    for( ii = 15; ii < 30; ++ii )
    {
        if( adAngleDiffs[ii] > kTRIA_FACTOR )
            adAngleDiffs[ii] *= pow( 2.0, ( adAngleDiffs[ii] - kTRIA_FACTOR ) );
    }
    //Sum up angle diffs for all possible quads and their corresponding tria.
    //The pair with the sum closest to zero is the one to use. BestIndex
    //points to the node left out of the quad.
    iBestIndex = -1;
    dBestSum = 1000.0;
    jj = 2;
    kk = 3;
    mm = 9;
    nn = 11;
    ita = 15;
    itb = 21;
    itc = 29;
    for( ii = 0; ii < 5; ++ii )
    {
        dSumDiffs = adAngleDiffs[jj] + adAngleDiffs[kk] + 
            adAngleDiffs[mm] + adAngleDiffs[nn] + 
            adAngleDiffs[ita] + adAngleDiffs[itb] + adAngleDiffs[itc];

        if( dSumDiffs < dBestSum )
        {
            dBestSum = dSumDiffs;
            iBestIndex = ii;
        }
        ++jj;
        if( jj == 5 ) jj = 0;
        ++kk;
        if( kk == 5 ) kk = 0;
        ++mm;
        if( mm == 10 ) mm = 5;
        ++nn;
        if( nn == 15 ) nn = 10;
        ++ita;
        ++itb;
        if( itb == 25 ) itb = 20;
        ++itc;
        if( itc == 30 ) itc = 25;
    }
    //If there is no new position, return
    if( iBestIndex == -1 ) return fmu_kCLEAN_NO_ACTION;
    //If the new position of the quad is the same as what is already there
    //and the single quad has been passed in, then simply return.
    if( iBestIndex == 4 && pCQuad ) return fmu_kCLEAN_NO_ACTION;

    if( pCQuad )
    {
        
		
		// Compute Jac ratio of the would-be elements
		dQ = test_quad_quality( CRingNodes.next()->node_x(),  CRingNodes.next()->node_y(),
                                           CRingNodes.next(2)->node_x(), CRingNodes.next(2)->node_y(),
                                           CRingNodes.next(3)->node_x(), CRingNodes.next(3)->node_y(),
                                           CRingNodes.next(4)->node_x(), CRingNodes.next(4)->node_y());
        if (dQ < 0.1)
           return fmu_kCLEAN_NO_ACTION;
		
		// Compute area of input Tria
		pCTria->nodes( apCTriaNodes );
		dTriaAreaOrig = get_tria_area (apCTriaNodes[0]->node_x(),  apCTriaNodes[0]->node_y(),
                                       apCTriaNodes[1]->node_x(),  apCTriaNodes[1]->node_y(),
                                       apCTriaNodes[2]->node_x(),  apCTriaNodes[2]->node_y());

		// Compute area of the would-be-tria. Is it negative ? Or near zero ?
		dTriaArea = get_tria_area (CRingNodes.get()->node_x(),  CRingNodes.get()->node_y(),
                                   CRingNodes.next()->node_x(), CRingNodes.next()->node_y(),
                                   CRingNodes.next(4)->node_x(), CRingNodes.next(4)->node_y());
		
		// The original tria and the would-be-tria must have the same sense. Their areas must not change signs.
		// Furthermore the tria must not be collapsed.
		if (fabs(dTriaArea) < dAreaTol)
			return fmu_kCLEAN_NO_ACTION;
		else if ((dTriaAreaOrig/dTriaArea) < dAreaTol)
			return fmu_kCLEAN_NO_ACTION;


        if( lDrawCleaning )
        {
            //draw_open();
            //pCTria->draw_me( dp_kCLR_CYAN );
            //pCQuad->draw_me( dp_kCLR_LTBLUE );
            //draw_close();
        }
        delete_this_tria( pCTria );
        delete_this_quad( pCQuad );
        delete_this_edge( pCEdge );
    }

    CRingNodes.find( CRingNodes[iBestIndex ] );

    pCQuad = new elmQClass( CRingNodes.next(), CRingNodes.next(2),
        CRingNodes.next(3), CRingNodes.next(4),
        pCMaster );

    pCTria = new elmTClass( CRingNodes.get(), CRingNodes.next(),
        CRingNodes.next(4), pCMaster );

    smooth_mesh( &CRingNodes );
    if( lDrawCleaning )
    {
        //draw_mesh();
        //pCTria->draw_me( dp_kCLR_CYAN );
        //pCQuad->draw_me( dp_kCLR_LTBLUE );
        //draw_close();
    }

    return fmu_kCLEAN_MESH_CHANGED;
}


void TriaCleanTool::delete_this_tria( elmTClass *pCTria )
{
    if(!pCTria)
        return; // Sanity check - NM
    int iPosition = pCTria->id();
    (*pCMaster->global_trias())[iPosition] = NULL;
    if(iPosition < CTriaList.length())
        CTriaList[iPosition] = NULL;
    delete pCTria;
}

void TriaCleanTool::delete_this_quad( elmQClass *pCQuad )
{
    if(!pCQuad)
        return; // Sanity check - NM
    int iPosition = pCQuad->id();
    (*pCMaster->global_quads())[iPosition] = NULL;
    delete pCQuad;
}

void TriaCleanTool::delete_this_edge( fmuEdge *pCEdge )
{
    if(!pCEdge)
        return; // Sanity check - NM
    int iPosition = pCEdge->id();
    (*pCMaster->global_edges())[iPosition] = NULL;
    delete pCEdge;
}

void TriaCleanTool::delete_this_node( nodClass *pCNode )
{
    if (!pCNode)
        return; // Sanity check - NM
    int iPosition = pCMaster->global_interior()->find( pCNode );
    if( iPosition > -1 )
        (*pCMaster->global_interior())[iPosition] = NULL;
    delete pCNode;
}


void TriaCleanTool::isolate_trias()
{
    int iDebugFlag = 0;
    int ii = 0, jj = 0, kk = 0, mm = 0, nn = 0, iIndex = 0;
    int nMatches = 0;

    bool lFound;
    elmTClassDynArray *pCTriaList = pCMaster->global_trias();
    elmTClass *pCTria = NULL;
    nodClassDynArrayDynArray *pCBoundaries = pCMaster->global_boundary();
    nodClassDynArray *pCNodes = pCMaster->global_interior();
    nodClassDynArray *pCBoundary = NULL; 

    typedef std::vector< nodClassDynArray*, SAM_STL_ALLOC( nodClassDynArray* ) > VEC_nodClassDynArray;
    VEC_nodClassDynArray vecMatchLoops( 200, NULL );
    VEC_int vecIndex( 200, 0 );

    nodClass *pCNode = NULL;
    nodClass *apCTriaNodes[3];
    fmuEdge *pCEdge = NULL;
    int iPosition1, iPosition2;

    for( ii = 0; ii < pCTriaList->length(); ++ii )
    {
        pCTria = (*pCTriaList)[ii];
        if( !pCTria )
            continue;

        pCTria->nodes( apCTriaNodes );

        //Make a boundary loop to separate the tria from the surrounding
        //quads. It is in reverse order as it is for the quads around it,
        //not the tria inside it.
        pCBoundary = new nodClassDynArray;
        (*pCBoundary)[0] = apCTriaNodes[2];
        (*pCBoundary)[1] = apCTriaNodes[1];
        (*pCBoundary)[2] = apCTriaNodes[0];
        pCBoundaries->append( pCBoundary );

        //A hardpoint or boundary node of the tria being analyzed should now
        //be a part of two boundary loops, one that describes an actual
        //surface boundary and one around the new loop. As this process
        //works around the tria, the two occurances may be on the same loop.

        for( jj = 0; jj < 3; ++jj )
        {
            pCNode = apCTriaNodes[jj];
            if( pCNode == NULL )
                continue;

            if( pCNode->boundary() )
            {
                lFound = false;
                nMatches = 0;
                for( kk = 0; kk < pCBoundaries->length() && !lFound; ++kk )
                {
                    pCBoundary = (*pCBoundaries)[kk];
                    if( pCBoundary->length() > vecMatchLoops.size() )
                    {
                        vecMatchLoops.resize( pCBoundary->length() );
                        vecIndex     .resize( pCBoundary->length() );
                    }
                    for( mm = 0; mm < pCBoundary->length() && !lFound; ++mm )
                    {
                        if( (*pCBoundary)[mm] == pCNode )
                        {
                            vecMatchLoops[nMatches] = pCBoundary;
                            vecIndex[nMatches] = mm;

                            ++nMatches;
                            if( nMatches > 200 ) 
                            {
                                lFound = true;
                            }
                        }
                    }
                }
                if(nMatches >= 2)
                {
                    // Get combination of loop locations which give an
                    // Interior Split/merge
                    lFound = false;
                    mm=0;
                    while(mm < nMatches-1 && !lFound)
                    {
                        nn=mm+1;
                        while(nn < nMatches && !lFound)
                        {
                            lFound = check_interior( vecMatchLoops[mm], vecIndex[mm],
                                                     vecMatchLoops[nn], vecIndex[nn] );
                            if(!lFound)
                                nn++;
                        }
                        if(!lFound)
                            mm++;
                    }

                    //If no interior matches don't split or merge
                    if(!lFound)
                        continue;

                    if( vecMatchLoops[mm] == vecMatchLoops[nn] )
                    {
                        // Node is in the same loop... SPLIT
                        pCBoundary = new nodClassDynArray;

                        const int indexEnd = vecIndex[nn];
                        for( kk = vecIndex[mm]; kk < indexEnd; ++kk )
                        {
                            pCBoundary->append( (*(vecMatchLoops[mm]))[kk] );
                        }
                        if( pCBoundary->length() <= 2)
                        {
                            delete pCBoundary;
                            pCBoundary = NULL;
                        }
                        else
                        {
                            compress_loop( pCBoundary );
                            if( pCBoundary->length() <= 2 )
                            {
                                delete pCBoundary;
                                pCBoundary = NULL;
                            }
                            else
                                pCBoundaries->append( pCBoundary );
                        }

                        vecMatchLoops[mm]->remove( vecIndex[mm], vecIndex[nn] );
                        if( vecMatchLoops[mm]->length() <= 2)
                        {
                            pCBoundaries->remove( pCBoundaries->find( vecMatchLoops[mm] ));
                            delete vecMatchLoops[mm];
                            vecMatchLoops[mm] = NULL;
                        }
                    }
                    else
                    {
                        // MERGE
                        // The node is in two separate loops, combine into one loop. 
                        // The node will be in that loop twice. Merge Loop nn into Loop mm
                        iPosition1 = vecIndex[mm];
                        iPosition2 = vecIndex[nn];
                        for( kk = 0; kk < vecMatchLoops[nn]->length(); ++kk )
                        {
                            pCNode = (*(vecMatchLoops[nn]))[iPosition2];
                            vecMatchLoops[mm]->insert( pCNode, iPosition1 );
                            iPosition1++;
                            iPosition2++;
                            if(iPosition2 >= vecMatchLoops[nn]->length())
                                iPosition2 = 0;
                        }
                        pCBoundaries->remove( pCBoundaries->find( vecMatchLoops[nn] ));
                        delete vecMatchLoops[nn];
                        vecMatchLoops[nn] = NULL;
                    }
                }
            }
            if( !apCTriaNodes[0]->boundary() )
            {
                apCTriaNodes[0]->promoted_to_boundary( true );
                apCTriaNodes[0]->boundary( true );
                apCTriaNodes[0]->hardpt( true );
                iIndex = pCNodes->find( apCTriaNodes[0] );
                if(iIndex > -1 )
                    (*pCNodes)[iIndex] = NULL;
            }
            if( !apCTriaNodes[1]->boundary() )
            {
                apCTriaNodes[1]->promoted_to_boundary( true );
                apCTriaNodes[1]->boundary( true );
                apCTriaNodes[1]->hardpt( true );
                iIndex = pCNodes->find( apCTriaNodes[1] );
                if(iIndex > -1 )
                    (*pCNodes)[iIndex] = NULL;
            }
            if( !apCTriaNodes[2]->boundary() )
            {
                apCTriaNodes[2]->promoted_to_boundary( true );
                apCTriaNodes[2]->boundary( true );
                apCTriaNodes[2]->hardpt( true );
                iIndex = pCNodes->find( apCTriaNodes[2] );
                if(iIndex > -1 )
                    (*pCNodes)[iIndex] = NULL;
            }
            pCEdge = apCTriaNodes[0]->shared_edge( apCTriaNodes[1] );
            pCEdge->frozen( true );
            pCEdge = apCTriaNodes[1]->shared_edge( apCTriaNodes[2] );
            pCEdge->frozen( true );
            pCEdge = apCTriaNodes[2]->shared_edge( apCTriaNodes[0] );
            pCEdge->frozen( true );

        }
    }

    if( iDebugFlag )
    {
        //draw_mesh();
        //for( kk = 0; kk < pCBoundaries->length(); ++kk )
        //{
        //pCBoundary = (*pCBoundaries)[kk];
        //for( mm = 0; mm < pCBoundary->length(); ++mm )
        //{
        //pCNode = pCBoundary->get();
        //pCNodeA = pCBoundary->next();
        //pCEdge = pCNode->shared_edge( pCNodeA );
        //if( pCEdge ) pCEdge->draw_me( dp_kCLR_CYAN );
        //}
        //}
        //draw_close();
    }
}

bool TriaCleanTool::check_interior(nodClassDynArray *pCLoop1,
                                   int iIndex1,
                                   nodClassDynArray *pCLoop2,
                                   int iIndex2)
{
    // Added sanity check - NM
    if( pCLoop1 == NULL )
        return false;
    if( pCLoop2 == NULL )
        return false;

    pCLoop1->set(iIndex1);
    fmuVector CPrevVec1(pCLoop1->get(), pCLoop1->prev());
    fmuVector CNextVec1(pCLoop1->get(), pCLoop1->next());
    pCLoop2->set(iIndex2);
    fmuVector CPrevVec2(pCLoop2->get(), pCLoop2->prev());
    fmuVector CNextVec2(pCLoop2->get(), pCLoop2->next());

    // dAngle2 should be "interior" to dAngle1
    // Check if pCLoop1 and iIndex1 has 
    // pCLoop2->next() in its interior path

    double dAngle1 = CNormal.angle( CNextVec1, CPrevVec1);
    double dAngle2 = CNormal.angle( CNextVec2, CPrevVec1);
    if(dAngle1 >= dAngle2)
    {
        //Check if pCLoop2 at iIndex2 has 
        //pCLoop1->next() in its interior path

        dAngle2 = CNormal.angle( CNextVec1, CPrevVec2);
        dAngle1 = CNormal.angle( CNextVec2, CPrevVec2);
        if( dAngle1 >= dAngle2)
            return true;
    }

    return false;
}

void TriaCleanTool::compress_loop(nodClassDynArray *pCLoop)
{
    for( int ii=0; ii < pCLoop->length(); ii++)
    {
        pCLoop->set(ii);
        while(pCLoop->prev()->id() == pCLoop->next()->id()) 
        {
            if( pCLoop->length() <= 2 ) return;
            if( ii == pCLoop->length()-1 )
            {
                pCLoop->remove(ii);
                pCLoop->remove(ii-1);
            }
            else
            {
                pCLoop->remove(ii);
                pCLoop->remove(ii);
            }
            if(ii > 0)ii--;
            pCLoop->set(ii);
        }
    }
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
double TriaCleanTool::test_quad_quality( double dx1, double dy1,
                                         double dx2, double dy2,
                                         double dx3, double dy3,
                                         double dx4, double dy4)
{
    double dQ=-1.0;
    double adQuad[8]={0.0};
    double adDetJac[4] = {0.0, 0.0, 0.0, 0.0}, dMinJac =0.0, dMaxJac=0.0, dArea=0.0;

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

double TriaCleanTool::get_tria_area (double dx1, double dy1,
                                     double dx2, double dy2,
                                     double dx3, double dy3)
{

	double adP0[2] = {0.0}, adP1[2] = {0.0}, adP2[2]={0.0};

	adP0[0] = dx1;
	adP0[1] = dy1;

	adP1[0] = dx2;
	adP1[1] = dy2;

	adP2[0] = dx3;
	adP2[1] = dy3;

	// Find the area
    double dArea = MUTAREA2D(adP0,adP1,adP2);

	return (dArea);
}
