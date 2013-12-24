/*

*/

#include "cleaner/fmuclean.hxx"
#include "cleaner/fmuedgar.hxx"
#include "cleaner/fmuedge.hxx"
#include "cleaner/nodclass.hxx"
#include "cleaner/fmumaster.hxx"
#include "cleaner/fmuvector.hxx"



///

CleanTool::CleanTool() :
    pCMaster( NULL ),
    CNormal(),
    pSMGradeFunc( NULL ),
    lDrawCleaning( false )
{
}


// SDRC FE quad mesh cleaning Member Function Descriptor Block -----
//
// CleanTool::destructor
// class - CleanTool
// date 04-23-97  SDRC/P. Kinney
//
// Description:
//
// destructor for the class.
//
// Access:
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
// ------------------------------------------
///

CleanTool::~CleanTool()
{
    this->pSMGradeFunc = NULL;
}

void CleanTool::SetSizeMap(PtrToSMGradeFunc pSMGradeFunc)
{
    this->pSMGradeFunc = pSMGradeFunc;
}


// SDRC FE quad mesh cleaning Member Function Descriptor Block -----
//
// CleanTool::smooth_mesh
// class - CleanTool
// date 04-23-97  SDRC/P. Kinney
//
// Description:
//
// constructor for the class.
//
// Access:
//
// smooth_mesh( pCInitialList )
//
// Input Parameters:
//
// pCInitialList         nodClassDynArray*      list of nodes to be smoothed
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
// If a node is declared to have moved enough to influence its neighbors,
// those nodes are in turn smoothed. This continues until the node movement
// is too small to matter or until the iteration limit is reached.
//
// Nodes do not have to be in the smooth list twice. Speed analysis
// determined that it is faster to query the node as to whether it is in
// the smooth list instead of querying the list to see if it contains
// the node. The node has a flag for this query. This is not ideal C++,
// but is done only for speed.
//
// ------------------------------------------
///

void CleanTool::smooth_mesh( nodClassDynArray *pCInitialList )
{
    //int iDrawSmooth = 0;
    int ii, jj, kk;
    int iMax = 20;
    bool lNodeMoved;
    nodClass *pCNode, *pCOtherNode;
    nodClassDynArray CSmoothList;
    nodClassDynArray CNextSmoothList;
    fmuEdge *pCEdge;
    fmuEdgeDynArray CEdges;

    CSmoothList = *pCInitialList;

    for( jj = 0; (jj < CSmoothList.length()) ; ++jj )
        if( CSmoothList[jj] )
            CSmoothList[jj]->in_smooth_list(true);

    for( jj = 0; jj < iMax && CSmoothList.length() > 0; ++jj )
    {
        for( ii = 0; ii < CSmoothList.length(); ++ii )
        {
            pCNode = CSmoothList[ii];
            CSmoothList[ii] = NULL;
            if( pCNode)
            {
                pCNode->in_smooth_list(false);

                //            if( iDrawSmooth )
                //{
                //   draw_mesh();
                //   draw_close();
                //   for( kk = ii+1; kk < CSmoothList.length(); ++kk )
                //   {
                //      pCOtherNode = CSmoothList[kk];
                //      if( pCOtherNode ) pCOtherNode->draw_me( dp_kCLR_YELLOW );
                //   }
                //   pCNode->draw_me( dp_kCLR_RED );
                //   for( kk = 0; kk < CNextSmoothList.length(); ++kk )
                //   {
                //      pCOtherNode = CNextSmoothList[kk];
                //      if( pCOtherNode ) pCOtherNode->draw_me( dp_kCLR_GREEN );
                //   }
                //   pCNode->draw_me( dp_kCLR_RED );
                //}

                lNodeMoved = pCNode->smooth_me( &CNormal );
                if( lNodeMoved )
                {
                    //The node moved enough to influence neighbors, resmooth them.
                    CEdges.clear();
                    pCNode->edge_list( CEdges );
                    for( kk = 0; kk < CEdges.length(); ++kk )
                    {
                        pCEdge = CEdges[kk];
                        pCOtherNode = pCEdge->other_node( pCNode );
                        if( !(pCOtherNode->in_smooth_list()) && !pCOtherNode->hardpt() )
                        {
                            CNextSmoothList.append( pCOtherNode );
                            pCOtherNode->in_smooth_list(true);
                        }
                    }
                }
            }
        }

        CSmoothList.absorb( CNextSmoothList );
    }

    for( ii = 0; ii < CSmoothList.length(); ++ii )
    {
        pCNode = CSmoothList[ii];
        if(pCNode)
            pCNode->in_smooth_list(false);
    }
}


// SDRC FE quad mesh cleaning Member Function Descriptor Block -----
//
// CleanTool::one_middle_node
// class - CleanTool
// date 04-23-97  SDRC/P. Kinney
//
// Description:
//
// Create a node that is half way between two existing nodes
//
// Access:
//
// one_middle_node( pCNode1, pCNode2 )
//
// Input Parameters:
//
//  pCNode1         nodClass*         First node
//  pCNode1         nodClass*         Second node
//
// Output Parameters:
//
// none
//
// Return Code:
//
// nodclass*        pointer to new node which is halfway between the
//                  original nodes
//
// Error conditions:
//
// none
//
// Additional comments:
//
// ------------------------------------------
///

nodClass* CleanTool::one_middle_node( nodClass* pCNode1, nodClass* pCNode2 )
{
    fmuVector CDirectionVec( pCNode1, pCNode2 );
    CDirectionVec /= 2.0;
    fmuVector CVec1( pCNode1 );
    fmuVector CNodeCoords = CVec1 + CDirectionVec;
    nodClass *pCNewNode = new nodClass( CNodeCoords, pCMaster );
    return pCNewNode;
}


// SDRC FE quad mesh cleaning Member Function Descriptor Block -----
//
// CleanTool::draw_mesh
// class - CleanTool
// date 04-23-97  SDRC/P. Kinney
//
// Description:
//
// Draw all the edges of the mesh in the pure graphics mode
//
// Access:
//
// draw_mesh
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
// This routine calls draw_open as that has to be after pure graphics
// erase. This routine does not call draw_close so that additional quad
// and node displays can be included in the same graphics update.
//
// ------------------------------------------
///

void CleanTool::draw_mesh()
{
    int  ii;
    fmuEdge *pCEdge;

    //Erase what is there
    //gg3pger();

    draw_open();

    //int iColor = dp_kCLR_WHITE;
    for( ii = 0; ii < pCMaster->global_edges()->length(); ++ii )
    {
        pCEdge = (*pCMaster->global_edges())[ii];
        if( pCEdge)
        {
            //if( pCEdge->frozen() )
            //pCEdge->draw_me( dp_kCLR_YELLOW );
            //else
            //pCEdge->draw_me( iColor );
        }
    }
}


// SDRC FE quad mesh cleaning Member Function Descriptor Block -----
//
// CleanTool::draw_open
// class - CleanTool
// date 04-23-97  SDRC/P. Kinney
//
// Description:
//
// Interface to the routine to open pure graphcs
//
// Access:
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
// Pure graphcs are set up as used by mesh cleaning
//
// ------------------------------------------
///

void CleanTool::draw_open()
{
    //int iConst_0 = 0;
    //int iConst_1 = 1;
    //int iMask = 1;
    //int iRet;
    //CALL_FORTRAN (sgropn) ( &iConst_0, &iConst_1, &iConst_0, &iConst_1, &iMask,
}


// SDRC FE quad mesh cleaning Member Function Descriptor Block -----
//
// CleanTool::draw_close
// class - CleanTool
// date 04-23-97  SDRC/P. Kinney
//
// Description:
//
// close the pure graphcs as used by mesh cleaning and dump the buffer
// to the screen.
//
// Access:
//
// draw_close
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
//
// ------------------------------------------
///

void CleanTool::draw_close()
{
    //int iConst_0 = 0;
    //CALL_FORTRAN (sgrcls) (&iConst_0, &iConst_0);
    //CALL_FORTRAN (lstout) ();   
}
