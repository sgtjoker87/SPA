/*

*/

#include "cleaner/fmucwrap.hxx"
#include "cleaner/fmuqclean.hxx"
#include "cleaner/fmutclean.hxx"
#include "cleaner/fmuedge.hxx"
#include "cleaner/elmqclass.hxx"
#include "cleaner/elmtclass.hxx"
#include "cleaner/nodclass.hxx"
#include "cleaner/fmumaster.hxx"

#include "errorhandler/ErrorReporter.hxx"


///

int fmuCleanInterface
(
    int                 idFace,
    int                 iNumLoops,
    fmu_TzContourLoop   *azContourLoops,
    bool                fAllowTrias,
    PtrToSMGradeFunc    pSMGradeFunc,
    sam_TzNodeList      *pzINodeList,
    sam_TzElemList      *pzIElemList,
    sam_TzNodeList      *pzONodeList,
    sam_TzElemList      *pzOElemList,
    int                 *nElem
)
{
    int ii, jj, iret = 0;
    int callstatus = fmu_kCLEAN_NO_ACTION;
    int status = fmu_kCLEAN_NO_ACTION;
    int iFirstElem;
    bool do_analyzer = false;

    // This flag controls acess to the 
    // analyzer. To use the analyzer, change
    // this to true and recompile or use the
    // debugger to switch it.
    static const char scFuncName[] = "cleaner/fmucwrap";

    int iRet = 0;

    try
    {
        fmuMaster master;               // Holds all nodes, edges, elements
        elmTClassDynArray *trias;
        elmQClassDynArray *quads;
        nodClassDynArray *nodes;
        fmuEdgeDynArray *edges;
        nodClassDynArrayDynArray *boundaries;
        elmQClass *pQuad;
        elmTClass *pTria;
        nodClass *pNode;
        fmuEdge *pEdge;
        fmuVector CNormal;

        // Initialize the Cleaning Abort Flag
        master.abort(false);

        // Set flag if quad mesh allows tria elements 
        master.triqua_mesh_flag(fAllowTrias);

        master.FaceId( idFace );    // debugging

        // Convert from c structures to C++ classes
        iret = fmuCleanUnpack(  iNumLoops, azContourLoops,
                                pzINodeList, pzIElemList,
                                &iFirstElem, 
                                &master );

master.dump( 1 );

#if JC_DEBUG_ALLOW_TRIAS
        // Check if this is a special case 
        if (fAllowTrias)
        {
            int nTotalInteriorNodes=0, nTotalBoundaryNodes=0, nBoundaries=0;
            int nNodesOnThisBoundary=0;
            nTotalInteriorNodes=master.global_interior()->length();

            //Boundary analysis
            nodClassDynArrayDynArray *boundaries = master.global_boundary();
            nBoundaries = boundaries->length();
            for( ii = 0; ii < nBoundaries; ++ii )
            {
                nodClassDynArray *boundary = (*boundaries)[ii];
                nNodesOnThisBoundary=boundary->length();

                nTotalBoundaryNodes+=nNodesOnThisBoundary;
            }

            if ( (0 == nTotalInteriorNodes) && (nTotalBoundaryNodes <= 5) )
            {
                // Nothing much to do here!
                // We should not be coming to quad cleaning with triangles 
                // All triangles should have been isolated prior to quad cleaning.

                // Copy all the boundary nodes to output
                fmuCopyBoundaryNodes(iNumLoops, azContourLoops, pzONodeList,&iret);

                // Pack it all back into the c structures if any cleaning was done.
                iret = fmuCleanPack( master, iFirstElem,
                    pzINodeList,
                    pzONodeList, pzOElemList,
                    nElem );

            master.abort(true);
            }
        }

#endif 
        if( iret == fmu_kCLEAN_BAD_INPUT ||
            master.abort() )
        {
            ErrorReporter::report_warning(  scFuncName,sam_eBAD_INPUT,
                                            0, 0);
        }
        else if( iret != fmu_kCLEAN_BAD_INPUT &&
                 !master.abort() )
        {
            // Compute the normal of the infinite plane so angles are computed
            // properly. Abort cleaning if a normal is not found.
            iret = fmuPlaneNormal( master, CNormal );

            if( iret != fmu_kCLEAN_ABORT )
            {
                // Construct and execute the mixed element cleaner. If it is not able
                // to elminate all the quads, it creates boundaries around the
                // remaining trias.
                if( master.global_trias()->length() > 0 )
                {
                    TriaCleanTool triaCleaner;
                    callstatus = triaCleaner.CleanMesh( &master, CNormal );
                    if( callstatus == fmu_kCLEAN_MESH_CHANGED )
                        status = fmu_kCLEAN_MESH_CHANGED;
                }

master.dump( 2 );

                // Construct and execute the quad cleaner if nothing fatal in tria clean
                QuadCleanTool quadCleaner;
                if( callstatus == fmu_kCLEAN_NO_ACTION ||
                    callstatus == fmu_kCLEAN_MESH_CHANGED )
                {
                    // Set the sizemap
                    quadCleaner.SetSizeMap(pSMGradeFunc);
                    callstatus = quadCleaner.CleanMesh( &master, CNormal );
                    if( callstatus == fmu_kCLEAN_MESH_CHANGED )
                        status = fmu_kCLEAN_MESH_CHANGED;
                }

master.dump( 3 );

                // Execute the analyzer
                if( do_analyzer )
                {
                    quadCleaner.AnalyzeMesh( &master );
                }

                // Copy all the boundary nodes to output
                fmuCopyBoundaryNodes(iNumLoops, azContourLoops, pzONodeList,&iret);

                // Pack it all back into the c structures if any cleaning was done.
                if( status == fmu_kCLEAN_MESH_CHANGED ||
                    status == fmu_kCLEAN_NO_ACTION)
                {
                    iret = fmuCleanPack( master, iFirstElem,
                                        pzINodeList,
                                        pzONodeList, pzOElemList,
                                        nElem );
                }
                else
                {
                    ErrorReporter::report_warning(scFuncName,sam_eMESH_NOT_CLEANED,0,0);
                }
            }
        }

        //Deallocate all that was allocated
        trias = master.global_trias();
        quads = master.global_quads();
    nodes = master.global_interior();

        //Deletion of quads, trias, edges is straighforward
        for( ii = 0; ii < quads->length(); ++ii)
        {
            pQuad = (*quads)[ii];
            if( pQuad )
                delete pQuad;
        }
        for( ii = 0; ii < trias->length(); ++ii)
        {
            pTria = (*trias)[ii];
            if( pTria )
                delete pTria;
        }
        edges = master.global_edges();
        for( ii = 0; ii < edges->length(); ++ii)
        {
            pEdge = (*edges)[ii];
            if( pEdge )
                delete pEdge;
        }

#if 0
        nodes = master.global_interior();
        //Deletion of interior nodes is also straightforward
        for( ii = 0; ii < nodes->length(); ++ii)
        {
            pNode = (*nodes)[ii];
            if( pNode ) 
            {
                delete pNode;
                pNode = NULL;
            }
        }
#endif

        //Delete nodes through the original boundary structure to avoid
        //duplicate nodes when deleting boundaries that have been modified
        //by trias. It is possible that nodes are in the boundary twice so
        //sort the nodes of the boundary for easy checking so they aren't
        //deleted twice. It is also possible that the same nodes can exist on
        //different loops as in compound section cases thus we need to restore
        //and sort the nodes across all loops to avoid duplicate deletion
        //(H. Xiao, 06/15/01, Q9-13645)

#if 0
        // 06-Oct-2011 <spglen> was leaking memory, don't know why...

        nodClassDynArray tempNodes;
        int nNodes = 0;

        boundaries = master.original_boundary();
        for( jj = 0; jj < boundaries->length(); ++jj )
        {
            nodes = (*boundaries)[jj];
            nNodes = nodes->length();
            for( ii = 0; ii < nNodes; ++ii)
            {         
                pNode = (*nodes)[ii];
                tempNodes.append(pNode);
            }
        }

        nNodes = tempNodes.length();
        if(nNodes > 0)
            tempNodes.sort();

        for( jj = 0; jj < nNodes; ++jj )
        {
            pNode = tempNodes[jj];
            if( jj < nNodes-1 )
                if( pNode == tempNodes[jj+1] ) continue;
            if( pNode ) 
            {
                delete pNode;
                pNode = NULL;
            }
        }

        for( jj = 0; jj < boundaries->length(); ++jj )
        {
            nodes = (*boundaries)[jj];
            delete nodes;
        }

        //Delete the boundary structures, nodes associated with
        //these loops have already been deleted.
        boundaries = master.global_boundary();
        for( jj = 0; jj < boundaries->length(); ++jj )
        {
            nodes = (*boundaries)[jj];
            delete nodes;
        }
#else

        typedef std::vector< nodClass*, SAM_STL_ALLOC( nodClass* ) > VEC_nodClassPtr;
        VEC_nodClassPtr vecN;

        // count all node references

        nodes = master.global_interior();
        int nNodes = nodes->length();

        boundaries = master.original_boundary();
        for( jj = 0; jj < boundaries->length(); ++jj )
        {
            nodes = (*boundaries)[jj];
            nNodes += nodes->length();
        }
        boundaries = master.global_boundary();
        for( jj = 0; jj < boundaries->length(); ++jj )
        {
            nodes = (*boundaries)[jj];
            nNodes += nodes->length();
        }

        vecN.reserve( nNodes );

        // get all node references
        nodes = master.global_interior();
        nNodes = nodes->length();
        for( ii = 0; ii < nNodes; ++ii )
            vecN.push_back( (*nodes)[ii] );

        boundaries = master.original_boundary();
        for( jj = 0; jj < boundaries->length(); ++jj )
        {
            nodes = (*boundaries)[jj];
            nNodes = nodes->length();
            for( ii = 0; ii < nNodes; ++ii )
                vecN.push_back( (*nodes)[ii] );
        }

        boundaries = master.global_boundary();
        for( jj = 0; jj < boundaries->length(); ++jj )
        {
            nodes = (*boundaries)[jj];
            nNodes = nodes->length();
            for( ii = 0; ii < nNodes; ++ii )
                vecN.push_back( (*nodes)[ii] );
        }

        // get a sorted, unique array of nodes
        size_t count1 = std::distance( vecN.begin(), vecN.end() );
        std::sort( vecN.begin(), vecN.end() );
        size_t count2 = std::distance( vecN.begin(), vecN.end() );
        VEC_nodClassPtr::iterator itrEndNew = std::unique( vecN.begin(), vecN.end() );
        size_t count3 = std::distance( vecN.begin(), itrEndNew );

        // delete the nodes
        VEC_nodClassPtr::iterator itr = vecN.begin();
        for( int ic=0; itr != itrEndNew; ++itr, ++ic )
        {
            pNode = *itr;
            delete pNode;
            pNode = *itr = NULL;
        }

        // now, delete the container structures
        boundaries = master.original_boundary();
        for( jj = 0; jj < boundaries->length(); ++jj )
        {
            nodes = (*boundaries)[jj];
            delete nodes;
            nodes = (*boundaries)[jj] = NULL;
        }

        boundaries = master.global_boundary();
        for( jj = 0; jj < boundaries->length(); ++jj )
        {
            nodes = (*boundaries)[jj];
            delete nodes;
            nodes = (*boundaries)[jj] = NULL;
        }

        // note, the vector will go away at scope close...

#endif

    }
    catch( std::bad_alloc ex )
    {
        iRet = -1;
    }
    catch( ... )
    {
        iRet = -2;
    }

    //Abort means cleaning wasn't done, mesher can proceed.
    if( iret == fmu_kCLEAN_ABORT )
        iret = 0;

    return iret;
}


int fmuCleanUnpack
(
    int                 iNumLoops,
    fmu_TzContourLoop   *azContourLoops,
    sam_TzNodeList      *pzINodeList,
    sam_TzElemList      *pzElemList,
    int                 *iFirstElem,
    fmuMaster           *master
)
{
    int iret = 0, iType = 0;
    int iNode = 0;
    double dX2d, dY2d;

    nodClass *pCNode;
    nodClass *nodelist[4];

    NodeTable zNodeLink; 

    fmu_TzLink zLink;

    nodClassDynArray boundary_nodes, interior_nodes;
    nodClassDynArray *nodes = master->global_interior();
    nodClassDynArrayDynArray *boundary = master->global_boundary();
    nodClassDynArrayDynArray *origBoundary = master->original_boundary();

    // Process nodes.
    int iHighLabel = 0;
    int nInterior = (int)(*pzINodeList).size();

    for( int ii = 0; ii < nInterior; ++ii )
    {
        dX2d = (*pzINodeList)[ii].adCoord[0];
        dY2d = (*pzINodeList)[ii].adCoord[1];

        iNode = (*pzINodeList)[ii].iLabel;

        if( iret != 0 )
        {
            printf(" Problems unpacking node %d.\n", iNode );
            return fmu_kCLEAN_BAD_INPUT;
        }

        pCNode = new nodClass( dX2d, dY2d, 0.0 );
        pCNode->id( iNode );
        (*nodes)[ii] = pCNode;
        zLink.iLabel = iNode;
        zLink.pCNode = pCNode;
        fhsNodePut(&zNodeLink, &zLink );

        if( iHighLabel < iNode )
            iHighLabel = iNode;
    }

    // The highest label currently in use
    master->highest_node_id( iHighLabel );

    // Process boundary
    iret = fmuProcessBoundary( iNumLoops, azContourLoops, &zNodeLink, master );
    if( iret == fmu_kCLEAN_BAD_INPUT )
        return iret;

    // Process elements
    *iFirstElem = 99999999;

    int nTotalElems = (int)(*pzElemList).size();
    int iTag = 0;
    for( int i = 0; i < nTotalElems; i++ )
    {
        iTag  = (*pzElemList)[i].iTag;
        iType = (*pzElemList)[i].nNode;

        if( *iFirstElem > iTag )
            *iFirstElem = iTag;

        if( iType == 3 )
        {
            // Process 3 nodes for a tria
            for( int ii = 0; ii < 3; ++ii )
            {
                zLink.iLabel = (*pzElemList)[i].aiNode[ii];

                iret = 0;
                iret = fhsNodeGet( &zNodeLink, &zLink );
                if (iret == 1)
                {
                    printf("Problem unpacking element with node %d\n", zLink.iLabel );
                    return fmu_kCLEAN_BAD_INPUT;
                }
                nodelist[ii] = zLink.pCNode;
            }
            new elmTClass( nodelist[0], nodelist[1], nodelist[2], master );
        }
        else
        {
            // Process 4 nodes for a quad
            for( int ii = 0; ii < 4; ++ii )
            {
                zLink.iLabel = (*pzElemList)[i].aiNode[ii];

                iret = 0;
                iret = fhsNodeGet( &zNodeLink, &zLink );
                if (iret == 1)
                {
                    printf("Problem unpacking element with node %d\n", zLink.iLabel );
                    return fmu_kCLEAN_BAD_INPUT;
                }

                nodelist[ii] = zLink.pCNode;
            }
            new elmQClass( nodelist[0], nodelist[1], nodelist[2], nodelist[3], master );
        }
    }

    return fmu_kCLEAN_NO_ACTION;
}


int fmuProcessBoundary( int             iNumLoops,
                       fmu_TzContourLoop   *azContourLoops,
                       NodeTable           *pzNodeLink,
                       fmuMaster           *pCMaster )
{
    int iDebugFlag = 0;
    int ii, jj, jk, kl;
    int iTag, iret = 0;
    int iLoopLength;

    double dX2d, dY2d;
    double dG;

    nodClassDynArray *pCLoop, *pCOLoop;
    nodClass *pCNode, *pCNode1;

    fmu_TzLink zLink;
    sam_TeNodeType eType = sam_eBDRY_NODE;

    fmuEdge *pCEdge;  

    nodClassDynArray *nodes = pCMaster->global_interior();
    nodClassDynArrayDynArray *boundary = pCMaster->global_boundary();
    nodClassDynArrayDynArray *origBoundary = pCMaster->original_boundary();

    //Make two copies of the boundary structure, one for general use that
    //tria processing can modify to go around triangles, one left unchanged
    //for accurate deletion.
    for( ii = 0; ii < iNumLoops; ++ii )
    {
        pCLoop = new nodClassDynArray;
        (*boundary)[ii] = pCLoop;

        if( iDebugFlag )
            printf("Loop %3d\n", ii );

        for( jj = 0; jj < azContourLoops[ii].nNodes; ++jj )
        {
            iTag = azContourLoops[ii].aiNodes[jj];
            dX2d = azContourLoops[ii].adCoordX[jj];
            dY2d = azContourLoops[ii].adCoordY[jj];
            dG   = azContourLoops[ii].adG[jj];

            if( iret != 0 )
            {
                printf(" Problems unpacking boundary %d.\n", iTag );
                return fmu_kCLEAN_BAD_INPUT;
            }

            zLink.iLabel = iTag;
            //Check to see if this node already exists e.g. weldline case
            iret = 0;
            iret = fhsNodeGet( pzNodeLink, &zLink );
            //if(!pzLink)
            if (iret == 1) // node non-existent
            {
                iret = 0;
                pCNode = new nodClass( dX2d, dY2d, 0.0 );
                pCNode->id( iTag );
                pCNode->hardpt( true );
                pCNode->gradient( dG );
                pCNode->boundary( true );
                pCNode->type( eType );
                zLink.pCNode = pCNode;
                /*iret =*/ fhsNodePut( pzNodeLink, &zLink );
                //if( iret < 0 )
                //   return fmu_kCLEAN_BAD_INPUT;
            }
            else
            {
                pCNode = zLink.pCNode;
            }

            iLoopLength = pCLoop->length();
            if( iLoopLength >= 2)
            {
                if(pCNode == (*pCLoop)[iLoopLength-2])
                {
                    pCNode1 = (*pCLoop)[iLoopLength-1];
                    nodes->append( pCNode1 );
                    pCNode1->boundary( false );
                    pCNode1->is_on_degenerate_loop( true );
                    pCLoop->remove( iLoopLength-1 );
                    if( iDebugFlag )
                    {
                        printf(" Node %d: delete\n", iTag );
                        kl = 0;
                        for( jk = 0; jk < pCLoop->length(); ++jk )
                        {
                            printf("%5d", (*pCLoop)[jk]->id() );
                            ++kl;
                            if( kl > 14 )
                            {
                                kl = 0;
                                printf("\n");
                            }
                        }
                        printf("\n");
                    }
                    continue;
                }
            }
            pCLoop->append( pCNode );
            iLoopLength = pCLoop->length();
            if( iLoopLength >= 2)
            {
                pCEdge = new fmuEdge( (*pCLoop)[iLoopLength-1],
                                      (*pCLoop)[iLoopLength-2],
                                      pCMaster );
                pCEdge->frozen( true );
            }
            if( iDebugFlag )
            {
                printf(" Node %d: add\n", iTag );
                kl = 0;
                for( jk = 0; jk < pCLoop->length(); ++jk )
                {
                    printf("%5d", (*pCLoop)[jk]->id() );
                    ++kl;
                    if( kl > 14 )
                    {
                        kl = 0;
                        printf("\n");
                    }
                }
                printf("\n");
            }
        }

        iLoopLength = pCLoop->length();

        //Check for degeneracies
        if(iLoopLength == 1) // Anchor Node
        {
            if( iDebugFlag )
                printf( "Anchor node\n" );

            pCNode = (*pCLoop)[0];
            nodes->append( pCNode );
            pCNode->boundary( false );
            pCNode->is_on_degenerate_loop( true );
            pCLoop->remove(iLoopLength-1);
        }
        else if(iLoopLength == 2) // This is what is left of a weldline
        {
            if( iDebugFlag )
                printf( "Weldline\n" );

            while( pCLoop->length() > 0 )
            {
                pCNode = (*pCLoop)[pCLoop->length()-1];
                nodes->append( pCNode );
                pCNode->boundary( false );
                pCNode->is_on_degenerate_loop( true );
                pCLoop->remove(pCLoop->length()-1);
            }
        }
        else if( (*pCLoop)[1] == (*pCLoop)[iLoopLength-1] ) //loop start in crack
        {
            if( iDebugFlag )
                printf( "Loop started in crack\n" );

            while( pCLoop->length() > 0 &&
                   (*pCLoop)[1] == (*pCLoop)[pCLoop->length()-1] )
            {
                pCNode = (*pCLoop)[0];
                nodes->append( pCNode );
                pCNode->boundary( false );
                pCNode->is_on_degenerate_loop( true );
                pCLoop->remove(pCLoop->length()-1);
                pCLoop->remove(0);
            }
        }
        else  // Not degenerate. Freeze the edge between last and first node
        {
            pCEdge = new fmuEdge( (*pCLoop)[iLoopLength-1], (*pCLoop)[0], pCMaster );
            pCEdge->frozen( true );
        }

        pCOLoop = new nodClassDynArray( *pCLoop );
        (*origBoundary)[ii] = pCOLoop;
    }

    //Get rid of NULL loops or zero length loops on the boundary
    for(ii = 0; ii < boundary->length(); /* don't increment ii here */ )
    {
        pCLoop = (*boundary)[ii];
        if( pCLoop->length() == 0 )
        {
            // remove compresses the container at current index, 
            // so that "next" becomes "current".  if we remove,
            // don't increment index.
            boundary->remove( boundary->find( pCLoop ));
            delete pCLoop;
        }
        else
            // increment only when we don't remove.
            ++ii;
    }

    return fmu_kCLEAN_NO_ACTION;
}


int fmuPlaneNormal( fmuMaster &pCMaster, fmuVector &CNormal )
{
    //Determine the normal of the infinite plane. Use one quad to do it. If
    //the quad happens to be a chevron or bowtie (both unlikely) try again
    //until an unambiguous normal is found.
    bool lFound = false;
    int ii, jj;
    elmQClass *pCQuad;
    elmTClass *pCTria;
    nodClass *pCThisNode, *pCNextNode, *PCPrevNode, *pCNodeList[4];
    fmuVector CCrossVec;

    for( ii = 0; ii < pCMaster.global_quads()->length() && !lFound; ++ii )
    {
        pCQuad = (*pCMaster.global_quads())[ii];
        pCQuad->nodes( pCNodeList );

        for( jj = 0; jj < 4; ++jj )
        {
            pCThisNode = pCNodeList[jj];
            pCNextNode = pCQuad->next_node( pCThisNode );
            PCPrevNode = pCQuad->prev_node( pCThisNode );
            fmuVector CNextVec( pCThisNode, pCNextNode );
            fmuVector CPrevVec( pCThisNode, PCPrevNode );
            CPrevVec.unitize();
            CNextVec.unitize();
            CCrossVec = CNextVec * CPrevVec;
            //Check for zero length cross, implies colienar vectors. Since
            //units are in meters, length tolerance needs to be pretty small.
            if( CCrossVec.length() < 1.0e-12 )
                break;
            CCrossVec.unitize();
            if( jj == 0 )
                CNormal = CCrossVec;
            else
            {
                if( CNormal % CCrossVec < 0.75 )
                    break;
            }
        }
        if( jj == 4 )
            lFound = true;
    }

    for( ii = 0; ii < pCMaster.global_trias()->length() && !lFound; ++ii )
    {
        pCTria = (*pCMaster.global_trias())[ii];
        pCTria->nodes( pCNodeList );

        pCThisNode = pCNodeList[0];
        pCNextNode = pCTria->next_node( pCThisNode );
        PCPrevNode = pCTria->prev_node( pCThisNode );
        fmuVector CNextVec( pCThisNode, pCNextNode );
        fmuVector CPrevVec( pCThisNode, PCPrevNode );
        CPrevVec.unitize();
        CNextVec.unitize();
        CCrossVec = CNextVec * CPrevVec;
        //Check for zero length cross, implies colienar vectors. Since
        //units are in meters, length squared needs to be pretty small.
        if( CCrossVec.length() < 1.0e-12 )
            break;
        CCrossVec.unitize();
        CNormal = CCrossVec;
        lFound = true;
    }

    //The infinite plane has a normal of 0,0,1 or 0,0,-1. The value from
    //computations may not be exactly either one. If a long way from either
    //(normal in X-Y plane) abort the whole cleaning process.
    if( CNormal.z() < 0.001 && CNormal.z() > -0.001 ) 
        return fmu_kCLEAN_ABORT;

    if( CNormal.z() > 0.0 )
        CNormal.z( 1.0 );
    else
        CNormal.z( -1.0 );

    CNormal.x( 0.0 );
    CNormal.y( 0.0 );

    return fmu_kCLEAN_NO_ACTION;
}

void fmuCopyBoundaryNodes(    int           iNumLoops,
                          fmu_TzContourLoop *azContourLoops,
                          sam_TzNodeList    *pzONodeList,
                          int           *piRet)
{

    int nSize = 0, i = 0, ii;
    sam_TzNode  zNode;


    *piRet = 0;

    for(i = 0; i < iNumLoops; i++)
    {
        nSize = azContourLoops[i].nNodes;

        for(ii=0; ii < nSize; ii++)
        {
            zNode.iLabel = azContourLoops[i].aiNodes[ii];
            zNode.eNodeType = sam_eBDRY_NODE;
            zNode.adCoord[0] = azContourLoops[i].adCoordX[ii];
            zNode.adCoord[1] = azContourLoops[i].adCoordY[ii];
            zNode.adCoord[2] = 0.0;
            (*pzONodeList).push_back(zNode);

        }
    }
}

int fmuCleanPack
(
    fmuMaster       &CMaster, 
    int             iFirstElem,
    sam_TzNodeList  *pzINodeList,
    sam_TzNodeList  *pzONodeList,
    sam_TzElemList  *pzOElemList,
    int             *nElem
)
{
    int ii, jj, kk, iNodeTag = 0, nInterior = 0, iElemTag = 0;
    int iret = 0;
    int iNewId;
    //double dZero = 0.0;
    double dX, dY;
    nodClass *pNode;
    elmQClass *pQuad;
    elmTClass *pTria;
    nodClassDynArray *pCLoop;
    nodClass *node_list[4];

    sam_TzNode zNode;
    sam_TzElem zElem;
    sam_TeNodeType eInsideType = sam_eINSIDE_NODE;

    nodClassDynArray *nodes = CMaster.global_interior();
    nodClassDynArrayDynArray *boundary = CMaster.global_boundary();
    elmTClassDynArray *trias = CMaster.global_trias();
    elmQClassDynArray *quads = CMaster.global_quads();

    if( trias->length() > 0 )
    {
        for( jj = 0; jj < boundary->length(); ++jj )
        {
            pCLoop = (*boundary)[jj];

            for( ii = 0; ii < pCLoop->length(); ++ii)
            {
                pNode = (*pCLoop)[ii];
                if( !pNode )
                    continue;

                // iNodeTag = pNode->id();

                //Check if this node is on a tria and "promoted" from interior
                //to boundary. If so, stick it back in the interior list so
                //that it gets processed correctly.
                if(pNode->promoted_to_boundary())
                {  
                    if(nodes->find(pNode) < 0)
                    {
                        nodes->append( pNode );
                    }
                }
            }
        }
    }

    jj = 0;

    // Find the highest boundary node label
    int iNumBdryNodes=static_cast<int>(pzONodeList->size());
    int iHighBdryNodeLabel=0;
    for (ii=0;ii<iNumBdryNodes;ii++)
    {
        iHighBdryNodeLabel=std::max(iHighBdryNodeLabel,(*pzONodeList)[ii].iLabel);
    }

    //Store interior nodes. Relabel nodes by first using the labels from
    //the input list of labels. When those are gone, use zero so that
    //tqmStorePoint will assign an appropriate label. That label is put
    //back into the node so that elements will refer to the right ones.
    int iHighNodeLabel=iHighBdryNodeLabel;
    for( ii = 0; ii < nodes->length(); ++ii)
    {
        pNode = (*nodes)[ii];
        if( !pNode ) continue;
        //Point is on a weldline or crack so is really a "boundary" node. Skip it
        if( pNode->is_on_degenerate_loop() ) 
        {
            //Before skipping, nullify the address to make sure any further
            //deletion attempt on this array does not lead to a "trying to
            //free already freed memory" situation - Nilanjan Mukherjee 02-23-99

            (*nodes)[ii] = NULL;
            continue;
        }

        if( jj < nInterior )
            iNewId = (*pzINodeList)[jj].iLabel;
        else
            iNewId = 0;

        dX = pNode->node_x();
        dY = pNode->node_y();

        iNodeTag = iHighNodeLabel + 1;
        zNode.eNodeType = eInsideType;

        zNode.adCoord[0] = dX;
        zNode.adCoord[1] = dY;
        zNode.adCoord[2] = 0.0;

        zNode.iLabel     = iNodeTag;  

        (*pzONodeList).push_back(zNode);
        // increment counter 
        iHighNodeLabel++;

        iNewId = iNodeTag;
        pNode->id( iNewId );

        if( iret != 0 )
        {
            printf(" Problems storing node %d.\n", pNode->id() );
            return iret;
        }

        ++jj;
    }

    nInterior = jj;

    //Store quads
    jj = 0;
    for( ii = 0; ii < quads->length(); ++ii )
    {
        pQuad = (*quads)[ii];
      if( !pQuad )
        continue;

        pQuad->nodes( node_list );

        zElem.iPavingZone = -1;

        for( kk = 0; kk < 4; ++kk )
        {
            zElem.aiNode[kk] = node_list[kk]->id();
        }
        iElemTag = (int)(*pzOElemList).size() + 1;
        zElem.iTag = iElemTag;
        zElem.nNode = 4;
        (*pzOElemList).push_back(zElem);
        ++jj;
    }

    //Store trias
    for( ii = 0; ii < trias->length(); ++ii )
    {
        pTria = (*trias)[ii];
        if( !pTria ) continue;

        pTria->nodes( node_list );

        zElem.iPavingZone = -1;

        for( kk = 0; kk < 3; ++kk )
        {
            zElem.aiNode[kk] = node_list[kk]->id();
        }
        iElemTag = (int)(*pzOElemList).size() + 1;
        zElem.iTag = iElemTag;
        zElem.nNode = 3;
        (*pzOElemList).push_back(zElem);
        ++jj;
    }

    *nElem = jj;
    return iret;
}
