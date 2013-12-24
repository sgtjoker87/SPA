/*
 * CSPAMeshEdgeSupport.cpp - implementation for mesh edge Support class
 *
*/

#include <cstdlib>
#include "CSPAMeshEdgeSupport.h"

// Initialize the static fields here
CSPAMeshEdgeSupport*   CSPAMeshEdgeSupport::m_pzInstance = 0;

IMeshEdgeSupport* CSPAMeshEdgeSupport::createInstance( SPA_TzMeshEdgeSupport& zMeshEdgeSupport )
{
    CSPAMeshEdgeSupport::deleteInstance();

    m_pzInstance = new CSPAMeshEdgeSupport(zMeshEdgeSupport);

    return m_pzInstance;
}
        
IMeshEdgeSupport* CSPAMeshEdgeSupport::getInstance()
{
    return m_pzInstance;
}

void CSPAMeshEdgeSupport::deleteInstance()
{
    if ( m_pzInstance )
        delete m_pzInstance;

    m_pzInstance = 0;
}

CSPAMeshEdgeSupport::CSPAMeshEdgeSupport(
                     SPA_TzMeshEdgeSupport&  zMeshEdgeSupport )
                     : m_zMeshEdgeSupport(zMeshEdgeSupport)
{
}

/* mesh definitions (return -1 if not defined) */
int CSPAMeshEdgeSupport::queryElementSize(int iEdgeId, double* pdElemSize)
{
	 *pdElemSize = -1.0; //default

     if ( m_zMeshEdgeSupport.pxQueryElementSize )
         return m_zMeshEdgeSupport.pxQueryElementSize( iEdgeId, pdElemSize);

	 return -1;
}

/* frozen status */
bool CSPAMeshEdgeSupport::IsFrozen(int iEdgeId)
{
    if ( m_zMeshEdgeSupport.pxIsFrozen )
        return m_zMeshEdgeSupport.pxIsFrozen(iEdgeId) ? true : false;

    return false;
}


/* local size definition */
int CSPAMeshEdgeSupport::queryNumLocalSizePts(int iEdgeId)
{
    if ( m_zMeshEdgeSupport.pxQueryNumLocalSizePts )
        return m_zMeshEdgeSupport.pxQueryNumLocalSizePts( iEdgeId);

	return 0;
}

int CSPAMeshEdgeSupport::queryLocalSizePts(int iEdgeId, double adXYZ[],
                                           double adLocalSize[])
{
    if ( m_zMeshEdgeSupport.pxQueryLocalSizePts )
        return m_zMeshEdgeSupport.pxQueryLocalSizePts( iEdgeId, adXYZ, adLocalSize );

	return -1;
}

/**hard points definition */
int CSPAMeshEdgeSupport::queryNumHardPoints(int iEdgeId)
{
    if ( m_zMeshEdgeSupport.pxQueryNumHardPoints )
        return m_zMeshEdgeSupport.pxQueryNumHardPoints( iEdgeId );

	return 0;
}

int CSPAMeshEdgeSupport::queryHardPoints(int iEdgeId, int aiHardPtIds[], double adXYZ[])
{
    if ( m_zMeshEdgeSupport.pxQueryHardPoints )
        return m_zMeshEdgeSupport.pxQueryHardPoints( iEdgeId, aiHardPtIds, adXYZ );

	return -1;
}


/* frozen nodes */
int CSPAMeshEdgeSupport::queryNumFrozenNodes(int iEdgeId)
{
    if ( m_zMeshEdgeSupport.pxQueryNumFrozenNodes )
        return m_zMeshEdgeSupport.pxQueryNumFrozenNodes(iEdgeId);

	return 0;
}

int CSPAMeshEdgeSupport::queryFrozenNodes(int iEdgeId, int aiNodeId[], double adXYZ[],
										  double adT[], int *piElemOrder)
{
    if ( m_zMeshEdgeSupport.pxQueryFrozenNodes )
        return m_zMeshEdgeSupport.pxQueryFrozenNodes( iEdgeId, aiNodeId, adXYZ, adT, piElemOrder );

	return 0;
}

int CSPAMeshEdgeSupport::queryElementOrder(int iEdgeId)
{
    return(m_zMeshEdgeSupport.pxQueryElementOrder(iEdgeId));
}

/*element count and biasing definition */
int CSPAMeshEdgeSupport::queryElementCountDef(int iEdgeId, int* piNumElem)
{
    *piNumElem = 0;

    if ( m_zMeshEdgeSupport.pxQueryElementCountDef )
        return m_zMeshEdgeSupport.pxQueryElementCountDef ( iEdgeId, piNumElem );

	return -1;
}

int CSPAMeshEdgeSupport::queryBiasingDef(int iEdgeId, double* pdBias)
{
    *pdBias = 0.0;

    if ( m_zMeshEdgeSupport.pxQueryBiasingDef )
        return m_zMeshEdgeSupport.pxQueryBiasingDef (iEdgeId, pdBias);

	return -1;
}

int CSPAMeshEdgeSupport::queryNumSlaves(int iEdgeId)
{
    if ( m_zMeshEdgeSupport.pxQueryNumSlaves )
        return m_zMeshEdgeSupport.pxQueryNumSlaves( iEdgeId );

    return 0;
}

int CSPAMeshEdgeSupport::querySlaves(int iEdgeId, int aiSlaveEdges[])
{
    if ( m_zMeshEdgeSupport.pxQuerySlaves )
        return m_zMeshEdgeSupport.pxQuerySlaves( iEdgeId, aiSlaveEdges );

    return 0;
}
bool CSPAMeshEdgeSupport::isSlave(int iEdgeId)
{
    if (m_zMeshEdgeSupport.pxIsSlave)
        return ( m_zMeshEdgeSupport.pxIsSlave(iEdgeId) != 0 );
    
    return false;
}

int CSPAMeshEdgeSupport::queryNumEdgeSections( int iEdgeId ) 
{
    if ( m_zMeshEdgeSupport.pxQueryNumEdgeSections )
        return m_zMeshEdgeSupport.pxQueryNumEdgeSections( iEdgeId );

    return 0;
}

int CSPAMeshEdgeSupport::queryEdgeSections( int                 iEdgeId,  
                                            SPA_TzEdgeSection*  apzEdgeSections[] )
{
    if ( m_zMeshEdgeSupport.pxQueryEdgeSections )
        return m_zMeshEdgeSupport.pxQueryEdgeSections( iEdgeId, apzEdgeSections );

    return 0;
}

int CSPAMeshEdgeSupport::queryEdgeDensity( int                 iEdgeId,  
                                           SPA_TzEdgeDensityData *pzEdgeDensityData )
{
    int iRet=EXIT_FAILURE;
    if ( m_zMeshEdgeSupport.pxQueryEdgeDensity )
        return m_zMeshEdgeSupport.pxQueryEdgeDensity( iEdgeId, pzEdgeDensityData );

    return iRet;
}

int CSPAMeshEdgeSupport::isCollapsible( int iEdgeId )
{
	if ( m_zMeshEdgeSupport.pxIsCollapsible )
        return m_zMeshEdgeSupport.pxIsCollapsible(iEdgeId);

	return 0;
}

bool CSPAMeshEdgeSupport::isWeldRow(int iEdgeId)
{
	if ( m_zMeshEdgeSupport.pxIsWeldRow )
        return m_zMeshEdgeSupport.pxIsWeldRow(iEdgeId);

	return 0;
}

int  CSPAMeshEdgeSupport::queryNumWeldRowFaces(int iEdgeId)
{
	if ( m_zMeshEdgeSupport.pxQueryNumWeldRowFaces )
        return m_zMeshEdgeSupport.pxQueryNumWeldRowFaces(iEdgeId);

	return 0;
}

int CSPAMeshEdgeSupport::queryWeldRowFaces(int iEdgeId, int aiFaces[])
{
	if ( m_zMeshEdgeSupport.pxQueryWeldRowFaces )
        return m_zMeshEdgeSupport.pxQueryWeldRowFaces(iEdgeId, aiFaces);

	return 0;
}

int CSPAMeshEdgeSupport::queryWeldRowFaceSides(int iEdgeId, int aiFaceSides[])
{

	if ( m_zMeshEdgeSupport.pxQueryWeldRowFaceSides )
        return m_zMeshEdgeSupport.pxQueryWeldRowFaceSides(iEdgeId, aiFaceSides);

	return 0;
}

int CSPAMeshEdgeSupport::queryWeldRowFaceNumLayers(int iEdgeId, int aiFaceNumLayers[])
{

	if ( m_zMeshEdgeSupport.pxQueryWeldRowFaceNumLayers )
        return m_zMeshEdgeSupport.pxQueryWeldRowFaceNumLayers(iEdgeId, aiFaceNumLayers);

	return 0;
}

int CSPAMeshEdgeSupport::queryWeldRowFaceOffsets (int iEdgeId, double adFaceOffsets[])
{

	if ( m_zMeshEdgeSupport.pxQueryWeldRowFaceOffsets )
        return m_zMeshEdgeSupport.pxQueryWeldRowFaceOffsets(iEdgeId, adFaceOffsets);

	return 0;
}

double CSPAMeshEdgeSupport::queryWeldRowSpacing(int iEdgeId)
{
	if ( m_zMeshEdgeSupport.pxQueryWeldRowSpacing )
        return ( m_zMeshEdgeSupport.pxQueryWeldRowSpacing(iEdgeId));

	return 0.0;
}


int CSPAMeshEdgeSupport::queryMappedLoopNumLayers (int iEdgeId)
{
	if ( m_zMeshEdgeSupport.pxQueryMappedLoopNumLayers )
        return ( m_zMeshEdgeSupport.pxQueryMappedLoopNumLayers(iEdgeId));
	return 0;

}

double CSPAMeshEdgeSupport::queryMappedLoopOffset(int iEdgeId)
{
	if ( m_zMeshEdgeSupport.pxQueryMappedLoopOffset )
        return ( m_zMeshEdgeSupport.pxQueryMappedLoopOffset(iEdgeId));
	return 0.0;
}
 
double CSPAMeshEdgeSupport::queryMappedLoopSpacing(int iEdgeId)
{
	if ( m_zMeshEdgeSupport.pxQueryMappedLoopSpacing )
        return ( m_zMeshEdgeSupport.pxQueryMappedLoopSpacing(iEdgeId));
	return 0.0;
}


bool CSPAMeshEdgeSupport::isRidgeEdge(int iEdgeId)
{

	if ( m_zMeshEdgeSupport.pxIsRidgeEdge )
        return ( m_zMeshEdgeSupport.pxIsRidgeEdge(iEdgeId));
	else
		return false;

}
	
int  CSPAMeshEdgeSupport::queryRidgeEdgeData (	int		iEdgeId, 
												int		*pnLayers, 
												double	*pdTotalLayerDepth, 
												double	*pdSizeAlongLength, 
												int		*pnCountAlongLength)
{
	// Initialize
	*pnLayers = 0;
	*pnCountAlongLength = 0;
	*pdTotalLayerDepth = 0.0;
	*pdSizeAlongLength = 0.0;

	if ( m_zMeshEdgeSupport.pxQueryRidgeEdgeData )
        return ( m_zMeshEdgeSupport.pxQueryRidgeEdgeData(iEdgeId, pnLayers, 
		                                                 pdTotalLayerDepth, 
														 pdSizeAlongLength, 
														 pnCountAlongLength));
	else
		return 0;
}



