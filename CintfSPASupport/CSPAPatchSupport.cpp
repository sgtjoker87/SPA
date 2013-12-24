/*===========================================================================
 * CSPAPatchSupport.cpp - class implementation for Patch mesh surface/edge
 *                                   Support class
 *
*/

#include "CSPAPatchSupport.h"
#include "mesh3d/SPASurfaceAdaptiveLocalRemesh.hxx"

// Initialize the static fields here
CSPAPatchSupport*   CSPAPatchSupport::m_pzInstance = 0;

IPatchSupport* CSPAPatchSupport::createInstance( SPA_TzPatchSupport& zPatchSupport )
{
    CSPAPatchSupport::deleteInstance();

    m_pzInstance = new CSPAPatchSupport(zPatchSupport);

    return m_pzInstance;
}
        
IPatchSupport* CSPAPatchSupport::getInstance()
{
    return m_pzInstance;
}

void CSPAPatchSupport::deleteInstance()
{
    if ( m_pzInstance )
        delete m_pzInstance;

    m_pzInstance = 0;
}

CSPAPatchSupport::CSPAPatchSupport( SPA_TzPatchSupport& zPatchSupport )
: m_zPatchSupport(zPatchSupport)
{
}


int CSPAPatchSupport::queryPatchElementSize( double* pdPatchElemSize, double *pdPatchElementSizeFactor)
{

    SPASurfaceAdaptiveLocalRemesh* pzAdaptiveLocalPatch = SPASurfaceAdaptiveLocalRemesh::getInstance();
    if( pzAdaptiveLocalPatch )
        return pzAdaptiveLocalPatch->queryPatchElementSize(pdPatchElemSize ,pdPatchElementSizeFactor);
    if ( m_zPatchSupport.pxQueryPatchElementSize )
        return (m_zPatchSupport.pxQueryPatchElementSize( pdPatchElemSize ,pdPatchElementSizeFactor));

	return -1;
}

int CSPAPatchSupport::queryPatchNumElements( int iFaceId )
{ 
    SPASurfaceAdaptiveLocalRemesh* pzAdaptiveLocalPatch = SPASurfaceAdaptiveLocalRemesh::getInstance();
    if( pzAdaptiveLocalPatch )
        return (pzAdaptiveLocalPatch->queryPatchNumElements(iFaceId));    
    else if ( m_zPatchSupport.pxQueryPatchNumElements )
        return (m_zPatchSupport.pxQueryPatchNumElements(iFaceId)); 

	return -1;
}

int CSPAPatchSupport::queryPatchElements( int iFaceId, int aiElems[])
{ 
	SPASurfaceAdaptiveLocalRemesh* pzAdaptiveLocalPatch = SPASurfaceAdaptiveLocalRemesh::getInstance();
    if( pzAdaptiveLocalPatch )
        return (pzAdaptiveLocalPatch->queryPatchElements(iFaceId, aiElems));    
    else if ( m_zPatchSupport.pxQueryPatchElements )
        return (m_zPatchSupport.pxQueryPatchElements(iFaceId, aiElems)); 

	return -1;
}


int CSPAPatchSupport::queryPatchTransitionFactor(double* pdTransitionFactor )
{
    SPASurfaceAdaptiveLocalRemesh* pzAdaptiveLocalPatch = SPASurfaceAdaptiveLocalRemesh::getInstance();
    if( pzAdaptiveLocalPatch )
    {
        return pzAdaptiveLocalPatch->queryPatchTransitionFactor(pdTransitionFactor);
    }   
	else if ( m_zPatchSupport.pxQueryPatchTransitionFactor )
        return (m_zPatchSupport.pxQueryPatchTransitionFactor( pdTransitionFactor)); 

	return -1;
}

/* frozen element edges */
int CSPAPatchSupport::queryPatchNumFrozenElementEdges(int iEdgeId)
{
    SPASurfaceAdaptiveLocalRemesh* pzAdaptiveLocalPatch = SPASurfaceAdaptiveLocalRemesh::getInstance();
    if( pzAdaptiveLocalPatch )
    {
        return pzAdaptiveLocalPatch->queryPatchNumFrozenElementEdges(iEdgeId);
    }
    else if ( m_zPatchSupport.pxQueryPatchNumFrozenElementEdges )
        return (m_zPatchSupport.pxQueryPatchNumFrozenElementEdges(iEdgeId));

	return 0;
}

int CSPAPatchSupport::queryPatchFrozenElementEdges(int iEdgeId, int aiElemEdges[])
{
    SPASurfaceAdaptiveLocalRemesh* pzAdaptiveLocalPatch = SPASurfaceAdaptiveLocalRemesh::getInstance();
    if( pzAdaptiveLocalPatch )
    {
        return pzAdaptiveLocalPatch->queryPatchFrozenElementEdges(iEdgeId, aiElemEdges );
    }
    if ( m_zPatchSupport.pxQueryPatchFrozenElementEdges )
        return (m_zPatchSupport.pxQueryPatchFrozenElementEdges( iEdgeId, aiElemEdges ));

	return 0;
}

