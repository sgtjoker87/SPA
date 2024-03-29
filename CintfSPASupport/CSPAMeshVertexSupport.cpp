/*
 * CSPAMeshVertexSupport.cpp - Implementation for Mesh Vertex Support class
 *
 */

/* Generated by Together */


#include "CintfSPASupport/CSPAMeshVertexSupport.h"

// Initialize the static fields here
CSPAMeshVertexSupport*   CSPAMeshVertexSupport::m_pzInstance = 0;

IMeshVertexSupport* CSPAMeshVertexSupport::createInstance( SPA_TzMeshVertexSupport& zMeshVertexSupport )
{
    CSPAMeshVertexSupport::deleteInstance();

    m_pzInstance = new CSPAMeshVertexSupport(zMeshVertexSupport);

    return m_pzInstance;
}
        
IMeshVertexSupport* CSPAMeshVertexSupport::getInstance()
{
    return m_pzInstance;
}

void CSPAMeshVertexSupport::deleteInstance()
{
    if ( m_pzInstance )
        delete m_pzInstance;

    m_pzInstance = 0;
}

CSPAMeshVertexSupport::CSPAMeshVertexSupport( SPA_TzMeshVertexSupport& zMeshVertexSupport )
 : m_zMeshVertexSupport( zMeshVertexSupport )
{
}

/* mesh definitions (return -1 if not defined) */
int CSPAMeshVertexSupport::queryElementSize(int iVertexId, double* pdElemSize)
{
	if ( m_zMeshVertexSupport.pxQueryElementSize )
        return m_zMeshVertexSupport.pxQueryElementSize(iVertexId, pdElemSize);
    else
        *pdElemSize = -1;

    return 0;
}

/* frozen status */
bool CSPAMeshVertexSupport::IsFrozen(int iVertexId)
{
    if ( m_zMeshVertexSupport.pxIsFrozen )
        return m_zMeshVertexSupport.pxIsFrozen(iVertexId) ? true : false;
    
    return false;
}

/* frozen node */
int CSPAMeshVertexSupport::queryFrozenNode(int iVertexId, int * piNodeId )
{
	if ( m_zMeshVertexSupport.pxQueryFrozenNode )
        return m_zMeshVertexSupport.pxQueryFrozenNode(
                                   iVertexId, piNodeId );

    return -1;
}
