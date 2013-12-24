/*
 * CSPAVertexSupport.cpp - implemetation for Vertex Support class
 *
 *
 */

#include "CSPAVertexSupport.h"

// Initialize the static fields here
CSPAVertexSupport*   CSPAVertexSupport::m_pzInstance = 0;

IVertexSupport* CSPAVertexSupport::createInstance( SPA_TzVertexSupport& zVertexSupport )
{
    CSPAVertexSupport::deleteInstance();

    m_pzInstance = new CSPAVertexSupport(zVertexSupport);

    return m_pzInstance;
}
        
IVertexSupport* CSPAVertexSupport::getInstance()
{
    return m_pzInstance;
}

void CSPAVertexSupport::deleteInstance()
{
    if ( m_pzInstance )
        delete m_pzInstance;

    m_pzInstance = 0;
}

CSPAVertexSupport::CSPAVertexSupport(SPA_TzVertexSupport& zVertexSupport)
                                     :m_zVertexSupport(zVertexSupport)
{
}

/* Get the location of the given vertex */
void CSPAVertexSupport::getLocation(int iVertexId, double adXYZ[3])
{
    if ( m_zVertexSupport.pxGetLocation )
        m_zVertexSupport.pxGetLocation(iVertexId, adXYZ);
}

int CSPAVertexSupport::queryAppId(int iVertexId)
{
	if ( m_zVertexSupport.pxQueryAppId )
		return (m_zVertexSupport.pxQueryAppId( iVertexId));

		return (iVertexId);

}
