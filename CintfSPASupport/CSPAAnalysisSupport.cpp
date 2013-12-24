/*===========================================================================
 * CSPAAnalysisSupport.cpp - class definition for Analysis (ie Adaptivity
 *                               
 *
 * ==========================================================================
*/

#include "CSPAAnalysisSupport.h"

// Initialize the static fields here
CSPAAnalysisSupport*   CSPAAnalysisSupport::m_pzInstance = 0;

IAnalysisSupport* CSPAAnalysisSupport::createInstance( SPA_TzAnalysisSupport& zAnalysisSupport )
{
    CSPAAnalysisSupport::deleteInstance();

    m_pzInstance = new CSPAAnalysisSupport(zAnalysisSupport);

    return m_pzInstance;
}
        
IAnalysisSupport* CSPAAnalysisSupport::getInstance()
{
    return m_pzInstance;
}

void CSPAAnalysisSupport::deleteInstance()
{
    if ( m_pzInstance )
        delete m_pzInstance;

    m_pzInstance = 0;
}

CSPAAnalysisSupport::CSPAAnalysisSupport( SPA_TzAnalysisSupport& zAnalysisSupport )
: m_zAnalysisSupport(zAnalysisSupport)
{
    return;
}


int CSPAAnalysisSupport::queryAnalysisNumIterations()
{ 
	if ( m_zAnalysisSupport.pxQueryAnalysisNumIterations )
        return (m_zAnalysisSupport.pxQueryAnalysisNumIterations()); 

	return -1;
}


int CSPAAnalysisSupport::queryAnalysisNumDataPoints( )
{ 
	if ( m_zAnalysisSupport.pxQueryAnalysisNumDataPoints )
        return (m_zAnalysisSupport.pxQueryAnalysisNumDataPoints()); 

	return -1;
}

int CSPAAnalysisSupport::queryAnalysisDataSizeMap(int aiNodeID[], double adTargetSize[])
{ 
	if ( m_zAnalysisSupport.pxQueryAnalysisDataSizeMap )
        return (m_zAnalysisSupport.pxQueryAnalysisDataSizeMap(aiNodeID, adTargetSize)); 

	return -1;
}

// Max. number of nodes - default -1 when not in use
int CSPAAnalysisSupport::queryAnalysisMaxNodes()
{ 
	if ( m_zAnalysisSupport.pxQueryAnalysisMaxNodes )
        return (m_zAnalysisSupport.pxQueryAnalysisMaxNodes()); 

	return -1;
}


// min elt size - default 0.0 when not in use 
double CSPAAnalysisSupport::queryAnalysisMinElementSize()
{ 
	if ( m_zAnalysisSupport.pxQueryAnalysisMinElementSize )
        return (m_zAnalysisSupport.pxQueryAnalysisMinElementSize()); 

	return 0.0;
}


// User specified adaptivity - local or global remesh desired - default is false
bool CSPAAnalysisSupport::queryAnalysisGlobalReanalyize()
{ 
	if ( m_zAnalysisSupport.pxQueryAnalysisGlobalReanalyize )
        return (m_zAnalysisSupport.pxQueryAnalysisGlobalReanalyize()); 

	return false;
}



