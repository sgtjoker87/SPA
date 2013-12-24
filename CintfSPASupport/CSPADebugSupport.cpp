/*
 * CSPADebugSupport.cpp - implementation for Edge Support class
 *
 * */

#include "CintfSPASupport/CSPADebugSupport.h"

// Initialize the static fields here
CSPADebugSupport*   CSPADebugSupport::m_pzInstance = 0;

IDebugSupport* CSPADebugSupport::createInstance( SPA_TzDebugSupport& zDebugSupport )
{
    CSPADebugSupport::deleteInstance();

    m_pzInstance = new CSPADebugSupport(zDebugSupport);

    return m_pzInstance;
}
        
IDebugSupport* CSPADebugSupport::getInstance()
{
    return m_pzInstance;
}

void CSPADebugSupport::deleteInstance()
{
    if ( m_pzInstance )
        delete m_pzInstance;

    m_pzInstance = 0;
}

CSPADebugSupport::CSPADebugSupport( SPA_TzDebugSupport& zDebugSupport )
                                    : m_zDebugSupport( zDebugSupport )
{
}

SPA_TeDebugLevel CSPADebugSupport::queryDebugLevel(void)
{
    if ( m_zDebugSupport.pxQueryDebugLevel )
        return m_zDebugSupport.pxQueryDebugLevel();

    return SPA_eDEBUG_OFF;
}

void CSPADebugSupport::setDebugLevel(SPA_TeDebugLevel eDebugLevel)
{
    if ( m_zDebugSupport.pxSetDebugLevel )
        m_zDebugSupport.pxSetDebugLevel(eDebugLevel);

    return;
}

bool CSPADebugSupport::debugSurfaceMesher(void)
{
    if ( m_zDebugSupport.pxDebugSurfaceMesher )
        return ( m_zDebugSupport.pxDebugSurfaceMesher() != 0 );

    return false;
}

bool CSPADebugSupport::debug2DMesher(void)
{
    if ( m_zDebugSupport.pxDebug2DMesher )
        return ( m_zDebugSupport.pxDebug2DMesher() != 0 );

    return false;
}

bool CSPADebugSupport::debugCleaner(void)
{
    if ( m_zDebugSupport.pxDebugCleaner )
        return ( m_zDebugSupport.pxDebugCleaner() != 0 );

    return false;
}

bool CSPADebugSupport::debugSmoother(void)
{
    if ( m_zDebugSupport.pxDebugSmoother )
        return ( m_zDebugSupport.pxDebugSmoother() != 0 ) ;

    return false;
}

SPA_TeTimerOption CSPADebugSupport::timerOption(void)
{
    if ( m_zDebugSupport.pxTimerOption )
        return m_zDebugSupport.pxTimerOption();

    return SPA_eTIMER_OFF;
}

int CSPADebugSupport::display3dNodes( int     iEntType,
                                      int     iEntId,
                                      int     iNumNodes, 
                                      int     aiNodeId[], 
                                      double  adXYZ[] )
{
    if ( m_zDebugSupport.pxDisplay3dNodes )
        return m_zDebugSupport.pxDisplay3dNodes( iEntType,
                                                 iEntId,
                                                 iNumNodes,
                                                 aiNodeId,
                                                 adXYZ );

    return 0;
}

int CSPADebugSupport::display3dMesh( int      iEntType,
                                     int      iEntId,
                                     int      iNumNodes, 
                                     int      aiNodeId[], 
                                     double   adXYZ[],
                                     int      iNumElem, 
                                     int      aiElemId[], 
                                     int      aiElemType[],
                                     int      aiElemConn[] )
{
    if ( m_zDebugSupport.pxDisplay3dMesh )
        return m_zDebugSupport.pxDisplay3dMesh( iEntType,
                                                iEntId,
                                                iNumNodes,
                                                aiNodeId,
                                                adXYZ,
                                                iNumElem,
                                                aiElemId,
                                                aiElemType,
                                                aiElemConn );

    return 0;
}

int CSPADebugSupport::display2dNodes( int     iEntType,
                                      int     iEntId,
                                      int     iNumNodes, 
                                      int     aiNodeId[], 
                                      double  adXY[] )
{
    if ( m_zDebugSupport.pxDisplay2dNodes )
        return m_zDebugSupport.pxDisplay2dNodes( iEntType,
                                                 iEntId,
                                                 iNumNodes,
                                                 aiNodeId,
                                                 adXY );

    return 0;
}

int CSPADebugSupport::display2dMesh( int      iEntType,
                                     int      iEntId,
                                     int      iNumNodes, 
                                     int      aiNodeId[], 
                                     double   adXY[],
                                     int      iNumElem, 
                                     int      aiElemId[], 
                                     int      aiElemType[],
                                     int      aiElemConn[] )
{
    if ( m_zDebugSupport.pxDisplay2dMesh )
        return m_zDebugSupport.pxDisplay2dMesh( iEntType,
                                                iEntId,
                                                iNumNodes,
                                                aiNodeId,
                                                adXY,
                                                iNumElem,
                                                aiElemId,
                                                aiElemType,
                                                aiElemConn );

    return 0;
}

int CSPADebugSupport::display2dContour( int      iEntType,
                                        int      iEntId,
                                        int      iNumNodes, 
                                        int      aiNodeId[], 
                                        double   adXY[] )
{
    if ( m_zDebugSupport.pxDisplay2dContour )
        return m_zDebugSupport.pxDisplay2dContour( iEntType,
                                                  iEntId,
                                                  iNumNodes,
                                                  aiNodeId,
                                                  adXY );

    return 0;
}

int CSPADebugSupport::remove2dDisplay()
{
    if ( m_zDebugSupport.pxRemove2dDisplay )
        return m_zDebugSupport.pxRemove2dDisplay();

    return 0;
}

int CSPADebugSupport::remove3dDisplay()
{
    if ( m_zDebugSupport.pxRemove3dDisplay )
        return m_zDebugSupport.pxRemove3dDisplay();

    return 0;
}

/** Meshing status . */
int CSPADebugSupport::status( const char *pcMessage )
{
    if ( m_zDebugSupport.pxStatus )
        return m_zDebugSupport.pxStatus( pcMessage );

    return 0;
}

int CSPADebugSupport::ProcessStatus ( int *piProcessType, 
		                              int *piEntityId, 
						              int *piEntityType,
						              int *piMessageId,
						              double *pdProcessComplete 
						             )
{
	if ( m_zDebugSupport.pxProcessStatus )
        return m_zDebugSupport.pxProcessStatus( piProcessType,
		                                        piEntityId,
												piEntityType,
												piMessageId,
												pdProcessComplete);

    return 0;

}

int CSPADebugSupport::message( const char *pcMessage )
{
    if ( m_zDebugSupport.pxMessage )
        return m_zDebugSupport.pxMessage( pcMessage );

    return 0;
}

int CSPADebugSupport::prompt( const char *pcMessage, int iOption )
{
    if ( m_zDebugSupport.pxPrompt )
        return m_zDebugSupport.pxPrompt( pcMessage );

    return 0;
}

// root file name string (iOpt: =0, file name root; !=0, path name)
const char* CSPADebugSupport::meshFileNameRoot( int iOpt )
{
    if ( m_zDebugSupport.pxMeshFileNameRoot )
        return m_zDebugSupport.pxMeshFileNameRoot( iOpt );

    return 0;
}

const char* CSPADebugSupport::meshSummaryFileName()
{
    if ( m_zDebugSupport.pxMeshSummaryFileName )
        return m_zDebugSupport.pxMeshSummaryFileName();

    return 0;
}

const char* CSPADebugSupport::mappedmeshSummaryFileName()
{
    if ( m_zDebugSupport.pxMappedMeshSummaryFileName )
        return m_zDebugSupport.pxMappedMeshSummaryFileName();

    return 0;
}
