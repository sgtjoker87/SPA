/*
* CSPAElementSupport.cpp - implementation for element Support class
*
*
*/
// 01-31-07     Nilanjan Mukherjee  SCS57344 - Added AspectRatio/TetCollapse appl. queries

#include "CintfSPASupport/CSPAElementSupport.h"

// Initialize the static fields here
CSPAElementSupport*   CSPAElementSupport::m_pzInstance = 0;

IElementSupport* CSPAElementSupport::createInstance( SPA_TzElementSupport   *pzElementSupport )
{
    CSPAElementSupport::deleteInstance();

    if ( pzElementSupport )
        m_pzInstance = new CSPAElementSupport( pzElementSupport );

    return m_pzInstance;
}
IElementSupport* CSPAElementSupport::getInstance()
{
    return m_pzInstance;
}

void CSPAElementSupport::deleteInstance()
{
    if ( m_pzInstance )
        delete m_pzInstance;

    m_pzInstance = 0;
}


CSPAElementSupport::CSPAElementSupport(SPA_TzElementSupport *pzElementSupport) :
    m_pzElementSupport(pzElementSupport)
{
    return;
}

bool CSPAElementSupport::isStraightenMidNodes( const int shape, const double xyz[][3] ) const
{
    return( m_pzElementSupport && m_pzElementSupport->pxIsStraightenMidNodes ) ?
        m_pzElementSupport->pxIsStraightenMidNodes( shape, xyz ) :
        false;
}

bool CSPAElementSupport::isSplitQuad( const double xyz[][3] ) const
{
    return( m_pzElementSupport && m_pzElementSupport->pxIsSplitQuad ) ?
        m_pzElementSupport->pxIsSplitQuad( xyz ) :
        false;
}

bool CSPAElementSupport::isTetraTooFlat( const double xyz[][3] ) const
{
    return( m_pzElementSupport && m_pzElementSupport->pxIsTetraTooFlat ) ?
        m_pzElementSupport->pxIsTetraTooFlat( xyz ) :
        false;
}

bool CSPAElementSupport::doIPassAllEQC( const int shape, const int count, const double xyz[][3] ) const
{
    return( m_pzElementSupport && m_pzElementSupport->pxDoIPassAllEQC ) ?
        m_pzElementSupport->pxDoIPassAllEQC( shape, count, xyz ) :
        false;
}

bool CSPAElementSupport::doesElementFailTetCollapse( int       shape,
                                                     int       count,
                                                     double    xyz[][3],
                                                     double    *pdTetCollapseRatio) const
{
    if( m_pzElementSupport && m_pzElementSupport->pxDoesElementFailTetCollapse)
        return m_pzElementSupport->pxDoesElementFailTetCollapse( shape, count, xyz, pdTetCollapseRatio);

    return false;
}

bool CSPAElementSupport::doesElementFailAspectRatio( int       shape,
                                                     int       count,
                                                     double    xyz[][3],
                                                     double    *pdAspectRatio) const
{
    if( m_pzElementSupport && m_pzElementSupport->pxDoesElementFailAspectRatio)
        return m_pzElementSupport->pxDoesElementFailAspectRatio( shape, count, xyz, pdAspectRatio);

    return false;
}
