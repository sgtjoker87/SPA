/*
 * CSPAEdgeSupport.cpp - implementation for Edge Support class
 *
 */
#include "CSPAEdgeSupport.h"

// Initialize the static fields here
CSPAEdgeSupport*   CSPAEdgeSupport::m_pzInstance = 0;

IEdgeSupport* CSPAEdgeSupport::createInstance( SPA_TzEdgeSupport& zEdgeSupport )
{
    CSPAEdgeSupport::deleteInstance();

    m_pzInstance = new CSPAEdgeSupport(zEdgeSupport);

    return m_pzInstance;
}
        
IEdgeSupport* CSPAEdgeSupport::getInstance()
{
    return m_pzInstance;
}

void CSPAEdgeSupport::deleteInstance()
{
    if ( m_pzInstance )
        delete m_pzInstance;

    m_pzInstance = 0;
}

CSPAEdgeSupport::CSPAEdgeSupport(SPA_TzEdgeSupport& zEdgeSupport)
                                 :m_zEdgeSupport(zEdgeSupport)
{
}

/* Query the start vertex of the input edge */
int CSPAEdgeSupport::queryStartVertex(int iEdgeId)
{
    if ( m_zEdgeSupport.pxQueryStartVertex )
        return m_zEdgeSupport.pxQueryStartVertex(iEdgeId);

    return -1;
}

/* Query the end vertex of the input edge .*/
int CSPAEdgeSupport::queryEndVertex(int iEdgeId)
{
    if ( m_zEdgeSupport.pxQueryEndVertex )
        return m_zEdgeSupport.pxQueryEndVertex(iEdgeId);

    return -1;
}

int CSPAEdgeSupport::queryVertices(int iEdgeId, int aiVertexId[])
{
    if ( m_zEdgeSupport.pxQueryVertices )
        return m_zEdgeSupport.pxQueryVertices(iEdgeId, aiVertexId);

    return -1;
}

double CSPAEdgeSupport::queryArcLength(int iEdgeId )
{
    if ( m_zEdgeSupport.pxQueryArcLength )
        return m_zEdgeSupport.pxQueryArcLength( iEdgeId );

    return -1.0;
}

int CSPAEdgeSupport::queryTangent(int iEdgeId, double dX, double adTang[])
{
    if ( m_zEdgeSupport.pxQueryTangent )
        return m_zEdgeSupport.pxQueryTangent( iEdgeId, dX, adTang );

    return -1;
}

int CSPAEdgeSupport::queryParamLocation(int iEdgeId, double adXYZ[],
                             double* pdParamLoc)
{
    if ( m_zEdgeSupport.pxQueryParamLocation )
        return m_zEdgeSupport.pxQueryParamLocation( iEdgeId, adXYZ, pdParamLoc);

    return -1;
}

int CSPAEdgeSupport::queryXYZLocation(int iEdgeId, double dParamLoc, double adXYZ[])
{
    if ( m_zEdgeSupport.pxQueryXYZLocation )
        return m_zEdgeSupport.pxQueryXYZLocation( iEdgeId, dParamLoc, adXYZ );

    return -1;
}

int CSPAEdgeSupport::queryCurvature(int iEdgeId, double dParamLoc, double* pdCurvature)
{
    if ( m_zEdgeSupport.pxQueryCurvature )
        return m_zEdgeSupport.pxQueryCurvature( iEdgeId, dParamLoc, pdCurvature);

    return -1;
}

int CSPAEdgeSupport::GetNumControlPoints(int iEdgeId)
{
    if ( m_zEdgeSupport.pxGetNumControlPoints )
        return m_zEdgeSupport.pxGetNumControlPoints( iEdgeId );

    return -1;
}

int CSPAEdgeSupport::GetXYZControlPoints(int iEdgeId, int aiPts[], double adXYZ[])
{
    if ( m_zEdgeSupport.pxGetXYZControlPoints )
        return m_zEdgeSupport.pxGetXYZControlPoints( iEdgeId, aiPts, adXYZ );

    return -1;
}

int CSPAEdgeSupport::GetControlPointsParam(int iEdgeId, double adParam[])
{
    if ( m_zEdgeSupport.pxGetControlPointsParam )
        return m_zEdgeSupport.pxGetControlPointsParam( iEdgeId, adParam );

    return -1;
}

int CSPAEdgeSupport::GetLineRadiusCurvatureAtControlPoints(int iEdgeId, double adCurvature[])
{
    if ( m_zEdgeSupport.pxGetLineRadiusCurvatureAtControlPoints )
        return m_zEdgeSupport.pxGetLineRadiusCurvatureAtControlPoints( iEdgeId, adCurvature );

    return -1;
}

int CSPAEdgeSupport::IsSeam(int iEdgeId)
{
	if (m_zEdgeSupport.pxIsSeam)
		return m_zEdgeSupport.pxIsSeam(iEdgeId);

	return -1;
}

double CSPAEdgeSupport::queryFacetEdgeTolerance( int iEdgeId )
{
	if ( m_zEdgeSupport.pxQueryFacetEdgeTolerance )
         return m_zEdgeSupport.pxQueryFacetEdgeTolerance(iEdgeId);

	return 0.0;
}

int CSPAEdgeSupport::queryAppId( int iEdgeId )
{
	if ( m_zEdgeSupport.pxQueryAppId )
         return (m_zEdgeSupport.pxQueryAppId(iEdgeId));

	return (iEdgeId);
}