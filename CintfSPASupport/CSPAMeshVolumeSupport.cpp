/*
 * CSAMMeshVolumeSupport.cpp - implementation for Volume Support class
 *
 *
 */

#include "CSAMMeshVolumeSupport.h"

// Initialize the static fields here
CSAMMeshVolumeSupport*   CSAMMeshVolumeSupport::m_pzInstance = 0;

IMeshVolumeSupport* CSAMMeshVolumeSupport::createInstance( sam_TzMeshVolumeSupport& zMeshVolumeSupport )
{
  CSAMMeshVolumeSupport::deleteInstance();

  m_pzInstance = new CSAMMeshVolumeSupport( zMeshVolumeSupport );

  return m_pzInstance;
}
        
IMeshVolumeSupport* CSAMMeshVolumeSupport::getInstance()
{
  return m_pzInstance;
}

void CSAMMeshVolumeSupport::deleteInstance()
{
  if ( m_pzInstance )
      delete m_pzInstance;

  m_pzInstance = 0;
}

CSAMMeshVolumeSupport::CSAMMeshVolumeSupport( sam_TzMeshVolumeSupport& zMeshVolumeSupport )
: m_zMeshVolumeSupport( zMeshVolumeSupport )
{
}



/** mesh definitions (return -1 if not defined) */

/** frozen status */
bool CSAMMeshVolumeSupport::IsFrozen(int iVolumeId)
{
	if ( m_zMeshVolumeSupport.pxIsFrozen )
        return m_zMeshVolumeSupport.pxIsFrozen(iVolumeId) ? true : false;

	return false;
}

/* frozen data access */
int CSAMMeshVolumeSupport::queryNumInternalFrozenNodes(int iVolumeId)
{
    if (m_zMeshVolumeSupport.pxQueryNumInternalFrozenNodes )
       return m_zMeshVolumeSupport.pxQueryNumInternalFrozenNodes(iVolumeId);
    
    return 0;
}
   
int  CSAMMeshVolumeSupport::queryInternalFrozenNodes(int iVolumeId, int aiNodes[], double adXYZ[])
{
    if (m_zMeshVolumeSupport.pxQueryInternalFrozenNodes )
       return m_zMeshVolumeSupport.pxQueryInternalFrozenNodes(iVolumeId, aiNodes, adXYZ);

    return 0;
};
 
int CSAMMeshVolumeSupport::queryNumFrozenElements(int iVolumeId)
{
    if (m_zMeshVolumeSupport.pxQueryNumFrozenElements )
       return m_zMeshVolumeSupport.pxQueryNumFrozenElements(iVolumeId);
    
    return 0;
}
 
int CSAMMeshVolumeSupport::queryFrozenElements(int iVolumeId, int aiElems[], int aiElemTypes[], 
                                     int aiElemConn[], int *piElemOrder)
{
    if (m_zMeshVolumeSupport.pxQueryFrozenElements )
       return m_zMeshVolumeSupport.pxQueryFrozenElements(iVolumeId, aiElems, aiElemTypes, 
                                                aiElemConn, piElemOrder);
    return 0;
}



/** mesh methods */
/** 1: free, 2: free-mapped, 3: mapped */
int CSAMMeshVolumeSupport::queryMeshMethod( int iVolumeId )
{
	if ( m_zMeshVolumeSupport.pxQueryMeshMethod )
		return m_zMeshVolumeSupport.pxQueryMeshMethod( iVolumeId );

	return -1;
} 

/* Space Strategy 1: PS, 2: MAP, 3: Flattening */
int CSAMMeshVolumeSupport::querySpaceStrategy( int iVolumeId )
{
  if ( m_zMeshVolumeSupport.pxQuerySpaceStrategy )
    return m_zMeshVolumeSupport.pxQuerySpaceStrategy( iVolumeId );

  return 3;
} 

/** 0: off, 1: on */
void
CSAMMeshVolumeSupport::
querySmoothCleanOptions( int iVolumeId, int* piSmoothOption, int* piCleanOption )
{
  if ( m_zMeshVolumeSupport.pxQuerySmoothCleanOptions )
    m_zMeshVolumeSupport.pxQuerySmoothCleanOptions( iVolumeId, piSmoothOption, piCleanOption );
  else
  {
    *piSmoothOption = 0;
    *piCleanOption = 0;
  }
}

int CSAMMeshVolumeSupport::queryElementSize( int iVolumeId, double* pdElemSize )
{
	if ( m_zMeshVolumeSupport.pxQueryElementSize )
		return m_zMeshVolumeSupport.pxQueryElementSize( iVolumeId, pdElemSize );

	return -1;
}

/** 4: Quad, 3: tria, 5:quad-dominant */
int CSAMMeshVolumeSupport::queryElementType( int iVolumeId )
{
	if ( m_zMeshVolumeSupport.pxQueryElementType )
		return m_zMeshVolumeSupport.pxQueryElementType( iVolumeId );

	return -1;
} 

// =1, linear; =2, parabolic
int CSAMMeshVolumeSupport::queryElementOrder( int iVolumeId )
{
	if ( m_zMeshVolumeSupport.pxQueryElementOrder )
		return m_zMeshVolumeSupport.pxQueryElementOrder( iVolumeId );

	return -1;
}

int  
CSAMMeshVolumeSupport::
queryStraightEdgesOption(int iVolumeId)
{
	if ( m_zMeshVolumeSupport.pxQueryStraightEdgesOption )
		return m_zMeshVolumeSupport.pxQueryStraightEdgesOption(iVolumeId);

	return -1;
}

void
CSAMMeshVolumeSupport::
queryMidNodeQualityChecks(  int iVolumeId,
                            int*    piCheckJacobian,
                            double* pdMinJacobian,
                            int*    piCheckDeviationRatio,
                            double* pdMaxDeviationRatio,
                            int*    piCheckNastranEPLR 
                         )
{
  if ( m_zMeshVolumeSupport.pxQueryMidNodeQualityChecks ) {
    m_zMeshVolumeSupport.pxQueryMidNodeQualityChecks( iVolumeId, piCheckJacobian, 
                                                      pdMinJacobian, piCheckDeviationRatio, 
                                                      pdMaxDeviationRatio, piCheckNastranEPLR );
  }
  else {
    *piCheckJacobian       = 0;
    *pdMinJacobian         = 0.0;
    *piCheckDeviationRatio = 0;
    *pdMaxDeviationRatio   = 0.0;
    *piCheckNastranEPLR = 0;
  }
}

double CSAMMeshVolumeSupport::queryGrowthFactor()
{
	return (m_zMeshVolumeSupport.pxQueryGrowthFactor());
}

