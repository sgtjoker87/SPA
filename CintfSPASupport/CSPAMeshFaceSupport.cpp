/*
 * CSPAMeshSurfaceSupport.cpp - implementation for mesh surface Support class
 *
 */


#include <cstdlib>
#include "CSPAMeshFaceSupport.h"

// Initialize the static fields here
CSPAMeshFaceSupport*   CSPAMeshFaceSupport::m_pzInstance = 0;

IMeshFaceSupport* CSPAMeshFaceSupport::createInstance( SPA_TzMeshFaceSupport& zMeshFaceSupport )
{
    CSPAMeshFaceSupport::deleteInstance();

    m_pzInstance = new CSPAMeshFaceSupport(zMeshFaceSupport);

    return m_pzInstance;
}
        
IMeshFaceSupport* CSPAMeshFaceSupport::getInstance()
{
    return m_pzInstance;
}

void CSPAMeshFaceSupport::deleteInstance()
{
    if ( m_pzInstance )
        delete m_pzInstance;

    m_pzInstance = 0;
}

CSPAMeshFaceSupport::CSPAMeshFaceSupport( SPA_TzMeshFaceSupport& zMeshFaceSupport )
: m_zMeshFaceSupport(zMeshFaceSupport)
{
}

/* mesh definitions (return -1 if not defined) */
int CSPAMeshFaceSupport::queryElementSize(int iFaceId, double* pdElemSize)
{
    if ( m_zMeshFaceSupport.pxQueryElementSize )
        return m_zMeshFaceSupport.pxQueryElementSize( iFaceId, pdElemSize );

	return -1;
}

/* frozen status */
bool CSPAMeshFaceSupport::IsFrozen(int iFaceId)
{
    if ( m_zMeshFaceSupport.pxIsFrozen )
        return m_zMeshFaceSupport.pxIsFrozen( iFaceId ) ? true : false;
	
    return false;
}

/* frozen data access */
int CSPAMeshFaceSupport::queryNumFrozenNodes(int iFaceId)
{
    if (m_zMeshFaceSupport.pxQueryNumFrozenNodes )
       return m_zMeshFaceSupport.pxQueryNumFrozenNodes(iFaceId);
    
    return 0;
}
   
int  CSPAMeshFaceSupport::queryFrozenNodes(int iFaceId, int aiNodes[], double adXYZ[])
{
    if (m_zMeshFaceSupport.pxQueryFrozenNodes )
       return m_zMeshFaceSupport.pxQueryFrozenNodes(iFaceId, aiNodes, adXYZ);

    return 0;
};
 
int CSPAMeshFaceSupport::queryNumFrozenElements(int iFaceId)
{
    if (m_zMeshFaceSupport.pxQueryNumFrozenElements )
       return m_zMeshFaceSupport.pxQueryNumFrozenElements(iFaceId);
    
    return 0;
}
 
int CSPAMeshFaceSupport::queryFrozenElements(int iFaceId, int aiElems[], int aiElemTypes[], 
                                     int aiElemConn[], int *piElemOrder)
{
    if (m_zMeshFaceSupport.pxQueryFrozenElements )
       return m_zMeshFaceSupport.pxQueryFrozenElements(iFaceId, aiElems, aiElemTypes, 
                                                aiElemConn, piElemOrder);
    return 0;
}

/* mesh methods 1: free, 2: free-mapped, 3: mapped */
int CSPAMeshFaceSupport::queryMeshMethod(int iFaceId)
{
    if (m_zMeshFaceSupport.pxQueryMeshMethod )
        return m_zMeshFaceSupport.pxQueryMeshMethod(iFaceId);

	return 1; //free mesh
}

/* Space Strategy 1: PS, 2: MAP, 3: Flattening */
int CSPAMeshFaceSupport::querySpaceStrategy(int iFaceId)
{
    if ( m_zMeshFaceSupport.pxQuerySpaceStrategy)
        return m_zMeshFaceSupport.pxQuerySpaceStrategy(iFaceId);

	return 3;
} 

/* Element Type 4: Quad, 3: tria  */
int CSPAMeshFaceSupport::queryElementType(int iFaceId)
{
    if ( m_zMeshFaceSupport.pxQueryElementType )
        return m_zMeshFaceSupport.pxQueryElementType(iFaceId);

	return -1;
} 

/* =1, linear; =2, parabolic  */
int CSPAMeshFaceSupport::queryElementOrder(int iFaceId)
{
    if ( m_zMeshFaceSupport.pxQueryElementOrder )
        return m_zMeshFaceSupport.pxQueryElementOrder(iFaceId);

	return -1;
} 

int CSPAMeshFaceSupport::queryStraightEdgesOption(int iFaceId)
{
    if ( m_zMeshFaceSupport.pxQueryStraightEdgesOption )
        return m_zMeshFaceSupport.pxQueryStraightEdgesOption(iFaceId);

	return -1;
} 

/* Query the smoother and cleaner options [ON(1) or OFF(0)] */
void CSPAMeshFaceSupport::querySmoothCleanOptions( int iFaceId, int* piSmoothOption, int* piCleanOption )
{
  if ( m_zMeshFaceSupport.pxQuerySmoothCleanOptions )
      m_zMeshFaceSupport.pxQuerySmoothCleanOptions( iFaceId, piSmoothOption, piCleanOption );
  else
  {
      *piSmoothOption = 0;
      *piCleanOption = 0;
  }

  return;
}

// query mid-node quality check options and data
void
CSPAMeshFaceSupport::
queryMidNodeQualityChecks(  int iFaceId, 
                            int*    piCheckJacobian,
                            double* pdMinJacobian,
                            int*    piCheckDeviationRatio,
                            double* pdMaxDeviationRatio,
                            int*    piCheckNastranEPLR 
                         )
{
  if ( m_zMeshFaceSupport.pxQueryMidNodeQualityChecks ) {
    m_zMeshFaceSupport.pxQueryMidNodeQualityChecks( iFaceId, piCheckJacobian, 
        pdMinJacobian, piCheckDeviationRatio, pdMaxDeviationRatio, piCheckNastranEPLR );
  }
  else {
    *piCheckJacobian       = 0;
    *pdMinJacobian         = 0.0;
    *piCheckDeviationRatio = 0;
    *pdMaxDeviationRatio   = 0.0;
    *piCheckNastranEPLR    = 0;
  }
}


/* element size */
int CSPAMeshFaceSupport::queryNumLocalSizePts(int iFaceId)
{
    if ( m_zMeshFaceSupport.pxQueryNumLocalSizePts )
        return m_zMeshFaceSupport.pxQueryNumLocalSizePts(iFaceId);

    return -1;
}

int CSPAMeshFaceSupport::queryLocalSizePts(int iFaceId, double adXYZ[],
                                           double adLocalSize[])
{
    if ( m_zMeshFaceSupport.pxQueryLocalSizePts )
        return m_zMeshFaceSupport.pxQueryLocalSizePts( iFaceId, adXYZ, adLocalSize);

	return -1;
}

/* hard points and hard lines */
int CSPAMeshFaceSupport::queryNumHardPoints(int iFaceId)
{
    if (m_zMeshFaceSupport.pxQueryNumHardPoints )
        return m_zMeshFaceSupport.pxQueryNumHardPoints(iFaceId);

    return -1;
}

int CSPAMeshFaceSupport::queryHardPoints(int iFaceId, int aiHardPtIds[], double adXYZ[])
{
    if ( m_zMeshFaceSupport.pxQueryHardPoints )
        return m_zMeshFaceSupport.pxQueryHardPoints( iFaceId, aiHardPtIds, adXYZ );

	return -1;
}

int CSPAMeshFaceSupport::queryNumHardLines(int iFaceId)
{
    if ( m_zMeshFaceSupport.pxQueryNumHardLines )
        return m_zMeshFaceSupport.pxQueryNumHardLines(iFaceId);

	return 0;
}

int CSPAMeshFaceSupport::queryNumPointsOnHardLines(int iFaceId, int aiHardLineId[], int aiNumPts[])
{
    if ( m_zMeshFaceSupport.pxQueryNumPointsOnHardLines )
        return m_zMeshFaceSupport.pxQueryNumPointsOnHardLines(
                                           iFaceId, aiHardLineId, aiNumPts);

	return 0;
}

int CSPAMeshFaceSupport::queryHardLinePoints(int iFaceId, int iHardLineId, double adXYZ[])
{
    if ( m_zMeshFaceSupport.pxQueryHardLinePoints )
        return m_zMeshFaceSupport.pxQueryHardLinePoints (
                                           iFaceId, iHardLineId, adXYZ );
	return 0;
}



int CSPAMeshFaceSupport::queryNumMinElmOnEdg(int iFaceId, int *piCylCircumType)
{
    int nElm = 0;
	if (piCylCircumType)
		*piCylCircumType = 0;

    if ( m_zMeshFaceSupport.pxQueryNumMinElmOnEdg )
        nElm = m_zMeshFaceSupport.pxQueryNumMinElmOnEdg(iFaceId, piCylCircumType);

    return(nElm);
}

int CSPAMeshFaceSupport::queryNumCorners(int iFaceId)
{
    int nElm = 0;

    if ( m_zMeshFaceSupport.pxQueryNumCorners )
        nElm = m_zMeshFaceSupport.pxQueryNumCorners(iFaceId);

    return(nElm);
}

int CSPAMeshFaceSupport::queryCorners(int iFaceId, double a3DCorners[])
{
    int nElm = 0;

    if ( m_zMeshFaceSupport.pxQueryCorners )
        nElm = m_zMeshFaceSupport.pxQueryCorners(iFaceId, a3DCorners);

    return(nElm);
}

int CSPAMeshFaceSupport::queryElemSizeAlongCylHeight( int iFaceId, double *pdSize )
{
	int nL = 0;
	*pdSize=0.0;

	if (m_zMeshFaceSupport.pxQueryElemSizeAlongCylHeight)
		nL = m_zMeshFaceSupport.pxQueryElemSizeAlongCylHeight(iFaceId, pdSize);

	return(nL);
}

double CSPAMeshFaceSupport::queryElemSizeAlongFilletAxis( int iFaceId )
{
	double dSize=0.0;

	if (m_zMeshFaceSupport.pxQueryElemSizeAlongFilletAxis)
		dSize = m_zMeshFaceSupport.pxQueryElemSizeAlongFilletAxis(iFaceId);

	return(dSize);
}

double CSPAMeshFaceSupport::queryElemSizeAlongCircumference( int iFaceId )
{
	double dSize=0.0;

	if (m_zMeshFaceSupport.pxQueryElemSizeAlongCircumference)
		dSize = m_zMeshFaceSupport.pxQueryElemSizeAlongCircumference(iFaceId);

	return(dSize);
}

int CSPAMeshFaceSupport::queryElemCountAlongCircumference( int iFaceId )
{
	int nCnt=0;

	if (m_zMeshFaceSupport.pxQueryElemCountAlongCircumference)
		nCnt = m_zMeshFaceSupport.pxQueryElemCountAlongCircumference(iFaceId);

	return(nCnt);
}

double CSPAMeshFaceSupport::queryMinimumElementSize (int iFaceId)
{
	double dSize = 0.0;

	if (m_zMeshFaceSupport.pxQueryMinimumElementSize)
		dSize = m_zMeshFaceSupport.pxQueryMinimumElementSize(iFaceId);

	return(dSize);
}

double CSPAMeshFaceSupport::queryAspectRatioLimit (int iFaceId)
{
	double dARLimit = 0.0;
	
	if (m_zMeshFaceSupport.pxQueryAspectRatioLimit)
		dARLimit = m_zMeshFaceSupport.pxQueryAspectRatioLimit(iFaceId);

	return(dARLimit);
}

int CSPAMeshFaceSupport::queryFilletSelectionParams(int iFaceId, double adFilletSelectionParam[])
{
    int iRet=EXIT_FAILURE;

    // Initialise 
    adFilletSelectionParam[0]=0.0;   // Min Radius
    adFilletSelectionParam[1]=1e+30; // Max radius
    adFilletSelectionParam[2]=0.0;   // Min angle
    adFilletSelectionParam[3]=360.0; // Max angle

    if ( m_zMeshFaceSupport.pxQueryFilletSelectionParams )
        iRet = m_zMeshFaceSupport.pxQueryFilletSelectionParams(iFaceId, adFilletSelectionParam);

    return(iRet);
}


/* is this a feature face? */
bool CSPAMeshFaceSupport::isFeatureFace(int iFaceId)
{
    if ( m_zMeshFaceSupport.pxIsFeatureFace )
        return m_zMeshFaceSupport.pxIsFeatureFace( iFaceId ) ? true : false;

    return false;
}

#if 0
int CSPAMeshFaceSupport::queryCylinder(int iFaceId)
{
    int fCylinder = 0;

    if ( m_zMeshFaceSupport.pxQueryCylinder )
        fCylinder = m_zMeshFaceSupport.pxQueryCylinder(iFaceId);

    return(fCylinder);
}
#endif 
