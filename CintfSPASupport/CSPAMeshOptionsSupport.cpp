/* 


/* Generated by Together */

#include "CintfSPASupport/CSPAMeshOptionsSupport.h"

static SPA_TzMapMeshingOptions m_zMapMeshingOptions;
static SPA_TzGlobalMeshingOptions m_zGlobalMeshingOptions;
static int s_iHighOrderInterpolation = static_cast<int>(SPA_eLINEAR_TRI);
static bool s_bAppHighOrderInterpolationIsRead=false;
static double s_dAbstactionTolerance=0.0;
static double s_dMinElementSizeReduction=0.0;
static double s_dGlobalMEL=0.0;
static double s_dFESF=1.0; // Feature Element Size Factor

static bool   s_bAppVolumePyramidTransitionOnReadFirstTime=true;
static int   s_fVolumePyramidTransitionOn=0;
static bool  s_bAppVolumeTwoElementsThruThicknessOnReadFirstTime=true;
static int   s_fVolumeTwoElementsThruThicknessOn=0;

static bool  s_bAppHighNodeLabelReadFirstTime=true;
static int   s_AppHighNodeLabel=0;

static bool  s_bAppHighElementLabelReadFirstTime=true;
static int   s_AppHighElementLabel=0;

static int   s_fQuadOnly=-1; 
static int   s_fMidNodeUpdateOption=-1;
static int   s_fSmoothTarget=-1; 
static int   s_nHexLayers = -1;
static int   s_iFluidDomainSurfaceMeshOption=0; 
static int   s_fTrgMinElmEdgLen = -1;
static int   s_fFeatureAbstractionOption = 0;

static bool  s_bSPAEnableLoadFacetDataAsFrozenReadFirstTime=true;

// Initialize the static fields here
CSPAMeshOptionsSupport*   CSPAMeshOptionsSupport::m_pzInstance = 0;

IMeshOptionsSupport* CSPAMeshOptionsSupport::createInstance( SPA_TzMeshOptionsSupport& zMeshOptionsSupport )
{
    CSPAMeshOptionsSupport::deleteInstance();

    m_pzInstance = new CSPAMeshOptionsSupport(zMeshOptionsSupport);

    return m_pzInstance;
}
        
IMeshOptionsSupport* CSPAMeshOptionsSupport::getInstance()
{
    return m_pzInstance;
}

void CSPAMeshOptionsSupport::deleteInstance()
{
    if ( m_pzInstance )
        delete m_pzInstance;

    m_pzInstance = 0;
}

// -----------------------------------------------------------------------------
CSPAMeshOptionsSupport::
CSPAMeshOptionsSupport( SPA_TzMeshOptionsSupport& zMeshOptionsSupport )
:
  m_zMeshOptionsSupport( zMeshOptionsSupport )
{
}

// -----------------------------------------------------------------------------
double CSPAMeshOptionsSupport::getElementSize(void)
{
  return ( m_zMeshOptionsSupport.pxGetElementSizeOption ) ? m_zMeshOptionsSupport.pxGetElementSizeOption() : -1.0;
}

// -----------------------------------------------------------------------------
SPA_eElementShape CSPAMeshOptionsSupport::getElementShape(void)
{
  SPA_TeElementShape eShape = ( m_zMeshOptionsSupport.pxGetElementShapeOption ) ? (SPA_eElementShape) m_zMeshOptionsSupport.pxGetElementShapeOption() : SPA_eTRIA;
  return eShape;
}


// -----------------------------------------------------------------------------
SPA_TeElementOrder CSPAMeshOptionsSupport::getElementOrder(void)
{
  SPA_TeElementOrder eOrder = ( m_zMeshOptionsSupport.pxGetElementOrderOption ) ? (SPA_TeElementOrder) m_zMeshOptionsSupport.pxGetElementOrderOption() : SPA_eLINEAR;
  return eOrder;
}


// -----------------------------------------------------------------------------
SPA_TeElementFamily CSPAMeshOptionsSupport::getElementFamily(void)
{
  SPA_TeElementFamily eFamily = ( m_zMeshOptionsSupport.pxGetElementFamilyOption ) ? (SPA_TeElementFamily) m_zMeshOptionsSupport.pxGetElementFamilyOption() : SPA_eTHINSHELL;
  return eFamily;
}

// -----------------------------------------------------------------------------
int CSPAMeshOptionsSupport::getDebugFlagOption(void)
{
  return ( m_zMeshOptionsSupport.pxGetDebugFlag ) ? m_zMeshOptionsSupport.pxGetDebugFlag() : -1;
}

// -----------------------------------------------------------------------------
int CSPAMeshOptionsSupport::getCleanFlagOption(void)
{
  return ( m_zMeshOptionsSupport.pxGetCleanFlagOption ) ? m_zMeshOptionsSupport.pxGetCleanFlagOption() : -1;
}

// -----------------------------------------------------------------------------
int CSPAMeshOptionsSupport::getSmoothFlagOption(void)
{
  return ( m_zMeshOptionsSupport.pxGetSmoothFlagOption ) ? m_zMeshOptionsSupport.pxGetSmoothFlagOption() : -1;
}

// -----------------------------------------------------------------------------
int CSPAMeshOptionsSupport::getStraightEdgesFlagOption(void)
{
  return ( m_zMeshOptionsSupport.pxGetStraightEdgesFlagOption ) ? m_zMeshOptionsSupport.pxGetStraightEdgesFlagOption() : -1;
}

// -----------------------------------------------------------------------------
int CSPAMeshOptionsSupport::get2DMesherOption(void)
{
  return ( m_zMeshOptionsSupport.pxGet2DMesherOption ) ? m_zMeshOptionsSupport.pxGet2DMesherOption() : -1;
}

// -----------------------------------------------------------------------------
int CSPAMeshOptionsSupport::getCleanerOption(void)
{
  return ( m_zMeshOptionsSupport.pxGetCleanerOption ) ? m_zMeshOptionsSupport.pxGetCleanerOption() : -1;
}

// -----------------------------------------------------------------------------
int CSPAMeshOptionsSupport::getMeshingStrategy(void)
{
  return ( m_zMeshOptionsSupport.pxGetMeshingStrategy ) ? m_zMeshOptionsSupport.pxGetMeshingStrategy() : -1;
}

// -----------------------------------------------------------------------------
int CSPAMeshOptionsSupport::getSmootherOption(void)
{
  return ( m_zMeshOptionsSupport.pxGetSmootherOption ) ? m_zMeshOptionsSupport.pxGetSmootherOption() : -1;
}

// -----------------------------------------------------------------------------
int CSPAMeshOptionsSupport::getSpaceDevelopmentOption(void)
{
  return ( m_zMeshOptionsSupport.pxGetSpaceDevelopmentOption ) ? m_zMeshOptionsSupport.pxGetSpaceDevelopmentOption() : -1;
}

// -----------------------------------------------------------------------------
double CSPAMeshOptionsSupport::getMinElementSizeReduction(void)
{
	 return ( m_zMeshOptionsSupport.pxGetMinElementSizeReduction ) ? m_zMeshOptionsSupport.pxGetMinElementSizeReduction() : s_dMinElementSizeReduction;
}

// -----------------------------------------------------------------------------
int CSPAMeshOptionsSupport::getSizeControlOption(void)
{
  return ( m_zMeshOptionsSupport.pxGetSizeControlOption ) ? m_zMeshOptionsSupport.pxGetSizeControlOption() : -1;
}

// -----------------------------------------------------------------------------
int CSPAMeshOptionsSupport::getCheckMinJacobianOption(void)
{
  return ( m_zMeshOptionsSupport.pxGetCheckMinJacobianOption ) ? m_zMeshOptionsSupport.pxGetCheckMinJacobianOption() : -1;
}

// -----------------------------------------------------------------------------
int CSPAMeshOptionsSupport::getCheckMaxDeviationOption(void)
{
  return ( m_zMeshOptionsSupport.pxGetCheckMaxDeviationOption ) ? m_zMeshOptionsSupport.pxGetCheckMaxDeviationOption() : -1;
}

// -----------------------------------------------------------------------------
double CSPAMeshOptionsSupport::getMaxDeviationRatioOption(void)
{
  return ( m_zMeshOptionsSupport.pxGetMaxDeviationRatioOption ) ? m_zMeshOptionsSupport.pxGetMaxDeviationRatioOption() : -1.0;
}

// -----------------------------------------------------------------------------
double CSPAMeshOptionsSupport::getMinJacobianRatioOption(void)
{
  return ( m_zMeshOptionsSupport.pxGetMinJacobianRatioOption ) ? m_zMeshOptionsSupport.pxGetMinJacobianRatioOption() : -1.0;
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
int
CSPAMeshOptionsSupport::getLocalRemeshOption()
{
  return ( m_zMeshOptionsSupport.pxGetLocalRemeshOption ) ? m_zMeshOptionsSupport.pxGetLocalRemeshOption() : -1;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void CSPAMeshOptionsSupport::getQualityMeasureOptions(
                                                       double *pdMinSkew,
                                                       double *pdMinWarp,
                                                       double *pdMinTaper, 
                                                       double *pdMinAspect
                                                     )
{
	if ( m_zMeshOptionsSupport.pxGetQualityMeasureOptions )
	{
		m_zMeshOptionsSupport.pxGetQualityMeasureOptions(
                                                      pdMinSkew,
                                                      pdMinWarp,
                                                      pdMinTaper, 
                                                      pdMinAspect
                                                    );
	}
}

// -----------------------------------------------------------------------------
bool CSPAMeshOptionsSupport::getQualityCheckOption(void)
{
    if ( m_zMeshOptionsSupport.pxGetQualityCheckOption )
        return m_zMeshOptionsSupport.pxGetQualityCheckOption()? true : false;

    return false;
}

// -----------------------------------------------------------------------------
void CSPAMeshOptionsSupport::initializeMapMeshingOptions()
{
    m_zMapMeshingOptions.dMaxCornerAngle        = 135.0;
    m_zMapMeshingOptions.dMinCornerAngle        = 45.0;
    m_zMapMeshingOptions.dLongShortSideVariation= 5.0;
    m_zMapMeshingOptions.fSmooth                = 0;
    m_zMapMeshingOptions.fAlternatingTriangles  = 0;
    m_zMapMeshingOptions.fClamp                 = 0;
    m_zMapMeshingOptions.fSubmap                = 0;
    m_zMapMeshingOptions.fCheckEdgeDistribution = 0;
    m_zMapMeshingOptions.fMapMeshCylinders      = 0;
    m_zMapMeshingOptions.fModifyElmSizOnShortLoops=0;
    m_zMapMeshingOptions.iCurvatureEffect       = 0;
}

void CSPAMeshOptionsSupport::setMapMeshingOptions( SPA_TzMapMeshingOptions zMMO)
{
    m_zMapMeshingOptions.dMaxCornerAngle        = zMMO.dMaxCornerAngle;
    m_zMapMeshingOptions.dMinCornerAngle        = zMMO.dMinCornerAngle;
    m_zMapMeshingOptions.dLongShortSideVariation    = zMMO.dLongShortSideVariation;
    m_zMapMeshingOptions.fSmooth                = zMMO.fSmooth;
    m_zMapMeshingOptions.fAlternatingTriangles  = zMMO.fAlternatingTriangles;
    m_zMapMeshingOptions.fClamp                 = zMMO.fClamp;
    m_zMapMeshingOptions.fSubmap                = zMMO.fSubmap;
    m_zMapMeshingOptions.fCheckEdgeDistribution = zMMO.fCheckEdgeDistribution;
    m_zMapMeshingOptions.fMapMeshCylinders      = zMMO.fMapMeshCylinders;
    m_zMapMeshingOptions.fModifyElmSizOnShortLoops=zMMO.fModifyElmSizOnShortLoops;
    m_zMapMeshingOptions.iCurvatureEffect       = zMMO.iCurvatureEffect;
}

SPA_TzMapMeshingOptions CSPAMeshOptionsSupport::queryMapMeshingOptions()
{
    return(m_zMapMeshingOptions);
}

void CSPAMeshOptionsSupport::mmTurnOnSmoothing()
{
}

void CSPAMeshOptionsSupport::mmTurnOffSmoothing()
{

}
   
void CSPAMeshOptionsSupport::mmTurnOnClamping()
{

}
  
void CSPAMeshOptionsSupport::mmTurnOffClamping()
{

}
    
void CSPAMeshOptionsSupport::mmTurnOnSubmapping()
{

}
   
void CSPAMeshOptionsSupport::mmTurnOffSubmapping()
{

}

void CSPAMeshOptionsSupport::setProcessNeed(SPA_TeMeshingProcessStage eProcessNeed)
{

}
   
SPA_TeMeshingProcessStage CSPAMeshOptionsSupport::getProcessNeed()
{
    return(SPA_eMESH_ALL);
}

int CSPAMeshOptionsSupport::getHighNodeLabel()
{
	// JC - Call back the application 
	//      Cache the data in a static variable -
	if (true == s_bAppHighNodeLabelReadFirstTime)
	{
		if ( m_zMeshOptionsSupport.pxGetHighNodeLabel )
		{
			s_AppHighNodeLabel=m_zMeshOptionsSupport.pxGetHighNodeLabel();
		}
		else
		{
			s_AppHighNodeLabel=0;
		}
		s_bAppHighNodeLabelReadFirstTime=false;
	}

	return s_AppHighNodeLabel;
}

int CSPAMeshOptionsSupport::getHighElementLabel()
{
	// JC - Call back the application 
	//      Cache the data in a static variable -
	if (true == s_bAppHighElementLabelReadFirstTime)
	{
		if ( m_zMeshOptionsSupport.pxGetHighElementLabel )
		{
			s_AppHighElementLabel=m_zMeshOptionsSupport.pxGetHighElementLabel();
		}
		else
		{
			s_AppHighElementLabel=0;
		}
		s_bAppHighElementLabelReadFirstTime=false;
	}

	return s_AppHighElementLabel;
}

void CSPAMeshOptionsSupport::getSplitQuadOption(int *pfSpliQuad, double *pdMaxWarp, double *pdMaxJacobian)
{
    if ( m_zMeshOptionsSupport.pxGetSplitQuadOption )
	{
		m_zMeshOptionsSupport.pxGetSplitQuadOption(pfSpliQuad,pdMaxWarp,pdMaxJacobian);
	}
    return;
}
int CSPAMeshOptionsSupport::getAttemptMapMeshFlagOption(void)
{
    return ( m_zMeshOptionsSupport.pxGetAttemptMapMeshFlagOption ) ? m_zMeshOptionsSupport.pxGetAttemptMapMeshFlagOption() : 0;
}

void CSPAMeshOptionsSupport::setAttemptMapMeshFlagOption(int fAMFlag)
{
}

int CSPAMeshOptionsSupport::getVertexMappingOption(void)
{
    return ( m_zMeshOptionsSupport.pxGetVertexMappingOption ) ? m_zMeshOptionsSupport.pxGetVertexMappingOption() : 0;
}

void CSPAMeshOptionsSupport::setVertexMappingOption(int fVMFlag)
{
	m_zGlobalMeshingOptions.fVertexMapping = fVMFlag;
}

int CSPAMeshOptionsSupport::getFlipDiagonalsOption(void)
{
    return ( m_zMeshOptionsSupport.pxGetFlipDiagonalsOption ) ? m_zMeshOptionsSupport.pxGetFlipDiagonalsOption() : 0;
}

void CSPAMeshOptionsSupport::setFlipDiagonalsOption(int fFDFlag)
{
	m_zGlobalMeshingOptions.fFlipDiagonals = fFDFlag;
}

int CSPAMeshOptionsSupport::previewBoundary()
{
 return ( m_zMeshOptionsSupport.pxPreviewBoundary ) ? m_zMeshOptionsSupport.pxPreviewBoundary() : 0;
}

void CSPAMeshOptionsSupport::setSmallFeatureFactor(double dSFF)
{

}
double CSPAMeshOptionsSupport::querySmallFeatureFactor()
{
	return(m_zMeshOptionsSupport.pxQuerySmallFeatureFactor());
}

int   CSPAMeshOptionsSupport::queryHigherOrderInterpolation()
{
    // JC - Call back the application ... 
    // return(m_zMeshOptionsSupport.pxQueryHigherOrderInterpolation());
    // int iHighOrderInterpolation = static_cast<int>(SPA_eWALTON_MEEK_TRI);

    if (false == s_bAppHighOrderInterpolationIsRead)
    {
        if (m_zMeshOptionsSupport.pxQueryHigherOrderInterpolation)
        {
            s_iHighOrderInterpolation=m_zMeshOptionsSupport.pxQueryHigherOrderInterpolation();
        }
        
        s_bAppHighOrderInterpolationIsRead=true;
    }

    return(s_iHighOrderInterpolation);
}

double CSPAMeshOptionsSupport::queryAbstractionTolerance(void)
{
	return(m_zMeshOptionsSupport.pxQueryAbstractionTolerance());
}

void CSPAMeshOptionsSupport::setAbstractionTolerance(double dAbstractionTolerance)
{
   s_dAbstactionTolerance=dAbstractionTolerance;
   return;
}

double CSPAMeshOptionsSupport::getAbstractionTolerance()
{
    return s_dAbstactionTolerance;
}

void CSPAMeshOptionsSupport::setMinElementSizeReduction(double dMinElementSizeReduction)
{
   s_dMinElementSizeReduction=dMinElementSizeReduction;
   return;
}

void CSPAMeshOptionsSupport::setNumLayers(int nLayers)
{
	s_nHexLayers = nLayers;
}

int CSPAMeshOptionsSupport::queryNumLayers()
{
	if (s_nHexLayers < 0)
		return (m_zMeshOptionsSupport.pxQueryNumHexLayers());
	else
		return(s_nHexLayers);
}

int CSPAMeshOptionsSupport::queryVolumeRemeshingOption()
{
	return(m_zMeshOptionsSupport.pxQueryVolumeRemeshingOption());
}

void CSPAMeshOptionsSupport::setVolumeRemeshingOption(int fVolRemeshing)
{
	m_zGlobalMeshingOptions.fVolRemeshing = fVolRemeshing;
}

int CSPAMeshOptionsSupport::getVolumeRemeshingOption()
{
	return (m_zGlobalMeshingOptions.fVolRemeshing);
}

int CSPAMeshOptionsSupport::queryVolumePyramidTransitionOn()
{
	// JC - Call back the application 
	//      Cache the data in a static variable - 
    if (true == s_bAppVolumePyramidTransitionOnReadFirstTime)
    {
        if (m_zMeshOptionsSupport.pxQueryVolumePyramidTransitionOn)
        {
            s_fVolumePyramidTransitionOn=m_zMeshOptionsSupport.pxQueryVolumePyramidTransitionOn();
        }
        
        s_bAppVolumePyramidTransitionOnReadFirstTime=false;
    }

    return(s_fVolumePyramidTransitionOn);
}

void CSPAMeshOptionsSupport::initAppsCachedMeshOptions(void)
{
	// Static variable that need to work as state machines 
	// Each time the application calls the tet mesher 
	// the init method needs to be called 
	s_bAppVolumePyramidTransitionOnReadFirstTime=true;
	s_bAppVolumeTwoElementsThruThicknessOnReadFirstTime=true;
	s_bAppHighNodeLabelReadFirstTime=true;
	s_bAppHighElementLabelReadFirstTime=true;
    s_bSPAEnableLoadFacetDataAsFrozenReadFirstTime=true;

	s_bAppHighOrderInterpolationIsRead=false;
}

void   CSPAMeshOptionsSupport::setVolumePyramidTransitionOn(void)
{
	// This is a call for the unitest when SPA runs stand alone
	return;
}

int CSPAMeshOptionsSupport::getAllQuadOption(void)
{
	if (s_fQuadOnly < 0)
		return(m_zMeshOptionsSupport.pxQueryAllQuadOption());
	else
		return(s_fQuadOnly);
}

void CSPAMeshOptionsSupport::setAllQuadOption(int fQuadOnly)
{
	s_fQuadOnly=fQuadOnly;
}

void CSPAMeshOptionsSupport::setSmoothTargetOption(int fSmoothTarget)
{
	s_fSmoothTarget=fSmoothTarget;
}

int CSPAMeshOptionsSupport::getSmoothTargetOption(void)
{
	if (s_fSmoothTarget < 0)
		return(m_zMeshOptionsSupport.pxQueryTargetFaceSmoothOption());
	else
		return(s_fSmoothTarget);
}


int CSPAMeshOptionsSupport::getMidNodeUpdateOption(void)
{
	if (s_fMidNodeUpdateOption < 0)
		return(m_zMeshOptionsSupport.pxQueryMidNodeUpdateOption());
	else
		return(s_fMidNodeUpdateOption);
}

void CSPAMeshOptionsSupport::setMidNodeUpdateOption(int fMidNodeUpdateOption)
{
	s_fMidNodeUpdateOption=fMidNodeUpdateOption;
}

int CSPAMeshOptionsSupport::queryVolumeTwoElementsThruThicknessOn()
{
	// JC - Call back the application 
	//      Cache the data in a static variable - 
	if (true == s_bAppVolumeTwoElementsThruThicknessOnReadFirstTime)
	{
        if (m_zMeshOptionsSupport.pxQueryVolumeTwoElementsThruThicknessOn)
        {
            s_fVolumeTwoElementsThruThicknessOn=m_zMeshOptionsSupport.pxQueryVolumeTwoElementsThruThicknessOn();
        }
		s_bAppVolumeTwoElementsThruThicknessOnReadFirstTime=false;
	}

    return(s_fVolumeTwoElementsThruThicknessOn);
}

void   CSPAMeshOptionsSupport::setVolumeTwoElementsThruThicknessOn(void)
{
	// This is a call for the unitest when SPA runs stand alone
	return;
}

void   CSPAMeshOptionsSupport::setFluidDomainSurfaceMeshOption(int iFluidDomainSurfaceMeshOption)
{
	// This is a call for the unitest when SPA runs stand alone
    s_iFluidDomainSurfaceMeshOption=iFluidDomainSurfaceMeshOption;
	return;
}


int CSPAMeshOptionsSupport::queryFluidDomainSurfaceMeshOption()
{
    //   =============================================================================
	//   Maya Surface Wrapping Fluid domain surface mesh option
    //
    //   NO GUARANTEES STAND ALONE WOULD WORK - NEEDS TO BE IN NX.
    // 
    //   Comments : 
    //      1. Trying to replicate the results obtained in NX stand alone might not work 
    //         SFMI_option_ask_HighNodeLabel in NX guarantees that the start node in 
    //         SPA will be higher than the facet node label ie the facet node data 
    //         can be used as frozen data.
    //
    //      2. In SPA stand alone (SPAtest) or GUI there is no such implementation 
    //         a. we need to get the total # facet nodes 
    //         b. With env var turned ON make sure that with start at this high water mark  
    //
    //     search SPA_SMART_LOAD_FACET_DATA_AS_FROZEN_DATA in SPA_MOS_GetHighNodeLabel
    //     search getSPAPartialLoadFacetDataAsFrozenData in MeshOptionsSupportSPA::getHighNodeLabel
    //     For now instead of getting the total # of facet points in the model 
    //     we are using a hard coded value of 10000000
    //
    //     QA results for this env var are @ E:\workdir\PR\surface_wrapping\mayaNoDecimationFrozenFacet
    //   
    //   ============================================================================= 
	if (true == s_bSPAEnableLoadFacetDataAsFrozenReadFirstTime)
	{ 
        if (m_zMeshOptionsSupport.pxQueryFluidDomainSurfaceMeshOption)
        {
            // Get the fluid domain surface mesh option from the UI.
            s_iFluidDomainSurfaceMeshOption=m_zMeshOptionsSupport.pxQueryFluidDomainSurfaceMeshOption();
        }
        s_bSPAEnableLoadFacetDataAsFrozenReadFirstTime=false;
	}

    return(s_iFluidDomainSurfaceMeshOption);
}

int CSPAMeshOptionsSupport::queryTrgMinElmEdgLenOption()
{
		if (m_zMeshOptionsSupport.pxQueryTrgMinElmEdgLenOption)
			return(m_zMeshOptionsSupport.pxQueryTrgMinElmEdgLenOption());
		else
			return (0);
	}

void CSPAMeshOptionsSupport::setTrgMinElmEdgLenOption(int fFlag)
{
	s_fTrgMinElmEdgLen = fFlag;
}

int  CSPAMeshOptionsSupport::queryFeatureAbstractionOption(void)
{
	if (s_fFeatureAbstractionOption <= 0)
		return(m_zMeshOptionsSupport.pxQueryFeatureAbstractionOption());
	else
		return(s_fFeatureAbstractionOption);
}

void  CSPAMeshOptionsSupport::setFeatureAbstractionOption(int fFlag)
{
	s_fFeatureAbstractionOption = fFlag;
}

void  CSPAMeshOptionsSupport::setGlobalMinimumElementLength(double dMEL)
{
	s_dGlobalMEL = dMEL;
}

double CSPAMeshOptionsSupport::queryGlobalMinimumElementLength(void)
{
	if (s_dGlobalMEL <= 0.0)
		return(m_zMeshOptionsSupport.pxQueryGlobalMinimumElementLength());
	else
		return(s_dGlobalMEL);
}

void  CSPAMeshOptionsSupport::setFeatureElemSizeFactor(double dFESF)
{
	s_dFESF = dFESF;
}

double CSPAMeshOptionsSupport::queryFeatureElemSizeFactor(void)
{
	if (s_dFESF <= 0.0)
		return(m_zMeshOptionsSupport.pxQueryFeatureElemSizeFactor());
	else
		return(s_dFESF);
}