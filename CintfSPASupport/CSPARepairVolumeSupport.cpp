/*
 * CSAMRepairVolumeSupport.cpp - implementation for Repair Volume Support class
 *
 *
 */

#include "CSAMRepairVolumeSupport.h"

// Initialize the static fields here
CSAMRepairVolumeSupport*   CSAMRepairVolumeSupport::m_pzInstance = 0;

IRepairVolumeSupport* CSAMRepairVolumeSupport::createInstance( sam_TzRepairVolumeSupport& zRepairVolumeSupport )
{
  CSAMRepairVolumeSupport::deleteInstance();

  m_pzInstance = new CSAMRepairVolumeSupport( zRepairVolumeSupport );

  return m_pzInstance;
}
        
IRepairVolumeSupport* CSAMRepairVolumeSupport::getInstance()
{
  return m_pzInstance;
}

void CSAMRepairVolumeSupport::deleteInstance()
{
  if ( m_pzInstance )
      delete m_pzInstance;

  m_pzInstance = 0;
}

CSAMRepairVolumeSupport::CSAMRepairVolumeSupport( sam_TzRepairVolumeSupport& zRepairVolumeSupport )
: m_zRepairVolumeSupport( zRepairVolumeSupport )
{
}


int CSAMRepairVolumeSupport::queryElementSize( int iRepVolId, double* pdElemSize )
{
	if ( m_zRepairVolumeSupport.pxQueryRepairVolumeElementSize )
		return m_zRepairVolumeSupport.pxQueryRepairVolumeElementSize( iRepVolId, pdElemSize );

	return -1;
}

/** 4: Quad, 3: tria, 5:quad-dominant */
int CSAMRepairVolumeSupport::queryElementType( int iRepVolId )
{
	if ( m_zRepairVolumeSupport.pxQueryRepairVolumeElementType )
		return m_zRepairVolumeSupport.pxQueryRepairVolumeElementType( iRepVolId );

	return -1;
} 

// =1, linear; =2, parabolic
int CSAMRepairVolumeSupport::queryElementOrder( int iRepVolId )
{
	if ( m_zRepairVolumeSupport.pxQueryRepairVolumeElementOrder )
		return m_zRepairVolumeSupport.pxQueryRepairVolumeElementOrder( iRepVolId );

	return -1;
}

// Queries the application for the number of layers
int  CSAMRepairVolumeSupport::queryNumLayers (int iRepVolId)
{
	if (m_zRepairVolumeSupport.pxQueryRepairVolumeNumLayers)
		return m_zRepairVolumeSupport.pxQueryRepairVolumeNumLayers (iRepVolId);
	
	return 0;
}


